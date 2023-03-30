


####----DRPPM PATH SURVEIOR Post-Processing LINCS L1000 Drug Prioritization----####


## Read in libraries
library(readr)
library(dplyr)
library(ggplot2)
library(stringr)


## Set "Step4_LINCSL1000_Drug_Prioritization" as working directory


##--Post-process drug prioritization for full KMT2A Subset--##

## File paths
LINCS_Annotation_File <- "Input/LINCS_compoundinfo_beta.txt"
CoxH_Ranked_Output_File <- "Input/AML_TARGET_KMT2A_Fusion_EFS_LINCS_DN_coxh_ranked.txt"


## Read in files
CoxH_KMT2A <- as.data.frame(read_delim(CoxH_Ranked_Output_File, delim = '\t', col_names = T))
anno <- as.data.frame(read_delim(LINCS_Annotation_File, delim = '\t', col_names = T))


## Subset pipeline results for significant and high-risk results
# The previous version of the pipeline this was run on associated hazard ratios < 1 to be High Risk
# The previous version also labels the hazard ratio column as "estimate"
CoxH_KMT2A <- CoxH_KMT2A %>%
  relocate(variable,estimate,ci,p.value)
CoxH_KMT2A_HR <- CoxH_KMT2A[which(CoxH_KMT2A$estimate <= 1),]
CoxH_KMT2A_HR$p.value <- as.numeric(CoxH_KMT2A_HR$p.value)
CoxH_KMT2A_HR_Sig <- CoxH_KMT2A_HR[which(is.na(CoxH_KMT2A_HR$p.value) | CoxH_KMT2A_HR$p.value <= 0.05),]
CoxH_KMT2A_HR_Sig$p.value <- as.character(CoxH_KMT2A_HR_Sig$p.value)
CoxH_KMT2A_HR_Sig[is.na(CoxH_KMT2A_HR_Sig$p.value),4] <- "<0.001"


## Annotate Gene Sets with LINCS L1000 Data
CoxH_KMT2A_HR_Sig_Anno <- CoxH_KMT2A_HR_Sig %>%
  select(variable,estimate,ci,p.value)
var_list <- CoxH_KMT2A_HR_Sig_Anno$variable
pert_id <- c()
for (i in var_list) {
  first <- strsplit(i, ".", fixed = T)[[1]][1]
  second <- strsplit(i, ".", fixed = T)[[1]][2]
  id <- paste(first,"-",second,sep = "")
  pert_id <- c(pert_id,id)
}
CoxH_KMT2A_HR_Sig_Anno$pert_id <- pert_id
CoxH_KMT2A_HR_Sig_Anno <- merge(CoxH_KMT2A_HR_Sig_Anno,anno[,c(1,2,4)], by = "pert_id",all.x = T)
CoxH_KMT2A_HR_Sig_Anno <- unique(CoxH_KMT2A_HR_Sig_Anno)

##--Frquency of MOA and Cmap names--##

## Get frequency of unique mechanisms of action (MOA)
moa_freq <- as.data.frame(table(CoxH_KMT2A_HR_Sig_Anno$moa))
moa_freq <- moa_freq[order(moa_freq$Freq,decreasing = T),]
rownames(moa_freq) <- 1:nrow(moa_freq)
colnames(moa_freq) <- c("MOA","Frequency")

## Get frequency of unique Cmap names
cmap_freq <- as.data.frame(table(CoxH_KMT2A_HR_Sig_Anno$cmap_name))
cmap_freq <- cmap_freq[order(cmap_freq$Freq,decreasing = T),]
rownames(cmap_freq) <- 1:nrow(cmap_freq)
colnames(cmap_freq) <- c("CMap_Name","Frequency")

## Get frequency of unique ATPase inhibitor Cmap names
CoxH_KMT2A_HR_Sig_Anno_ATPase <- CoxH_KMT2A_HR_Sig_Anno[which(CoxH_KMT2A_HR_Sig_Anno$moa == "ATPase inhibitor"),]
cmapFreq_ATPaseInhib <- as.data.frame(table(CoxH_KMT2A_HR_Sig_Anno_ATPase$cmap_name))
cmapFreq_ATPaseInhib <- cmapFreq_ATPaseInhib[order(-cmapFreq_ATPaseInhib$Freq),]
colnames(cmapFreq_ATPaseInhib) <- c("CMap_Name","Frequency")
write_delim(cmapFreq_ATPaseInhib, "Output/TARGET_AML_KMT2A_EFS_LINCS_ATPaseHits_Frequency.txt", delim = '\t')


## Plot Histograms of MOA and Cmap frequency
MOA_Freq <- ggplot(data = moa_freq[c(1:15),], aes(x = reorder(MOA, Frequency), y = Frequency)) +
  geom_bar(stat = "identity", color='skyblue',fill='steelblue') +
  theme_minimal() +
  coord_flip() +
  xlab("Mechanism of Action (MOA)") +
  ggtitle("MOA Frquency of Appearance in LINCS L1000 High Risk Gene Sets")
MOA_Freq

Cmap_Freq <- ggplot(data = cmap_freq[c(1:15),], aes(x = reorder(CMap_Name, Frequency), y = Frequency)) +
  geom_bar(stat = "identity", color='skyblue',fill='steelblue') +
  theme_minimal() +
  coord_flip() +
  xlab("CMap Name") +
ggtitle("CMap Name Frquency of Appearance in LINCS L1000 High Risk Gene Sets")
Cmap_Freq

ATPase_Freq <- ggplot(data = cmapFreq_ATPaseInhib, aes(x = reorder(CMap_Name, Frequency), y = Frequency)) +
  geom_bar(stat = "identity", color='skyblue',fill='steelblue') +
  theme_minimal() +
  coord_flip() +
  xlab("ATPase Inhibitor")
ATPase_Freq


##--Odds Ratio of MOA and Cmap Names--##

## Annotated full CoxH results with LINCS L1000 MOA Cmap names
CoxH_KMT2A_Anno <- CoxH_KMT2A
var_list_full <- CoxH_KMT2A_Anno$variable
pert_id_full <- c()
for (i in var_list_full) {
  first <- strsplit(i, ".", fixed = T)[[1]][1]
  second <- strsplit(i, ".", fixed = T)[[1]][2]
  id <- paste(first,"-",second,sep = "")
  pert_id_full <- c(pert_id_full,id)
}
CoxH_KMT2A_Anno$pert_id <- pert_id_full
CoxH_KMT2A_Anno <- merge(CoxH_KMT2A_Anno,anno[,c(1,2,4)], by = "pert_id",all.x = T)
CoxH_KMT2A_Anno <- unique(CoxH_KMT2A_Anno)

## Perform Fishers Exact Test to get odds ratio for each MOA and Cmap Name
## MOA
moa_top15 <- as.vector(moa_freq[c(1:15),1])
MOA_Fishers_Res <- list()
moa_pval <- c()
moa_or <- c()
for (i in moa_top15) {
  
  moa <- i
  
  moa_sig <- length(which(CoxH_KMT2A_HR_Sig_Anno$moa == moa))
  not_moa_sig <- length(which(CoxH_KMT2A_HR_Sig_Anno$moa != moa | is.na(CoxH_KMT2A_HR_Sig_Anno$moa)))
  moa_not_sig <- length(which(CoxH_KMT2A_Anno$moa == moa))
  not_moa_not_sig <- length(which(CoxH_KMT2A_Anno$moa != moa | is.na(CoxH_KMT2A_Anno$moa)))
  
  fish_moa_mat <- matrix(data = c(moa_sig,not_moa_sig,moa_not_sig,not_moa_not_sig), nrow = 2)
  fish_moa_res <- fisher.test(fish_moa_mat, alternative = "two.sided")
  
  moa_p <- fish_moa_res[[1]]
  moa_o <- as.numeric(fish_moa_res[[3]])
  
  MOA_Fishers_Res[[moa]] <- fish_moa_res
  moa_pval <- c(moa_pval,moa_p)
  moa_or <- c(moa_or,moa_o)
  
}
moa_fisher_tab <- data.frame(MOA = moa_top15,Pval = moa_pval,Odds_Ratio = moa_or)

## Cmap
cmap_top15 <- as.vector(cmap_freq[c(1:15),1])
CMAP_Fishers_Res <- list()
cmap_pval <- c()
cmap_or <- c()
for (i in cmap_top15) {
  
  cmap <- i
  
  cmap_sig <- length(which(CoxH_KMT2A_HR_Sig_Anno$cmap_name == cmap))
  not_cmap_sig <- length(which(CoxH_KMT2A_HR_Sig_Anno$cmap_name != cmap | is.na(CoxH_KMT2A_HR_Sig_Anno$cmap_name)))
  cmap_not_sig <- length(which(CoxH_KMT2A_Anno$cmap_name == cmap))
  not_cmap_not_sig <- length(which(CoxH_KMT2A_Anno$cmap_name != cmap | is.na(CoxH_KMT2A_Anno$cmap_name)))
  
  fish_cmap_mat <- matrix(data = c(cmap_sig,not_cmap_sig,cmap_not_sig,not_cmap_not_sig), nrow = 2)
  fish_cmap_res <- fisher.test(fish_cmap_mat, alternative = "two.sided")
  
  cmap_p <- fish_cmap_res[[1]]
  cmap_o <- as.numeric(fish_cmap_res[[3]])
  
  CMAP_Fishers_Res[[cmap]] <- fish_cmap_res
  cmap_pval <- c(cmap_pval,cmap_p)
  cmap_or <- c(cmap_or,cmap_o)
  
}
cmap_fisher_tab <- data.frame(CMap_Name = cmap_top15,Pval = cmap_pval,Odds_Ratio = cmap_or)


moa_tab <- merge(moa_fisher_tab,moa_freq, by = "MOA", all.x = T)
cmap_tab <- merge(cmap_fisher_tab,cmap_freq, by = "CMap_Name", all.x = T)

moa_tab2 <- moa_tab[which(moa_tab$Odds_Ratio > 1),]
cmap_tab2 <- cmap_tab[which(cmap_tab$Odds_Ratio > 1),]

moa_tab2 <- moa_tab2 %>%
  mutate(pval_star = case_when(
    Pval < 0.001 ~ "***",
    Pval < 0.01 ~ "**",
    Pval < 0.05 ~ "*"
  ))
cmap_tab2 <- cmap_tab2 %>%
  mutate(pval_star = case_when(
    Pval < 0.001 ~ "***",
    Pval < 0.01 ~ "**",
    Pval < 0.05 ~ "*"
  ))


MOA_OR <- ggplot(data = moa_tab2, aes(x = reorder(MOA, Odds_Ratio), y = Odds_Ratio)) +
  geom_bar(stat = "identity", color='skyblue',fill='steelblue') +
  theme_minimal() +
  coord_flip() +
  xlab("Mechanism of Action (MOA)") +
  ylab("Odds Ratio") +
  ggtitle("Odds Ratios of Top Mechanisms of Action for\nLINCS L1000 High Risk Gene Sets") +
  geom_text(aes(label = pval_star), hjust = -1, color = "black")
MOA_OR
ggsave(MOA_OR,file = "Output/TARGET_AML_KMT2A_EFS_TopMOA_OR_BarPlot.svg", height = 6, width = 10)

Cmap_OR <- ggplot(data = cmap_tab2, aes(x = reorder(CMap_Name, Odds_Ratio), y = Odds_Ratio)) +
  geom_bar(stat = "identity", color='skyblue',fill='steelblue') +
  theme_minimal() +
  coord_flip() +
  xlab("CMap Name") +
  ylab("Odds Ratio") +
  ggtitle("Odds Ratios of Top CMap Names for\nLINCS L1000 High Risk Gene Sets") +
  geom_text(aes(label = pval_star), hjust = -1, color = "black")
Cmap_OR
ggsave(Cmap_OR,file = "Output/TARGET_AML_KMT2A_EFS_TopCMap_OR_BarPlot.svg", height = 6, width = 10)



##--Post-process drug prioritization for 75% bootstrapped KMT2A Subset--##

## File paths
LINCS_Annotation_File <- "Input/LINCS_compoundinfo_beta.txt"
CoxH_Boot_Ranked_Output_File <- "Input/AML_TARGET_KMT2A_EFS_75Bootstrap_LINCS.DN_coxh_ranked.txt"


## Read in files
CoxH_Boot <- as.data.frame(read_delim(CoxH_Boot_Ranked_Output_File, delim = '\t', col_names = T, comment = "#"))
anno <- as.data.frame(read_delim(LINCS_Annotation_File, delim = '\t', col_names = T))


## Subset pipeline results for significant and high-risk results
# The bootstrap was run on the newer version of the DRPPM PATH SURVEIOR Pipeline, where Hazard Ratio > 1 is associated with High Risk
# The p.value column used to filter in the previous version is not labeled "High_AboveMedian_Pval"
# The estimate column used to filter hazard ratio in the previous version is now labeled "Hazard_Ratio"
CoxH_Boot_HR <- CoxH_Boot[which(CoxH_Boot$Hazard_Ratio > 1),]
CoxH_Boot_HR$High_AboveMedian_Pval <- as.numeric(CoxH_Boot_HR$High_AboveMedian_Pval)
CoxH_Boot_HR_Sig <- CoxH_Boot_HR[which(is.na(CoxH_Boot_HR$High_AboveMedian_Pval) | CoxH_Boot_HR$High_AboveMedian_Pval <= 0.05),]
CoxH_Boot_HR_Sig$High_AboveMedian_Pval <- as.character(CoxH_Boot_HR_Sig$High_AboveMedian_Pval)
CoxH_Boot_HR_Sig[is.na(CoxH_Boot_HR_Sig$High_AboveMedian_Pval),4] <- "<0.001"


## Annotate Gene Sets with LINCS L1000 Data
CoxH_Boot_HR_Sig_Anno <- CoxH_Boot_HR_Sig %>%
  select(variable,Hazard_Ratio,ci,High_AboveMedian_Pval)
var_list <- CoxH_Boot_HR_Sig_Anno$variable
pert_id <- c()
for (i in var_list) {
  first <- strsplit(i, ".", fixed = T)[[1]][1]
  second <- strsplit(i, ".", fixed = T)[[1]][2]
  id <- paste(first,"-",second,sep = "")
  pert_id <- c(pert_id,id)
}
CoxH_Boot_HR_Sig_Anno$pert_id <- pert_id
CoxH_Boot_HR_Sig_Anno <- merge(CoxH_Boot_HR_Sig_Anno,anno[,c(1,2,4)], by = "pert_id",all.x = T)
CoxH_Boot_HR_Sig_Anno <- unique(CoxH_Boot_HR_Sig_Anno)

##--Frquency of MOA and Cmap names--##

## Get frequency of unique mechanisms of action (MOA)
moa_freq_boot <- as.data.frame(table(CoxH_Boot_HR_Sig_Anno$moa))
moa_freq_boot <- moa_freq_boot[order(moa_freq_boot$Freq,decreasing = T),]
rownames(moa_freq_boot) <- 1:nrow(moa_freq_boot)
colnames(moa_freq_boot) <- c("MOA","Frequency")

## Get frequency of unique Cmap names
cmap_freq_boot <- as.data.frame(table(CoxH_Boot_HR_Sig_Anno$cmap_name))
cmap_freq_boot <- cmap_freq_boot[order(cmap_freq_boot$Freq,decreasing = T),]
rownames(cmap_freq_boot) <- 1:nrow(cmap_freq_boot)
colnames(cmap_freq_boot) <- c("CMap_Name","Frequency")

## Plot Histograms of MOA and Cmap frequency
MOA_Freq_boot <- ggplot(data = moa_freq_boot[c(1:15),], aes(x = reorder(MOA, Frequency), y = Frequency)) +
  geom_bar(stat = "identity", color='skyblue',fill='steelblue') +
  theme_minimal() +
  coord_flip() +
  xlab("Mechanism of Action (MOA)") +
  ggtitle("MOA Frquency of Appearance in\nLINCS L1000 High Risk Gene Sets\nFrom 75% Bootstrap")
MOA_Freq_boot

Cmap_Freq_boot <- ggplot(data = cmap_freq_boot[c(1:15),], aes(x = reorder(CMap_Name, Frequency), y = Frequency)) +
  geom_bar(stat = "identity", color='skyblue',fill='steelblue') +
  theme_minimal() +
  coord_flip() +
  xlab("CMap Name") +
  ggtitle("CMap Name Frquency of Appearance in\nLINCS L1000 High Risk Gene Sets\nFrom 75% Bootstrap")
Cmap_Freq_boot


##--Odds Ratio of MOA and Cmap Names--##

## Annotated full CoxH results with LINCS L1000 MOA Cmap names
CoxH_Boot_Anno <- CoxH_Boot
var_list_boot <- CoxH_Boot_Anno$variable
pert_id_boot <- c()
for (i in var_list_boot) {
  first <- strsplit(i, ".", fixed = T)[[1]][1]
  second <- strsplit(i, ".", fixed = T)[[1]][2]
  id <- paste(first,"-",second,sep = "")
  pert_id_boot <- c(pert_id_boot,id)
}
CoxH_Boot_Anno$pert_id <- pert_id_boot
CoxH_Boot_Anno <- merge(CoxH_Boot_Anno,anno[,c(1,2,4)], by = "pert_id",all.x = T)
CoxH_Boot_Anno <- unique(CoxH_Boot_Anno)

## Perform Fishers Exact Test to get odds ratio for each MOA and Cmap Name
## MOA
moa_top15_boot <- as.vector(moa_freq_boot[c(1:15),1])
MOA_Fishers_Res_boot <- list()
moa_pval_boot <- c()
moa_or_boot <- c()
for (i in moa_top15_boot) {
  
  moa <- i
  
  moa_sig <- length(which(CoxH_Boot_HR_Sig_Anno$moa == moa))
  not_moa_sig <- length(which(CoxH_Boot_HR_Sig_Anno$moa != moa | is.na(CoxH_Boot_HR_Sig_Anno$moa)))
  moa_not_sig <- length(which(CoxH_Boot_Anno$moa == moa))
  not_moa_not_sig <- length(which(CoxH_Boot_Anno$moa != moa | is.na(CoxH_Boot_Anno$moa)))
  
  fish_moa_mat <- matrix(data = c(moa_sig,not_moa_sig,moa_not_sig,not_moa_not_sig), nrow = 2)
  fish_moa_res <- fisher.test(fish_moa_mat, alternative = "two.sided")
  
  moa_p <- fish_moa_res[[1]]
  moa_o <- as.numeric(fish_moa_res[[3]])
  
  MOA_Fishers_Res_boot[[moa]] <- fish_moa_res
  moa_pval_boot <- c(moa_pval_boot,moa_p)
  moa_or_boot <- c(moa_or_boot,moa_o)
  
}
moa_fisher_tab_boot <- data.frame(MOA = moa_top15_boot,Pval = moa_pval_boot,Odds_Ratio = moa_or_boot)

## Cmap
cmap_top15_boot <- as.vector(cmap_freq_boot[c(1:15),1])
CMAP_Fishers_Res_boot <- list()
cmap_pval_boot <- c()
cmap_or_boot <- c()
for (i in cmap_top15_boot) {
  
  cmap <- i
  
  cmap_sig <- length(which(CoxH_Boot_HR_Sig_Anno$cmap_name == cmap))
  not_cmap_sig <- length(which(CoxH_Boot_HR_Sig_Anno$cmap_name != cmap | is.na(CoxH_Boot_HR_Sig_Anno$cmap_name)))
  cmap_not_sig <- length(which(CoxH_Boot_Anno$cmap_name == cmap))
  not_cmap_not_sig <- length(which(CoxH_Boot_Anno$cmap_name != cmap | is.na(CoxH_Boot_Anno$cmap_name)))
  
  fish_cmap_mat <- matrix(data = c(cmap_sig,not_cmap_sig,cmap_not_sig,not_cmap_not_sig), nrow = 2)
  fish_cmap_res <- fisher.test(fish_cmap_mat, alternative = "two.sided")
  
  cmap_p <- fish_cmap_res[[1]]
  cmap_o <- as.numeric(fish_cmap_res[[3]])
  
  CMAP_Fishers_Res_boot[[cmap]] <- fish_cmap_res
  cmap_pval_boot <- c(cmap_pval_boot,cmap_p)
  cmap_or_boot <- c(cmap_or_boot,cmap_o)
  
}
cmap_fisher_tab_boot <- data.frame(CMap_Name = cmap_top15_boot,Pval = cmap_pval_boot,Odds_Ratio = cmap_or_boot)


moa_tab_boot <- merge(moa_fisher_tab_boot,moa_freq_boot, by = "MOA", all.x = T)
cmap_tab_boot <- merge(cmap_fisher_tab_boot,cmap_freq_boot, by = "CMap_Name", all.x = T)

moa_tab2_boot <- moa_tab_boot[which(moa_tab_boot$Odds_Ratio > 1),]
cmap_tab2_boot <- cmap_tab_boot[which(cmap_tab_boot$Odds_Ratio > 1),]

moa_tab2_boot <- moa_tab2_boot %>%
  mutate(pval_star = case_when(
    Pval < 0.001 ~ "***",
    Pval < 0.01 ~ "**",
    Pval < 0.05 ~ "*"
  ))
cmap_tab2_boot <- cmap_tab2_boot %>%
  mutate(pval_star = case_when(
    Pval < 0.001 ~ "***",
    Pval < 0.01 ~ "**",
    Pval < 0.05 ~ "*"
  ))


MOA_OR_boot <- ggplot(data = moa_tab2_boot, aes(x = reorder(MOA, Odds_Ratio), y = Odds_Ratio)) +
  geom_bar(stat = "identity", color='skyblue',fill='steelblue') +
  theme_minimal() +
  coord_flip() +
  xlab("Mechanism of Action (MOA)") +
  ylab("Odds Ratio") +
  ggtitle("Odds Ratios of Top Mechanisms of Action for\nLINCS L1000 High Risk Gene Sets\nfrom 75% Bootstrap") +
  geom_text(aes(label = pval_star), hjust = -1, color = "black")
MOA_OR_boot
ggsave(MOA_OR_boot,file = "Output/TARGET_AML_KMT2A_EFS_75Bootstrap_TopMOA_OR_BarPlot.svg", height = 6, width = 10)

Cmap_OR_boot <- ggplot(data = cmap_tab2_boot, aes(x = reorder(CMap_Name, Odds_Ratio), y = Odds_Ratio)) +
  geom_bar(stat = "identity", color='skyblue',fill='steelblue') +
  theme_minimal() +
  coord_flip() +
  xlab("CMap Name") +
  ylab("Odds Ratio") +
  ggtitle("Odds Ratios of Top CMap Names for\nLINCS L1000 High Risk Gene Sets\nfrom 75% Bootstrap") +
  geom_text(aes(label = pval_star), hjust = -1, color = "black")
Cmap_OR_boot
ggsave(Cmap_OR_boot,file = "Output/TARGET_AML_KMT2A_EFS_75Bootstrap_TopCMap_OR_BarPlot.svg", height = 6, width = 10)




