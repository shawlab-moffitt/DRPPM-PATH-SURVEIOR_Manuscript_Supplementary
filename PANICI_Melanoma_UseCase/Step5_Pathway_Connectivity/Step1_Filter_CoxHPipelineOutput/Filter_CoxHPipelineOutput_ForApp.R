


####----Filter DRPPM PATH SURVEIOR CoxH Ranked Immune Signature OR Output for Pathway Connectivity App----###



## Set "Step1_Filter_CoxHPipelineOutput_ForApp" as working directory


## Load in libraries
library(readr)
library(dplyr)


## File Paths
ImmSig_CoxH_OS_Ranked_File <- "Input/ICI_iAtlas_Skin_Pre_ImmuneSig_OS_coxh_ranked.txt"
ImmSig_CoxH_PFI_Ranked_File <- "Input/ICI_iAtlas_Skin_Pre_ImmuneSig_PFI_coxh_ranked.txt"


## Read in files
ImmSig_OS_CoxH <- as.data.frame(read_delim(ImmSig_CoxH_OS_Ranked_File, delim = '\t', col_names = T))
ImmSig_PFI_CoxH <- as.data.frame(read_delim(ImmSig_CoxH_PFI_Ranked_File, delim = '\t', col_names = T))


## Filter for low risk features
ImmSig_OS_CoxH_LR <- ImmSig_OS_CoxH[which(ImmSig_OS_CoxH$Hazard_Ratio < 1),]
ImmSig_PFI_CoxH_LR <- ImmSig_PFI_CoxH[which(ImmSig_PFI_CoxH$Hazard_Ratio < 1),]


## Filter for adj P val (BH) < 0.05
ImmSig_OS_CoxH_LR_Sig <- ImmSig_OS_CoxH_LR[which(ImmSig_OS_CoxH_LR$Likelihood_Ratio_AdjPval_BH < 0.05),]
ImmSig_PFI_CoxH_LR_Sig <- ImmSig_PFI_CoxH_LR[which(ImmSig_PFI_CoxH_LR$Likelihood_Ratio_AdjPval_BH < 0.05),]


## Write out to file
write_delim(ImmSig_OS_CoxH_LR_Sig,"Output/ICI_iAtlas_Skin_Pre_ImmuneSig_OS_coxh_ranked_Filtered.txt", delim = '\t')
write_delim(ImmSig_PFI_CoxH_LR_Sig,"Output/ICI_iAtlas_Skin_Pre_ImmuneSig_PFI_coxh_ranked_Filtered.txt", delim = '\t')
