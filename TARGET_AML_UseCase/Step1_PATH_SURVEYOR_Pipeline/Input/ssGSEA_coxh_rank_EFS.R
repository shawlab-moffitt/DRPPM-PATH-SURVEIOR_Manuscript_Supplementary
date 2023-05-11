
####----Command Line Arguments----####

args=commandArgs(trailingOnly = T)

if (length(args) < 6) {
  stop("Incorrect Number of Arguments. Argument order: [Project Name] [Expression Matrix File] [Meta Data File] [Gene set File] [Output File Path] [Optional: Annotation File] [Optional: Number of Threads]", call.=FALSE)
}
if (length(args) >= 6) {
  
  # Project Name
  Project_Name <- args[1]
  # Expression File
  Expression_Matrix_file <- args[2]
  # Meta File
  Meta_Data_File <- args[3]
  # Gene Set File
  Gene_Set_File <- args[4]
  # Output File Path
  Output_File_Path <- args[5]
  # Add forward slash if user forgets
  #if (str_sub(Output_File_Path,-1) != '/') {
  #  Output_File_Path <- paste(Output_File_Path,'/',sep = "")
  #}
  # Annotation File
  Annotation_File <- args[6]
  # Number of Threads
  Number_of_Threads <- args[7]
  
}


####----Packages----####

# Base Packages
packages <- c("gtsummary","tidyr","dplyr","DT","ggpubr","tibble",
              "survival","readr","survminer","gridExtra")
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
invisible(lapply(packages, library, character.only = TRUE))
# Bioconductor packages
bioCpacks <- c("GSVA","clusterProfiler")
installed_packages_BIOC <- bioCpacks %in% rownames(installed.packages())
if (any(installed_packages_BIOC == FALSE)) {
  BiocManager::install(bioCpacks[!installed_packages_BIOC], ask = F)
}
invisible(lapply(bioCpacks, library, character.only = TRUE))


####----User Input----####

## Project Name
#Project_Name <- "AML_TARGET_NUP98_fusion"
## Expression File
#Expression_Matrix_file <- "~/R/data/AML_TARGET_data/AML_TARGET_NUP98_expr.txt"
## Meta File
#Meta_Data_File <- "~/R/data/AML_TARGET_data/AML_TARGET_NUP98_meta.txt"
## Gene Set File
#Gene_Set_File <- "~/R/data/GeneSetData/LINCS_L1000_gs_DN_HS.RData"
## Output File Path
#Output_File_Path <- "~/R/data/AML_TARGET_data/"
## Number of threads - OPTIONAL
##Number_of_Threads <- NA
## ssGSEA File - OPTIONAL
##ssgsea_file <- ""
## Annotation data
#Annotation_File <- "~/R/data/GeneSetData/LINCS_compoundinfo_beta.txt"



####----Read in Files----####

print("####----Loading in Expression Matirx----####")

##--Expression--##
expr <- as.data.frame(read_delim(Expression_Matrix_file,delim = '\t', col_names = T))
rownames(expr) <- expr[,1]
expr <- as.matrix(expr[,-1])

##--Meta--##
meta <- as.data.frame(read_delim(Meta_Data_File,delim = '\t',col_names = T))

##--Annotation--##
anno <- as.data.frame(read_delim(Annotation_File,delim = '\t'))

print("####----Loading in Gene Set----####")

##--Gene Set--##
# If user loads GMT File
if (tools::file_ext(Gene_Set_File) == "gmt") {
  GeneSet <- read.gmt(Gene_Set_File)
  colnames(GeneSet) <- c("term","gene")
  GeneSetList <- list()
  for (i in unique(GeneSet[,1])){
    GeneSetList[[i]] <- GeneSet[GeneSet[,1] == i,]$gene
  }
}
# If user loads TSV/TXT file
if (tools::file_ext(Gene_Set_File) == "tsv" || tools::file_ext(Gene_Set_File) == "txt") {
  GeneSet <- read.delim(Gene_Set_File, header = T, sep = '\t')
  colnames(GeneSet) <- c("term","gene")
  GeneSetList <- list()
  for (i in unique(GeneSet[,1])){
    GeneSetList[[i]] <- GeneSet[GeneSet[,1] == i,]$gene
  }
}
# If user loads RData list file
if (tools::file_ext(Gene_Set_File) == "RData") {
  loadRData <- function(fileName){
    #loads an RData file, and returns it
    load(fileName)
    get(ls()[ls() != "fileName"])
  }
  GeneSetList <- loadRData(Gene_Set_File)
}

## High-Low
highlow = function(mat) {
  new_mat = mat;
  new_mat[mat > quantile(mat)[3]] = "High_AboveMedian";
  new_mat[mat <= quantile(mat)[3]] = "Low_BelowMedian";
  return (new_mat)
}


####----Perform Functions----####

if (is.na(Number_of_Threads) == TRUE) {
  
  ##--ssGSEA--##
  print("####----Performing ssGSEA Analysis----####")
  ssgsea <- gsva(expr,GeneSetList,method = "ssgsea")
  ssgsea <- as.data.frame(t(ssgsea))
  ssgsea$SampleName <- rownames(ssgsea)
  ssgsea <- ssgsea %>%
    relocate(SampleName)
  write_delim(as.data.frame(ssgsea),paste(Output_File_Path,gsub(" ","_",Project_Name),"_ssGSEA_score.txt", sep = ""), delim = '\t')
  
  ###--GSVA--##
  #print("####----Performing GSVA Analysis----####")
  #gsva <- gsva(expr,GeneSetList,method = "gsva")
  #gsva <- as.data.frame(gsva)
  #gsva$GeneSet <- rownames(gsva)
  #gsva <- gsva %>%
  #  relocate(GeneSet)
  #write_delim(as.data.frame(gsva),paste(Output_File_Path,gsub(" ","_",Project_Name),"_GSVA_score.txt", sep = ""), delim = '\t')
  #
  ###--zScore--##
  #print("####----Performing zScore Analysis----####")
  #zscore <- gsva(expr,GeneSetList,method = "zscore")
  #zscore <- as.data.frame(zscore)
  #zscore$GeneSet <- rownames(zscore)
  #zscore <- zscore %>%
  #  relocate(GeneSet)
  #write_delim(as.data.frame(zscore),paste(Output_File_Path,gsub(" ","_",Project_Name),"_zScore_score.txt", sep = ""), delim = '\t')
  #
  ###--Plage--##
  #print("####----Performing Plage Analysis----####")
  #plage <- gsva(expr,GeneSetList,method = "plage")
  #plage <- as.data.frame(plage)
  #plage$GeneSet <- rownames(plage)
  #plage <- plage %>%
  #  relocate(GeneSet)
  #write_delim(as.data.frame(plage),paste(Output_File_Path,gsub(" ","_",Project_Name),"_Plage_score.txt", sep = ""), delim = '\t')
  
}

if (is.na(Number_of_Threads) == FALSE) {
  
  ##--ssGSEA--##
  print("####----Performing ssGSEA Analysis----####")
  ssgsea <- gsva(expr,GeneSetList,method = "ssgsea", parallel.sz = Number_of_Threads)
  ssgsea <- as.data.frame(ssgsea)
  ssgsea$GeneSet <- rownames(ssgsea)
  ssgsea <- ssgsea %>%
    relocate(GeneSet)
  write_delim(as.data.frame(ssgsea),paste(Output_File_Path,gsub(" ","_",Project_Name),"_ssGSEA_score.txt", sep = ""), delim = '\t')
  
  ###--GSVA--##
  #print("####----Performing GSVA Analysis----####")
  #gsva <- gsva(expr,GeneSetList,method = "gsva", parallel.sz = Number_of_Threads)
  #gsva <- as.data.frame(gsva)
  #gsva$GeneSet <- rownames(gsva)
  #gsva <- gsva %>%
  #  relocate(GeneSet)
  #write_delim(as.data.frame(gsva),paste(Output_File_Path,gsub(" ","_",Project_Name),"_GSVA_score.txt", sep = ""), delim = '\t')
  #
  ###--zScore--##
  #print("####----Performing zScore Analysis----####")
  #zscore <- gsva(expr,GeneSetList,method = "zscore", parallel.sz = Number_of_Threads)
  #zscore <- as.data.frame(zscore)
  #zscore$GeneSet <- rownames(zscore)
  #zscore <- zscore %>%
  #  relocate(GeneSet)
  #write_delim(as.data.frame(zscore),paste(Output_File_Path,gsub(" ","_",Project_Name),"_zScore_score.txt", sep = ""), delim = '\t')
  #
  ###--Plage--##
  #print("####----Performing Plage Analysis----####")
  #plage <- gsva(expr,GeneSetList,method = "plage", parallel.sz = Number_of_Threads)
  #plage <- as.data.frame(plage)
  #plage$GeneSet <- rownames(plage)
  #plage <- plage %>%
  #  relocate(GeneSet)
  #write_delim(as.data.frame(plage),paste(Output_File_Path,gsub(" ","_",Project_Name),"_Plage_score.txt", sep = ""), delim = '\t')
  
}



####----Binary Analysis----####

print("####----Performing Binary Analysis----####")

rownames(ssgsea) <- ssgsea[,1]
ssgsea <- ssgsea[,-1]

ssgsea <- as.data.frame(t(ssgsea))

## Perform High/Low function
new_col_list <- c()
for (i in colnames(ssgsea)) {
  
  ssgsea[,paste(i,'_BIN',sep = "")] <- highlow(ssgsea[, which(colnames(ssgsea) == i)])
  new_col_list <- c(new_col_list,paste(i,'_BIN',sep = ""))
  
}

## Reformat
ssGSEA_BIN <- ssgsea[,new_col_list]
ssGSEA_BIN$SampleName <- rownames(ssGSEA_BIN)
ssGSEA_BIN <- ssGSEA_BIN %>%
  relocate(SampleName)

write_delim(ssGSEA_BIN,paste(Output_File_Path,gsub(" ","_",Project_Name),"_BIN.txt", sep = ""), delim = '\t')


####----Coxh Analysis----####

print("####----Performing Coxh Analysis----####")

## Clean Gene Set names
colnames(ssGSEA_BIN) <- gsub("[[:punct:]]",".",colnames(ssGSEA_BIN))
colnames(ssGSEA_BIN) <- gsub(" ",".",colnames(ssGSEA_BIN))

## Merge with Meta data
ssGSEA_BIN_meta <- merge(meta,ssGSEA_BIN,by = "SampleName",all = T)

## Binary columns to loop through
calc_cols <- grep("..BIN$", colnames(ssGSEA_BIN_meta), value = T)

## Write out header
header <- c("variable","var_label","var_type","reference_row","row_type","header_row",
            "N_obs","N_event","N","coefficients_type","coefficients_label","label","term",
            "var_class","var_nlevels","contrasts","contrasts_type","n_obs","n_event",
            "exposure","estimate","std.error","statistic","nevent","conf.low",
            "conf.high","ci","p.value")
write(header, file = paste(Output_File_Path,gsub(" ","_",Project_Name),"_coxh.txt", sep = ""),
      append = T, sep = '\t', ncolumns = 28)

## Loop through to perform coxh function
tab_df_list <- list()
for (i in calc_cols) {
  
  meta_ssgsea_sub <- ssGSEA_BIN_meta[,c("EFS.time","Event.ID",i)]
  temp_tab <- coxph(as.formula(paste("Surv(EFS.time,Event.ID) ~ ",i,sep = "")), data = meta_ssgsea_sub) %>% 
    gtsummary::tbl_regression(exp = TRUE) %>%
    as_gt()
  temp_tab_df <- as.data.frame(temp_tab)
  temp_tab_vect <- as.character(temp_tab_df[3,]) #for binary
  tab_df_list[[i]] <- temp_tab_df[3,]
  write(temp_tab_vect, file = paste(Output_File_Path,gsub(" ","_",Project_Name),"_coxh.txt", sep = ""),
        append = T, sep = '\t', ncolumns = 28)
  
}

print("####----Ranking Coxh Analysis----####")


## Combine each gene set to data frame
tab_df <- do.call("rbind", tab_df_list)

## Rank by P.Value and write to file
tab_df$variable <- gsub(".BIN$","",tab_df$variable)
tab_df$p.value <- gsub(">0.9","0.9",tab_df$p.value)
tab_df$p.value <- as.numeric(tab_df$p.value)
tab_df_ordered <- tab_df[order(tab_df$p.value, decreasing = F, na.last = F),]
tab_df_ordered[which(is.na(tab_df_ordered$p.value)),"p.value"] <- "<0.001"
write_delim(tab_df_ordered, paste(Output_File_Path,gsub(" ","_",Project_Name),"_coxh_ranked.txt", sep = ""), delim = '\t')


####----Annotate Data----####

## Get pert_id

var_list <- tab_df_ordered$variable
pert_id <- c()
for (i in var_list) {
  first <- strsplit(i, ".", fixed = T)[[1]][1]
  second <- strsplit(i, ".", fixed = T)[[1]][2]
  id <- paste(first,"-",second,sep = "")
  pert_id <- c(pert_id,id)
}
tab_df_ordered$pert_id <- pert_id


##--Merge Data--##

tab_df_ordered_anno <- merge(tab_df_ordered,anno, by = "pert_id", all.x = T)

##--Reformatting--##

## Relocate columns of interest to the front
tab_df_ordered_anno <- tab_df_ordered_anno %>%
  relocate(pert_id,variable,target,moa,cmap_name,estimate,ci,p.value)
## Reorder by ranking again
tab_df_ordered_anno$p.value <- as.numeric(tab_df_ordered_anno$p.value)
tab_df_ordered_anno <- tab_df_ordered_anno[order(tab_df_ordered_anno$p.value, decreasing = F, na.last = F),]
tab_df_ordered_anno[which(is.na(tab_df_ordered_anno$p.value)),"p.value"] <- "<0.001"

## Write Ranked and annotated data to file
write_delim(tab_df_ordered_anno, paste(Output_File_Path,gsub(" ","_",Project_Name),"_coxh_ranked_anno.txt", sep = ""), delim = '\t')



