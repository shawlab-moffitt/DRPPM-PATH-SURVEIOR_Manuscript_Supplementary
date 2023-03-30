

####----Command Line Arguments----####

args=commandArgs(trailingOnly = T)

Param_File <- args[1]

#Param_File <- "~/R/data/PenileCancer/PSCC_HPV_CoxH_Ranking/PSCC_HPV_Positive_OS_Genes_CoxH_Params_local.txt"


####----Packages----####

# Base Packages
packages <- c("gtsummary","tidyr","dplyr","DT","ggpubr","tibble","stringr","tools",
              "survival","readr","survminer","gridExtra","BiocManager","matrixStats")
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


####----Read in Parameter File----####

params <- as.data.frame(read_delim(Param_File,delim = "\t", col_names = F))

## Project Name
Project_Name <- params[which(params[,1] == "Project_Name"),2]
Project_Name <- gsub(" ","_",Project_Name)
Project_Name <- gsub("[[:punct:]]","_",Project_Name)

## Expression File
Expression_Matrix_file <- params[which(params[,1] == "Expression_Matrix_file"),2]

## Meta File
Meta_Data_File <- params[which(params[,1] == "Meta_Data_File"),2]

## Gene Set File
Gene_Set_File <- params[which(params[,1] == "Gene_Set_File"),2]

## Gene Set Name
Gene_Set_Name <- params[which(params[,1] == "Gene_Set_Name"),2]
Gene_Set_Name <- gsub(" ",".",Gene_Set_Name)
Gene_Set_Name <- gsub("[[:punct:]]",".",Gene_Set_Name)
if (!exists("Gene_Set_Name")) {
  Gene_Set_Name <- NA
}

## Output File Path
Output_File_Path <- params[which(params[,1] == "Output_File_Path"),2]
if (!exists("Output_File_Path")) {
  Output_File_Path <- NA
}
if (is.na(Output_File_Path)) {
  Output_File_Path <- getwd()
}
# Add a forward slash if missing from the end of path
last_char <- str_sub(Output_File_Path,-1,-1)
if (last_char != "/") {
  Output_File_Path <- paste(Output_File_Path,"/",sep = "")
}

## Survival Time ID
Survival_Time <- params[which(params[,1] == "Survival_Time_Label"),2]

## Survival ID ID
Survival_ID <- params[which(params[,1] == "Survival_ID_Label"),2]

## Rank By Genes Option
Rank_Genes_Choice <- params[which(params[,1] == "Rank_Genes"),2]

## Covariate Column Name
Covariate_Column_Label <- params[which(params[,1] == "Covariate_Column_Label"),2]

## Covariate Reference
Covariate_Reference <- params[which(params[,1] == "Covariate_Reference"),2]



####----Read in Files----####

print("####----Loading in Files----####")

##--Expression--##
expr <- as.data.frame(read_delim(Expression_Matrix_file,delim = '\t', col_names = T))
colnames(expr)[1] <- "Symbol"
if (TRUE %in% duplicated(expr[,1])) {
  expr <- expr %>%
    group_by(Symbol) %>%
    summarise_all(max)
}
isChar <- unname(which(sapply(expr, function(x) is.character(x))))
isChar <-  isChar[-1]
if (length(isChar) > 0) {
  expr[isChar] <- sapply(expr[isChar],as.numeric)
}
expr <- as.data.frame(expr)
rownames(expr) <- expr[,1]
expr <- expr[,-1]
colnames(expr) <- gsub("[[:punct:]]",".",colnames(expr))
expr <- mutate_all(expr, function(x) as.numeric(as.character(x)))
expr <- as.matrix(expr)

##--Meta--##
meta <- as.data.frame(read_delim(Meta_Data_File,delim = '\t',col_names = T))
colnames(meta)[1] <- "SampleName"
meta$SampleName <- gsub("[[:punct:]]",".",meta$SampleName)


## Check that all samples have survival data and covariate data if needed
meta <- meta[!is.na(meta[,Survival_Time]),]
meta <- meta[!is.na(meta[,Survival_ID]),]

## Check if user provided covariate
if (!exists("Covariate_Column_Label")) {
  Covariate_Column_Label <- NA
  Covariate_Reference <- NA
}
if (!is.na(Covariate_Column_Label)) {
  meta <- meta[!is.na(meta[,Covariate_Column_Label]),]
}


##--Gene Set--##

## List to add gene sets - user may not provide any gene set
GeneSet_list <- list()
## If a gene set file is provided
if (!is.na(Gene_Set_File)) {
  ## If the user provides a .lst file of multiple gene sets
  if (tools::file_ext(Gene_Set_File) == "lst") {
    ## Read in as data frame, first column gene set name and second column path
    GeneSetFiles <- as.data.frame(read_delim(Gene_Set_File, delim = '\t',col_names = F))
    colnames(GeneSetFiles) <- c("GeneSetName","GeneSetFilePath")
    ## Read in each gene set file and add to gene set list
    for (i in 1:nrow(GeneSetFiles)) {
      gs_name <- GeneSetFiles[i,1]
      gs_path <- GeneSetFiles[i,2]
      # If user loads GMT File
      if (tools::file_ext(gs_path) == "gmt") {
        GeneSet <- read.gmt(gs_path)
        colnames(GeneSet) <- c("term","gene")
        GeneSetList <- list()
        for (j in unique(GeneSet[,1])){
          GeneSetList[[j]] <- GeneSet[GeneSet[,1] == j,]$gene
        }
        GeneSet_list[[gs_name]] <- GeneSetList
      }
      # If user loads TSV/TXT file
      if (tools::file_ext(gs_path) == "tsv" || tools::file_ext(gs_path) == "txt") {
        GeneSet <- read.delim(gs_path, header = T, sep = '\t')
        colnames(GeneSet) <- c("term","gene")
        GeneSetList <- list()
        for (j in unique(GeneSet[,1])){
          GeneSetList[[j]] <- GeneSet[GeneSet[,1] == j,]$gene
        }
        GeneSet_list[[gs_name]] <- GeneSetList
      }
      # If user loads RData list file
      if (tools::file_ext(gs_path) == "RData") {
        loadRData <- function(fileName){
          #loads an RData file, and returns it
          load(fileName)
          get(ls()[ls() != "fileName"])
        }
        GeneSetList <- loadRData(gs_path)
        GeneSet_list[[gs_name]] <- GeneSetList
      }
    }
  }
  ## User upload only 1 gene set
  if (tools::file_ext(Gene_Set_File) != "lst") {
    if (is.na(Gene_Set_Name)) {
      Gene_Set_Name <- "UserGeneSet"
    }
    # If user loads GMT File
    if (tools::file_ext(Gene_Set_File) == "gmt") {
      GeneSet <- read.gmt(Gene_Set_File)
      colnames(GeneSet) <- c("term","gene")
      GeneSetList <- list()
      for (i in unique(GeneSet[,1])){
        GeneSetList[[i]] <- GeneSet[GeneSet[,1] == i,]$gene
      }
      GeneSet_list[[Gene_Set_Name]] <- GeneSetList
    }
    # If user loads TSV/TXT file
    if (tools::file_ext(Gene_Set_File) == "tsv" || tools::file_ext(Gene_Set_File) == "txt") {
      GeneSet <- read.delim(Gene_Set_File, header = T, sep = '\t')
      colnames(GeneSet) <- c("term","gene")
      GeneSetList <- list()
      for (i in unique(GeneSet[,1])){
        GeneSetList[[i]] <- GeneSet[GeneSet[,1] == i,]$gene
      }
      GeneSet_list[[Gene_Set_Name]] <- GeneSetList
    }
    # If user loads RData list file
    if (tools::file_ext(Gene_Set_File) == "RData") {
      loadRData <- function(fileName){
        #loads an RData file, and returns it
        load(fileName)
        get(ls()[ls() != "fileName"])
      }
      GeneSetList <- loadRData(Gene_Set_File)
      GeneSet_list[[Gene_Set_Name]] <- GeneSetList
    }
  }
}


## High-Low function
highlow = function(mat) {
  new_mat = mat;
  new_mat[mat > quantile(as.numeric(mat),na.rm = T)[3]] = "High_AboveMedian";
  new_mat[mat <= quantile(as.numeric(mat),na.rm = T)[3]] = "Low_BelowMedian";
  return (new_mat)
}



####----Perform Functions----####

##--ssGSEA--##
print("####----Performing ssGSEA Analysis----####")
## List to add ssGSEA df results to
ssGSEA_tabs <- list()
ssGSEA_summ_tabs <- list()
## If gene set provided
if (length(GeneSet_list) > 0) {
  ## For each gene set provided
  for (i in 1:length(GeneSet_list)){
    gs_name <- names(GeneSet_list)[i]           #Gene set name
    gs <- GeneSet_list[[i]]                     #Gene Set contents
    ssgsea <- gsva(expr,gs,method = "ssgsea")   #Outputs sample names as column names
    ssgsea <- as.data.frame(t(ssgsea))          #Transpose to make column names gs names
    ssgsea$SampleName <- rownames(ssgsea)       #Change gs names to row names
    ssgsea <- ssgsea %>%                        #Make rownames first column to output to file
      relocate(SampleName)
    ssGSEA_tabs[[paste("ssgsea",gs_name,sep = "_")]] <- ssgsea
    ## Write to file
    # If onle one gene set
    if (gs_name == "UserGeneSet") {
      write_delim(as.data.frame(ssgsea),paste(Output_File_Path,Project_Name,"_ssGSEA_score.txt", sep = ""), delim = '\t')
    }
    # If more than one gene set and gene set names are provided
    if (gs_name != "UserGeneSet") {
      write_delim(as.data.frame(ssgsea),paste(Output_File_Path,Project_Name,"_",gs_name,"_ssGSEA_score.txt", sep = ""), delim = '\t')
    }
    
    ssgsea2 <- ssgsea
    
    ssgsea2 <- as.data.frame(t(ssgsea2))
    ssgsea2 <- ssgsea2[-1,]
    ssgsea2[] <- lapply(ssgsea2, as.numeric)
    col_to_calc <- ncol(ssgsea2)
    ssgsea2$Median_Score <- rowMedians(as.matrix(ssgsea2[,c(1:col_to_calc)]))
    ssgsea2$Mean_Score <- rowMeans(as.matrix(ssgsea2[,c(1:col_to_calc)]))
    ssgsea2$Min_Score <- apply(ssgsea2[,c(1:col_to_calc)], 1, FUN = min)
    ssgsea2$Max_Score <- apply(ssgsea2[,c(1:col_to_calc)], 1, FUN = max)
    
    ssgsea3 <- ssgsea2 %>%                        #Make rownames first column to output to file
      select(Median_Score,Mean_Score,Min_Score,Max_Score)
    rownames(ssgsea3) <- gsub("[[:punct:]]",".",rownames(ssgsea3))
    ssgsea3$Feature <- rownames(ssgsea3)       #Change gs names to row names
    ssgsea3 <- ssgsea3 %>%                        #Make rownames first column to output to file
      relocate(Feature)
    
    ssGSEA_summ_tabs[[paste("ssgseaSumm",gs_name,sep = "_")]] <- ssgsea3
    
    
  }
}




####----Binary Analysis----####

print("####----Performing Binary Analysis----####")

## List to add ssGSEA binary df results to
ssGSEA_BIN_tabs <- list()
ssGSEA_GeneSumm_tabs <- list()
## if the user chooses to rank genes by raw gene expression
if (Rank_Genes_Choice == TRUE) {
  
  expr2 <- as.data.frame(expr)
  expr2$Median_Score <- rowMedians(as.matrix(expr2[,c(1:ncol(expr))]), na.rm = T)
  expr2$Mean_Score <- rowMeans(as.matrix(expr2[,c(1:ncol(expr))]), na.rm = T)
  expr2$Min_Score <- apply(expr2[,c(1:ncol(expr))], 1, FUN = min)
  expr2$Max_Score <- apply(expr2[,c(1:ncol(expr))], 1, FUN = max)
  
  expr3 <- expr2 %>%
    select(Median_Score,Mean_Score,Min_Score,Max_Score)
  rownames(expr3) <- gsub("[[:punct:]]",".",rownames(expr3))
  expr3$Feature <- rownames(expr3)
  expr3 <- expr3 %>%
    relocate(Feature)
  
  ssGSEA_GeneSumm_tabs[[paste("ssgseaSumm","Genes",sep = "_")]] <- expr3 
  
  ## Transpose expression data
  ssgsea <- as.data.frame(t(expr))
  ## loop through each gene set to get BIN
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
  ## Add ssGSEA_BIN df result to list
  ssGSEA_BIN_tabs[[paste("ssGSEA.BIN","Genes",sep = "_")]] <- ssGSEA_BIN
  ## Write to file
  write_delim(ssGSEA_BIN,paste(Output_File_Path,Project_Name,"_Genes_BIN.txt", sep = ""), delim = '\t')
}
## If a gene set is provided
if (length(ssGSEA_tabs) > 0) {
  for (i in 1:length(ssGSEA_tabs)) {
    ## Get variables
    ssgsea <- ssGSEA_tabs[[i]]
    table_name <- names(ssGSEA_tabs)[i]
    gs_name <- unlist(str_split(table_name,"_",))[2]
    # Place sample names as row names
    rownames(ssgsea) <- ssgsea[,1]
    ssgsea <- ssgsea[,-1]
    ## loop through each gene set to get BIN
    ## Perform High/Low function
    new_col_list <- c()
    for (j in colnames(ssgsea)) {
      ssgsea[,paste(j,'_BIN',sep = "")] <- highlow(ssgsea[, which(colnames(ssgsea) == j)])
      new_col_list <- c(new_col_list,paste(j,'_BIN',sep = ""))
    }
    ## Reformat
    ssGSEA_BIN <- ssgsea[,new_col_list]
    ssGSEA_BIN$SampleName <- rownames(ssGSEA_BIN)
    ssGSEA_BIN <- ssGSEA_BIN %>%
      relocate(SampleName)
    ssGSEA_BIN_tabs[[paste("ssGSEA.BIN",gs_name,sep = "_")]] <- ssGSEA_BIN
    if (gs_name == "UserGeneSet") {
      write_delim(ssGSEA_BIN,paste(Output_File_Path,Project_Name,"_BIN.txt", sep = ""), delim = '\t')
    }
    if (gs_name != "UserGeneSet") {
      write_delim(ssGSEA_BIN,paste(Output_File_Path,Project_Name,"_",gs_name,"_BIN.txt", sep = ""), delim = '\t')
    }
  }
}



####----Coxh Analysis----####

print("####----Performing Coxh Analysis----####")

ssGSEA_BIN_meta_tabs <- list()
for (i in 1:length(ssGSEA_BIN_tabs)) {
  ## Prep and merge tables
  ssGSEA_BIN <- ssGSEA_BIN_tabs[[i]]
  table_name <- names(ssGSEA_BIN_tabs)[i]
  gs_name <- unlist(str_split(table_name,"_",))[2]
  ## Clean Gene Set names
  colnames(ssGSEA_BIN) <- gsub("[[:punct:]]",".",colnames(ssGSEA_BIN))
  colnames(ssGSEA_BIN) <- gsub(" ",".",colnames(ssGSEA_BIN))
  ## Merge with Meta data
  ssGSEA_BIN_meta <- merge(meta,ssGSEA_BIN,by = "SampleName",all.y = T)
  ssGSEA_BIN_meta_tabs[[paste("ssGSEA.BIN.meta",gs_name,sep = "_")]] <- ssGSEA_BIN_meta
  ## Binary columns to loop through
  calc_cols <- grep("..BIN$", colnames(ssGSEA_BIN_meta), value = T)
  
  ## Put parameters and header together
  ProjName_Line <- paste("## Project Name:",Project_Name)
  Survival_Line <- paste("## Survival/Event Variables:",Survival_Time,"///",Survival_ID)
  Samples_Line <- paste("## Number of Samples:",nrow(meta))
  expr_file <- paste("## Expression Matrix File Name:",basename(Expression_Matrix_file))
  meta_file <- paste("## Meta Data File Name:",basename(Meta_Data_File))
  header <- c("variable","var_label","var_type","reference_row","row_type","header_row",
              "N_obs","N_event","N","coefficients_type","coefficients_label","Criteria","term",
              "var_class","var_nlevels","contrasts","contrasts_type","n_obs","n_event",
              "exposure","Hazard_Ratio","std.error","statistic","nevent","conf.low",
              "conf.high","ci","High_AboveMedian_Pval","Concordance","Likelihood_Ratio_Pval","Wald_Test_Pval",
              "Logrank_Test_Pval","Median_Score","Mean_Score","Min_Score","Max_Score")
  ## If no covariate being analyzed
  if (is.na(Covariate_Column_Label)) {
    ## If performing coxh on gene express
    if (gs_name == "Genes") {
      line_analysis <- paste("## Samples were grouped to Above or Below median based on raw gene expression then ran through Cox Proportional Hazard Analysis")
      line_analysis2 <- paste("## Genes were ranked by Coxph Likelihood_Ratio_Pval")
      line_analysis3 <- paste("## Hazard Ratios > 1 represent high (above median) raw gene expression associated with High-Risk")
      param_lines <- paste(ProjName_Line,expr_file,meta_file,Survival_Line,Samples_Line,line_analysis,line_analysis2,line_analysis3,sep = "\n")
      file_made <- paste(Output_File_Path,Project_Name,"_",gs_name,"_coxh.txt", sep = "")
      write(param_lines,file = file_made, append = T, sep = '\t', ncolumns = 36)
    }
    ## If performing coxh on pathways
    if (gs_name != "Genes") {
      geneset_file <- paste("## Gene Set File Name:",basename(Gene_Set_File))
      if (tools::file_ext(geneset_file) == "lst") {
        geneset_file <- basename(GeneSetFiles[which(GeneSetFiles$GeneSetName == gs_name),2])
      }
      line_geneset <- paste("## Gene Set:",gs_name,"Of",length(GeneSet_list[[gs_name]]),"Gene Sets")
      line_analysis <- paste("## ssGSEA was used to calculate pathway scores and group samples Above or Below median ssGSEA score per gene set")
      line_analysis2 <- paste("## Pathways were ranked by Coxph Likelihood_Ratio_Pval")
      line_analysis3 <- paste("## Hazard Ratios > 1 represent pathway scores associated with High-Risk")
      param_lines <- paste(ProjName_Line,expr_file,meta_file,geneset_file,Survival_Line,Samples_Line,line_geneset,line_analysis,line_analysis2,line_analysis3,sep = "\n")
      if (gs_name == "UserGeneSet") {
        file_made <- paste(Output_File_Path,Project_Name,"Pathway_coxh.txt", sep = "")
        write(param_lines,file = file_made, append = T, sep = '\t', ncolumns = 36)
      }
      if (gs_name != "UserGeneSet") {
        file_made <- paste(Output_File_Path,Project_Name,"_",gs_name,"_coxh.txt", sep = "")
        write(param_lines,file = file_made, append = T, sep = '\t', ncolumns = 36)
      }
    }
  }
  ## If there is a covariate to be analyzed
  if (!is.na(Covariate_Column_Label)) {
    line_covar <- paste("## Covariate and Reference Variable:",Covariate_Column_Label,"Referencing",Covariate_Reference)
    covar_brk <- as.data.frame(table(ssGSEA_BIN_meta[,Covariate_Column_Label]))
    covars <- c()
    for (k in 1:nrow(covar_brk)) {
      var <- as.character(covar_brk[k,1])
      num <- as.character(covar_brk[k,2])
      varnum <- c("### Variable:",var,"(","N =",num,")")
      covars <- c(covars,paste(varnum,collapse = " "))
    }
    line_covar2 <- paste("## Covariate Breakdown:\n",paste(covars,collapse = "\n"),sep = "")
    ## If gene expression and covariate being analyzed
    if (gs_name == "Genes") {
      line_analysis_int <- paste("## Samples were grouped to Above or Below median based on raw gene expression then ran through interactive Cox Proportional Hazard Analysis with the covariate")
      line_analysis_add <- paste("## Samples were grouped to Above or Below median based on raw gene expression then ran through additive Cox Proportional Hazard Analysis with the covariate")
      line_analysis2 <- paste("## Genes were ranked by Coxph Likelihood_Ratio_Pval")
      line_analysis3 <- paste("## Hazard Ratios > 1 represent high (above median) raw gene expression associated with High-Risk")
      param_lines_int <- paste(ProjName_Line,expr_file,meta_file,Survival_Line,Samples_Line,line_analysis_int,line_analysis2,line_analysis3,line_covar,line_covar2,sep = "\n")
      param_lines_add <- paste(ProjName_Line,expr_file,meta_file,Survival_Line,Samples_Line,line_analysis_add,line_analysis2,line_analysis3,line_covar,line_covar2,sep = "\n")
      file_made_int <- paste(Output_File_Path,Project_Name,"_",gs_name,"_",Covariate_Column_Label,"_Interactive_coxh.txt", sep = "")
      file_made_add <- paste(Output_File_Path,Project_Name,"_",gs_name,"_",Covariate_Column_Label,"_Additive_coxh.txt", sep = "")
      write(param_lines_int,file = file_made_int, append = T, sep = '\t', ncolumns = 36)
      write(param_lines_add,file = file_made_add, append = T, sep = '\t', ncolumns = 36)
    }
    ## If pathway and gene expression being analyzed
    if (gs_name != "Genes") {
      geneset_file <- paste("## Gene Set File Name:",basename(Gene_Set_File))
      if (tools::file_ext(geneset_file) == "lst") {
        geneset_file <- basename(GeneSetFiles[which(GeneSetFiles$GeneSetName == gs_name),2])
      }
      line_geneset <- paste("## Gene Set:",gs_name,"Of",length(GeneSet_list[[gs_name]]),"Gene Sets")
      line_analysis_int <- paste("## ssGSEA was used to calculate pathway scores and group samples Above or Below median ssGSEA score per gene set then ran through interactive Cox Proportional Hazard Analysis with the covariate")
      line_analysis_add <- paste("## ssGSEA was used to calculate pathway scores and group samples Above or Below median ssGSEA score per gene set then ran through additive Cox Proportional Hazard Analysis with the covariate")
      line_analysis2 <- paste("## Pathways were ranked by Coxph Likelihood_Ratio_Pval")
      line_analysis3 <- paste("## Hazard Ratios > 1 represent pathway scores associated with High-Risk")
      param_lines_int <- paste(ProjName_Line,expr_file,meta_file,geneset_file,Survival_Line,Samples_Line,line_geneset,line_analysis_int,line_analysis2,line_analysis3,line_covar,line_covar2,sep = "\n")
      param_lines_add <- paste(ProjName_Line,expr_file,meta_file,geneset_file,Survival_Line,Samples_Line,line_geneset,line_analysis_add,line_analysis2,line_analysis3,line_covar,line_covar2,sep = "\n")
      if (gs_name == "UserGeneSet") {
        file_made_int <- paste(Output_File_Path,Project_Name,"_",Covariate_Column_Label,"_Interactive_coxh.txt", sep = "")
        file_made_add <- paste(Output_File_Path,Project_Name,"_",Covariate_Column_Label,"_Additive_coxh.txt", sep = "")
        write(param_lines_int,file = file_made_int, append = T, sep = '\t', ncolumns = 36)
        write(param_lines_add,file = file_made_add, append = T, sep = '\t', ncolumns = 36)
      }
      if (gs_name != "UserGeneSet") {
        file_made_int <- paste(Output_File_Path,Project_Name,"_",gs_name,"_",Covariate_Column_Label,"_Interactive_coxh.txt", sep = "")
        file_made_add <- paste(Output_File_Path,Project_Name,"_",gs_name,"_",Covariate_Column_Label,"_Additive_coxh.txt", sep = "")
        write(param_lines_int,file = file_made_int, append = T, sep = '\t', ncolumns = 36)
        write(param_lines_add,file = file_made_add, append = T, sep = '\t', ncolumns = 36)
      }
    }
  }
  
  ## For non-interactive pathway analysis
  if (is.na(Covariate_Column_Label)) {
    out_df <- as.data.frame(read_delim(file_made,delim = '\t', col_names = F))
    out_df[c(2:36)] <- NA
    out_df <- rbind(out_df,header)
    write(header,file = file_made, append = T, sep = '\t', ncolumns = 36)
    for (j in calc_cols) {
      meta_ssgsea_sub <- ssGSEA_BIN_meta[,c(Survival_Time,Survival_ID,j)]
      meta_ssgsea_sub[,j] <- as.factor(meta_ssgsea_sub[,j])
      meta_ssgsea_sub[,j] <- relevel(meta_ssgsea_sub[,j], ref = "Low_BelowMedian")
      colnames(meta_ssgsea_sub)[which(colnames(meta_ssgsea_sub) == Survival_Time)] <- "time"
      colnames(meta_ssgsea_sub)[which(colnames(meta_ssgsea_sub) == Survival_ID)] <- "ID"
      if (length(levels(as.factor(meta_ssgsea_sub[,j]))) > 1) {
        if (!is.na(as.numeric(str_sub(j,1,1)))) {
          colnames(meta_ssgsea_sub)[3] <- paste("n",j,sep = "")
          temp_tab <- coxph(as.formula(paste("Surv(time,ID) ~ ",paste("n",j,sep = ""),sep = "")),
                            data = meta_ssgsea_sub) %>% 
            gtsummary::tbl_regression(exp = TRUE) %>%
            as_gt()
          temp_tab_df <- as.data.frame(temp_tab)
          tab <- coxph(as.formula(paste("Surv(time,ID) ~ ",paste("n",j,sep = ""),sep = "")),
                       data = meta_ssgsea_sub)
          out <- capture.output(summary(tab))
          con_line <- grep("^Concordance=",out,value = T)
          con_v <- str_split(con_line," ")[[1]][2]
          lik_line <- grep("^Likelihood ratio test=",out,value = T)
          lik_p <- str_split(lik_line,"=")[[1]][3]
          wal_line <- grep("^Wald test",out,value = T)
          wal_p <- str_split(wal_line,"=")[[1]][3]
          sco_line <- grep("^Score ",out,value = T)
          sco_p <- str_split(sco_line,"=")[[1]][3]
          temp_tab_df[3,c(1,2,13)] <- sub(".","",temp_tab_df[3,c(1,2,13)])
          temp_tab_vect <- as.character(c(temp_tab_df[3,]))
          temp_tab_vect <- c(temp_tab_vect,con_v,lik_p,wal_p,sco_p)
          
          variable <- file_path_sans_ext(temp_tab_vect[1])
          if (gs_name == "Genes") {
            score_df <- ssGSEA_GeneSumm_tabs[[1]]
          }
          if (gs_name != "Genes") {
            score_df <- ssGSEA_summ_tabs[[i]]
          }
          var_scores <- unname(unlist(c(score_df[variable,c(2:5)])))
          
          temp_tab_vect <- c(temp_tab_vect,var_scores)
          
          out_df <- rbind(out_df,temp_tab_vect)
          write(temp_tab_vect,file = file_made, append = T, sep = '\t', ncolumns = 36)
        }
        if (is.na(as.numeric(str_sub(j,1,1)))) {
          tab <- coxph(as.formula(paste("Surv(time,ID) ~ ",j,sep = "")),
                       data = meta_ssgsea_sub)
          temp_tab <- tab %>% 
            gtsummary::tbl_regression(exp = TRUE) %>%
            as_gt()
          temp_tab_df <- as.data.frame(temp_tab)
          out <- capture.output(summary(tab))
          con_line <- grep("^Concordance=",out,value = T)
          con_v <- str_split(con_line," ")[[1]][2]
          lik_line <- grep("^Likelihood ratio test=",out,value = T)
          lik_p <- str_split(lik_line,"=")[[1]][3]
          wal_line <- grep("^Wald test",out,value = T)
          wal_p <- str_split(wal_line,"=")[[1]][3]
          sco_line <- grep("^Score ",out,value = T)
          sco_p <- str_split(sco_line,"=")[[1]][3]
          temp_tab_vect <- as.character(c(temp_tab_df[3,]))
          temp_tab_vect <- c(temp_tab_vect,con_v,lik_p,wal_p,sco_p)
          
          variable <- file_path_sans_ext(temp_tab_vect[1])
          if (gs_name == "Genes") {
            score_df <- ssGSEA_GeneSumm_tabs[[1]]
          }
          if (gs_name != "Genes") {
            score_df <- ssGSEA_summ_tabs[[i]]
          }
          var_scores <- unname(unlist(c(score_df[variable,c(2:5)])))
          
          temp_tab_vect <- c(temp_tab_vect,var_scores)
          
          out_df <- rbind(out_df,temp_tab_vect)
          write(temp_tab_vect,file = file_made, append = T, sep = '\t', ncolumns = 36)
        }
      }
    }
    out_df_top <- out_df[grep("##",out_df[,1]),]
    out_df_top[,c(37:39)] <- NA
    tab_df3 <- out_df[grep("##",out_df[,1],invert = T),]
    colnames(tab_df3) <- header
    tab_df3 <- tab_df3[-1,]
    Likelihood_Ratio_adjPval_BH <- p.adjust(as.numeric(tab_df3$Likelihood_Ratio_Pval), method = "BH")
    Wald_Test_adjPval_BH <- p.adjust(as.numeric(tab_df3$Wald_Test_Pval), method = "BH")
    Logrank_Test_adjPval_BH <- p.adjust(as.numeric(tab_df3$Logrank_Test_Pval), method = "BH")
    tab_df3 <- cbind(tab_df3,Likelihood_Ratio_adjPval_BH,Wald_Test_adjPval_BH,Logrank_Test_adjPval_BH)
    tab_df3$variable <- gsub(".BIN$","",tab_df3$variable)
    tab_df3$High_AboveMedian_Pval <- gsub(">0.9","0.9",tab_df3$High_AboveMedian_Pval)
    tab_df3$High_AboveMedian_Pval <- as.numeric(tab_df3$High_AboveMedian_Pval)
    tab_df3$Likelihood_Ratio_Pval <- as.numeric(tab_df3$Likelihood_Ratio_Pval)
    tab_df3_ordered <- tab_df3[order(tab_df3$Likelihood_Ratio_Pval, decreasing = F),]
    tab_df3_ordered$Likelihood_Ratio_Pval <- formatC(tab_df3_ordered$Likelihood_Ratio_Pval)
    tab_df3_ordered$Wald_Test_Pval <- formatC(as.numeric(tab_df3_ordered$Wald_Test_Pval))
    tab_df3_ordered$Logrank_Test_Pval <- formatC(as.numeric(tab_df3_ordered$Logrank_Test_Pval))
    #tab_df3_ordered <- tab_df3[order(tab_df3$High_AboveMedian_Pval, decreasing = F, na.last = F),]
    tab_df3_ordered[which(is.na(tab_df3_ordered$High_AboveMedian_Pval)),"High_AboveMedian_Pval"] <- "<0.001"
    tab_df3_ordered <- tab_df3_ordered %>%
      relocate(variable,Hazard_Ratio,ci,High_AboveMedian_Pval,Concordance,Likelihood_Ratio_Pval,Wald_Test_Pval,
               Logrank_Test_Pval,Likelihood_Ratio_adjPval_BH,Wald_Test_adjPval_BH,Logrank_Test_adjPval_BH,Criteria,
               Median_Score,Mean_Score,Min_Score,Max_Score)
    new_header <- colnames(tab_df3_ordered)
    colnames(tab_df3_ordered) <- colnames(out_df_top)
    tab_df3_ordered <- rbind(out_df_top,new_header,tab_df3_ordered)
    new_file <- gsub(".txt","_ranked.txt",file_made)
    write_delim(tab_df3_ordered, new_file, delim = '\t', col_names = F, na = "")
    
    
  }
  if (!is.na(Covariate_Column_Label)) {
    out_df2 <- as.data.frame(read_delim(file_made_int,delim = '\t', col_names = F))
    out_df2[c(2:36)] <- NA
    out_df2 <- rbind(out_df2,header)
    out_df3 <- as.data.frame(read_delim(file_made_add,delim = '\t', col_names = F))
    out_df3[c(2:36)] <- NA
    out_df3 <- rbind(out_df3,header)
    write(header,file = file_made_int, append = T, sep = '\t', ncolumns = 36)
    write(header,file = file_made_add, append = T, sep = '\t', ncolumns = 36)
    for (j in calc_cols) {
      meta_ssgsea_sub <- ssGSEA_BIN_meta[,c(Survival_Time,Survival_ID,j,Covariate_Column_Label)]
      # factor hi/lo column
      meta_ssgsea_sub[,j] <- as.factor(meta_ssgsea_sub[,j])
      meta_ssgsea_sub[,j] <- relevel(meta_ssgsea_sub[,j], ref = "Low_BelowMedian")
      # factor reference column
      meta_ssgsea_sub[,Covariate_Column_Label] <- as.factor(meta_ssgsea_sub[,Covariate_Column_Label])
      if (is.na(Covariate_Reference)) {
        current_ref <- levels(meta_ssgsea_sub[,Covariate_Column_Label])[1]
        print(paste("No Covariate Reference Given. Using ",current_ref," as Coxh reference variable.",sep = ""))
      }
      else {
        meta_ssgsea_sub[,Covariate_Column_Label] <- relevel(meta_ssgsea_sub[,Covariate_Column_Label], ref = Covariate_Reference)
      }
      if (length(levels(as.factor(meta_ssgsea_sub[,j]))) > 1) {
        if (!is.na(as.numeric(str_sub(j,1,1)))) {
          colnames(meta_ssgsea_sub)[3] <- paste("n",j,sep = "")
          ## Additive
          temp_tab_add <- coxph(as.formula(paste("Surv(",Survival_Time,",",Survival_ID,") ~ ",paste("n",j,sep = "")," + ",Covariate_Column_Label,sep = "")),
                                data = meta_ssgsea_sub) %>% 
            gtsummary::tbl_regression(exp = TRUE) %>%
            as_gt()
          temp_tab_df_add <- as.data.frame(temp_tab_add)
          tab_add <- coxph(as.formula(paste("Surv(",Survival_Time,",",Survival_ID,") ~ ",paste("n",j,sep = "")," + ",Covariate_Column_Label,sep = "")),
                           data = meta_ssgsea_sub)
          out_add <- capture.output(summary(tab_add))
          con_line_add <- grep("^Concordance=",out_add,value = T)
          con_v_add <- str_split(con_line_add," ")[[1]][2]
          lik_line_add <- grep("^Likelihood ratio test=",out_add,value = T)
          lik_p_add <- str_split(lik_line_add,"=")[[1]][3]
          wal_line_add <- grep("^Wald test",out_add,value = T)
          wal_p_add <- str_split(wal_line_add,"=")[[1]][3]
          sco_line_add <- grep("^Score ",out_add,value = T)
          sco_p_add <- str_split(sco_line_add,"=")[[1]][3]
          temp_tab_df_add[3,c(1,2,13)] <- sub(".","",temp_tab_df_add[3,c(1,2,13)])
          temp_tab_vect_add <- as.character(c(temp_tab_df_add[3,]))
          temp_tab_vect_add <- c(temp_tab_vect_add,con_v_add,lik_p_add,wal_p_add,sco_p_add)
          
          variable <- file_path_sans_ext(temp_tab_vect_add[1])
          if (gs_name == "Genes") {
            score_df <- ssGSEA_GeneSumm_tabs[[1]]
          }
          if (gs_name != "Genes") {
            score_df <- ssGSEA_summ_tabs[[i]]
          }
          var_scores <- unname(unlist(c(score_df[variable,c(2:5)])))
          
          temp_tab_vect_add <- c(temp_tab_vect_add,var_scores)
          
          out_df2 <- rbind(out_df2,temp_tab_vect_add)
          write(temp_tab_vect_add,file = file_made_add,append = T, sep = '\t', ncolumns = 36)
          ## Interactive
          temp_tab_int <- coxph(as.formula(paste("Surv(",Survival_Time,",",Survival_ID,") ~ ",paste("n",j,sep = "")," * ",Covariate_Column_Label,sep = "")),
                                data = meta_ssgsea_sub) %>% 
            gtsummary::tbl_regression(exp = TRUE) %>%
            as_gt()
          temp_tab_df_int <- as.data.frame(temp_tab_int)
          tab_int <- coxph(as.formula(paste("Surv(",Survival_Time,",",Survival_ID,") ~ ",paste("n",j,sep = "")," * ",Covariate_Column_Label,sep = "")),
                           data = meta_ssgsea_sub)
          out_int <- capture.output(summary(tab_int))
          con_line_int <- grep("^Concordance=",out_int,value = T)
          con_v_int <- str_split(con_line_int," ")[[1]][2]
          lik_line_int <- grep("^Likelihood ratio test=",out_int,value = T)
          lik_p_int <- str_split(lik_line_int,"=")[[1]][3]
          wal_line_int <- grep("^Wald test",out_int,value = T)
          wal_p_int <- str_split(wal_line_int,"=")[[1]][3]
          sco_line_int <- grep("^Score ",out_int,value = T)
          sco_p_int <- str_split(sco_line_int,"=")[[1]][3]
          temp_tab_df_int[3,c(1,2,13)] <- sub(".","",temp_tab_df_int[3,c(1,2,13)])
          temp_tab_vect_int <- as.character(c(temp_tab_df_int[3,]))
          temp_tab_vect_int <- c(temp_tab_vect_int,con_v_int,lik_p_int,wal_p_int,sco_p_int)
          
          variable <- file_path_sans_ext(temp_tab_vect_int[1])
          if (gs_name == "Genes") {
            score_df <- ssGSEA_GeneSumm_tabs[[1]]
          }
          if (gs_name != "Genes") {
            score_df <- ssGSEA_summ_tabs[[i]]
          }
          var_scores <- unname(unlist(c(score_df[variable,c(2:5)])))
          
          temp_tab_vect_int <- c(temp_tab_vect_int,var_scores)
          
          out_df3 <- rbind(out_df3,temp_tab_vect_int)
          write(temp_tab_vect_int,file = file_made_int, append = T, sep = '\t', ncolumns = 36)
        }
        if (is.na(as.numeric(str_sub(j,1,1)))) {
          ## Additive
          temp_tab_add <- coxph(as.formula(paste("Surv(",Survival_Time,",",Survival_ID,") ~ ",j," + ",Covariate_Column_Label,sep = "")),
                                data = meta_ssgsea_sub) %>% 
            gtsummary::tbl_regression(exp = TRUE) %>%
            as_gt()
          temp_tab_df_add <- as.data.frame(temp_tab_add)
          tab_add <- coxph(as.formula(paste("Surv(",Survival_Time,",",Survival_ID,") ~ ",j," + ",Covariate_Column_Label,sep = "")),
                           data = meta_ssgsea_sub)
          out_add <- capture.output(summary(tab_add))
          con_line_add <- grep("^Concordance=",out_add,value = T)
          con_v_add <- str_split(con_line_add," ")[[1]][2]
          lik_line_add <- grep("^Likelihood ratio test=",out_add,value = T)
          lik_p_add <- str_split(lik_line_add,"=")[[1]][3]
          wal_line_add <- grep("^Wald test",out_add,value = T)
          wal_p_add <- str_split(wal_line_add,"=")[[1]][3]
          sco_line_add <- grep("^Score ",out_add,value = T)
          sco_p_add <- str_split(sco_line_add,"=")[[1]][3]
          temp_tab_vect_add <- as.character(c(temp_tab_df_add[3,]))
          temp_tab_vect_add <- c(temp_tab_vect_add,con_v_add,lik_p_add,wal_p_add,sco_p_add)
          
          variable <- file_path_sans_ext(temp_tab_vect_add[1])
          if (gs_name == "Genes") {
            score_df <- ssGSEA_GeneSumm_tabs[[1]]
          }
          if (gs_name != "Genes") {
            score_df <- ssGSEA_summ_tabs[[i]]
          }
          var_scores <- unname(unlist(c(score_df[variable,c(2:5)])))
          
          temp_tab_vect_add <- c(temp_tab_vect_add,var_scores)
          
          out_df2 <- rbind(out_df2,temp_tab_vect_add)
          write(temp_tab_vect_add,file = file_made_add, append = T, sep = '\t', ncolumns = 36)
          ## Interactive
          temp_tab_int <- coxph(as.formula(paste("Surv(",Survival_Time,",",Survival_ID,") ~ ",j," * ",Covariate_Column_Label,sep = "")),
                                data = meta_ssgsea_sub) %>% 
            gtsummary::tbl_regression(exp = TRUE) %>%
            as_gt()
          temp_tab_df_int <- as.data.frame(temp_tab_int)
          tab_int <- coxph(as.formula(paste("Surv(",Survival_Time,",",Survival_ID,") ~ ",j," * ",Covariate_Column_Label,sep = "")),
                           data = meta_ssgsea_sub)
          out_int <- capture.output(summary(tab_int))
          con_line_int <- grep("^Concordance=",out_int,value = T)
          con_v_int <- str_split(con_line_int," ")[[1]][2]
          lik_line_int <- grep("^Likelihood ratio test=",out_int,value = T)
          lik_p_int <- str_split(lik_line_int,"=")[[1]][3]
          wal_line_int <- grep("^Wald test",out_int,value = T)
          wal_p_int <- str_split(wal_line_int,"=")[[1]][3]
          sco_line_int <- grep("^Score ",out_int,value = T)
          sco_p_int <- str_split(sco_line_int,"=")[[1]][3]
          temp_tab_vect_int <- as.character(c(temp_tab_df_int[3,]))
          temp_tab_vect_int <- c(temp_tab_vect_int,con_v_int,lik_p_int,wal_p_int,sco_p_int)
          
          variable <- file_path_sans_ext(temp_tab_vect_int[1])
          if (gs_name == "Genes") {
            score_df <- ssGSEA_GeneSumm_tabs[[1]]
          }
          if (gs_name != "Genes") {
            score_df <- ssGSEA_summ_tabs[[i]]
          }
          var_scores <- unname(unlist(c(score_df[variable,c(2:5)])))
          
          temp_tab_vect_int <- c(temp_tab_vect_int,var_scores)
          
          out_df3 <- rbind(out_df3,temp_tab_vect_int)
          write(temp_tab_vect_int,file = file_made_int, append = T, sep = '\t', ncolumns = 36)
        }
      }
    }
    ## Additive
    out_df2_top <- out_df2[grep("##",out_df2[,1]),]
    out_df2_top[,c(37:39)] <- NA
    tab_df1 <- out_df2[grep("##",out_df2[,1],invert = T),]
    colnames(tab_df1) <- header
    tab_df1 <- tab_df1[-1,]
    Likelihood_Ratio_adjPval_BH <- p.adjust(as.numeric(tab_df1$Likelihood_Ratio_Pval), method = "BH")
    Wald_Test_adjPval_BH <- p.adjust(as.numeric(tab_df1$Wald_Test_Pval), method = "BH")
    Logrank_Test_adjPval_BH <- p.adjust(as.numeric(tab_df1$Logrank_Test_Pval), method = "BH")
    tab_df1 <- cbind(tab_df1,Likelihood_Ratio_adjPval_BH,Wald_Test_adjPval_BH,Logrank_Test_adjPval_BH)
    tab_df1$variable <- gsub(".BIN$","",tab_df1$variable)
    tab_df1$High_AboveMedian_Pval <- gsub(">0.9","0.9",tab_df1$High_AboveMedian_Pval)
    tab_df1$High_AboveMedian_Pval <- as.numeric(tab_df1$High_AboveMedian_Pval)
    tab_df1$Likelihood_Ratio_Pval <- as.numeric(tab_df1$Likelihood_Ratio_Pval)
    tab_df1_ordered <- tab_df1[order(tab_df1$Likelihood_Ratio_Pval, decreasing = F),]
    tab_df1_ordered$Likelihood_Ratio_Pval <- formatC(tab_df1_ordered$Likelihood_Ratio_Pval)
    tab_df1_ordered$Wald_Test_Pval <- formatC(as.numeric(tab_df1_ordered$Wald_Test_Pval))
    tab_df1_ordered$Logrank_Test_Pval <- formatC(as.numeric(tab_df1_ordered$Logrank_Test_Pval))
    tab_df1_ordered[which(is.na(tab_df1_ordered$High_AboveMedian_Pval)),"High_AboveMedian_Pval"] <- "<0.001"
    tab_df1_ordered <- tab_df1_ordered %>%
      relocate(variable,Hazard_Ratio,ci,High_AboveMedian_Pval,Concordance,Likelihood_Ratio_Pval,Wald_Test_Pval,
               Logrank_Test_Pval,Likelihood_Ratio_adjPval_BH,Wald_Test_adjPval_BH,Logrank_Test_adjPval_BH,Criteria,
               Median_Score,Mean_Score,Min_Score,Max_Score)
    new_header <- colnames(tab_df1_ordered)
    colnames(tab_df1_ordered) <- colnames(out_df2_top)
    tab_df1_ordered <- rbind(out_df2_top,new_header,tab_df1_ordered)
    new_file_add <- gsub(".txt","_ranked.txt",file_made_add)
    write_delim(tab_df1_ordered, new_file_add, delim = '\t', col_names = F, na = "")
    
    
    
    ## Interactive
    out_df3_top <- out_df3[grep("##",out_df3[,1]),]
    out_df3_top[,c(37:39)] <- NA
    tab_df2 <- out_df3[grep("##",out_df3[,1],invert = T),]
    colnames(tab_df2) <- header
    tab_df2 <- tab_df2[-1,]
    Likelihood_Ratio_adjPval_BH <- p.adjust(as.numeric(tab_df2$Likelihood_Ratio_Pval), method = "BH")
    Wald_Test_adjPval_BH <- p.adjust(as.numeric(tab_df2$Wald_Test_Pval), method = "BH")
    Logrank_Test_adjPval_BH <- p.adjust(as.numeric(tab_df2$Logrank_Test_Pval), method = "BH")
    tab_df2 <- cbind(tab_df2,Likelihood_Ratio_adjPval_BH,Wald_Test_adjPval_BH,Logrank_Test_adjPval_BH)
    tab_df2$variable <- gsub(".BIN$","",tab_df2$variable)
    tab_df2$High_AboveMedian_Pval <- gsub(">0.9","0.9",tab_df2$High_AboveMedian_Pval)
    tab_df2$High_AboveMedian_Pval <- as.numeric(tab_df2$High_AboveMedian_Pval)
    tab_df2$Likelihood_Ratio_Pval <- as.numeric(tab_df2$Likelihood_Ratio_Pval)
    tab_df2_ordered <- tab_df2[order(tab_df2$Likelihood_Ratio_Pval, decreasing = F),]
    tab_df2_ordered$Likelihood_Ratio_Pval <- formatC(tab_df2_ordered$Likelihood_Ratio_Pval)
    tab_df2_ordered$Wald_Test_Pval <- formatC(as.numeric(tab_df2_ordered$Wald_Test_Pval))
    tab_df2_ordered$Logrank_Test_Pval <- formatC(as.numeric(tab_df2_ordered$Logrank_Test_Pval))
    tab_df2_ordered[which(is.na(tab_df2_ordered$High_AboveMedian_Pval)),"High_AboveMedian_Pval"] <- "<0.001"
    tab_df2_ordered <- tab_df2_ordered %>%
      relocate(variable,Hazard_Ratio,ci,High_AboveMedian_Pval,Concordance,Likelihood_Ratio_Pval,Wald_Test_Pval,
               Logrank_Test_Pval,Likelihood_Ratio_adjPval_BH,Wald_Test_adjPval_BH,Logrank_Test_adjPval_BH,Criteria,
               Median_Score,Mean_Score,Min_Score,Max_Score)
    new_header <- colnames(tab_df2_ordered)
    colnames(tab_df2_ordered) <- colnames(out_df3_top)
    tab_df2_ordered <- rbind(out_df3_top,new_header,tab_df2_ordered)
    new_file_int <- gsub(".txt","_ranked.txt",file_made_int)
    write_delim(tab_df2_ordered, new_file_int, delim = '\t', col_names = F, na = "")
    
  }
}







