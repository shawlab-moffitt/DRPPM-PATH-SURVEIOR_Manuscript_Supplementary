

####----Command Line Arguments----####

args=commandArgs(trailingOnly = T)

Param_File <- args[1]




####----Packages----####

# Base Packages
packages <- c("gtsummary","tidyr","dplyr","DT","ggpubr","tibble","stringr",
              "survival","readr","survminer","gridExtra","BiocManager")
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

####----Read in files----####

params <- as.data.frame(read_delim(Param_File,delim = "\t", col_names = F))

# Project Name
Project_Name <- params[which(params[,1] == "Project_Name"),2]
Project_Name <- gsub(" ","_",Project_Name)
Project_Name <- gsub("[[:punct:]]","_",Project_Name)
# Expression File
Expression_Matrix_file <- params[which(params[,1] == "Expression_Matrix_file"),2]
# Meta File
Meta_Data_File <- params[which(params[,1] == "Meta_Data_File"),2]
# Gene Set File
Gene_Set_File <- params[which(params[,1] == "Gene_Set_File"),2]
# Output File Path
Output_File_Path <- params[which(params[,1] == "Output_File_Path"),2]
last_char <- str_sub(Output_File_Path,-1,-1)
if (last_char != "/") {
  Output_File_Path <- paste(Output_File_Path,"/",sep = "")
}
# Survival Time ID
Survival_Time_Label <- params[which(params[,1] == "Survival_Time_Label"),2]
surv_time_labs <- unlist(str_split(Survival_Time_Label,","))
# Survival ID ID
Survival_ID_Label <- params[which(params[,1] == "Survival_ID_Label"),2]
surv_ID_labs <- unlist(str_split(Survival_ID_Label,","))
# Rank By Genes Option
Rank_Genes_Choice <- params[which(params[,1] == "Rank_Genes"),2]
# Covariate Column Name
Covariate_Column_Label <- params[which(params[,1] == "Covariate_Column_Label"),2]
# covariate reference
Covariate_Reference <- params[which(params[,1] == "Covariate_Reference"),2]


####----Read in Files----####

print("####----Loading in Expression Matirx----####")

##--Expression--##
expr <- as.data.frame(read_delim(Expression_Matrix_file,delim = '\t', col_names = T))
rownames(expr) <- expr[,1]
expr <- as.matrix(expr[,-1])

##--Meta--##
meta <- as.data.frame(read_delim(Meta_Data_File,delim = '\t',col_names = T))
colnames(meta)[1] <- "SampleName"

if (length(Covariate_Column_Label) > 0) {
  Covariate_Column_Label <- NA
  Covariate_Reference <- NA
}

## Remove samples with NA in the survival data columns
meta_tabs <- list()
expr_tabs <- list()
surv_pairs <- list()
for (j in 1:length(surv_time_labs)) {
  
  Survival_Time <- surv_time_labs[j]
  Survival_ID <- surv_ID_labs[j]
  surv_ID <- unlist(str_split(Survival_Time,"[.]",))[1]
  surv_pairs[[surv_ID]] <- c(Survival_Time,Survival_ID)
  
  meta2 <- meta[!is.na(meta[,Survival_Time]),]
  meta2 <- meta2[!is.na(meta2[,Survival_ID]),]
  if (!is.na(Covariate_Column_Label)) {
    meta2 <- meta2[!is.na(meta2[,Covariate_Column_Label]),]
  }
  expr2 <- expr[,which(colnames(expr) %in% meta2[,1])]
  
  meta_tabs[[paste("meta",surv_ID,sep = "_")]] <- meta2
  expr_tabs[[paste("expr",surv_ID,sep = "_")]] <- expr2
}



print("####----Loading in Gene Set----####")

##--Gene Set--##

GeneSet_list <- list()

#if (Rank_Genes_Choice == TRUE) {
#  
#  gmt <- data.frame(term = rownames(expr),gene = rownames(expr))
#  GeneSetList <- list()
#  for (i in unique(gmt[,1])){
#    GeneSetList[[i]] <- gmt[gmt[,1] == i,]$gene
#  }
#  GeneSet_list[["Genes"]] <- GeneSetList
#  
#  
#}

if (!is.na(Gene_Set_File)) {
  
  if (tools::file_ext(Gene_Set_File) == "lst") {
    
    GeneSetFiles <- as.data.frame(read_delim(Gene_Set_File, delim = '\t',col_names = F))
    colnames(GeneSetFiles) <- c("GeneSetName","GeneSetFilePath")
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
    
    # If user loads GMT File
    if (tools::file_ext(Gene_Set_File) == "gmt") {
      GeneSet <- read.gmt(Gene_Set_File)
      colnames(GeneSet) <- c("term","gene")
      GeneSetList <- list()
      for (i in unique(GeneSet[,1])){
        GeneSetList[[i]] <- GeneSet[GeneSet[,1] == i,]$gene
      }
      GeneSet_list[["UserGeneSet"]] <- GeneSetList
    }
    # If user loads TSV/TXT file
    if (tools::file_ext(Gene_Set_File) == "tsv" || tools::file_ext(Gene_Set_File) == "txt") {
      GeneSet <- read.delim(Gene_Set_File, header = T, sep = '\t')
      colnames(GeneSet) <- c("term","gene")
      GeneSetList <- list()
      for (i in unique(GeneSet[,1])){
        GeneSetList[[i]] <- GeneSet[GeneSet[,1] == i,]$gene
      }
      GeneSet_list[["UserGeneSet"]] <- GeneSetList
    }
    # If user loads RData list file
    if (tools::file_ext(Gene_Set_File) == "RData") {
      loadRData <- function(fileName){
        #loads an RData file, and returns it
        load(fileName)
        get(ls()[ls() != "fileName"])
      }
      GeneSetList <- loadRData(Gene_Set_File)
      GeneSet_list[["UserGeneSet"]] <- GeneSetList
    }
    
  }
  
}



## High-Low
highlow = function(mat) {
  new_mat = mat;
  new_mat[mat > quantile(mat)[3]] = "High_AboveMedian";
  new_mat[mat <= quantile(mat)[3]] = "Low_BelowMedian";
  return (new_mat)
}

####----Perform Functions----####


##--ssGSEA--##
print("####----Performing ssGSEA Analysis----####")
ssGSEA_tabs <- list()
if (length(GeneSet_list) > 0) {
  
  for (i in 1:length(GeneSet_list)){
    
    gs_name <- names(GeneSet_list)[i]
    gs <- GeneSet_list[[i]]
    
    for (j in 1:length(surv_pairs)) {
      
      Survival_Time <- surv_pairs[[j]][1]
      Survival_ID <- surv_pairs[[j]][2]
      surv_ID <- names(surv_pairs)[j]
      
      ssgsea <- gsva(expr_tabs[[j]],gs,method = "ssgsea")  #Outputs sample names as column names
      ssgsea <- as.data.frame(t(ssgsea))                  #Transpose to make column names gs names
      ssgsea$SampleName <- rownames(ssgsea)               #Change gs names to row names
      ssgsea <- ssgsea %>%                                #Make rownames first column to output to file
        relocate(SampleName)
      ssGSEA_tabs[[paste("ssgsea",gs_name,surv_ID,sep = "_")]] <- ssgsea
      if (gs_name == "UserGeneSet") {
        write_delim(as.data.frame(ssgsea),paste(Output_File_Path,Project_Name,"_",surv_ID,"_ssGSEA_score.txt", sep = ""), delim = '\t')
      }
      if (gs_name != "UserGeneSet") {
        write_delim(as.data.frame(ssgsea),paste(Output_File_Path,Project_Name,"_",surv_ID,"_",gs_name,"_ssGSEA_score.txt", sep = ""), delim = '\t')
      }
    }
    
  }
  
}





####----Binary Analysis----####

print("####----Performing Binary Analysis----####")

ssGSEA_BIN_tabs <- list()
if (Rank_Genes_Choice == TRUE) {
  
  for (j in 1:length(surv_pairs)) {
    
    Survival_Time <- surv_pairs[[j]][1]
    Survival_ID <- surv_pairs[[j]][2]
    surv_ID <- names(surv_pairs)[j]
    
    expr_temp <- expr_tabs[[j]]
    expr_temp2 <- as.data.frame(t(expr_temp))
    ssgsea <- expr_temp2
    
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
    ssGSEA_BIN_tabs[[paste("ssGSEA.BIN","Genes",surv_ID,sep = "_")]] <- ssGSEA_BIN
    write_delim(ssGSEA_BIN,paste(Output_File_Path,Project_Name,"_",surv_ID,"_Genes","_BIN.txt", sep = ""), delim = '\t')
    
  }
  
}
if (length(ssGSEA_tabs) > 0) {
  
  for (i in 1:length(ssGSEA_tabs)) {
    
    ssgsea <- ssGSEA_tabs[[i]]
    table_name <- names(ssGSEA_tabs)[i]
    gs_name <- unlist(str_split(table_name,"_",))[2]
    surv_ID <- unlist(str_split(table_name,"_",))[3]
    
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
    ssGSEA_BIN_tabs[[paste("ssGSEA.BIN",gs_name,surv_ID,sep = "_")]] <- ssGSEA_BIN
    if (gs_name == "UserGeneSet") {
      write_delim(ssGSEA_BIN,paste(Output_File_Path,Project_Name,"_",surv_ID,"_BIN.txt", sep = ""), delim = '\t')
    }
    if (gs_name != "UserGeneSet") {
      write_delim(ssGSEA_BIN,paste(Output_File_Path,Project_Name,"_",surv_ID,"_",gs_name,"_BIN.txt", sep = ""), delim = '\t')
    }
    
  }
  
}


####----Coxh Analysis----####

print("####----Performing Coxh Analysis----####")

ssGSEA_BIN_meta_tabs <- list()

for (i in 1:length(ssGSEA_BIN_tabs)) {
  
  ssGSEA_BIN <- ssGSEA_BIN_tabs[[i]]
  table_name <- names(ssGSEA_BIN_tabs)[i]
  gs_name <- unlist(str_split(table_name,"_",))[2]
  surv_ID <- unlist(str_split(table_name,"_",))[3]
  Survival_Time <- surv_pairs[[grep(surv_ID,names(surv_pairs),value = T)]][1]
  Survival_ID <- surv_pairs[[grep(surv_ID,names(surv_pairs),value = T)]][2]
  
  ## Clean Gene Set names
  colnames(ssGSEA_BIN) <- gsub("[[:punct:]]",".",colnames(ssGSEA_BIN))
  colnames(ssGSEA_BIN) <- gsub(" ",".",colnames(ssGSEA_BIN))
  
  ## Merge with Meta data
  ssGSEA_BIN_meta <- merge(meta,ssGSEA_BIN,by = "SampleName",all.y = T)
  ssGSEA_BIN_meta_tabs[[paste("ssGSEA.BIN.meta",gs_name,surv_ID,sep = "_")]] <- ssGSEA_BIN_meta
  
  ## Binary columns to loop through
  calc_cols <- grep("..BIN$", colnames(ssGSEA_BIN_meta), value = T)
  
  ## Write out header
  header <- c("variable","var_label","var_type","reference_row","row_type","header_row",
              "N_obs","N_event","N","coefficients_type","coefficients_label","Criteria","term",
              "var_class","var_nlevels","contrasts","contrasts_type","n_obs","n_event",
              "exposure","Hazard_Ratio","std.error","statistic","nevent","conf.low",
              "conf.high","ci","p.value","Concordance","Likelihood_Ratio_Pval","Wald_Test_Pval",
              "Logrank_Test_Pval","Likelihood_Ratio_AdjPval_BH","Wald_Test_AdjPval_BH","Logrank_Test_AdjPval_BH")
  if (gs_name == "UserGeneSet") {
    write(header,file = paste(Output_File_Path,Project_Name,"_",surv_ID,"_coxh.txt", sep = ""),
          append = T, sep = '\t', ncolumns = 35)
  }
  if (gs_name != "UserGeneSet") {
    write(header,file = paste(Output_File_Path,Project_Name,"_",surv_ID,"_",gs_name,"_coxh.txt", sep = ""),
          append = T, sep = '\t', ncolumns = 35)
  }
  
  out_df <- data.frame(matrix(header,1))
  colnames(out_df) <- out_df[1,]
  out_df2 <- data.frame(matrix(header,1))
  colnames(out_df2) <- out_df2[1,]
  out_df3 <- data.frame(matrix(header,1))
  colnames(out_df3) <- out_df3[1,]
  
  
  ## For non-interactive pathway analysis
  if (is.na(Covariate_Column_Label)) {
    
    ## Loop through to perform coxh function
    for (j in calc_cols) {
      
      meta_ssgsea_sub <- ssGSEA_BIN_meta[,c(Survival_Time,Survival_ID,j)]
      meta_ssgsea_sub[,j] <- as.factor(meta_ssgsea_sub[,j])
      meta_ssgsea_sub[,j] <- relevel(meta_ssgsea_sub[,j], ref = "Low_BelowMedian")
      
      if (length(levels(as.factor(meta_ssgsea_sub[,j]))) > 1) {
        
        if (!is.na(as.numeric(str_sub(j,1,1)))) {
          j <- paste("n",j,sep = "")
          colnames(meta_ssgsea_sub)[3] <- paste("n",j,sep = "")
        }
        
        temp_tab <- coxph(as.formula(paste("Surv(",Survival_Time,",",Survival_ID,") ~ ",j,sep = "")),
                          data = meta_ssgsea_sub) %>% 
          gtsummary::tbl_regression(exp = TRUE) %>%
          as_gt()
        temp_tab_df <- as.data.frame(temp_tab)
        tab <- coxph(as.formula(paste("Surv(",Survival_Time,",",Survival_ID,") ~ ",j,sep = "")),
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
        adj.p <- p.adjust(as.numeric(c(lik_p,wal_p,sco_p)),method = "BH")
        if (!is.na(as.numeric(str_sub(j,1,1)))) {
          temp_tab_df[3,c(1,2,13)] <- sub(".","",temp_tab_df[3,c(1,2,13)])
        }
        temp_tab_vect <- as.character(c(temp_tab_df[3,]))
        temp_tab_vect <- c(temp_tab_vect,con_v,lik_p,wal_p,sco_p,adj.p)
        out_df <- rbind(out_df,temp_tab_vect)
        if (gs_name == "UserGeneSet") {
          write(temp_tab_vect,file = paste(Output_File_Path,Project_Name,"_",surv_ID,"_coxh.txt", sep = ""),
                append = T, sep = '\t', ncolumns = 35)
        }
        if (gs_name != "UserGeneSet") {
          write(temp_tab_vect,file = paste(Output_File_Path,Project_Name,"_",surv_ID,"_",gs_name,"_coxh.txt", sep = ""),
                append = T, sep = '\t', ncolumns = 35)
        }
        
      }
      
    }
    
    out_df <- out_df[-1,]
    tab_df3 <- out_df
    tab_df3$variable <- gsub(".BIN$","",tab_df3$variable)
    tab_df3$p.value <- gsub(">0.9","0.9",tab_df3$p.value)
    tab_df3$p.value <- as.numeric(tab_df3$p.value)
    tab_df3_ordered <- tab_df3[order(tab_df3$p.value, decreasing = F, na.last = F),]
    tab_df3_ordered[which(is.na(tab_df3_ordered$p.value)),"p.value"] <- "<0.001"
    tab_df3_ordered <- tab_df3_ordered %>%
      relocate(variable,Hazard_Ratio,ci,p.value,Concordance,Likelihood_Ratio_Pval,Wald_Test_Pval,Logrank_Test_Pval,Likelihood_Ratio_AdjPval_BH,Wald_Test_AdjPval_BH,Logrank_Test_AdjPval_BH,Criteria)
    if (gs_name == "UserGeneSet") {
      write_delim(tab_df3_ordered, paste(Output_File_Path,Project_Name,"_",surv_ID,"_coxh_ranked.txt", sep = ""), delim = '\t')
    }
    if (gs_name != "UserGeneSet") {
      write_delim(tab_df3_ordered, paste(Output_File_Path,Project_Name,"_",surv_ID,"_",gs_name,"_coxh_ranked.txt", sep = ""), delim = '\t')
    }
    
  }
  
    
  if (!is.na(Covariate_Column_Label)) {
    
    ## Loop through to perform coxh function
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
        meta_ssgsea_sub[,j] <- relevel(meta_ssgsea_sub[,j], ref = Covariate_Reference)
      }
      
      if (length(levels(as.factor(meta_ssgsea_sub[,j]))) > 1) {
        
        if (!is.na(as.numeric(str_sub(j,1,1)))) {
          j <- paste("n",j,sep = "")
          colnames(meta_ssgsea_sub)[3] <- paste("n",j,sep = "")
        }
        
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
        adj.p_add <- p.adjust(as.numeric(c(lik_p_add,wal_p_add,sco_p_add)),method = "BH")
        if (!is.na(as.numeric(str_sub(j,1,1)))) {
          temp_tab_df_add[3,c(1,2,13)] <- sub(".","",temp_tab_df_add[3,c(1,2,13)])
        }
        temp_tab_vect_add <- as.character(c(temp_tab_df_add[3,]))
        temp_tab_vect_add <- c(temp_tab_vect_add,con_v_add,lik_p_add,wal_p_add,sco_p_add,adj.p_add)
        out_df2 <- rbind(out_df2,temp_tab_vect_add)
        if (gs_name == "UserGeneSet") {
          write(temp_tab_vect_add,file = paste(Output_File_Path,Project_Name,"_",surv_ID,"_Additive_coxh.txt", sep = ""),
                append = T, sep = '\t', ncolumns = 35)
        }
        if (gs_name != "UserGeneSet") {
          write(temp_tab_vect_add,file = paste(Output_File_Path,Project_Name,"_",surv_ID,"_",gs_name,"_Additive_coxh.txt", sep = ""),
                append = T, sep = '\t', ncolumns = 35)
        }
        
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
        adj.p_int <- p.adjust(as.numeric(c(lik_p_int,wal_p_int,sco_p_int)),method = "BH")
        if (!is.na(as.numeric(str_sub(j,1,1)))) {
          temp_tab_df_int[3,c(1,2,13)] <- sub(".","",temp_tab_df_int[3,c(1,2,13)])
        }
        temp_tab_vect_int <- as.character(c(temp_tab_df_int[3,]))
        temp_tab_vect_int <- c(temp_tab_vect_int,con_v_int,lik_p_int,wal_p_int,sco_p_int,adj.p_int)
        out_df3 <- rbind(out_df3,temp_tab_vect_int)
        if (gs_name == "UserGeneSet") {
          write(temp_tab_vect_int,file = paste(Output_File_Path,Project_Name,"_",surv_ID,"_Interactive_coxh.txt", sep = ""),
                append = T, sep = '\t', ncolumns = 35)
        }
        if (gs_name != "UserGeneSet") {
          write(temp_tab_vect_int,file = paste(Output_File_Path,Project_Name,"_",surv_ID,"_",gs_name,"_Interactive_coxh.txt", sep = ""),
                append = T, sep = '\t', ncolumns = 35)
        }
        
      }
      
    }
    
    out_df2 <- out_df2[-1,]
    tab_df1 <- out_df2
    ## Rank by P.Value and write to file
    tab_df1$variable <- gsub(".BIN$","",tab_df1$variable)
    tab_df1$p.value <- gsub(">0.9","0.9",tab_df1$p.value)
    tab_df1$p.value <- as.numeric(tab_df1$p.value)
    tab_df1_ordered <- tab_df1[order(tab_df1$p.value, decreasing = F, na.last = F),]
    tab_df1_ordered[which(is.na(tab_df1_ordered$p.value)),"p.value"] <- "<0.001"
    tab_df1_ordered <- tab_df1_ordered %>%
      relocate(variable,Hazard_Ratio,ci,p.value,Concordance,Likelihood_Ratio_Pval,Wald_Test_Pval,Logrank_Test_Pval,Likelihood_Ratio_AdjPval_BH,Wald_Test_AdjPval_BH,Logrank_Test_AdjPval_BH,Criteria)
    if (gs_name == "UserGeneSet") {
      write(tab_df1_ordered,file = paste(Output_File_Path,Project_Name,"_",surv_ID,"_Additive_coxh_ranked.txt", sep = ""),
            append = T, sep = '\t', ncolumns = 35)
    }
    if (gs_name != "UserGeneSet") {
      write(tab_df1_ordered,file = paste(Output_File_Path,Project_Name,"_",surv_ID,"_",gs_name,"_Additive_coxh_ranked.txt", sep = ""),
            append = T, sep = '\t', ncolumns = 35)
    }

    out_df3 <- out_df3[-1,]
    tab_df2 <- out_df3
    ## Rank by P.Value and write to file
    tab_df2$variable <- gsub(".BIN$","",tab_df2$variable)
    tab_df2$p.value <- gsub(">0.9","0.9",tab_df2$p.value)
    tab_df2$p.value <- as.numeric(tab_df2$p.value)
    tab_df2_ordered <- tab_df2[order(tab_df2$p.value, decreasing = F, na.last = F),]
    tab_df2_ordered[which(is.na(tab_df2_ordered$p.value)),"p.value"] <- "<0.001"
    tab_df2_ordered <- tab_df2_ordered %>%
      relocate(variable,Hazard_Ratio,ci,p.value,Concordance,Likelihood_Ratio_Pval,Wald_Test_Pval,Logrank_Test_Pval,Likelihood_Ratio_AdjPval_BH,Wald_Test_AdjPval_BH,Logrank_Test_AdjPval_BH,Criteria)
    if (gs_name == "UserGeneSet") {
      write(tab_df2_ordered,file = paste(Output_File_Path,Project_Name,"_",surv_ID,"_Interactive_coxh_ranked.txt", sep = ""),
            append = T, sep = '\t', ncolumns = 35)
      }
    if (gs_name != "UserGeneSet") {
      write(tab_df2_ordered,file = paste(Output_File_Path,Project_Name,"_",surv_ID,"_",gs_name,"_Interactive_coxh_ranked.txt", sep = ""),
            append = T, sep = '\t', ncolumns = 35)
      }
    }
}

