
####----Subsetting Samples For Bootstrap Analysis----####

## Read in  libraries
library(readr)
library(dplyr)

## Set "Step2_Subset_Data_For_Bootstrap" as working directory


## File paths
Expression_File <- "Input/AML_TARGET_KMT2A_fusion_expr.txt"

Meta_File <- "Input/AML_TARGET_KMT2A_fusion_meta.txt"

Meta_Param <- "Input/AML_TARGET_KMT2A_fusion_meta_parameters.txt"


## Read in files
expr <- as.data.frame(read_delim(Expression_File,delim = '\t', col_names = T))
meta <- as.data.frame(read_delim(Meta_File,delim = '\t', col_names = T))
metaP<- as.data.frame(read_delim(Meta_Param,delim = '\t', col_names = F))


## Subset 75% of samples to bootstrap

Bootstrap_Sample_Num <- round(length(meta$SampleName) * (75/100))
Bootstrap_Samples <- sample(meta$SampleName,Bootstrap_Sample_Num)

Bootstrap_Expr <- expr[,c("GeneSymbol",Bootstrap_Samples)]
Bootstrap_Meta <- meta[which(meta$SampleName %in% Bootstrap_Samples),]


## Write bootstrap subset data to files
write_delim(Bootstrap_Expr,"Output/AML_TARGET_KMT2A75_expr.txt", delim = '\t')
write_delim(Bootstrap_Meta,"Output/AML_TARGET_KMT2A75_meta.txt", delim = '\t')


## Add bootstrap annotation to meta data
meta$BootstrapSamples_75 <- ifelse((meta$SampleName %in% Bootstrap_Samples), "TRUE", "FALSE")

metaP <- rbind(metaP,
               data.frame(X1 = "BootstrapSamples_75",
                          X2 = "Feature"))


write_delim(meta,"Output/AML_TARGET_KMT2A_fusion_meta_BootstrapAnno.txt",delim = '\t')
write_delim(metaP,"Output/AML_TARGET_KMT2A_fusion_meta_param_BootstrapAnno.txt",delim = '\t', col_names = F)


