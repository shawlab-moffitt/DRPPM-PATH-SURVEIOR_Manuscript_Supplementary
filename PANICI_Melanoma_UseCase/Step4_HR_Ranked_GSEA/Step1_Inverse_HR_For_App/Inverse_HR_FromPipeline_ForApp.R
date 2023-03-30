


####----Pre-Process Input data for Hazard Ratio ranked GSEA----####

## Set "Hazard_Ratio_Ranked_GSEA" to working directory

## Read in files
OS_Ranked_Genes <- read.delim("Input/ICI_iAtlas_Skin_Pre_OS_Genes_coxh_ranked.txt",
                              sep = '\t')
PFI_Ranked_Genes <- read.delim("Input/ICI_iAtlas_Skin_Pre_PFI_Genes_coxh_ranked.txt",
                               sep = '\t')

## Get inverse hazard ratio
OS_Ranked_Genes$Hazard_Ratio <- round(1/OS_Ranked_Genes$Hazard_Ratio,2)
PFI_Ranked_Genes$Hazard_Ratio <- round(1/PFI_Ranked_Genes$Hazard_Ratio,2)


## Write out new files
write.table(OS_Ranked_Genes,
            file = "Output/ICI_iAtlas_Skin_Pre_OS_Genes_coxh_ranked_Inverse.txt",
            sep = '\t',
            row.names = F)
write.table(PFI_Ranked_Genes,
            file = "Output/ICI_iAtlas_Skin_Pre_PFI_Genes_coxh_ranked_Inverse.txt",
            sep = '\t',
            row.names = F)
