####----Install and load packages----####

packages <- c("shiny","shinythemes","shinyjqui","pheatmap","RColorBrewer","umap",
              "ggdendro","factoextra","dplyr","DT","viridis","readr","tidyverse","ggrepel",
              "shinycssloaders","stringr","tools","plotly","reshape2","ggpubr","gridExtra")

installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
invisible(lapply(packages, library, character.only = TRUE))
#bioconductor packages
bioCpacks <- c("clusterProfiler")
installed_packages_BIOC <- bioCpacks %in% rownames(installed.packages())
if (any(installed_packages_BIOC == FALSE)) {
  BiocManager::install(bioCpacks[!installed_packages_BIOC], ask = F)
}
invisible(lapply(bioCpacks, library, character.only = TRUE))


####---- Read in Files----####

GeneSet_File <- "GeneSet_List_HS_v5.RData"

#GeneSetCat_File <- "GeneSet_CatTable_v5.zip"

PreSet_CoxH_Ranked_Feature_file <- "ICI_iAtlas_Skin_Pre_ImmuneSig_OS_coxh_ranked_Filtered.txt"


# R Data list load function for naming
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
gs <- loadRData(GeneSet_File)
names(gs) <- gsub("[[:punct:]]",".",names(gs))

#gs_cat <- as.data.frame(read_delim(GeneSetCat_File, delim = '\t', col_names = T))
#genesetCats <- unique(gs_cat[,1])
#gs_cat2 <- gs_cat[,-2]
#colnames(gs_cat2)[c(2,3)] <- c("Gene Set Category","Gene Set Name")

PreSet_CoxH_Ranked_Features <- as.data.frame(read_delim(PreSet_CoxH_Ranked_Feature_file, delim = '\t', col_names = T, comment = "#"))


#learn Jaccard index function
jaccard <- function(a, b) {  
  a = unique(a)
  b = unique(b)
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (1-(intersection/union))
}

#CV function for variance
cv <- function(x){
  (sd(x)/mean(x))*100
}

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

## Initial Pathway Selected
#Path_Selec <- NULL

shinytheme("sandstone")

#increase file upload size
options(shiny.maxRequestSize=5000*1024^2)



server <- function(input, output, session) {
  
  ####----Jaccard Connectivity----####
  
  ####----Render UI----####
  
  output$rendTopFeatureSelect <- renderUI({
    
    #req(input$UserPathwayFile)
    if (is.null(input$UserPathwayFile)) {
      ranked_file <- PreSet_CoxH_Ranked_Features
      if (ncol(ranked_file) > 2) {
        numericInput("TopFeatureSelect","Top Number of Features to Process:",
                     value = 50, min = 3, width = "250px")
      }
    } else {
      header_check <- input$HeaderCheckIntra
      gs.u <- input$UserPathwayFile
      ext <- tools::file_ext(gs.u$datapath)
      ranked_file <- as.data.frame(read_delim(gs.u$datapath, delim = '\t', col_names = header_check, comment = "#"))
      if (ext == "txt" | ext == "tsv") {
        if (ncol(ranked_file) > 2) {
          numericInput("TopFeatureSelect","Top Number of Features to Process:",
                       value = 50, min = 3, width = "250px")
        }
      }
    }
    
    
    
  })
  
  output$rendClustTabIntra <- renderUI({
    
    if (input$ViewClustTabIntra == TRUE) {
      
      div(DT::dataTableOutput("ClustTabIntra"), style = "font-size:10px; height:400px; overflow-X: scroll")
      
    }
    
  })
  
  output$rendSIFPreview <- renderUI({
    
    if (input$PrevSIF == TRUE) {
      
      div(DT::dataTableOutput("SIFPreview"), style = "font-size:10px; height:400px; overflow-X: scroll")
      
    }
    
  })
  
  output$renddnldJaccTabIntra <- renderUI({
    
    #req(input$UserPathwayFile)
    downloadButton("dnldJaccTabIntra","Download Jaccard Connectivity Table")
    
  })
  
  output$renddnldClusterTabAnno <- renderUI({
    
    req(input$UserAnnotationFile)
    downloadButton("dnldClusterTabAnno","Download Cluster Annotation Table")
    
  })
  
  output$rendJaccDendo <- renderUI({
    
    VisType <- input$ConnView
    if (VisType == "rectangle" | VisType == "phylogenic") {
      
      withSpinner(jqui_resizable(plotlyOutput('JaccDendoIntra', width = "100%", height = "900px")), type = 6)
      
    }
    else if (VisType == "circular") {
      
      withSpinner(jqui_resizable(plotOutput('JaccDendoIntraCirc', width = "100%", height = "900px")), type = 6)
      
    }
    
  })
  
  output$rendPhyloLayout <- renderUI({
    
    if (input$ConnView == "phylogenic") {
      
      selectInput("PhyloLayout","Choose Phylogenic Tree Layout:",
                  choices = c("Auto Layout" = "layout.auto",
                              "Force Directed (DrL) Layout" = "layout_with_drl",
                              "Tree Layout" = "layout_as_tree",
                              "GEM Layout" = "layout.gem",
                              "Mulitdimensional Scaling Layout" = "layout.mds",
                              "Large Graph Layout" = "layout_with_lgl"))
      
    }
    
  })
  
  ####----Reactives----####
  
  user_file_loaded <- reactive({
    
    if (is.null(input$UserPathwayFile)) {
      ranked_file <- PreSet_CoxH_Ranked_Features
      ranked_file
    } else {
      header_check <- input$HeaderCheckIntra
      gs.u <- input$UserPathwayFile
      ext <- tools::file_ext(gs.u$datapath)
      req(gs.u)
      validate(need(ext == c("tsv","txt","gmt"), "Please upload .tsv, .txt, or .gmt file"))
      
      if (ext == "gmt") {
        gmt <- read.gmt(gs.u$datapath)
        ranked_file <- data.frame(GeneSet = unique(gmt[,1]))
      }
      else {
        ranked_file <- as.data.frame(read_delim(gs.u$datapath, delim = '\t', col_names = header_check, comment = "#"))
      }
      
      ranked_file[,1] <- gsub("[[:punct:]]",".",ranked_file[,1])
      
      ranked_file
    }
    
  })
  
  ## User Upload gene set list reactive
  user_gs <- reactive({
    
    #req(input$UserPathwayFile)
    
    if (is.null(input$UserPathwayFile)) {
      ranked_file <- PreSet_CoxH_Ranked_Features
      if (ncol(ranked_file) == 2) {
        ranked_file <- ranked_file[which(ranked_file[,2] != "NA"),]
        gs_u2 <- list()
        for (i in unique(ranked_file[,1])){
          gs_u2[[i]] <- ranked_file[ranked_file[,1] == i,][,2]
        }
      }
      else {
        TopFeat <- input$TopFeatureSelect
        paths <- as.vector(ranked_file[,1])
        paths <- gsub("[[:punct:]]",".",paths)
        gs_u <- gs[unique(paths)]
        gs_u2 <- gs_u[c(1:TopFeat)]
        gs_u2 <- Filter(Negate(is.null), gs_u2)
      }
      gs_u2
    } else {
      header_check <- input$HeaderCheckIntra
      gs.u <- input$UserPathwayFile
      ext <- tools::file_ext(gs.u$datapath)
      
      ranked_file <- user_file_loaded()
      
      if (ext == "gmt") {
        gmt <- read.gmt(gs.u$datapath)
        gmt <- gmt[which(gmt[,2] != "NA"),]
        gs_u2 <- list()
        for (i in unique(gmt[,1])){
          gs_u2[[i]] <- gmt[gmt[,1] == i,]$gene
        }
      }
      else if (ncol(ranked_file) == 2) {
        ranked_file <- ranked_file[which(ranked_file[,2] != "NA"),]
        gs_u2 <- list()
        for (i in unique(ranked_file[,1])){
          gs_u2[[i]] <- ranked_file[ranked_file[,1] == i,][,2]
        }
      }
      else {
        TopFeat <- input$TopFeatureSelect
        paths <- as.vector(ranked_file[,1])
        paths <- gsub("[[:punct:]]",".",paths)
        gs_u <- gs[unique(paths)]
        gs_u2 <- gs_u[c(1:TopFeat)]
        gs_u2 <- Filter(Negate(is.null), gs_u2)
      }
      
      gs_u2
    }
    
  })
  
  user_anno <- reactive({
    
    gs.u <- input$UserAnnotationFile
    ext <- tools::file_ext(gs.u$datapath)
    req(gs.u)
    validate(need(ext == c("tsv","txt","csv"), "Please upload .tsv, .txt, or .csv file"))
    if (ext == "txt" | ext == "tsv") {
      
      df <- as.data.frame(read_delim(gs.u$datapath, delim = '\t', col_names = T, comment = "##"))
      
    }
    else if (ext == "csv") {
      
      df <- as.data.frame(read_delim(gs.u$datapath, delim = ',', col_names = T))
      
    }
    
    df
    
  })
  
  ## Jaccard Matrix Generation
  jacc_react_Intra <- reactive({
    
    if (length(user_gs()) > 0) {
      
      gs1 <- user_gs() #User Gene Sets Input
      gs2 <- user_gs() #User Gene Sets Input
      
      jac_list <- list()
      for (i in 1:length(gs1)) {
        
        jac_vect <- c()
        
        for (j in 1:length(gs2)) {
          
          gs1_genes <- gs1[[i]]
          gs2_genes <- gs2[[j]]
          jac_index <- jaccard(gs1_genes,gs2_genes)
          jac_vect <- c(jac_vect,jac_index)
          
        }
        
        jac_list[[paste(names(gs1)[i],sep = "")]] <- jac_vect
        
      }
      
      jac_df_Intra <- do.call(rbind,jac_list)
      colnames(jac_df_Intra) <- names(gs2)
      
      jac_df_Intra
      
    }
    
  })
  
  ClusterTabAnno_react <- reactive({
    
    if (is.null(input$UserAnnotationFile) == F) {
      
      ## user annotation table
      anno_df <- user_anno()
      
    }
    
    
    ## user gs list
    gs_u <- user_gs()
    
    ## cluster table
    ## Data matrix
    jacc_df <- jacc_react_Intra()
    jacc_mat <- as.matrix(jacc_df)
    clust_me <- input$ClustMethodIntra
    hm_res <- pheatmap::pheatmap(jacc_mat,
                                 clustering_method = clust_me,
                                 silent = T)
    jacc_df <- jacc_react_Intra()
    cut_k <- input$NumClusters
    clustTab <- data.frame(Pathways = rownames(jacc_df),
                           Cluster = cutree(hm_res$tree_row,k = cut_k))
    
    ## table of gene set with genes
    nameVector <- unlist(mapply(function(x,y){ rep(y, length(x)) }, gs_u, names(gs_u)))
    resultDF <- cbind.data.frame(Pathways = nameVector,Genes = unlist(gs_u))
    rownames(resultDF) <- 1:nrow(resultDF)
    
    ## double check punct matches
    clustTab$Pathways <- gsub("[[:punct:]]",".",clustTab$Pathways)
    resultDF$Pathways <- gsub("[[:punct:]]",".",resultDF$Pathways)
    
    ## Merge geneset and cluster table
    df2 <- merge(resultDF,clustTab, by = "Pathways")
    
    if (is.null(input$UserAnnotationFile) == F) {
      
      colnames(anno_df)[1] <- "Genes"
      df3 <- merge(df2,anno_df,by = "Genes", all.x = T)
      
    }
    else if (is.null(input$UserAnnotationFile) == T) {
      
      df3 <- df2
      
    }
    
    
    df4 <- df3 %>%
      relocate(Pathways,Cluster,Genes)
    
    df5 <- df4[order(df4$Pathways),]
    df5
    
  })
  
  ####----Data Tables----####
  
  ## Jaccard Matrix Table - INTRA pathway
  output$JaccTableIntra <- DT::renderDataTable({
    
    ## Jaccard Table
    jacc_df <- as.data.frame(jacc_react_Intra())
    colnames(jacc_df) <- gsub("[[:punct:]]", " ",colnames(jacc_df))
    
    DT::datatable(jacc_df,
                  extensions = "FixedColumns",
                  options = list(lengthMenu = c(10,20,50,100,1000,5000,10000),
                                 pageLength = 20,
                                 scrollX = TRUE,
                                 autoWidth = TRUE,
                                 fixedColumns = list(leftColumns = 1)),
                  selection=list(mode = "multiple")) %>%
      formatRound(columns = c(1:ncol(jacc_df)), digits = 4)
    
  })
  
  ## Cluster Table
  output$ClustTabIntra <- DT::renderDataTable({
    
    ## Data matrix
    jacc_df <- jacc_react_Intra()
    jacc_mat <- as.matrix(jacc_df)
    clust_me <- input$ClustMethodIntra
    
    hm_res <- pheatmap::pheatmap(jacc_mat,
                                 clustering_method = clust_me,
                                 silent = T)
    
    jacc_df <- jacc_react_Intra()
    cut_k <- input$NumClusters
    clustTab <- data.frame(Pathways = rownames(jacc_df),
                           Cluster = cutree(hm_res$tree_row,k = cut_k))
    
    DT::datatable(clustTab,
                  options = list(keys = TRUE,
                                 searchHighlight = TRUE,
                                 pageLength = 10,
                                 lengthMenu = c("10", "25", "50", "100")),
                  rownames = F)
    
  })
  
  output$ClusterTabAnno <- DT::renderDataTable({
    
    df5 <- ClusterTabAnno_react()
    
    DT::datatable(df5,
                  options = list(lengthMenu = c(10,20,50,100,1000,5000,10000),
                                 pageLength = 20,
                                 scrollX = TRUE),
                  rownames = F,
                  selection=list(mode = "multiple"))
    
    
  })
  
  output$SIFPreview <- DT::renderDataTable({
    
    jacc_df <- user_matrix_loaded()
    jacc_df <- as.data.frame(jacc_df)
    jacc_df_melt <- melt(jacc_df)
    jacc_cols <- colnames(jacc_df)
    jacc_df_melt$Pathway_B <- rep_len(jacc_cols,length.out = nrow(jacc_df_melt))
    colnames(jacc_df_melt)[c(1,2)] <- c("Pathway_A","Jaccard_Distance")
    jacc_df_melt <- jacc_df_melt[which(round(jacc_df_melt$Jaccard_Distance,4) <= dist_cutoff),]
    jacc_df_melt <- jacc_df_melt[which(jacc_df_melt$Pathway_A != jacc_df_melt$Pathway_B),]
    jacc_df_melt <- jacc_df_melt[!duplicated(t(apply(jacc_df_melt,1,sort))),]
    DT::datatable(jacc_df_melt,
                  options = list(keys = TRUE,
                                 searchHighlight = TRUE,
                                 pageLength = 10,
                                 lengthMenu = c("10", "25", "50", "100")),
                  rownames = F)
    
  })
  
  ####----Plots----####
  
  ## Jaccard Heatmap - INTRA
  JaccHeatmapIntra_react <- reactive({
    
    ## User inputs
    row_font <- input$HeatRowFontIntra
    col_font <- input$HeatColFontIntra
    row_name <- input$HeatRowNamesIntra
    col_name <- input$HeatColNamesIntra
    col_pall <- input$ColorPaletteIntra
    clust_me <- input$ClustMethodIntra
    row_dend <- input$HeatRowDendHeight
    col_dend <- input$HeatColDendHeight
    
    ## Data matrix
    jacc_df <- jacc_react_Intra()
    jacc_mat <- as.matrix(jacc_df)
    
    ## Color choice
    minimum = min(jacc_mat)
    maximum = max(jacc_mat)
    bk = c(seq(minimum,minimum/2, length=100), seq(minimum/2,maximum/2,length=100),seq(maximum/2,maximum,length=100))
    #Heatmap color
    col_sets <- c("OrRd","PuBu","Greens","YlGnBu")
    if (col_pall == "original") {
      HeatMap_Colors <- c("dark blue","blue","white","red", "dark red")
      hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
    }
    else if (col_pall %in% col_sets) {
      HeatMap_Colors <- brewer.pal(n = 5, col_pall)
      hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
    }
    else if (col_pall == "Inferno") {
      hmcols <- inferno(500)
    }
    else if (col_pall == "Viridis") {
      hmcols <- viridis(500)
    }
    else if (col_pall == "Plasma") {
      hmcols <- plasma(500)
    }
    else if (col_pall == "OmniBlueRed") {
      hmcols<- colorRampPalette(c("#1984C5", "#22A7F0", "#63BFF0", "#A7D5ED", "#E2E2E2", "#E1A692", "#DE6E56", "#E14B31", "#C23728"))(length(bk)-1)
    }
    else if (col_pall == "LightBlueBlackRed") {
      hmcols<- colorRampPalette(c("#34C5FD","black","red"))(length(bk)-1)
    }
    else if (col_pall == "GreenBlackRed") {
      hmcols<- colorRampPalette(c("green","black","red"))(length(bk)-1)
    }
    heat <- pheatmap::pheatmap(jacc_mat,
                               fontsize_row = row_font,
                               fontsize_col = col_font,
                               show_rownames = row_name,
                               show_colnames = col_name,
                               treeheight_col = col_dend,
                               treeheight_row = row_dend,
                               angle_col = 90,
                               border_color = NA,
                               clustering_method = clust_me,
                               color = hmcols)
    heat
    
  })
  
  output$JaccHeatmapIntra <- renderPlot({
    
    heat <- JaccHeatmapIntra_react()
    heat
    
  })
  
  JaccDendoIntra_react <- reactive({
    
    VisType <- input$ConnView
    
    if (VisType == "rectangle" | VisType == "phylogenic") {
      
      jacc_mat <- jacc_react_Intra()
      clust_me <- input$ClustMethodIntra
      nclust <- input$NumClusters
      LableChoice <- input$ShowConnLabels
      FontSize <- input$ConnFontSize
      
      
      hc <- hclust(dist(jacc_mat), method = clust_me)
      hc_dend <- as.dendrogram(hc)
      dend_data <- dendro_data(hc_dend)
      lab_order <- label(dend_data)[3]
      lab_vect <- lab_order$label
      clust_df <- data.frame(Pathways = rownames(jacc_mat),cluster = cutree(hc,k = nclust))
      clust_df2 <- clust_df[match(lab_vect,clust_df$Pathways),]
      colnames(clust_df2)[1] <- "label"
      clust_df2$cluster <- as.factor(clust_df2$cluster)
      
      dend_data[["labels"]] <- merge(dend_data[["labels"]],clust_df2, by="label")
      
      if (VisType == "rectangle") {
        
        FontSize <- FontSize * 6
        g <- ggdendrogram(dend_data,
                          rotate = TRUE,
                          theme_dendro = TRUE,
                          leaf_labels = F,
                          segments = T,
                          labels = F) + 
          labs(title = "") +
          scale_y_continuous(expand = c(0.4, 0))
        if (LableChoice == TRUE) {
          g <- g + geom_text(data=label(dend_data), aes(x, y, label=label, hjust=1, color=cluster),size=FontSize)
        }
        
        gp <- ggplotly(g)
        
      }
      else if (VisType == "phylogenic") {
        
        if (is.null(input$PhyloLayout) == TRUE) {
          PhyloLayout <- "layout.auto"
        }
        else if (is.null(input$PhyloLayout) == FALSE) {
          PhyloLayout <- input$PhyloLayout
        }
        
        
        if (LableChoice == TRUE) {
          k = nclust
          gp <- ggplotly(
            fviz_dend(
              hc_dend,
              cex = FontSize,
              k = nclust,
              show_labels = TRUE,
              type = "phylogenic",
              phylo_layout = PhyloLayout)
          )
          for(i in seq(k, 1)){
            gp$x$data[[i+1]]$text <- gp$x$data[[i+1+k]]$text
          }
          
        }
        else if (LableChoice == FALSE) {
          k = nclust
          gp <- ggplotly(
            fviz_dend(
              hc,
              cex = FontSize,
              k = k,
              show_labels = TRUE,
              type = "phylogenic",
              phylo_layout = PhyloLayout)
          )
          gp$x$data[[1]]$text <- NA 
          for(i in seq(k, 1)){
            gp$x$data[[i+1]]$text <- gp$x$data[[i+1+k]]$text
            gp$x$data[[i + 1 + k]] <- NULL
          }
          
        }
        
      }
      gp
    }
    
  })
  
  output$JaccDendoIntra <- renderPlotly({
    
    gp <- JaccDendoIntra_react()
    gp
    
  })
  
  JaccDendoIntraCirc_react <- reactive({
    
    VisType <- input$ConnView
    
    if (VisType == "circular") {
      
      jacc_mat <- jacc_react_Intra()
      clust_me <- input$ClustMethodIntra
      nclust <- input$NumClusters
      LableChoice <- input$ShowConnLabels
      FontSize <- input$ConnFontSize
      
      hc <- hclust(dist(jacc_mat), method = clust_me)
      hc_dend <- as.dendrogram(hc)
      f <- fviz_dend(hc_dend, k = nclust, cex = FontSize, color_labels_by_k = T, type = "circular", show_labels = LableChoice) +
        theme(legend.position="none") +
        labs(title = "") +
        theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.title.x = element_blank(),
              axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank())
      
      
      f
      
    }
    
  })
  
  output$JaccDendoIntraCirc <- renderPlot({
    
    f <- JaccDendoIntraCirc_react()
    f
    
  })
  
  
  ####----Matrix Clustering----####
  
  ####----Render UI----####
  
  output$rendClustTabMat <- renderUI({
    
    if (input$ViewClustTabMat == TRUE) {
      
      div(DT::dataTableOutput("ClustTabMat"), style = "font-size:10px; height:400px; overflow-X: scroll")
      
    }
    
  })
  
  output$rendSIFPreviewMat <- renderUI({
    
    if (input$PrevSIFMat == TRUE) {
      
      div(DT::dataTableOutput("SIFPreviewMat"), style = "font-size:10px; height:400px; overflow-X: scroll")
      
    }
    
  })
  
  output$rendUserMatixCLustering <- renderUI({
    
    VisType <- input$ConnViewMat
    if (VisType == "rectangle" | VisType == "phylogenic") {
      
      withSpinner(jqui_resizable(plotlyOutput('JaccDendoMat', width = "100%", height = "900px")), type = 6)
      
    }
    else if (VisType == "circular") {
      
      withSpinner(jqui_resizable(plotOutput('JaccDendoMatCirc', width = "100%", height = "900px")), type = 6)
      
    }
    
  })
  
  output$rendPhyloLayoutMat <- renderUI({
    
    if (input$ConnViewMat == "phylogenic") {
      
      selectInput("PhyloLayoutMat","Choose Phylogenic Tree Layout:",
                  choices = c("Auto Layout" = "layout.auto",
                              "Force Directed (DrL) Layout" = "layout_with_drl",
                              "Tree Layout" = "layout_as_tree",
                              "GEM Layout" = "layout.gem",
                              "Mulitdimensional Scaling Layout" = "layout.mds",
                              "Large Graph Layout" = "layout_with_lgl"))
      
    }
    
  })
  
  ####----Reactives----####
  
  user_matrix_loaded <- reactive({
    
    gs.u <- input$UserMatrixFile
    ext <- tools::file_ext(gs.u$datapath)
    req(gs.u)
    validate(need(ext == c("tsv","txt","csv"), "Please upload .tsv, .txt, or .csv file"))
    
    if (ext == "csv") {
      matrix.u <- as.data.frame(read_delim(gs.u$datapath, delim = ',', col_names = T))
    }
    else {
      matrix.u <- as.data.frame(read_delim(gs.u$datapath, delim = '\t', col_names = T))
    }
    
    
    matrix.u
    
  })
  
  
  ####----Data Tables----####
  
  output$UserMatrixPreview <- DT::renderDataTable({
    
    df <- user_matrix_loaded()
    
    DT::datatable(df,
                  extensions = "FixedColumns",
                  options = list(lengthMenu = c(10,20,50,100,1000,5000,10000),
                                 pageLength = 20,
                                 scrollX = TRUE,
                                 autoWidth = TRUE,
                                 fixedColumns = list(leftColumns = 1)),
                  rownames = F,
                  selection=list(mode = "multiple")) %>%
      formatRound(columns = c(2:ncol(df)), digits = 4)
    
    
  })
  
  ## Cluster Table
  output$ClustTabMat <- DT::renderDataTable({
    
    ## Data matrix
    df <- user_matrix_loaded()
    rownames(df) <- df[,1]
    df <- df[,-1]
    df_mat <- as.matrix(df)
    clust_me <- input$ClustMethodMat
    
    hm_res <- pheatmap::pheatmap(df_mat,
                                 clustering_method = clust_me,
                                 silent = T)
    
    cut_k <- input$NumClustersMat
    clustTab <- data.frame(Pathways = rownames(df),
                           Cluster = cutree(hm_res$tree_row,k = cut_k))
    
    DT::datatable(clustTab,
                  options = list(keys = TRUE,
                                 searchHighlight = TRUE,
                                 pageLength = 10,
                                 lengthMenu = c("10", "25", "50", "100")),
                  rownames = F)
    
  })
  
  output$SIFPreviewMat <- DT::renderDataTable({
    
    df <- user_matrix_loaded()
    rownames(df) <- df[,1]
    df <- df[,-1]
    df <- as.data.frame(df)
    df_melt <- melt(df)
    df_cols <- colnames(df)
    df_melt$Sample_B <- rep_len(df_cols,length.out = nrow(df_melt))
    colnames(df_melt)[c(1,2)] <- c("Sample_A","Distance")
    df_melt <- df_melt[which(df_melt$Sample_A != df_melt$Sample_B),]
    df_melt <- df_melt[!duplicated(t(apply(df_melt,1,sort))),]
    DT::datatable(df_melt,
                  options = list(keys = TRUE,
                                 searchHighlight = TRUE,
                                 pageLength = 10,
                                 lengthMenu = c("10", "25", "50", "100")),
                  rownames = F)
    
  })
  
  ####----Plots----####
  
  UserMatrixHeatmap_react <- reactive({
    
    ## User inputs
    row_font <- input$HeatRowFontMat
    col_font <- input$HeatColFontMat
    row_name <- input$HeatRowNamesMat
    col_name <- input$HeatColNamesMat
    col_pall <- input$ColorPaletteMat
    clust_me <- input$ClustMethodMat
    row_dend <- input$HeatRowDendHeightMat
    col_dend <- input$HeatColDendHeightMat
    
    ## Data matrix
    jacc_df <- user_matrix_loaded()
    rownames(jacc_df) <- jacc_df[,1]
    jacc_df <- jacc_df[,-1]
    jacc_mat <- as.matrix(jacc_df)
    
    ## Color choice
    minimum = min(jacc_mat)
    maximum = max(jacc_mat)
    bk = c(seq(minimum,minimum/2, length=100), seq(minimum/2,maximum/2,length=100),seq(maximum/2,maximum,length=100))
    #Heatmap color
    col_sets <- c("OrRd","PuBu","Greens","YlGnBu")
    if (col_pall == "original") {
      HeatMap_Colors <- c("dark blue","blue","white","red", "dark red")
      hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
    }
    else if (col_pall %in% col_sets) {
      HeatMap_Colors <- brewer.pal(n = 5, col_pall)
      hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
    }
    else if (col_pall == "Inferno") {
      hmcols <- inferno(500)
    }
    else if (col_pall == "Viridis") {
      hmcols <- viridis(500)
    }
    else if (col_pall == "Plasma") {
      hmcols <- plasma(500)
    }
    else if (col_pall == "OmniBlueRed") {
      hmcols<- colorRampPalette(c("#1984C5", "#22A7F0", "#63BFF0", "#A7D5ED", "#E2E2E2", "#E1A692", "#DE6E56", "#E14B31", "#C23728"))(length(bk)-1)
    }
    else if (col_pall == "LightBlueBlackRed") {
      hmcols<- colorRampPalette(c("#34C5FD","black","red"))(length(bk)-1)
    }
    else if (col_pall == "GreenBlackRed") {
      hmcols<- colorRampPalette(c("green","black","red"))(length(bk)-1)
    }
    heat <- pheatmap::pheatmap(jacc_mat,
                               fontsize_row = row_font,
                               fontsize_col = col_font,
                               show_rownames = row_name,
                               show_colnames = col_name,
                               treeheight_col = col_dend,
                               treeheight_row = row_dend,
                               angle_col = 90,
                               border_color = NA,
                               clustering_method = clust_me,
                               color = hmcols)
    heat
    
  })
  
  ## Jaccard Heatmap - INTRA
  output$UserMatrixHeatmap <- renderPlot({
    
    heat <- UserMatrixHeatmap_react()
    heat
    
  })
  
  JaccDendoMat_react <- reactive({
    
    VisType <- input$ConnViewMat
    
    if (VisType == "rectangle" | VisType == "phylogenic") {
      
      jacc_mat <- user_matrix_loaded()
      rownames(jacc_mat) <- jacc_mat[,1]
      jacc_mat <- jacc_mat[,-1]
      clust_me <- input$ClustMethodMat
      nclust <- input$NumClustersMat
      LableChoice <- input$ShowConnLabelsMat
      FontSize <- input$ConnFontSizeMat
      
      
      hc <- hclust(dist(jacc_mat), method = clust_me)
      hc_dend <- as.dendrogram(hc)
      dend_data <- dendro_data(hc_dend)
      lab_order <- label(dend_data)[3]
      lab_vect <- lab_order$label
      clust_df <- data.frame(Pathways = rownames(jacc_mat),cluster = cutree(hc,k = nclust))
      clust_df2 <- clust_df[match(lab_vect,clust_df$Pathways),]
      colnames(clust_df2)[1] <- "label"
      clust_df2$cluster <- as.factor(clust_df2$cluster)
      
      dend_data[["labels"]] <- merge(dend_data[["labels"]],clust_df2, by="label")
      
      if (VisType == "rectangle") {
        
        FontSize <- FontSize * 6
        g <- ggdendrogram(dend_data,
                          rotate = TRUE,
                          theme_dendro = TRUE,
                          leaf_labels = F,
                          segments = T,
                          labels = F) + 
          labs(title = "") +
          scale_y_continuous(expand = c(0.4, 0))
        if (LableChoice == TRUE) {
          g <- g + geom_text(data=label(dend_data), aes(x, y, label=label, hjust=1, color=cluster),size=FontSize)
        }
        
        gp <- ggplotly(g)
        
      }
      else if (VisType == "phylogenic") {
        
        if (is.null(input$PhyloLayout) == TRUE) {
          PhyloLayout <- "layout.auto"
        }
        else if (is.null(input$PhyloLayout) == FALSE) {
          PhyloLayout <- input$PhyloLayout
        }
        
        
        if (LableChoice == TRUE) {
          k = nclust
          gp <- ggplotly(
            fviz_dend(
              hc_dend,
              cex = FontSize,
              k = nclust,
              show_labels = TRUE,
              type = "phylogenic",
              phylo_layout = PhyloLayout)
          )
          for(i in seq(k, 1)){
            gp$x$data[[i+1]]$text <- gp$x$data[[i+1+k]]$text
          }
          
        }
        else if (LableChoice == FALSE) {
          k = nclust
          gp <- ggplotly(
            fviz_dend(
              hc,
              cex = FontSize,
              k = k,
              show_labels = TRUE,
              type = "phylogenic",
              phylo_layout = PhyloLayout)
          )
          gp$x$data[[1]]$text <- NA 
          for(i in seq(k, 1)){
            gp$x$data[[i+1]]$text <- gp$x$data[[i+1+k]]$text
            gp$x$data[[i + 1 + k]] <- NULL
          }
          
        }
        
      }
      gp
    }
    
  })
  
  output$JaccDendoMat <- renderPlotly({
    
    gp <- JaccDendoMat_react()
    gp
    
  })
  
  JaccDendoMatCirc_react <- reactive({
    
    VisType <- input$ConnViewMat
    
    if (VisType == "circular") {
      
      jacc_mat <- user_matrix_loaded()
      clust_me <- input$ClustMethodMat
      nclust <- input$NumClustersMat
      LableChoice <- input$ShowConnLabelsMat
      FontSize <- input$ConnFontSizeMat
      
      hc <- hclust(dist(jacc_mat), method = clust_me)
      hc_dend <- as.dendrogram(hc)
      f <- fviz_dend(hc_dend, k = nclust, cex = FontSize, color_labels_by_k = T, type = "circular", show_labels = LableChoice) +
        theme(legend.position="none") +
        labs(title = "") +
        theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.title.x = element_blank(),
              axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank())
      
      
      f
      
    }
    
  })
  
  output$JaccDendoMatCirc <- renderPlot({
    
    f <- JaccDendoMatCirc_react()
    f
    
  })
  
  
  #####----UMAP Cluster----####
  #
  #####----Render UI----####
  #
  #output$rendUMAPsampSelect <- renderUI({
  #  
  #  if (!is.null(input$UMAPmatrixFile) & !is.null(input$UMAPmetaFile)) {
  #    expr <- umap_matrix_loaded()
  #    meta <- umap_meta_loaded()
  #    samples <- c("Select All Samples",union(colnames(expr),meta[,1]))
  #    selectizeInput(
  #      "UMAPsampSelect", 
  #      label = "Select Samples:",
  #      choices = samples, 
  #      multiple = T,
  #      selected = "",
  #      options = list(delimiter = " ", create = T)
  #    )
  #  }
  #  else if (!is.null(input$UMAPmatrixFile) & is.null(input$UMAPmetaFile)) {
  #    expr <- umap_matrix_loaded()
  #    samples <- c("Select All Samples",colnames(expr))
  #    selectizeInput(
  #      "UMAPsampSelect", 
  #      label = "Select Samples:",
  #      choices = samples, 
  #      multiple = T,
  #      selected = "",
  #      options = list(delimiter = " ", create = T)
  #    )
  #  }
  #  else if (is.null(input$UMAPmatrixFile) & !is.null(input$UMAPmetaFile)) {
  #    meta <- umap_meta_loaded()
  #    samples <- c("Select All Samples",meta[,1])
  #    selectizeInput(
  #      "UMAPsampSelect", 
  #      label = "Select Samples:",
  #      choices = samples, 
  #      multiple = T,
  #      selected = "",
  #      options = list(delimiter = " ", create = T)
  #    )
  #  }
  #})
  #
  #output$rendUMAPannotateSamps <- renderUI({
  #  
  #  req(input$UMAPmetaFile)
  #  meta <- umap_meta_loaded()
  #  anno_options <- colnames(meta)[2:ncol(meta)]
  #  anno_options <- c(" ",anno_options)
  #  selectInput("UMAPannotateSamps","Annotate Samples By:", choices = anno_options, multiple = F)
  #  
  #})
  #
  #output$rendUMAPannoContCheck <- renderUI({
  #  
  #  req(input$UMAPmetaFile)
  #  checkboxInput("UMAPannoContCheck","Continuous Variable")
  #  
  #})
  #
  #output$rendGeneSelection <- renderUI({
  #  
  #  expr <- umap_matrix_loaded()
  #  genes <- rownames(expr)
  #  #selectInput("GeneSelection","Gene:", choices = genes)
  #  selectizeInput("GeneSelection","Gene:", choices = genes)
  #  
  #  
  #})
  #
  #output$rendGeneExprRange <- renderUI({
  #  
  #  geneSelec <- input$GeneSelection
  #  expr <- umap_matrix_loaded_og()
  #  LogChoice <- input$LogGeneSelection
  #  if (LogChoice == TRUE) {
  #    expr <- log2(as.matrix(expr) + 0.00001)
  #  }
  #  exprRange <- round(range(expr[geneSelec,]),2)
  #  exprRangeText <- paste(exprRange[1],",",exprRange[2],sep = "")
  #  textInput("GeneExprRange","Expression Range:",value = exprRangeText, placeholder = "min,max")
  #  
  #})
  #
  #output$rendUMAPtypeHeader <- renderUI({
  #  
  #  if (input$MainUMAPpan == 3) {
  #    h4("Pathway Selection")
  #  }
  #  else if (input$MainUMAPpan == 4) {
  #    h4("Most Variable Genes")
  #  }
  #  
  #})
  #
  ##selec_input <- reactiveValues(seleced_opt = 0)
  ##observeEvent(input$RadioPathSelect, {
  ##  selec_input$seleced_opt <- input$RadioPathSelect
  ##})
  #output$rendRadioPathSelect <- renderUI({
  #  
  #  if (input$MainUMAPpan == 3) {
  #    #if (selec_input$seleced_opt == 0) {
  #      radioButtons("RadioPathSelect","",choices = c("Provided Genesets","Upload Genesets"), inline = T)
  #    #}
  #    #if (selec_input$seleced_opt != 0) {
  #    #  radioButtons("RadioPathSelect","",choices = c("Provided Genesets","Upload Genesets"), inline = T, selected = selec_input$seleced_opt)
  #    #}
  #  }
  #  
  #})
  #
  ##observeEvent(input$RadioPathSelect, {
  ##  
  ##  #selec_input <- as.character(input$RadioPathSelect)
  ##  
  ##  updateRadioButtons(session, "RadioPathSelect",
  ##                     label = "",
  ##                     choices = c("Provided Genesets","Upload Genesets"),
  ##                     inline = T,
  ##                     selected = input$RadioPathSelect)
  ##  })
  #
  #output$rendUserGSupload <- renderUI({
  #  
  #  if (input$MainUMAPpan == 3) {
  #    if (input$RadioPathSelect == "Upload Genesets") {
  #      fileInput("UserGSupload","Upload Genesets", accept = c(".txt",".tsv",".gmt"))
  #    }
  #  }
  #  
  #})
  #
  #output$rendUserGSheaderCheck <- renderUI({
  #  
  #  if (input$MainUMAPpan == 3) {
  #    if (input$RadioPathSelect == "Upload Genesets") {
  #      checkboxInput("UserGSheaderCheck","Header")
  #    }
  #  }
  #  
  #})
  #
  #GS_CatSelec <- reactiveValues(GS_Cat = genesetCats[1])
  #observeEvent(input$GeneSetCatSelect, {
  #  GS_CatSelec$GS_Cat <- input$GeneSetCatSelect
  #})
  #output$rendGeneSetCatSelect <- renderUI({
  #  
  #  if (input$MainUMAPpan == 3) {
  #    if (input$RadioPathSelect == "Provided Genesets"){
  #      selectInput("GeneSetCatSelect","Select Gene Set Database:",choices = genesetCats, selected = GS_CatSelec$GS_Cat)
  #    }
  #  }
  #  
  #})
  #
  #output$rendPathSelectTable <- renderUI({
  #  
  #  if (input$MainUMAPpan == 3) {
  #    div(DT::dataTableOutput("PathSelectTable"), style = "font-size:10px")
  #  }
  #  
  #})
  #
  #MVG_TopNum <- reactiveValues(MVG_num = 1000)
  #observeEvent(input$TopNumMVG, {
  #  MVG_TopNum$MVG_num <- input$TopNumMVG
  #})
  #output$rendTopNumMVG <- renderUI({
  #  
  #  if (input$MainUMAPpan == 4) {
  #    numericInput("TopNumMVG","Number of Top Genes", value = MVG_TopNum$MVG_num, step = 1)
  #  }
  #    
  #})
  #
  #MVG_VarMethod <- reactiveValues(MVG_Var = "MAD")
  #observeEvent(input$VarMethodMVG, {
  #  MVG_VarMethod$MVG_Var <- input$VarMethodMVG
  #})
  #output$rendVarMethodMVG <- renderUI({
  #  
  #  if (input$MainUMAPpan == 4) {
  #    selectInput("VarMethodMVG","Variance Method", choices = c("MAD","VAR","CV"), selected = MVG_VarMethod$MVG_Var)
  #  }
  #  
  #})
  #
  #output$rendMVGlist <- renderUI({
  #  
  #  if (input$MainUMAPpan == 4) {
  #    if (input$VeiwMVGTable == T) {
  #        div(DT::dataTableOutput("MVGlist"), style = "font-size:12px; height:350px; overflow-Y: scroll")
  #    }
  #  }
  #  
  #})
  #
  #output$renddnldMVGtab <- renderUI({
  #  
  #  if (input$MainUMAPpan == 4) {
  #    downloadButton("dnldMVGtab","Download Most Variable Genes List")
  #  }
  #  
  #})
  #
  #output$rendSelectPreCalc1 <- renderUI({
  #  
  #  meta <- umap_meta_loaded()
  #  selectInput("SelectPreCalc1","Select X-Axis Coordinate Column:",
  #              choices = colnames(meta)[2:ncol(meta)], selected = colnames(meta)[2])
  #  
  #})
  #
  #output$rendSelectPreCalc2 <- renderUI({
  #  
  #  meta <- umap_meta_loaded()
  #  selectInput("SelectPreCalc2","Select Y-Axis Coordinate Column:",
  #              choices = colnames(meta)[2:ncol(meta)], selected = colnames(meta)[3])
  #  
  #})
  #
  #output$rendVeiwMVGTable <- renderUI({
  #  
  #  if (input$MainUMAPpan == 4) {
  #    checkboxInput("VeiwMVGTable","View Most Variable Genes List", value = F)
  #  }
  #  
  #})
  #
  #
  ##output$renddnldUMAP_SVG_PreC_clin <- renderUI({
  ##  req(input$UMAPmetaFile)
  ##  downloadButton("dnldUMAP_SVG_PreC_clin","Download SVG")
  ##})
  ##output$renddnldUMAP_PDF_PreC_clin <- renderUI({
  ##  req(input$UMAPmetaFile)
  ##  downloadButton("dnldUMAP_PDF_PreC_clin","Download PDF")
  ##})
  ##
  ##output$rendUMAPplot_PreC_expr <- renderUI({
  ##  req(input$UMAPmetaFile)
  ##  withSpinner(jqui_resizable(plotOutput('UMAPplot_PreC_expr', width = "100%", height = "400px")), type = 6)
  ##})
  ##output$renddnldUMAP_SVG_PreC_expr <- renderUI({
  ##  req(input$UMAPmatrixFile)
  ##  downloadButton("dnldUMAP_SVG_PreC_expr","Download SVG")
  ##})
  ##output$renddnldUMAP_PDF_PreC_expr <- renderUI({
  ##  req(input$UMAPmatrixFile)
  ##  downloadButton("dnldUMAP_PDF_PreC_expr","Download PDF")
  ##})
  ##
  ##output$rendUMAPplot_PreC_kmean <- renderUI({
  ##  req(input$UMAPmetaFile)
  ##  withSpinner(jqui_resizable(plotOutput('UMAPplot_PreC_kmean', width = "100%", height = "400px")), type = 6)
  ##})
  ##output$renddnldUMAP_SVG_PreC_kmean <- renderUI({
  ##  req(input$UMAPmatrixFile)
  ##  downloadButton("dnldUMAP_SVG_PreC_kmean","Download SVG")
  ##})
  ##output$renddnldUMAP_PDF_PreC_kmean <- renderUI({
  ##  req(input$UMAPmatrixFile)
  ##  downloadButton("dnldUMAP_PDF_PreC_kmean","Download PDF")
  ##})
  #
  #output$rendUMAPplot_MVG_ALL <- renderUI({
  #  
  #  ph <- input$UMAPplotHeight
  #  pw <- input$UMAPplotWidth
  #  withSpinner(jqui_resizable(plotlyOutput('UMAPplot_MVG_ALL', width = pw, height = ph)), type = 6)
  #  
  #})
  #output$rendUMAPplot_PATH_ALL <- renderUI({
  #  
  #  ph <- input$UMAPplotHeight
  #  pw <- input$UMAPplotWidth
  #  withSpinner(jqui_resizable(plotlyOutput('UMAPplot_PATH_ALL', width = pw, height = ph)), type = 6)
#
  #})
  #output$rendUMAPplot_ALL <- renderUI({
  #  
  #  ph <- input$UMAPplotHeight
  #  pw <- input$UMAPplotWidth
  #  withSpinner(jqui_resizable(plotlyOutput('UMAPplot_ALL', width = pw, height = ph)), type = 6)
  #  #withSpinner(jqui_resizable(plotOutput('UMAPplot_ALL', width = pw, height = ph)), type = 6)
  #  
  #})
  #output$rendUMAPplot_PreC_ALL <- renderUI({
  #  
  #  ph <- input$UMAPplotHeight
  #  pw <- input$UMAPplotWidth
  #  withSpinner(jqui_resizable(plotlyOutput('UMAPplot_PreC_ALL', width = pw, height = ph)), type = 6)
  #  #withSpinner(jqui_resizable(plotOutput('UMAPplot_PreC_ALL', width = pw, height = ph)), type = 6)
  #  
  #})
  #
  #output$rendUMAPmvgMainTitle <- renderUI({
  #  
  #  if (!is.null(input$TopNumMVG)) {
  #    topMVGnum <- input$TopNumMVG
  #  }
  #  if (is.null(input$TopNumMVG)) {
  #    topMVGnum <- 1000
  #  }
  #  if (!is.null(input$VarMethodMVG)) {
  #    topMVGMethod <- input$VarMethodMVG
  #  }
  #  if (is.null(input$VarMethodMVG)) {
  #    topMVGMethod <- "MAD"
  #  }
  #  
  #  h3(paste("UMAP of Top",topMVGnum,"Most Variable Genes by",topMVGMethod))
  #  
  #})
  #
  #output$rendUMAPpathMainTitle <- renderUI({
  #  
  #  if (input$RadioPathSelect == "Provided Genesets") {
  #    if (input$PathSelectTable_rows_selected > 0) {
  #      rowsSelected <- input$PathSelectTable_rows_selected
  #      GeneSet <- gs_cat2[rowsSelected,3]
  #    }
  #  }
  #  else if (input$RadioPathSelect == "Upload Genesets") {
  #    req(input$UserGSupload)
  #    gs_tab <- user_gs_upload()
  #    gs_tab2 <- data.frame(GeneSet = unique(gs_tab[,1]))
  #    gs_list <- user_gs_upload_asList()
  #    if (input$PathSelectTable_rows_selected > 0) {
  #      rowsSelected <- input$PathSelectTable_rows_selected
  #      GeneSet <- gs_tab2[rowsSelected,1]
  #    }
  #  }
  #  
  #  h3(paste("UMAP of",GeneSet,"Genes"))
  #  
  #})
  #
  #####----Reactives----####
  #
  ### User Matrix Upload
  #umap_matrix_loaded_og <- reactive({
  #  
  #  gs.u <- input$UMAPmatrixFile
  #  ext <- tools::file_ext(gs.u$datapath)
  #  req(gs.u)
  #  validate(need(ext == c("tsv","txt","csv","zip"), "Please upload .tsv, .txt, or .csv file"))
  #  
  #  if (ext == "csv") {
  #    matrix.u <- as.data.frame(read_delim(gs.u$datapath, delim = ',', col_names = T))
  #  }
  #  else {
  #    matrix.u <- as.data.frame(read_delim(gs.u$datapath, delim = '\t', col_names = T))
  #  }
  #  num_test <- apply(matrix.u[,-1],2, is.numeric)
  #  if (all(num_test) == FALSE) {
  #    matrix.u[,-1] <- mutate_all(matrix.u[,-1], function(x) as.numeric(as.character(x)))
  #    
  #  }
  #  colnames(matrix.u)[1] <- "Symbol"
  #  if (TRUE %in% duplicated(matrix.u[,1])) {
  #    matrix.u <- matrix.u %>%
  #      group_by(Symbol) %>%
  #      summarise_all(max)
  #  }
  #  rownames(matrix.u) <- matrix.u[,1]
  #  matrix.u <- matrix.u[,-1]
  #  
  #  matrix.u
  #  
  #})
  #
  #umap_matrix_loaded <- reactive({
  #  
  #  matrix.u <- umap_matrix_loaded_og()
  #  matrix.u2 <- matrix.u
  #  if (input$LogUMAPmatrix == T) {
  #    matrix.u <- log2(matrix.u + 0.00001)
  #  }
  #  if (input$NormUMAPmatrix == T) {
  #    matrix.u = apply(matrix.u, 1, scale)
  #    matrix.u = apply(matrix.u, 1, rev)
  #    colnames(matrix.u) <- colnames(matrix.u2)
  #  }
  #  matrix.u <- matrix.u[sort(rownames(matrix.u)),]
  #  
  #  matrix.u
  #  
  #})
  #
  ### User Meta Upload
  #umap_meta_loaded <- reactive({
  #  
  #  gs.u <- input$UMAPmetaFile
  #  ext <- tools::file_ext(gs.u$datapath)
  #  req(gs.u)
  #  validate(need(ext == c("tsv","txt","csv"), "Please upload .tsv, .txt, or .csv file"))
  #  
  #  if (ext == "csv") {
  #    meta.u <- as.data.frame(read_delim(gs.u$datapath, delim = ',', col_names = T))
  #  }
  #  else {
  #    meta.u <- as.data.frame(read_delim(gs.u$datapath, delim = '\t', col_names = T))
  #  }
  #  colnames(meta.u)[1] <- "SampleName"
  #  meta.u
  #  
  #})
  #
  ### User Geneset Upload
  #user_gs_upload <- reactive({
  #  
  #  header_check <- input$UserGSheaderCheck
  #  gs.u <- input$UserGSupload
  #  ext <- tools::file_ext(gs.u$datapath)
  #  req(gs.u)
  #  validate(need(ext == c("tsv","txt","gmt"), "Please upload .tsv, .txt, or .gmt file"))
  #  if (ext == "gmt") {
  #    ranked_file <- read.gmt(gs.u$datapath)
  #  }
  #  else {
  #    ranked_file <- as.data.frame(read_delim(gs.u$datapath, delim = '\t', col_names = header_check, comment = "#"))
  #  }
  #  
  #  colnames(ranked_file) <- c("term","gene")
  #  ranked_file
  #  
  #})
  #
  ### Convert Uploaded GS to list
  #user_gs_upload_asList <- reactive({
  #  
  #  user_gs <- user_gs_upload()
  #  gsDataList <- list()
  #  for (i in unique(user_gs[,1])){
  #    gsDataList[[i]] <- user_gs[user_gs[,1] == i,]$gene
  #  }
  #  gsDataList
  #  
  #})
  #
  ### Get MVGs from matrix
  #umap_mvg <- reactive({
  #  
  #  expr <- umap_matrix_loaded()
  #  top_probes <- input$TopNumMVG
  #  var_type <- input$VarMethodMVG
  #  
  #  exp <- expr
  #  mad <- NULL
  #  var <- NULL
  #  cv <- NULL
  #  
  #  if (var_type == "MAD"){
  #    mad <- apply(exp, 1, mad)
  #    mad <- sort(mad, decreasing = T)
  #    mad <- head(mad, n = (top_probes +1))
  #    out <- cbind(names(mad), mad[names(mad)], exp[names(mad),])
  #    colnames(out) <- c("Gene", "MAD", colnames(exp))
  #    dataset <- exp[names(mad),]
  #    variable_gene_list <- names(mad)
  #  }
  #  else if (var_type == "VAR"){
  #    var <- apply(exp, 1, var)
  #    var <- sort(var, decreasing = T)
  #    var <- head(var, n = (top_probes +1))
  #    out <- cbind(names(var), var[names(var)], exp[names(var),])
  #    colnames(out) <- c("Gene", "VAR", colnames(exp))
  #    dataset <- exp[names(var),]
  #    variable_gene_list <- names(var)
  #  }
  #  else if (var_type == "CV"){
  #    cv <- apply(exp, 1, cv)
  #    cv <- sort(cv, decreasing = T)
  #    cv <- head(cv, n = (top_probes +1))
  #    out <- cbind(names(cv), cv[names(cv)], exp[names(cv),])
  #    colnames(out) <- c("Gene", "CV", colnames(exp))
  #    dataset <- exp[names(cv),]
  #    variable_gene_list <- names(cv)
  #  }
  #  
  #  variable_gene_list
  #  
  #})
  #
  ### Get Genes from chosen pathway
  #path_specific_umap_genes <- reactive({
  #  
  #  if (input$RadioPathSelect == "Provided Genesets") {
  #    
  #    if (input$PathSelectTable_rows_selected > 0) {
  #      
  #      rowsSelected <- input$PathSelectTable_rows_selected
  #      GeneSet <- gs_cat2[rowsSelected,3]
  #      gs_sub <- gs[GeneSet]
  #      genesSelected <- unique(unname(unlist(gs_sub)))
  #      
  #    }
  #    
  #  }
  #  else if (input$RadioPathSelect == "Upload Genesets") {
  #    
  #    req(input$UserGSupload)
  #    gs_tab <- user_gs_upload()
  #    gs_tab2 <- data.frame(GeneSet = unique(gs_tab[,1]))
  #    gs_list <- user_gs_upload_asList()
  #    if (input$PathSelectTable_rows_selected > 0) {
  #      
  #      rowsSelected <- input$PathSelectTable_rows_selected
  #      GeneSet <- gs_tab2[rowsSelected,1]
  #      gs_sub <- gs_list[GeneSet]
  #      genesSelected <- unique(unname(unlist(gs_sub)))
  #      
  #    }
  #    
  #  }
  #  
  #  genesSelected
  #  
  #})
  #
  ### Get Cluster information on matrix
  #umapClusterTable_react <- reactive({
  #  
  #  if (input$MainUMAPpan == 2 | input$MainUMAPpan == 6) {
  #    df <- umap_matrix_loaded()
  #  }
  #  else if (input$MainUMAPpan == 4) {
  #    df <- umap_matrix_loaded()
  #    genesSelected <- umap_mvg()
  #    df <- df[which(rownames(df) %in% genesSelected),]
  #  }
  #  else if (input$MainUMAPpan == 3) {
  #    df <- umap_matrix_loaded()
  #    if (length(input$PathSelectTable_rows_selected) > 0) {
  #      genesSelected <- path_specific_umap_genes()
  #      df <- df[which(rownames(df) %in% genesSelected),]
  #    }
  #  }
  #  
  #  clust_me <- input$ClusterMethod
  #  cut_k <- input$ClusterNumber
  #  
  #  results2 = hclust(dist(as.matrix(t(df))), method = clust_me)
  #  m = sort(cutree(results2, k=cut_k))
  #  output = as.data.frame(cbind(colnames(m), as.matrix(m)))
  #  output$SampleName <- rownames(output)
  #  output2 <- output[,c(2,1)]
  #  colnames(output2)[2] <- "Cluster"
  #  output2$Cluster <- as.factor(output2$Cluster)
  #  
  #  output2
  #  
  #})
  #
  ### Base Matrix
  #umap_matrix <- reactive({
  #  
  #  df <- umap_matrix_loaded()
  #  df <- as.data.frame(df)
  #  isexpr <- rowSums(df>1) >= 1
  #  isexpr_t <- names(isexpr)[which(isexpr == T)]
  #  df2 <- df[rownames(df) %in% isexpr_t,]
  #  df3 = as.matrix(df2)
  #  df4 = apply(df3, 1, rev)
  #  
  #  df4
  #  
  #})
  #
  ### UMAP Coordinate Table  - Base
  #umap_plot_table_react <- reactive({
  #  
  #  tdata <- umap_matrix()
  #  set.seed(input$UserSetSeed)
  #  
  #  tdata_fit <- tdata %>%
  #    umap()
  #  
  #  tdata_fit_df <- tdata_fit$layout %>%
  #    as.data.frame()%>%
  #    rename(UMAP1="V1",
  #           UMAP2="V2") %>%
  #    mutate(ID=row_number())
  #  
  #  tdata_fit_df$SampleName <- rownames(tdata_fit_df)
  #  tdata_fit_df <- tdata_fit_df %>%
  #    relocate(SampleName)
  #  
  #  
  #  tdata_fit_df
  #  
  #})
  #
  ### Matrix containing only genes in chosen pathway
  #umap_matrix_pathway <- reactive({
  #  
  #  df <- umap_matrix_loaded()
  #  df <- as.data.frame(df)
  #  isexpr <- rowSums(df>1) >= 1
  #  isexpr_t <- names(isexpr)[which(isexpr == T)]
  #  df <- df[rownames(df) %in% isexpr_t,]
  #  genesSelected <- path_specific_umap_genes()
  #  df <- df[which(rownames(df) %in% genesSelected),]
  #  df3 = as.matrix(df)
  #  df4 = apply(df3, 1, rev)
  #  df4
  #  
  #})
  #
  ### UMAP Coordinate Table  - Pathway
  #umap_plot_table_PATH_react <- reactive({
  #  
  #  tdata <- umap_matrix_pathway()
  #  set.seed(input$UserSetSeed)
  #  
  #  tdata_fit <- tdata %>%
  #    umap()
  #  
  #  tdata_fit_df <- tdata_fit$layout %>%
  #    as.data.frame()%>%
  #    rename(UMAP1="V1",
  #           UMAP2="V2") %>%
  #    mutate(ID=row_number())
  #  
  #  tdata_fit_df$SampleName <- rownames(tdata_fit_df)
  #  tdata_fit_df <- tdata_fit_df %>%
  #    relocate(SampleName)
  #  
  #  
  #  tdata_fit_df
  #  
  #})
  #
  ### Matrix containing only MVGs
  #umap_matrix_mvg <- reactive({
  #  
  #  df <- umap_matrix_loaded()
  #  df <- as.data.frame(df)
  #  isexpr <- rowSums(df>1) >= 1
  #  isexpr_t <- names(isexpr)[which(isexpr == T)]
  #  df <- df[rownames(df) %in% isexpr_t,]
  #  genesSelected <- umap_mvg()
  #  df <- df[which(rownames(df) %in% genesSelected),]
  #  df3 = as.matrix(df)
  #  df4 = apply(df3, 1, rev)
  #  
  #  df4
  #  
  #})
  #
  ### UMAP Coordinate Table  - MVG
  #umap_plot_table_MVG_react <- reactive({
  #  
  #  tdata <- umap_matrix_mvg()
  #  set.seed(input$UserSetSeed)
  #  
  #  tdata_fit <- tdata %>%
  #    umap()
  #  
  #  tdata_fit_df <- tdata_fit$layout %>%
  #    as.data.frame()%>%
  #    rename(UMAP1="V1",
  #           UMAP2="V2") %>%
  #    mutate(ID=row_number())
  #  
  #  tdata_fit_df$SampleName <- rownames(tdata_fit_df)
  #  tdata_fit_df <- tdata_fit_df %>%
  #    relocate(SampleName)
  #  
  #  
  #  tdata_fit_df
  #  
  #})
  #
  ### UMAP Coordinate Table  - Pre-Calculated
  #umap_plot_table_PreC_react <- reactive({
  #  
  #  #req(input$UMAPmetaFile)
  #  meta <- umap_meta_loaded()
  #  umap1 <- input$SelectPreCalc1
  #  umap2 <- input$SelectPreCalc2
  #  
  #  tdata_fit_df <- meta[,c("SampleName",umap1,umap2)]
  #  if (ncol(tdata_fit_df) >= 3) {
  #    colnames(tdata_fit_df)[c(2,3)] <- c("UMAP1","UMAP2")
  #  }
  #  
  #  tdata_fit_df
  #  
  #})
  #
  #####----Data Tables----####
  #
  #output$UMAPMatrixPreview <- DT::renderDataTable({
  #  
  #  df <- umap_matrix_loaded()
  #  
  #  DT::datatable(df,
  #                extensions = "FixedColumns",
  #                options = list(lengthMenu = c(10,20,50,100,1000,5000,10000),
  #                               pageLength = 20,
  #                               scrollX = TRUE,
  #                               autoWidth = TRUE,
  #                               fixedColumns = list(leftColumns = 1)),
  #                selection=list(mode = "multiple")) %>%
  #    formatRound(columns = c(1:ncol(df)), digits = 4)
  #  
  #})
  #
  #output$UMAPMetaPreview <- DT::renderDataTable({
  #  
  #  df <- umap_meta_loaded()
  #  
  #  DT::datatable(df,
  #                extensions = "FixedColumns",
  #                options = list(lengthMenu = c(10,20,50,100,1000,5000,10000),
  #                               pageLength = 20,
  #                               scrollX = TRUE,
  #                               autoWidth = TRUE,
  #                               fixedColumns = list(leftColumns = 1)),
  #                rownames = F,
  #                selection=list(mode = "multiple"))
  #  
  #  
  #})
  #
  #output$umapClusterTable <- DT::renderDataTable({
  #  
  #  if (input$ViewUMAPclusterTab == T) {
  #    
  #    output2 <- umapClusterTable_react()
  #    DT::datatable(output2,
  #                  options = list(keys = TRUE,
  #                                 searchHighlight = TRUE,
  #                                 pageLength = 5,
  #                                 lengthMenu = c("5","10", "25", "50", "100")),
  #                  rownames = F)
  #  }
  #  
  #})
  #
  #Path_Selec <- reactiveValues(row_selec = 0)
  #observeEvent(input$PathSelectTable_rows_selected, {
  #  Path_Selec$row_selec <- input$PathSelectTable_rows_selected
  #})
  #output$PathSelectTable <- DT::renderDataTable({
  #  
  #  if (input$RadioPathSelect == "Provided Genesets") {
  #    if (Path_Selec$row_selec == 0) {
  #      chosen_cat <- input$GeneSetCatSelect
  #      gs_cat_sub <- gs_cat2[which(gs_cat2[,1] == chosen_cat),-1]
  #      DT::datatable(gs_cat_sub,
  #                    selection = list(mode = 'single', selected = 1),
  #                    options = list(keys = TRUE,
  #                                   searchHighlight = TRUE,
  #                                   pageLength = 10,
  #                                   scrollX = TRUE,
  #                                   lengthMenu = c("10", "25", "50", "100")),
  #                    rownames = F)
  #    }
  #    else if (Path_Selec$row_selec > 0) {
  #      chosen_cat <- input$GeneSetCatSelect
  #      gs_cat_sub <- gs_cat2[which(gs_cat2[,1] == chosen_cat),-1]
  #      DT::datatable(gs_cat_sub,
  #                    selection = list(mode = 'single', selected = Path_Selec$row_selec),
  #                    options = list(keys = TRUE,
  #                                   searchHighlight = TRUE,
  #                                   pageLength = 10,
  #                                   scrollX = TRUE,
  #                                   lengthMenu = c("10", "25", "50", "100")),
  #                    rownames = F)
  #    }
  #  }
  #  else if (input$RadioPathSelect == "Upload Genesets") {
  #    req(input$UserGSupload)
  #    if (Path_Selec$row_selec == 0) {
  #      gs_tab <- user_gs_upload()
  #      GeneSets <- unique(gs_tab[,1])
  #      gs_tab2 <- data.frame(GeneSet = GeneSets)
  #      DT::datatable(gs_tab2,
  #                    selection = list(mode = 'single', selected = 1),
  #                    options = list(keys = TRUE,
  #                                   searchHighlight = TRUE,
  #                                   pageLength = 10,
  #                                   scrollX = TRUE,
  #                                   lengthMenu = c("10", "25", "50", "100")),
  #                    rownames = F)
  #    }
  #    else if (Path_Selec$row_selec > 0) {
  #      gs_tab <- user_gs_upload()
  #      GeneSets <- unique(gs_tab[,1])
  #      gs_tab2 <- data.frame(GeneSet = GeneSets)
  #      DT::datatable(gs_tab2,
  #                    selection = list(mode = 'single', selected = Path_Selec$row_selec),
  #                    options = list(keys = TRUE,
  #                                   searchHighlight = TRUE,
  #                                   pageLength = 10,
  #                                   scrollX = TRUE,
  #                                   lengthMenu = c("10", "25", "50", "100")),
  #                    rownames = F)
  #    }
  #  }
  #  
  #})
  #
  #output$MVGlist <- renderDataTable({
  #  
  #  if (input$MainUMAPpan == 4) {
  #    
  #    if (input$VeiwMVGTable == T) {
  #      
  #      gene_list <- umap_mvg()
  #      top_num <- input$TopNumMVG
  #      gene_list_df <- data.frame(Feature = gene_list)
  #      DT::datatable(gene_list_df, options = list(paging = F,scrollY = TRUE), rownames = F)
  #      
  #    }
  #    
  #  }
  #  
  #})
  #
  #UMAP_MVG_CoordTable_react <- reactive({
  #  
  #  tdata_fit_df <- umap_plot_table_MVG_react()
  #  tdata_fit_df <- tdata_fit_df[,-4]
  #  
  #  ## Add Annotation column
  #  metaColanno <- input$UMAPannotateSamps
  #  if (is.null(input$UMAPmetaFile) == F) {
  #    if (input$UMAPannotateSamps != " ") {
  #      meta <- umap_meta_loaded()
  #      metaColanno <- input$UMAPannotateSamps
  #      tdata_fit_df <- merge(tdata_fit_df,meta[,c("SampleName",metaColanno)], by = "SampleName")
  #    }
  #  }
  #  
  #  ## Add Gene Expression Column
  #  expr <- umap_matrix_loaded_og()
  #  LogChoice <- input$LogGeneSelection
  #  if (LogChoice == TRUE) {
  #    expr <- log2(as.matrix(expr) + 0.00001)
  #  }
  #  if (!is.null(input$GeneSelection)) {
  #    metaColgene <- input$GeneSelection
  #    exprRange <- round(range(expr[metaColgene,]),2)
  #    expr2 <- as.data.frame(expr)
  #    expr_g <- expr2[metaColgene,]
  #    expr_g_t <- as.data.frame(t(expr_g))
  #    expr_g_t$SampleName <- rownames(expr_g_t)
  #    tdata_fit_df <- merge(tdata_fit_df,expr_g_t, by = "SampleName")
  #  }
  #  else if (is.null(input$GeneSelection)) {
  #    metaColgene <- rownames(expr)[1]
  #    exprRange <- round(range(expr[metaColgene,]),2)
  #    expr2 <- as.data.frame(expr)
  #    expr_g <- expr2[metaColgene,]
  #    expr_g_t <- as.data.frame(t(expr_g))
  #    expr_g_t$SampleName <- rownames(expr_g_t)
  #    tdata_fit_df <- merge(tdata_fit_df,expr_g_t, by = "SampleName")
  #  }
  #  
  #  ## Add Cluster Annotation Column
  #  cluster_tab <- umapClusterTable_react()
  #  metaColkmean <- "Cluster"
  #  tdata_fit_df <- merge(tdata_fit_df,cluster_tab, by = "SampleName")
  #  
  #  tdata_fit_df
  #  
  #  
  #})
  #
  #output$UMAP_MVG_CoordTable <- DT::renderDataTable({
  #  
  #  tdata_fit_df <- UMAP_MVG_CoordTable_react()
  #  
  #  #table output,
  #  DT::datatable(tdata_fit_df,
  #                options = list(keys = TRUE,
  #                               searchHighlight = TRUE,
  #                               pageLength = 20,
  #                               scrollY = T,
  #                               lengthMenu = c("10", "20", "50", "100")
  #                ),
  #                rownames = F,
  #                selection=list(mode = "multiple"))
  #  
  #})
  #
  #UMAP_PATH_CoordTable_react <- reactive({
  #  
  #  tdata_fit_df <- umap_plot_table_PATH_react()
  #  tdata_fit_df <- tdata_fit_df[,-4]
  #  
  #  ## Add Annotation column
  #  metaColanno <- input$UMAPannotateSamps
  #  if (is.null(input$UMAPmetaFile) == F) {
  #    if (input$UMAPannotateSamps != " ") {
  #      meta <- umap_meta_loaded()
  #      metaColanno <- input$UMAPannotateSamps
  #      tdata_fit_df <- merge(tdata_fit_df,meta[,c("SampleName",metaColanno)], by = "SampleName")
  #    }
  #  }
  #  
  #  ## Add Gene Expression Column
  #  expr <- umap_matrix_loaded_og()
  #  LogChoice <- input$LogGeneSelection
  #  if (LogChoice == TRUE) {
  #    expr <- log2(as.matrix(expr) + 0.00001)
  #  }
  #  if (!is.null(input$GeneSelection)) {
  #    metaColgene <- input$GeneSelection
  #    exprRange <- round(range(expr[metaColgene,]),2)
  #    expr2 <- as.data.frame(expr)
  #    expr_g <- expr2[metaColgene,]
  #    expr_g_t <- as.data.frame(t(expr_g))
  #    expr_g_t$SampleName <- rownames(expr_g_t)
  #    tdata_fit_df <- merge(tdata_fit_df,expr_g_t, by = "SampleName")
  #  }
  #  else if (is.null(input$GeneSelection)) {
  #    metaColgene <- rownames(expr)[1]
  #    exprRange <- round(range(expr[metaColgene,]),2)
  #    expr2 <- as.data.frame(expr)
  #    expr_g <- expr2[metaColgene,]
  #    expr_g_t <- as.data.frame(t(expr_g))
  #    expr_g_t$SampleName <- rownames(expr_g_t)
  #    tdata_fit_df <- merge(tdata_fit_df,expr_g_t, by = "SampleName")
  #  }
  #  
  #  ## Add Cluster Annotation Column
  #  cluter_tab <- umapClusterTable_react()
  #  metaColkmean <- "Cluster"
  #  tdata_fit_df <- merge(tdata_fit_df,cluter_tab, by = "SampleName")
  #  
  #  tdata_fit_df
  #  
  #})
  #
  #output$UMAP_PATH_CoordTable <- DT::renderDataTable({
  #  
  #  tdata_fit_df <- UMAP_PATH_CoordTable_react()
  #  #table output,
  #  DT::datatable(tdata_fit_df,
  #                options = list(keys = TRUE,
  #                               searchHighlight = TRUE,
  #                               pageLength = 20,
  #                               scrollY = T,
  #                               lengthMenu = c("10", "20", "50", "100")
  #                ),
  #                rownames = F,
  #                selection=list(mode = "multiple"))
  #  
  #})
  #
  #UMAP_CoordTable_react <- reactive({
  #  
  #  tdata_fit_df <- umap_plot_table_react()
  #  tdata_fit_df <- tdata_fit_df[,-4]
  #  
  #  ## Add Annotation column
  #  metaColanno <- input$UMAPannotateSamps
  #  if (is.null(input$UMAPmetaFile) == F) {
  #    if (input$UMAPannotateSamps != " ") {
  #      meta <- umap_meta_loaded()
  #      metaColanno <- input$UMAPannotateSamps
  #      tdata_fit_df <- merge(tdata_fit_df,meta[,c("SampleName",metaColanno)], by = "SampleName")
  #    }
  #  }
  #  
  #  ## Add Gene Expression Column
  #  expr <- umap_matrix_loaded_og()
  #  LogChoice <- input$LogGeneSelection
  #  if (LogChoice == TRUE) {
  #    expr <- log2(as.matrix(expr) + 0.00001)
  #  }
  #  if (!is.null(input$GeneSelection)) {
  #    metaColgene <- input$GeneSelection
  #    exprRange <- round(range(expr[metaColgene,]),2)
  #    expr2 <- as.data.frame(expr)
  #    expr_g <- expr2[metaColgene,]
  #    expr_g_t <- as.data.frame(t(expr_g))
  #    expr_g_t$SampleName <- rownames(expr_g_t)
  #    tdata_fit_df <- merge(tdata_fit_df,expr_g_t, by = "SampleName")
  #  }
  #  else if (is.null(input$GeneSelection)) {
  #    metaColgene <- rownames(expr)[1]
  #    exprRange <- round(range(expr[metaColgene,]),2)
  #    expr2 <- as.data.frame(expr)
  #    expr_g <- expr2[metaColgene,]
  #    expr_g_t <- as.data.frame(t(expr_g))
  #    expr_g_t$SampleName <- rownames(expr_g_t)
  #    tdata_fit_df <- merge(tdata_fit_df,expr_g_t, by = "SampleName")
  #  }
  #  
  #  ## Add Cluster Annotation Column
  #  cluter_tab <- umapClusterTable_react()
  #  metaColkmean <- "Cluster"
  #  tdata_fit_df <- merge(tdata_fit_df,cluter_tab, by = "SampleName")
  #  
  #  tdata_fit_df
  #  
  #  
  #})
  #
  #output$UMAP_CoordTable <- DT::renderDataTable({
  #  
  #  tdata_fit_df <- UMAP_CoordTable_react()
  #  #table output,
  #  DT::datatable(tdata_fit_df,
  #                options = list(keys = TRUE,
  #                               searchHighlight = TRUE,
  #                               pageLength = 20,
  #                               scrollY = T,
  #                               lengthMenu = c("10", "20", "50", "100")
  #                ),
  #                rownames = F,
  #                selection=list(mode = "multiple"))
  #  
  #})
  #
  #UMAP_PreC_CoordTable_react <- reactive({
  #  
  #  
  #  #meta <- umap_meta_loaded()
  #  #umap1 <- input$SelectPreCalc1
  #  #umap2 <- input$SelectPreCalc2
  #  #
  #  #tdata_fit_df <- meta[,c("SampleName",umap1,umap2)]
  #  #colnames(tdata_fit_df)[c(2,3)] <- c("UMAP1","UMAP2")
  #  ##if (ncol(tdata_fit_df) >= 3) {
  #  ##  colnames(tdata_fit_df)[c(2,3)] <- c("UMAP1","UMAP2")
  #  ##}
  #  
  #  tdata_fit_df <- umap_plot_table_PreC_react()
  #  
  #  ## Add Annotation column
  #  if (is.null(input$UMAPmetaFile) == F) {
  #    if (input$UMAPannotateSamps != " ") {
  #      meta <- umap_meta_loaded()
  #      metaColanno <- input$UMAPannotateSamps
  #      tdata_fit_df <- merge(tdata_fit_df,meta[,c("SampleName",metaColanno)], by = "SampleName")
  #    }
  #  }
  #  
  #  ## Add Gene Expression Column
  #  if (!is.null(input$UMAPmatrixFile)) {
  #    expr <- umap_matrix_loaded_og()
  #    LogChoice <- input$LogGeneSelection
  #    if (LogChoice == TRUE) {
  #      expr <- log2(as.matrix(expr) + 0.00001)
  #    }
  #    if (!is.null(input$GeneSelection)) {
  #      metaColgene <- input$GeneSelection
  #      exprRange <- round(range(expr[metaColgene,]),2)
  #      expr2 <- as.data.frame(expr)
  #      expr_g <- expr2[metaColgene,]
  #      expr_g_t <- as.data.frame(t(expr_g))
  #      expr_g_t$SampleName <- rownames(expr_g_t)
  #      tdata_fit_df <- merge(tdata_fit_df,expr_g_t, by = "SampleName")
  #    }
  #    else if (is.null(input$GeneSelection)) {
  #      metaColgene <- rownames(expr)[1]
  #      exprRange <- round(range(expr[metaColgene,]),2)
  #      expr2 <- as.data.frame(expr)
  #      expr_g <- expr2[metaColgene,]
  #      expr_g_t <- as.data.frame(t(expr_g))
  #      expr_g_t$SampleName <- rownames(expr_g_t)
  #      tdata_fit_df <- merge(tdata_fit_df,expr_g_t, by = "SampleName")
  #    }
  #    
  #    ## Add Cluster Annotation Column
  #    cluter_tab <- umapClusterTable_react()
  #    metaColkmean <- "Cluster"
  #    tdata_fit_df <- merge(tdata_fit_df,cluter_tab, by = "SampleName")
  #  }
  #  
  #  tdata_fit_df
  #  
  #  
  #})
  #
  #output$UMAP_PreC_CoordTable <- DT::renderDataTable({
  #  
  #  tdata_fit_df <- UMAP_PreC_CoordTable_react()
  #  #table output,
  #  DT::datatable(tdata_fit_df,
  #                options = list(keys = TRUE,
  #                               searchHighlight = TRUE,
  #                               pageLength = 20,
  #                               scrollY = T,
  #                               lengthMenu = c("10", "20", "50", "100")
  #                ),
  #                rownames = F,
  #                selection=list(mode = "multiple"))
  #  
  #})
  #
  #####----Plots----####
  #
  #####----MVG UMAP----####
  #
  ### UMAP Plot  - MVG - Clin
  #umap_plot_MVG_react_clin_base <- reactive({
  #  
  #  ## Variables
  #  sampSelected <- input$UMAPsampSelect         # Sample Names to annotate
  #  UMAPtitleText <- input$UMAPtitleTextSize     # Title text size
  #  UMAPaxisText <- input$UMAPaxisTextSize       # Axis text size
  #  UMAPlegendText <- input$UMAPlegendTextSize   # Legend text size
  #  UMAPdotSize <- input$UMAPdotSize             # Dot size
  #  metaColanno <- NULL
  #  
  #  plot_df <- UMAP_MVG_CoordTable_react()
  #  rownames(plot_df) <- plot_df[,1]
  #  if (is.null(input$UMAPmetaFile) == F) {
  #    if (input$UMAPannotateSamps != " ") {
  #      metaColanno <- input$UMAPannotateSamps
  #      if (input$UMAPannoContCheck != T) {
  #        plot_df[,metaColanno] <- as.factor(plot_df[,metaColanno])
  #      }
  #    }
  #  }
  #  expr <- umap_matrix_loaded_og()
  #  if (!is.null(input$GeneSelection)) {
  #    metaColgene <- input$GeneSelection
  #  }
  #  if (is.null(input$GeneSelection)) {
  #    metaColgene <- rownames(expr)[1]
  #  }
  #  metaColkmean <- "Cluster"
  #  
  #  ## Meta Annotation Plot
  #  if (!is.null(metaColanno)) {
  #    colnames(plot_df)[c(4,5)] <- c("AnnoName","GeneName")
  #    k <- plot_df %>%
  #      ggplot(aes(UMAP1, UMAP2, colour=AnnoName,
  #                 text = paste("</br> <b>Sample Name:</b> ", SampleName,
  #                              "</br> <b>",metaColanno,":</b> ", AnnoName,
  #                              "</br> <b>",metaColgene," Gene Expression:</b> ", round(GeneName,4),
  #                              "</br> <b>Cluster:</b> ", Cluster,
  #                              sep = "")))
  #  }
  #  if (is.null(metaColanno)) {
  #    colnames(plot_df)[4] <- "GeneName"
  #    k <- plot_df %>%
  #      ggplot(aes(UMAP1, UMAP2,
  #                 text = paste("</br> <b>Sample Name:</b> ", SampleName,
  #                              "</br> <b>",metaColgene," Gene Expression:</b> ", round(GeneName,4),
  #                              "</br> <b>Cluster:</b> ", Cluster,
  #                              sep = "")))
  #  }
  #  
  #  k <- k + geom_point(shape = 19,
  #                      size = UMAPdotSize) +
  #    
  #    theme_minimal()
  #  
  #  if (is.null(input$UMAPmetaFile) == F) {
  #    if (input$UMAPannotateSamps != " ") {
  #      metaColanno <- input$UMAPannotateSamps
  #      k <- k + labs(x = "UMAP1",
  #                    y = "UMAP2",
  #                    color = metaColanno)
  #    }
  #    if (input$UMAPannotateSamps == " ") {
  #      k <- k + labs(x = "UMAP1",
  #                    y = "UMAP2")
  #    }
  #  }
  #  if (is.null(input$UMAPmetaFile) == T) {
  #    k <- k + labs(x = "UMAP1",
  #                  y = "UMAP2")
  #  }
  #  
  #  if (is.null(input$UMAPmetaFile) == F) {
  #    if (input$UMAPannotateSamps != " ") {
  #      if (input$UMAPannoContCheck == T) {
  #        myPalette <- colorRampPalette(rev(brewer.pal(9, input$UMAPcolors)))
  #        k <- k + scale_colour_gradientn(colours = rev(myPalette(100)))
  #      }
  #    }
  #  }
  #  
  #  k <- k + theme(axis.text = element_text(size = UMAPaxisText),
  #                 axis.title = element_text(size = UMAPaxisText),
  #                 plot.title = element_text(size = UMAPtitleText),
  #                 legend.text=element_text(size=UMAPlegendText))
  #  
  #  if (length(sampSelected) > 0) {
  #    plotdata <- plot_df[sampSelected,]
  #    plotlabel <- rownames(plot_df[sampSelected,])
  #    if (!is.null(metaColanno)) {
  #      colnames(plotdata)[c(4,5)] <- c("AnnoName","GeneName")
  #    }
  #    if (is.null(metaColanno)) {
  #      colnames(plotdata)[4] <- "GeneName"
  #    }
  #    k <- k + geom_point(data = plotdata,
  #                        aes(x = UMAP1, y = UMAP2),
  #                        pch = 1,
  #                        color = "black",
  #                        size = UMAPdotSize,
  #                        stroke = .3)
  #  }
  #  k
  #  
  #})
  #
  #umap_plot_MVG_react_clin <- reactive({
  #  
  #  ## Variables
  #  sampSelected <- input$UMAPsampSelect         # Sample Names to annotate
  #  metaColanno <- NULL
#
  #  plot_df <- UMAP_MVG_CoordTable_react()
  #  rownames(plot_df) <- plot_df[,1]
  #  if (is.null(input$UMAPmetaFile) == F) {
  #      if (input$UMAPannotateSamps != " ") {
  #        metaColanno <- input$UMAPannotateSamps
  #        if (input$UMAPannoContCheck != T) {
  #          plot_df[,metaColanno] <- as.factor(plot_df[,metaColanno])
  #        }
  #      }
  #  }
  #  expr <- umap_matrix_loaded_og()
  #  if (!is.null(input$GeneSelection)) {
  #    metaColgene <- input$GeneSelection
  #  }
  #  if (is.null(input$GeneSelection)) {
  #    metaColgene <- rownames(expr)[1]
  #  }
  #  metaColkmean <- "Cluster"
  #  
  #  k <- umap_plot_MVG_react_clin_base()
  #  
  #  if (length(sampSelected) > 0) {
  #      plotdata <- plot_df[sampSelected,]
  #      plotlabel <- rownames(plot_df[sampSelected,])
  #      if (!is.null(metaColanno)) {
  #        colnames(plotdata)[c(4,5)] <- c("AnnoName","GeneName")
  #      }
  #      if (is.null(metaColanno)) {
  #        colnames(plotdata)[4] <- "GeneName"
  #      }
  #  }
  #  
  #  ## Make it plotly
  #  k2 <- ggplotly(k,
  #                 tooltip = "text") %>% 
  #    
  #    config(displayModeBar = F)  %>% 
  #    
  #    layout(font=list(color="#black"),
  #           xaxis=list(title="UMAP1",zeroline=F),
  #           yaxis=list(title="UMAP2",zeroline=F)) 
  #    
  #  k2 <- k2 %>%
  #    hide_legend() %>%
  #    hide_colorbar()
  #    
  #    if (length(sampSelected) > 0) {
  #      k2 <- k2 %>%
  #        add_annotations(x = plotdata$UMAP1,
  #                        y = plotdata$UMAP2,
  #                        text = rownames(plotdata),
  #                        showarrow = TRUE,
  #                        arrowhead = 4,
  #                        arrowsize = .5)
  #    }
  #    
  #  k2
  #  
  #})
  #
  #umap_plot_MVG_react_expr_base <- reactive({
  #  
  #  ## Variables
  #  sampSelected <- input$UMAPsampSelect         # Sample Names to annotate
  #  UMAPtitleText <- input$UMAPtitleTextSize     # Title text size
  #  UMAPaxisText <- input$UMAPaxisTextSize       # Axis text size
  #  UMAPlegendText <- input$UMAPlegendTextSize   # Legend text size
  #  UMAPdotSize <- input$UMAPdotSize             # Dot size
  #  metaColanno <- NULL
  #  
  #  plot_df <- UMAP_MVG_CoordTable_react()
  #  rownames(plot_df) <- plot_df[,1]
  #  if (is.null(input$UMAPmetaFile) == F) {
  #    if (input$UMAPannotateSamps != " ") {
  #      metaColanno <- input$UMAPannotateSamps
  #      if (input$UMAPannoContCheck != T) {
  #        plot_df[,metaColanno] <- as.factor(plot_df[,metaColanno])
  #      }
  #    }
  #  }
  #  expr <- umap_matrix_loaded_og()
  #  if (!is.null(input$GeneSelection)) {
  #    metaColgene <- input$GeneSelection
  #  }
  #  if (is.null(input$GeneSelection)) {
  #    metaColgene <- rownames(expr)[1]
  #  }
  #  LogChoice <- input$LogGeneSelection
  #  expr2 <- expr
  #  if (LogChoice == TRUE) {
  #    expr2 <- log2(as.matrix(expr) + 0.00001)
  #  }
  #  exprRange <- round(range(expr2[metaColgene,]),2)
  #  metaColkmean <- "Cluster"
  #  
  #  if (!is.null(input$GeneExprRange)) {
  #    exprMin <- as.numeric(gsub(" ","",strsplit(as.character(input$GeneExprRange),",")[[1]][1]))
  #    exprMax <- as.numeric(gsub(" ","",strsplit(as.character(input$GeneExprRange),",")[[1]][2]))
  #  }
  #  else if (is.null(input$GeneExprRange)) {
  #    exprMin <- as.numeric(gsub(" ","",strsplit(as.character(exprRange),",")[[1]][1]))
  #    exprMax <- as.numeric(gsub(" ","",strsplit(as.character(exprRange),",")[[1]][2]))
  #  }
  #  myPalette <- colorRampPalette(rev(brewer.pal(9, input$UMAPcolors)))
  #  
  #  ## Meta Annotation Plot
  #  if (!is.null(metaColanno)) {
  #    colnames(plot_df)[c(4,5)] <- c("AnnoName","GeneName")
  #    m <- plot_df %>%
  #      ggplot(aes(UMAP1, UMAP2, colour=GeneName,
  #                 text = paste("</br> <b>Sample Name:</b> ", SampleName,
  #                              "</br> <b>",metaColanno,":</b> ", AnnoName,
  #                              "</br> <b>",metaColgene," Gene Expression:</b> ", round(GeneName,4),
  #                              "</br> <b>Cluster:</b> ", Cluster,
  #                              sep = "")))
  #  }
  #  if (is.null(metaColanno)) {
  #    colnames(plot_df)[4] <- "GeneName"
  #    m <- plot_df %>%
  #      ggplot(aes(UMAP1, UMAP2, colour=GeneName,
  #                 text = paste("</br> <b>Sample Name:</b> ", SampleName,
  #                              "</br> <b>",metaColgene," Gene Expression:</b> ", round(GeneName,4),
  #                              "</br> <b>Cluster:</b> ", Cluster,
  #                              sep = "")))
  #  }
  #  m <- m + 
  #    geom_point(shape = 19,
  #               size = UMAPdotSize) +
  #    
  #    theme_minimal()
  #  
  #  m <- m + labs(x = "UMAP1",
  #                y = "UMAP2",
  #                color = metaColgene)
  #  
  #  m <-  m + scale_colour_gradientn(colours = rev(myPalette(100)), limits=c(exprMin, exprMax), oob = scales::squish)
  #  
  #  m <- m + theme(axis.text = element_text(size = UMAPaxisText),
  #                 axis.title = element_text(size = UMAPaxisText),
  #                 plot.title = element_text(size = UMAPtitleText),
  #                 legend.text=element_text(size=UMAPlegendText))
  #  
  #  if (length(sampSelected) > 0) {
  #    plotdata <- plot_df[sampSelected,]
  #    plotlabel <- rownames(plot_df[sampSelected,])
  #    if (!is.null(metaColanno)) {
  #      colnames(plotdata)[c(4,5)] <- c("AnnoName","GeneName")
  #    }
  #    if (is.null(metaColanno)) {
  #      colnames(plotdata)[4] <- "GeneName"
  #    }
  #    m <- m + geom_point(data = plotdata,
  #                        aes(x = UMAP1, y = UMAP2),
  #                        pch = 1,
  #                        color = "black",
  #                        size = UMAPdotSize,
  #                        stroke = .3)
  #  }
  #  m
  #  
  #})
  #
  ### UMAP Plot  - Base - expr
  #umap_plot_MVG_react_expr <- reactive({
  #  
  #  
  #  ## Variables}## Variables
  #  sampSelected <- input$UMAPsampSelect         # Sample Names to annotate
  #  metaColanno <- NULL
  #  
  #  plot_df <- UMAP_MVG_CoordTable_react()
  #  rownames(plot_df) <- plot_df[,1]
  #  if (is.null(input$UMAPmetaFile) == F) {
  #    if (input$UMAPannotateSamps != " ") {
  #      metaColanno <- input$UMAPannotateSamps
  #      if (input$UMAPannoContCheck != T) {
  #        plot_df[,metaColanno] <- as.factor(plot_df[,metaColanno])
  #      }
  #    }
  #  }
  #  expr <- umap_matrix_loaded_og()
  #  if (!is.null(input$GeneSelection)) {
  #    metaColgene <- input$GeneSelection
  #  }
  #  if (is.null(input$GeneSelection)) {
  #    metaColgene <- rownames(expr)[1]
  #  }
  #  LogChoice <- input$LogGeneSelection
  #  expr2 <- expr
  #  if (LogChoice == TRUE) {
  #    expr2 <- log2(as.matrix(expr) + 0.00001)
  #  }
  #  exprRange <- round(range(expr2[metaColgene,]),2)
  #  metaColkmean <- "Cluster"
  #  
  #  if (!is.null(input$GeneExprRange)) {
  #      exprMin <- as.numeric(gsub(" ","",strsplit(as.character(input$GeneExprRange),",")[[1]][1]))
  #      exprMax <- as.numeric(gsub(" ","",strsplit(as.character(input$GeneExprRange),",")[[1]][2]))
  #  }
  #  else if (is.null(input$GeneExprRange)) {
  #    exprMin <- as.numeric(gsub(" ","",strsplit(as.character(exprRange),",")[[1]][1]))
  #    exprMax <- as.numeric(gsub(" ","",strsplit(as.character(exprRange),",")[[1]][2]))
  #  }
  #  myPalette <- colorRampPalette(rev(brewer.pal(9, input$UMAPcolors)))
  #  
  #  m <- umap_plot_MVG_react_expr_base()
  #  
  #  if (length(sampSelected) > 0) {
  #    plotdata <- plot_df[sampSelected,]
  #    plotlabel <- rownames(plot_df[sampSelected,])
  #    if (!is.null(metaColanno)) {
  #      colnames(plotdata)[c(4,5)] <- c("AnnoName","GeneName")
  #    }
  #    if (is.null(metaColanno)) {
  #      colnames(plotdata)[4] <- "GeneName"
  #    }
  #  }
  #  
  #  ## Make it plotly
  #  m2 <- ggplotly(m,
  #                 tooltip = "text") %>% 
  #    
  #    config(displayModeBar = F)  %>% 
  #    
  #    layout(font=list(color="#black"),
  #           xaxis=list(title="UMAP1",zeroline=F),
  #           yaxis=list(title="UMAP2",zeroline=F)) %>% 
  #    
  #    hide_legend() %>%
  #    hide_colorbar()
  #  
  #  if (length(sampSelected) > 0) {
  #    m2 <- m2 %>%
  #      add_annotations(x = plotdata$UMAP1,
  #                      y = plotdata$UMAP2,
  #                      text = rownames(plotdata),
  #                      showarrow = TRUE,
  #                      arrowhead = 4,
  #                      arrowsize = .5)
  #  }
  #    
  #    
  #  m2
  #  
  #})
  #
  #umap_plot_MVG_react_kmean_base <- reactive({
  #  
  #  ## Variables
  #  sampSelected <- input$UMAPsampSelect         # Sample Names to annotate
  #  UMAPtitleText <- input$UMAPtitleTextSize     # Title text size
  #  UMAPaxisText <- input$UMAPaxisTextSize       # Axis text size
  #  UMAPlegendText <- input$UMAPlegendTextSize   # Legend text size
  #  UMAPdotSize <- input$UMAPdotSize             # Dot size
  #  metaColanno <- NULL
  #  
  #  plot_df <- UMAP_MVG_CoordTable_react()
  #  rownames(plot_df) <- plot_df[,1]
  #  if (is.null(input$UMAPmetaFile) == F) {
  #    if (input$UMAPannotateSamps != " ") {
  #      metaColanno <- input$UMAPannotateSamps
  #      if (input$UMAPannoContCheck != T) {
  #        plot_df[,metaColanno] <- as.factor(plot_df[,metaColanno])
  #      }
  #    }
  #  }
  #  expr <- umap_matrix_loaded_og()
  #  if (!is.null(input$GeneSelection)) {
  #    metaColgene <- input$GeneSelection
  #  }
  #  if (is.null(input$GeneSelection)) {
  #    metaColgene <- rownames(expr)[1]
  #  }
  #  metaColkmean <- "Cluster"
  #  
  #  ## Meta Annotation Plot
  #  if (!is.null(metaColanno)) {
  #    colnames(plot_df)[c(4,5)] <- c("AnnoName","GeneName")
  #    n <- plot_df %>%
  #      ggplot(aes(UMAP1, UMAP2, colour=Cluster,
  #                 text = paste("</br> <b>Sample Name:</b> ", SampleName,
  #                              "</br> <b>",metaColanno,":</b> ", AnnoName,
  #                              "</br> <b>",metaColgene," Gene Expression:</b> ", round(GeneName,4),
  #                              "</br> <b>Cluster:</b> ", Cluster,
  #                              sep = "")))
  #  }
  #  if (is.null(metaColanno)) {
  #    colnames(plot_df)[4] <- "GeneName"
  #    n <- plot_df %>%
  #      ggplot(aes(UMAP1, UMAP2, colour=Cluster,
  #                 text = paste("</br> <b>Sample Name:</b> ", SampleName,
  #                              "</br> <b>",metaColgene," Gene Expression:</b> ", round(GeneName,4),
  #                              "</br> <b>Cluster:</b> ", Cluster,
  #                              sep = "")))
  #  }
  #  n <- n + geom_point(shape = 19,
  #                      size = UMAPdotSize)+
  #    
  #    labs(x = "UMAP1",
  #         y = "UMAP2",
  #         color = "Cluster")  +
  #    
  #    theme_minimal()
  #  
  #  n <- n + theme(axis.text = element_text(size = UMAPaxisText),
  #                 axis.title = element_text(size = UMAPaxisText),
  #                 plot.title = element_text(size = UMAPtitleText),
  #                 legend.text=element_text(size=UMAPlegendText))
  #  
  #  if (length(sampSelected) > 0) {
  #    plotdata <- plot_df[sampSelected,]
  #    plotlabel <- rownames(plot_df[sampSelected,])
  #    if (!is.null(metaColanno)) {
  #      colnames(plotdata)[c(4,5)] <- c("AnnoName","GeneName")
  #    }
  #    if (is.null(metaColanno)) {
  #      colnames(plotdata)[4] <- "GeneName"
  #    }
  #    n <- n + geom_point(data = plotdata,
  #                        aes(x = UMAP1, y = UMAP2),
  #                        pch = 1,
  #                        color = "black",
  #                        size = UMAPdotSize,
  #                        stroke = .3)
  #  }
  #  n
  #  
  #})
  #
  ### UMAP Plot  - Base - kmean
  #umap_plot_MVG_react_kmean <- reactive({
  #  
  #  ## Variables
  #  sampSelected <- input$UMAPsampSelect         # Sample Names to annotate
  #  metaColanno <- NULL
  #  
  #  plot_df <- UMAP_MVG_CoordTable_react()
  #  rownames(plot_df) <- plot_df[,1]
  #  if (is.null(input$UMAPmetaFile) == F) {
  #    if (input$UMAPannotateSamps != " ") {
  #      metaColanno <- input$UMAPannotateSamps
  #      if (input$UMAPannoContCheck != T) {
  #        plot_df[,metaColanno] <- as.factor(plot_df[,metaColanno])
  #      }
  #    }
  #  }
  #  expr <- umap_matrix_loaded_og()
  #  if (!is.null(input$GeneSelection)) {
  #    metaColgene <- input$GeneSelection
  #  }
  #  if (is.null(input$GeneSelection)) {
  #    metaColgene <- rownames(expr)[1]
  #  }
  #  metaColkmean <- "Cluster"
  #  
  #  n <- umap_plot_MVG_react_kmean_base()
  #  
  #  if (length(sampSelected) > 0) {
  #    plotdata <- plot_df[sampSelected,]
  #    plotlabel <- rownames(plot_df[sampSelected,])
  #    if (!is.null(metaColanno)) {
  #      colnames(plotdata)[c(4,5)] <- c("AnnoName","GeneName")
  #    }
  #    if (is.null(metaColanno)) {
  #      colnames(plotdata)[4] <- "GeneName"
  #    }
  #  }
  #  
  #  ## Make it plotly
  #  n2 <- ggplotly(n,
  #                 tooltip = "text") %>% 
  #    
  #    config(displayModeBar = F)  %>% 
  #    
  #    layout(font=list(color="#black"),
  #           xaxis=list(title="UMAP1",zeroline=F),
  #           yaxis=list(title="UMAP2",zeroline=F)) %>% 
  #    
  #    hide_legend() %>%
  #    hide_colorbar()
  #  if (length(sampSelected) > 0) {
  #    n2 <- n2 %>%
  #      add_annotations(x = plotdata$UMAP1,
  #                      y = plotdata$UMAP2,
  #                      text = rownames(plotdata),
  #                      showarrow = TRUE,
  #                      arrowhead = 4,
  #                      arrowsize = .5)
  #  }
  #    
  #    
  #  n2
  #  
  #})
  #
  #output$UMAPmvgLegend <- renderPlot({
  #  
  #  k <- umap_plot_MVG_react_clin_base()
  #  m <- umap_plot_MVG_react_expr_base()
  #  n <- umap_plot_MVG_react_kmean_base()
  #  
  #  if (is.null(input$UMAPmetaFile) == F) {
  #    if (input$UMAPannotateSamps != " ") {
  #      legendk <- g_legend(k)
  #      legendm <- g_legend(m)
  #      legendn <- g_legend(n)
  #      legend_grid <- grid.arrange(legendk,legendm,legendn,nrow = 1)
  #    }
  #    else if (input$UMAPannotateSamps == " ") {
  #      legendm <- g_legend(m)
  #      legendn <- g_legend(n)
  #      legend_grid <- grid.arrange(legendm,legendn,nrow = 1)
  #    }
  #  }
  #  if (is.null(input$UMAPmetaFile) == T) {
  #      legendm <- g_legend(m)
  #      legendn <- g_legend(n)
  #      legend_grid <- grid.arrange(legendm,legendn,nrow = 1)
  #  }
  #  
  #  legend_grid
  #  
  #})
  #
  #UMAPplot_MVG_ALL_react <- reactive({
  #  
  #  k2 <- umap_plot_MVG_react_clin()
  #  m2 <- umap_plot_MVG_react_expr()
  #  n2 <- umap_plot_MVG_react_kmean()
  #  
  #  if (is.null(input$UMAPmetaFile) == F) {
  #    if (input$UMAPannotateSamps != " ") {
  #      UMAPtitlek <- paste("Annotated by",input$UMAPannotateSamps)
  #    }
  #    else if (input$UMAPannotateSamps == " ") {
  #      UMAPtitlek <- ""
  #    }
  #  }
  #  if (is.null(input$UMAPmetaFile) == T) {
  #    UMAPtitlek <- ""
  #  }
  #  expr <- umap_matrix_loaded_og()
  #  if (!is.null(input$GeneSelection)) {
  #    GeneSelec <- input$GeneSelection
  #  }
  #  if (is.null(input$GeneSelection)) {
  #    GeneSelec <- rownames(expr)[1]
  #  }
  #  
  #  if (input$LogGeneSelection == T) {
  #    GeneSelecExpr <- paste(GeneSelec,"Expression (Log2)")
  #  }
  #  if (input$LogGeneSelection == F) {
  #    GeneSelecExpr <- paste(GeneSelec,"Expression")
  #  }
  #  
  #  UMAPtitlem <- paste("Annotated by",GeneSelecExpr)
  #  UMAPtitlen <- paste("Annotated by",input$ClusterMethod,input$ClusterNumber,"Clusters")
  #  
  #  
  #  if (input$UMAPorientation == "Side-by-Side") {
  #    subplot_all_sbs <- subplot(k2, m2, n2, nrows = 1, shareY = TRUE, shareX = TRUE)
  #    plot_titles = list( 
  #      list( 
  #        x = 0,
  #        y = 1,
  #        text = UMAPtitlek,  
  #        xref = "x",  
  #        yref = "paper",  
  #        xanchor = "center",  
  #        yanchor = "bottom",  
  #        showarrow = FALSE 
  #      ),  
  #      list( 
  #        x = 0,
  #        y = 1,
  #        text = UMAPtitlem,  
  #        xref = "x2",
  #        yref = "paper",
  #        xanchor = "center",  
  #        yanchor = "bottom",    
  #        showarrow = FALSE 
  #      ),  
  #      list( 
  #        x = 0,
  #        y = 1,
  #        text = UMAPtitlen,  
  #        xref = "x3",
  #        yref = "paper",
  #        xanchor = "center",  
  #        yanchor = "bottom",  
  #        showarrow = FALSE 
  #      ))
  #  }
  #  else if (input$UMAPorientation == "Stacked") {
  #    subplot_all_sbs <- subplot(k2, m2, n2, nrows = 3, shareY = TRUE, shareX = TRUE, margin = 0.04)
  #    plot_titles = list( 
  #      list( 
  #        x = 0,
  #        y = 1,
  #        text = UMAPtitlek,  
  #        xref = "x",  
  #        yref = "paper",  
  #        xanchor = "center",  
  #        yanchor = "bottom",  
  #        showarrow = FALSE 
  #      ),  
  #      list( 
  #        x = 0,
  #        y = 0.65,
  #        text = UMAPtitlem,  
  #        xref = "x",
  #        yref = "paper",
  #        xanchor = "center",  
  #        yanchor = "bottom",    
  #        showarrow = FALSE 
  #      ),  
  #      list( 
  #        x = 0,
  #        y = 0.3,
  #        text = UMAPtitlen,  
  #        xref = "x",
  #        yref = "paper",
  #        xanchor = "center",  
  #        yanchor = "bottom",  
  #        showarrow = FALSE 
  #      ))
  #  }
  #  
  #  subplot_all_sbs <- subplot_all_sbs %>% layout(annotations = plot_titles)
  #  
  #  subplot_all_sbs
  #})
  #
  #output$UMAPplot_MVG_ALL <- renderPlotly({
  #  p_all <- UMAPplot_MVG_ALL_react()
  #  p_all
  #})
  #
  #####----PATH UMAP----####
  #
  #
  ### UMAP Plot  - PATH - Clin
  #umap_plot_PATH_react_clin_base <- reactive({
  #  
  #  ## Variables
  #  sampSelected <- input$UMAPsampSelect         # Sample Names to annotate
  #  UMAPtitleText <- input$UMAPtitleTextSize     # Title text size
  #  UMAPaxisText <- input$UMAPaxisTextSize       # Axis text size
  #  UMAPlegendText <- input$UMAPlegendTextSize   # Legend text size
  #  UMAPdotSize <- input$UMAPdotSize             # Dot size
  #  metaColanno <- NULL
  #  
  #  plot_df <- UMAP_PATH_CoordTable_react()
  #  rownames(plot_df) <- plot_df[,1]
  #  if (is.null(input$UMAPmetaFile) == F) {
  #    if (input$UMAPannotateSamps != " ") {
  #      metaColanno <- input$UMAPannotateSamps
  #      if (input$UMAPannoContCheck != T) {
  #        plot_df[,metaColanno] <- as.factor(plot_df[,metaColanno])
  #      }
  #    }
  #  }
  #  expr <- umap_matrix_loaded_og()
  #  if (!is.null(input$GeneSelection)) {
  #    metaColgene <- input$GeneSelection
  #  }
  #  if (is.null(input$GeneSelection)) {
  #    metaColgene <- rownames(expr)[1]
  #  }
  #  metaColkmean <- "Cluster"
  #  
  #  ## Meta Annotation Plot
  #  if (!is.null(metaColanno)) {
  #    colnames(plot_df)[c(4,5)] <- c("AnnoName","GeneName")
  #    k <- plot_df %>%
  #      ggplot(aes(UMAP1, UMAP2, colour=AnnoName,
  #                 text = paste("</br> <b>Sample Name:</b> ", SampleName,
  #                              "</br> <b>",metaColanno,":</b> ", AnnoName,
  #                              "</br> <b>",metaColgene," Gene Expression:</b> ", round(GeneName,4),
  #                              "</br> <b>Cluster:</b> ", Cluster,
  #                              sep = "")))
  #  }
  #  if (is.null(metaColanno)) {
  #    colnames(plot_df)[4] <- "GeneName"
  #    k <- plot_df %>%
  #      ggplot(aes(UMAP1, UMAP2,
  #                 text = paste("</br> <b>Sample Name:</b> ", SampleName,
  #                              "</br> <b>",metaColgene," Gene Expression:</b> ", round(GeneName,4),
  #                              "</br> <b>Cluster:</b> ", Cluster,
  #                              sep = "")))
  #  }
  #  
  #  k <- k + geom_point(shape = 19,
  #                      size = UMAPdotSize) +
  #    
  #    theme_minimal()
  #  
  #  if (is.null(input$UMAPmetaFile) == F) {
  #    if (input$UMAPannotateSamps != " ") {
  #      metaColanno <- input$UMAPannotateSamps
  #      k <- k + labs(x = "UMAP1",
  #                    y = "UMAP2",
  #                    color = metaColanno)
  #    }
  #    if (input$UMAPannotateSamps == " ") {
  #      k <- k + labs(x = "UMAP1",
  #                    y = "UMAP2")
  #    }
  #  }
  #  if (is.null(input$UMAPmetaFile) == T) {
  #    k <- k + labs(x = "UMAP1",
  #                  y = "UMAP2")
  #  }
  #  
  #  if (is.null(input$UMAPmetaFile) == F) {
  #    if (input$UMAPannotateSamps != " ") {
  #      if (input$UMAPannoContCheck == T) {
  #        myPalette <- colorRampPalette(rev(brewer.pal(9, input$UMAPcolors)))
  #        k <- k + scale_colour_gradientn(colours = rev(myPalette(100)))
  #      }
  #    }
  #  }
  #  
  #  k <- k + theme(axis.text = element_text(size = UMAPaxisText),
  #                 axis.title = element_text(size = UMAPaxisText),
  #                 plot.title = element_text(size = UMAPtitleText),
  #                 legend.text=element_text(size=UMAPlegendText))
  #  
  #  if (length(sampSelected) > 0) {
  #    plotdata <- plot_df[sampSelected,]
  #    plotlabel <- rownames(plot_df[sampSelected,])
  #    if (!is.null(metaColanno)) {
  #      colnames(plotdata)[c(4,5)] <- c("AnnoName","GeneName")
  #    }
  #    if (is.null(metaColanno)) {
  #      colnames(plotdata)[4] <- "GeneName"
  #    }
  #    k <- k + geom_point(data = plotdata,
  #                        aes(x = UMAP1, y = UMAP2),
  #                        pch = 1,
  #                        color = "black",
  #                        size = UMAPdotSize,
  #                        stroke = .3)
  #  }
  #  k
  #  
  #})
  #
  #umap_plot_PATH_react_clin <- reactive({
  #  
  #  ## Variables
  #  sampSelected <- input$UMAPsampSelect         # Sample Names to annotate
  #  metaColanno <- NULL
  #  
  #  plot_df <- UMAP_PATH_CoordTable_react()
  #  rownames(plot_df) <- plot_df[,1]
  #  if (is.null(input$UMAPmetaFile) == F) {
  #    if (input$UMAPannotateSamps != " ") {
  #      metaColanno <- input$UMAPannotateSamps
  #      if (input$UMAPannoContCheck != T) {
  #        plot_df[,metaColanno] <- as.factor(plot_df[,metaColanno])
  #      }
  #    }
  #  }
  #  expr <- umap_matrix_loaded_og()
  #  if (!is.null(input$GeneSelection)) {
  #    metaColgene <- input$GeneSelection
  #  }
  #  if (is.null(input$GeneSelection)) {
  #    metaColgene <- rownames(expr)[1]
  #  }
  #  metaColkmean <- "Cluster"
  #  
  #  k <- umap_plot_PATH_react_clin_base()
  #  
  #  if (length(sampSelected) > 0) {
  #    plotdata <- plot_df[sampSelected,]
  #    plotlabel <- rownames(plot_df[sampSelected,])
  #    if (!is.null(metaColanno)) {
  #      colnames(plotdata)[c(4,5)] <- c("AnnoName","GeneName")
  #    }
  #    if (is.null(metaColanno)) {
  #      colnames(plotdata)[4] <- "GeneName"
  #    }
  #  }
  #  
  #  ## Make it plotly
  #  k2 <- ggplotly(k,
  #                 tooltip = "text") %>% 
  #    
  #    config(displayModeBar = F)  %>% 
  #    
  #    layout(font=list(color="#black"),
  #           xaxis=list(title="UMAP1",zeroline=F),
  #           yaxis=list(title="UMAP2",zeroline=F)) 
  #  
  #  k2 <- k2 %>%
  #    hide_legend() %>%
  #    hide_colorbar()
  #  
  #  if (length(sampSelected) > 0) {
  #    k2 <- k2 %>%
  #      add_annotations(x = plotdata$UMAP1,
  #                      y = plotdata$UMAP2,
  #                      text = rownames(plotdata),
  #                      showarrow = TRUE,
  #                      arrowhead = 4,
  #                      arrowsize = .5)
  #  }
  #  
  #  k2
  #  
  #})
  #
  #umap_plot_PATH_react_expr_base <- reactive({
  #  
  #  ## Variables
  #  sampSelected <- input$UMAPsampSelect         # Sample Names to annotate
  #  UMAPtitleText <- input$UMAPtitleTextSize     # Title text size
  #  UMAPaxisText <- input$UMAPaxisTextSize       # Axis text size
  #  UMAPlegendText <- input$UMAPlegendTextSize   # Legend text size
  #  UMAPdotSize <- input$UMAPdotSize             # Dot size
  #  metaColanno <- NULL
  #  
  #  plot_df <- UMAP_PATH_CoordTable_react()
  #  rownames(plot_df) <- plot_df[,1]
  #  if (is.null(input$UMAPmetaFile) == F) {
  #    if (input$UMAPannotateSamps != " ") {
  #      metaColanno <- input$UMAPannotateSamps
  #      if (input$UMAPannoContCheck != T) {
  #        plot_df[,metaColanno] <- as.factor(plot_df[,metaColanno])
  #      }
  #    }
  #  }
  #  expr <- umap_matrix_loaded_og()
  #  if (!is.null(input$GeneSelection)) {
  #    metaColgene <- input$GeneSelection
  #  }
  #  if (is.null(input$GeneSelection)) {
  #    metaColgene <- rownames(expr)[1]
  #  }
  #  LogChoice <- input$LogGeneSelection
  #  expr2 <- expr
  #  if (LogChoice == TRUE) {
  #    expr2 <- log2(as.matrix(expr) + 0.00001)
  #  }
  #  exprRange <- round(range(expr2[metaColgene,]),2)
  #  metaColkmean <- "Cluster"
  #  
  #  if (!is.null(input$GeneExprRange)) {
  #    exprMin <- as.numeric(gsub(" ","",strsplit(as.character(input$GeneExprRange),",")[[1]][1]))
  #    exprMax <- as.numeric(gsub(" ","",strsplit(as.character(input$GeneExprRange),",")[[1]][2]))
  #  }
  #  else if (is.null(input$GeneExprRange)) {
  #    exprMin <- as.numeric(gsub(" ","",strsplit(as.character(exprRange),",")[[1]][1]))
  #    exprMax <- as.numeric(gsub(" ","",strsplit(as.character(exprRange),",")[[1]][2]))
  #  }
  #  myPalette <- colorRampPalette(rev(brewer.pal(9, input$UMAPcolors)))
  #  
  #  ## Meta Annotation Plot
  #  if (!is.null(metaColanno)) {
  #    colnames(plot_df)[c(4,5)] <- c("AnnoName","GeneName")
  #    m <- plot_df %>%
  #      ggplot(aes(UMAP1, UMAP2, colour=GeneName,
  #                 text = paste("</br> <b>Sample Name:</b> ", SampleName,
  #                              "</br> <b>",metaColanno,":</b> ", AnnoName,
  #                              "</br> <b>",metaColgene," Gene Expression:</b> ", round(GeneName,4),
  #                              "</br> <b>Cluster:</b> ", Cluster,
  #                              sep = "")))
  #  }
  #  if (is.null(metaColanno)) {
  #    colnames(plot_df)[4] <- "GeneName"
  #    m <- plot_df %>%
  #      ggplot(aes(UMAP1, UMAP2, colour=GeneName,
  #                 text = paste("</br> <b>Sample Name:</b> ", SampleName,
  #                              "</br> <b>",metaColgene," Gene Expression:</b> ", round(GeneName,4),
  #                              "</br> <b>Cluster:</b> ", Cluster,
  #                              sep = "")))
  #  }
  #  m <- m + 
  #    geom_point(shape = 19,
  #               size = UMAPdotSize) +
  #    
  #    theme_minimal()
  #  
  #  m <- m + labs(x = "UMAP1",
  #                y = "UMAP2",
  #                color = metaColgene)
  #  
  #  m <-  m + scale_colour_gradientn(colours = rev(myPalette(100)), limits=c(exprMin, exprMax), oob = scales::squish)
  #  
  #  m <- m + theme(axis.text = element_text(size = UMAPaxisText),
  #                 axis.title = element_text(size = UMAPaxisText),
  #                 plot.title = element_text(size = UMAPtitleText),
  #                 legend.text=element_text(size=UMAPlegendText))
  #  
  #  if (length(sampSelected) > 0) {
  #    plotdata <- plot_df[sampSelected,]
  #    plotlabel <- rownames(plot_df[sampSelected,])
  #    if (!is.null(metaColanno)) {
  #      colnames(plotdata)[c(4,5)] <- c("AnnoName","GeneName")
  #    }
  #    if (is.null(metaColanno)) {
  #      colnames(plotdata)[4] <- "GeneName"
  #    }
  #    m <- m + geom_point(data = plotdata,
  #                        aes(x = UMAP1, y = UMAP2),
  #                        pch = 1,
  #                        color = "black",
  #                        size = UMAPdotSize,
  #                        stroke = .3)
  #  }
  #  m
  #  
  #})
  #
  ### UMAP Plot  - Base - expr
  #umap_plot_PATH_react_expr <- reactive({
  #  
  #  
  #  ## Variables}## Variables
  #  sampSelected <- input$UMAPsampSelect         # Sample Names to annotate
  #  metaColanno <- NULL
  #  
  #  plot_df <- UMAP_PATH_CoordTable_react()
  #  rownames(plot_df) <- plot_df[,1]
  #  if (is.null(input$UMAPmetaFile) == F) {
  #    if (input$UMAPannotateSamps != " ") {
  #      metaColanno <- input$UMAPannotateSamps
  #      if (input$UMAPannoContCheck != T) {
  #        plot_df[,metaColanno] <- as.factor(plot_df[,metaColanno])
  #      }
  #    }
  #  }
  #  expr <- umap_matrix_loaded_og()
  #  if (!is.null(input$GeneSelection)) {
  #    metaColgene <- input$GeneSelection
  #  }
  #  if (is.null(input$GeneSelection)) {
  #    metaColgene <- rownames(expr)[1]
  #  }
  #  LogChoice <- input$LogGeneSelection
  #  expr2 <- expr
  #  if (LogChoice == TRUE) {
  #    expr2 <- log2(as.matrix(expr) + 0.00001)
  #  }
  #  exprRange <- round(range(expr2[metaColgene,]),2)
  #  metaColkmean <- "Cluster"
  #  
  #  if (!is.null(input$GeneExprRange)) {
  #    exprMin <- as.numeric(gsub(" ","",strsplit(as.character(input$GeneExprRange),",")[[1]][1]))
  #    exprMax <- as.numeric(gsub(" ","",strsplit(as.character(input$GeneExprRange),",")[[1]][2]))
  #  }
  #  else if (is.null(input$GeneExprRange)) {
  #    exprMin <- as.numeric(gsub(" ","",strsplit(as.character(exprRange),",")[[1]][1]))
  #    exprMax <- as.numeric(gsub(" ","",strsplit(as.character(exprRange),",")[[1]][2]))
  #  }
  #  myPalette <- colorRampPalette(rev(brewer.pal(9, input$UMAPcolors)))
  #  
  #  m <- umap_plot_PATH_react_expr_base()
  #  
  #  if (length(sampSelected) > 0) {
  #    plotdata <- plot_df[sampSelected,]
  #    plotlabel <- rownames(plot_df[sampSelected,])
  #    if (!is.null(metaColanno)) {
  #      colnames(plotdata)[c(4,5)] <- c("AnnoName","GeneName")
  #    }
  #    if (is.null(metaColanno)) {
  #      colnames(plotdata)[4] <- "GeneName"
  #    }
  #  }
  #  
  #  ## Make it plotly
  #  m2 <- ggplotly(m,
  #                 tooltip = "text") %>% 
  #    
  #    config(displayModeBar = F)  %>% 
  #    
  #    layout(font=list(color="#black"),
  #           xaxis=list(title="UMAP1",zeroline=F),
  #           yaxis=list(title="UMAP2",zeroline=F)) %>% 
  #    
  #    hide_legend() %>%
  #    hide_colorbar()
  #  
  #  if (length(sampSelected) > 0) {
  #    m2 <- m2 %>%
  #      add_annotations(x = plotdata$UMAP1,
  #                      y = plotdata$UMAP2,
  #                      text = rownames(plotdata),
  #                      showarrow = TRUE,
  #                      arrowhead = 4,
  #                      arrowsize = .5)
  #  }
  #  
  #  
  #  m2
  #  
  #})
  #
  #umap_plot_PATH_react_kmean_base <- reactive({
  #  
  #  ## Variables
  #  sampSelected <- input$UMAPsampSelect         # Sample Names to annotate
  #  UMAPtitleText <- input$UMAPtitleTextSize     # Title text size
  #  UMAPaxisText <- input$UMAPaxisTextSize       # Axis text size
  #  UMAPlegendText <- input$UMAPlegendTextSize   # Legend text size
  #  UMAPdotSize <- input$UMAPdotSize             # Dot size
  #  metaColanno <- NULL
  #  
  #  plot_df <- UMAP_PATH_CoordTable_react()
  #  rownames(plot_df) <- plot_df[,1]
  #  if (is.null(input$UMAPmetaFile) == F) {
  #    if (input$UMAPannotateSamps != " ") {
  #      metaColanno <- input$UMAPannotateSamps
  #      if (input$UMAPannoContCheck != T) {
  #        plot_df[,metaColanno] <- as.factor(plot_df[,metaColanno])
  #      }
  #    }
  #  }
  #  expr <- umap_matrix_loaded_og()
  #  if (!is.null(input$GeneSelection)) {
  #    metaColgene <- input$GeneSelection
  #  }
  #  if (is.null(input$GeneSelection)) {
  #    metaColgene <- rownames(expr)[1]
  #  }
  #  metaColkmean <- "Cluster"
  #  
  #  ## Meta Annotation Plot
  #  if (!is.null(metaColanno)) {
  #    colnames(plot_df)[c(4,5)] <- c("AnnoName","GeneName")
  #    n <- plot_df %>%
  #      ggplot(aes(UMAP1, UMAP2, colour=Cluster,
  #                 text = paste("</br> <b>Sample Name:</b> ", SampleName,
  #                              "</br> <b>",metaColanno,":</b> ", AnnoName,
  #                              "</br> <b>",metaColgene," Gene Expression:</b> ", round(GeneName,4),
  #                              "</br> <b>Cluster:</b> ", Cluster,
  #                              sep = "")))
  #  }
  #  if (is.null(metaColanno)) {
  #    colnames(plot_df)[4] <- "GeneName"
  #    n <- plot_df %>%
  #      ggplot(aes(UMAP1, UMAP2, colour=Cluster,
  #                 text = paste("</br> <b>Sample Name:</b> ", SampleName,
  #                              "</br> <b>",metaColgene," Gene Expression:</b> ", round(GeneName,4),
  #                              "</br> <b>Cluster:</b> ", Cluster,
  #                              sep = "")))
  #  }
  #  n <- n + geom_point(shape = 19,
  #                      size = UMAPdotSize)+
  #    
  #    labs(x = "UMAP1",
  #         y = "UMAP2",
  #         color = "Cluster")  +
  #    
  #    theme_minimal()
  #  
  #  n <- n + theme(axis.text = element_text(size = UMAPaxisText),
  #                 axis.title = element_text(size = UMAPaxisText),
  #                 plot.title = element_text(size = UMAPtitleText),
  #                 legend.text=element_text(size=UMAPlegendText))
  #  
  #  if (length(sampSelected) > 0) {
  #    plotdata <- plot_df[sampSelected,]
  #    plotlabel <- rownames(plot_df[sampSelected,])
  #    if (!is.null(metaColanno)) {
  #      colnames(plotdata)[c(4,5)] <- c("AnnoName","GeneName")
  #    }
  #    if (is.null(metaColanno)) {
  #      colnames(plotdata)[4] <- "GeneName"
  #    }
  #    n <- n + geom_point(data = plotdata,
  #                        aes(x = UMAP1, y = UMAP2),
  #                        pch = 1,
  #                        color = "black",
  #                        size = UMAPdotSize,
  #                        stroke = .3)
  #  }
  #  n
  #  
  #})
  #
  ### UMAP Plot  - Base - kmean
  #umap_plot_PATH_react_kmean <- reactive({
  #  
  #  ## Variables
  #  sampSelected <- input$UMAPsampSelect         # Sample Names to annotate
  #  metaColanno <- NULL
  #  
  #  plot_df <- UMAP_PATH_CoordTable_react()
  #  rownames(plot_df) <- plot_df[,1]
  #  if (is.null(input$UMAPmetaFile) == F) {
  #    if (input$UMAPannotateSamps != " ") {
  #      metaColanno <- input$UMAPannotateSamps
  #      if (input$UMAPannoContCheck != T) {
  #        plot_df[,metaColanno] <- as.factor(plot_df[,metaColanno])
  #      }
  #    }
  #  }
  #  expr <- umap_matrix_loaded_og()
  #  if (!is.null(input$GeneSelection)) {
  #    metaColgene <- input$GeneSelection
  #  }
  #  if (is.null(input$GeneSelection)) {
  #    metaColgene <- rownames(expr)[1]
  #  }
  #  metaColkmean <- "Cluster"
  #  
  #  n <- umap_plot_PATH_react_kmean_base()
  #  
  #  if (length(sampSelected) > 0) {
  #    plotdata <- plot_df[sampSelected,]
  #    plotlabel <- rownames(plot_df[sampSelected,])
  #    if (!is.null(metaColanno)) {
  #      colnames(plotdata)[c(4,5)] <- c("AnnoName","GeneName")
  #    }
  #    if (is.null(metaColanno)) {
  #      colnames(plotdata)[4] <- "GeneName"
  #    }
  #  }
  #  
  #  ## Make it plotly
  #  n2 <- ggplotly(n,
  #                 tooltip = "text") %>% 
  #    
  #    config(displayModeBar = F)  %>% 
  #    
  #    layout(font=list(color="#black"),
  #           xaxis=list(title="UMAP1",zeroline=F),
  #           yaxis=list(title="UMAP2",zeroline=F)) %>% 
  #    
  #    hide_legend() %>%
  #    hide_colorbar()
  #  if (length(sampSelected) > 0) {
  #    n2 <- n2 %>%
  #      add_annotations(x = plotdata$UMAP1,
  #                      y = plotdata$UMAP2,
  #                      text = rownames(plotdata),
  #                      showarrow = TRUE,
  #                      arrowhead = 4,
  #                      arrowsize = .5)
  #  }
  #  
  #  
  #  n2
  #  
  #})
  #
  #output$UMAPpathLegend <- renderPlot({
  #  
  #  k <- umap_plot_PATH_react_clin_base()
  #  m <- umap_plot_PATH_react_expr_base()
  #  n <- umap_plot_PATH_react_kmean_base()
  #  
  #  if (is.null(input$UMAPmetaFile) == F) {
  #    if (input$UMAPannotateSamps != " ") {
  #      legendk <- g_legend(k)
  #      legendm <- g_legend(m)
  #      legendn <- g_legend(n)
  #      legend_grid <- grid.arrange(legendk,legendm,legendn,nrow = 1)
  #    }
  #    else if (input$UMAPannotateSamps == " ") {
  #      legendm <- g_legend(m)
  #      legendn <- g_legend(n)
  #      legend_grid <- grid.arrange(legendm,legendn,nrow = 1)
  #    }
  #  }
  #  if (is.null(input$UMAPmetaFile) == T) {
  #    legendm <- g_legend(m)
  #    legendn <- g_legend(n)
  #    legend_grid <- grid.arrange(legendm,legendn,nrow = 1)
  #  }
  #  
  #  legend_grid
  #  
  #})
  #
  #UMAPplot_PATH_ALL_react <- reactive({
  #  
  #  k2 <- umap_plot_PATH_react_clin()
  #  m2 <- umap_plot_PATH_react_expr()
  #  n2 <- umap_plot_PATH_react_kmean()
  #  
  #  if (is.null(input$UMAPmetaFile) == F) {
  #    if (input$UMAPannotateSamps != " ") {
  #      UMAPtitlek <- paste("Annotated by",input$UMAPannotateSamps)
  #    }
  #    else if (input$UMAPannotateSamps == " ") {
  #      UMAPtitlek <- ""
  #    }
  #  }
  #  if (is.null(input$UMAPmetaFile) == T) {
  #    UMAPtitlek <- ""
  #  }
  #  expr <- umap_matrix_loaded_og()
  #  if (!is.null(input$GeneSelection)) {
  #    GeneSelec <- input$GeneSelection
  #  }
  #  if (is.null(input$GeneSelection)) {
  #    GeneSelec <- rownames(expr)[1]
  #  }
  #  
  #  if (input$LogGeneSelection == T) {
  #    GeneSelecExpr <- paste(GeneSelec,"Expression (Log2)")
  #  }
  #  if (input$LogGeneSelection == F) {
  #    GeneSelecExpr <- paste(GeneSelec,"Expression")
  #  }
  #  
  #  UMAPtitlem <- paste("Annotated by",GeneSelecExpr)
  #  UMAPtitlen <- paste("Annotated by",input$ClusterMethod,input$ClusterNumber,"Clusters")
  #  
  #  
  #  if (input$UMAPorientation == "Side-by-Side") {
  #    subplot_all_sbs <- subplot(k2, m2, n2, nrows = 1, shareY = TRUE, shareX = TRUE)
  #    plot_titles = list( 
  #      list( 
  #        x = 0,
  #        y = 1,
  #        text = UMAPtitlek,  
  #        xref = "x",  
  #        yref = "paper",  
  #        xanchor = "center",  
  #        yanchor = "bottom",  
  #        showarrow = FALSE 
  #      ),  
  #      list( 
  #        x = 0,
  #        y = 1,
  #        text = UMAPtitlem,  
  #        xref = "x2",
  #        yref = "paper",
  #        xanchor = "center",  
  #        yanchor = "bottom",    
  #        showarrow = FALSE 
  #      ),  
  #      list( 
  #        x = 0,
  #        y = 1,
  #        text = UMAPtitlen,  
  #        xref = "x3",
  #        yref = "paper",
  #        xanchor = "center",  
  #        yanchor = "bottom",  
  #        showarrow = FALSE 
  #      ))
  #  }
  #  else if (input$UMAPorientation == "Stacked") {
  #    subplot_all_sbs <- subplot(k2, m2, n2, nrows = 3, shareY = TRUE, shareX = TRUE, margin = 0.04)
  #    plot_titles = list( 
  #      list( 
  #        x = 0,
  #        y = 1,
  #        text = UMAPtitlek,  
  #        xref = "x",  
  #        yref = "paper",  
  #        xanchor = "center",  
  #        yanchor = "bottom",  
  #        showarrow = FALSE 
  #      ),  
  #      list( 
  #        x = 0,
  #        y = 0.65,
  #        text = UMAPtitlem,  
  #        xref = "x",
  #        yref = "paper",
  #        xanchor = "center",  
  #        yanchor = "bottom",    
  #        showarrow = FALSE 
  #      ),  
  #      list( 
  #        x = 0,
  #        y = 0.3,
  #        text = UMAPtitlen,  
  #        xref = "x",
  #        yref = "paper",
  #        xanchor = "center",  
  #        yanchor = "bottom",  
  #        showarrow = FALSE 
  #      ))
  #  }
  #  
  #  subplot_all_sbs <- subplot_all_sbs %>% layout(annotations = plot_titles)
  #  
  #  subplot_all_sbs
  #})
  #
  #output$UMAPplot_PATH_ALL <- renderPlotly({
  #  p_all <- UMAPplot_PATH_ALL_react()
  #  p_all
  #})
  #
  #####----BASE UMAP----####
  #
  ### UMAP Plot  - PATH - Clin
  #umap_plot_clin_react_base <- reactive({
  #  
  #  ## Variables
  #  sampSelected <- input$UMAPsampSelect         # Sample Names to annotate
  #  UMAPtitleText <- input$UMAPtitleTextSize     # Title text size
  #  UMAPaxisText <- input$UMAPaxisTextSize       # Axis text size
  #  UMAPlegendText <- input$UMAPlegendTextSize   # Legend text size
  #  UMAPdotSize <- input$UMAPdotSize             # Dot size
  #  metaColanno <- NULL
  #  
  #  plot_df <- UMAP_CoordTable_react()
  #  rownames(plot_df) <- plot_df[,1]
  #  if (is.null(input$UMAPmetaFile) == F) {
  #    if (input$UMAPannotateSamps != " ") {
  #      metaColanno <- input$UMAPannotateSamps
  #      if (input$UMAPannoContCheck != T) {
  #        plot_df[,metaColanno] <- as.factor(plot_df[,metaColanno])
  #      }
  #    }
  #  }
  #  expr <- umap_matrix_loaded_og()
  #  if (!is.null(input$GeneSelection)) {
  #    metaColgene <- input$GeneSelection
  #  }
  #  if (is.null(input$GeneSelection)) {
  #    metaColgene <- rownames(expr)[1]
  #  }
  #  metaColkmean <- "Cluster"
  #  
  #  ## Meta Annotation Plot
  #  if (!is.null(metaColanno)) {
  #    colnames(plot_df)[c(4,5)] <- c("AnnoName","GeneName")
  #    k <- plot_df %>%
  #      ggplot(aes(UMAP1, UMAP2, colour=AnnoName,
  #                 text = paste("</br> <b>Sample Name:</b> ", SampleName,
  #                              "</br> <b>",metaColanno,":</b> ", AnnoName,
  #                              "</br> <b>",metaColgene," Gene Expression:</b> ", round(GeneName,4),
  #                              "</br> <b>Cluster:</b> ", Cluster,
  #                              sep = "")))
  #  }
  #  if (is.null(metaColanno)) {
  #    colnames(plot_df)[4] <- "GeneName"
  #    k <- plot_df %>%
  #      ggplot(aes(UMAP1, UMAP2,
  #                 text = paste("</br> <b>Sample Name:</b> ", SampleName,
  #                              "</br> <b>",metaColgene," Gene Expression:</b> ", round(GeneName,4),
  #                              "</br> <b>Cluster:</b> ", Cluster,
  #                              sep = "")))
  #  }
  #  
  #  k <- k + geom_point(shape = 19,
  #                      size = UMAPdotSize) +
  #    
  #    theme_minimal()
  #  
  #  if (is.null(input$UMAPmetaFile) == F) {
  #    if (input$UMAPannotateSamps != " ") {
  #      metaColanno <- input$UMAPannotateSamps
  #      k <- k + labs(x = "UMAP1",
  #                    y = "UMAP2",
  #                    color = metaColanno)
  #    }
  #    if (input$UMAPannotateSamps == " ") {
  #      k <- k + labs(x = "UMAP1",
  #                    y = "UMAP2")
  #    }
  #  }
  #  if (is.null(input$UMAPmetaFile) == T) {
  #    k <- k + labs(x = "UMAP1",
  #                  y = "UMAP2")
  #  }
  #  
  #  if (is.null(input$UMAPmetaFile) == F) {
  #    if (input$UMAPannotateSamps != " ") {
  #      if (input$UMAPannoContCheck == T) {
  #        myPalette <- colorRampPalette(rev(brewer.pal(9, input$UMAPcolors)))
  #        k <- k + scale_colour_gradientn(colours = rev(myPalette(100)))
  #      }
  #    }
  #  }
  #  
  #  k <- k + theme(axis.text = element_text(size = UMAPaxisText),
  #                 axis.title = element_text(size = UMAPaxisText),
  #                 plot.title = element_text(size = UMAPtitleText),
  #                 legend.text=element_text(size=UMAPlegendText))
  #  
  #  if (length(sampSelected) > 0) {
  #    plotdata <- plot_df[sampSelected,]
  #    plotlabel <- rownames(plot_df[sampSelected,])
  #    if (!is.null(metaColanno)) {
  #      colnames(plotdata)[c(4,5)] <- c("AnnoName","GeneName")
  #    }
  #    if (is.null(metaColanno)) {
  #      colnames(plotdata)[4] <- "GeneName"
  #    }
  #    k <- k + geom_point(data = plotdata,
  #                        aes(x = UMAP1, y = UMAP2),
  #                        pch = 1,
  #                        color = "black",
  #                        size = UMAPdotSize,
  #                        stroke = .3)
  #  }
  #  k
  #  
  #})
  #
  #umap_plot_clin_react <- reactive({
  #  
  #  ## Variables
  #  sampSelected <- input$UMAPsampSelect         # Sample Names to annotate
  #  metaColanno <- NULL
  #  
  #  plot_df <- UMAP_CoordTable_react()
  #  rownames(plot_df) <- plot_df[,1]
  #  if (is.null(input$UMAPmetaFile) == F) {
  #    if (input$UMAPannotateSamps != " ") {
  #      metaColanno <- input$UMAPannotateSamps
  #      if (input$UMAPannoContCheck != T) {
  #        plot_df[,metaColanno] <- as.factor(plot_df[,metaColanno])
  #      }
  #    }
  #  }
  #  expr <- umap_matrix_loaded_og()
  #  if (!is.null(input$GeneSelection)) {
  #    metaColgene <- input$GeneSelection
  #  }
  #  if (is.null(input$GeneSelection)) {
  #    metaColgene <- rownames(expr)[1]
  #  }
  #  metaColkmean <- "Cluster"
  #  
  #  k <- umap_plot_clin_react_base()
  #  
  #  if (length(sampSelected) > 0) {
  #    plotdata <- plot_df[sampSelected,]
  #    plotlabel <- rownames(plot_df[sampSelected,])
  #    if (!is.null(metaColanno)) {
  #      colnames(plotdata)[c(4,5)] <- c("AnnoName","GeneName")
  #    }
  #    if (is.null(metaColanno)) {
  #      colnames(plotdata)[4] <- "GeneName"
  #    }
  #  }
  #  
  #  ## Make it plotly
  #  k2 <- ggplotly(k,
  #                 tooltip = "text") %>% 
  #    
  #    config(displayModeBar = F)  %>% 
  #    
  #    layout(font=list(color="#black"),
  #           xaxis=list(title="UMAP1",zeroline=F),
  #           yaxis=list(title="UMAP2",zeroline=F)) 
  #  
  #  k2 <- k2 %>%
  #    hide_legend() %>%
  #    hide_colorbar()
  #  
  #  if (length(sampSelected) > 0) {
  #    k2 <- k2 %>%
  #      add_annotations(x = plotdata$UMAP1,
  #                      y = plotdata$UMAP2,
  #                      text = rownames(plotdata),
  #                      showarrow = TRUE,
  #                      arrowhead = 4,
  #                      arrowsize = .5)
  #  }
  #  
  #  k2
  #  
  #})
  #
  #umap_plot_expr_react_base <- reactive({
  #  
  #  ## Variables
  #  sampSelected <- input$UMAPsampSelect         # Sample Names to annotate
  #  UMAPtitleText <- input$UMAPtitleTextSize     # Title text size
  #  UMAPaxisText <- input$UMAPaxisTextSize       # Axis text size
  #  UMAPlegendText <- input$UMAPlegendTextSize   # Legend text size
  #  UMAPdotSize <- input$UMAPdotSize             # Dot size
  #  metaColanno <- NULL
  #  
  #  plot_df <- UMAP_CoordTable_react()
  #  rownames(plot_df) <- plot_df[,1]
  #  if (is.null(input$UMAPmetaFile) == F) {
  #    if (input$UMAPannotateSamps != " ") {
  #      metaColanno <- input$UMAPannotateSamps
  #      if (input$UMAPannoContCheck != T) {
  #        plot_df[,metaColanno] <- as.factor(plot_df[,metaColanno])
  #      }
  #    }
  #  }
  #  expr <- umap_matrix_loaded_og()
  #  if (!is.null(input$GeneSelection)) {
  #    metaColgene <- input$GeneSelection
  #  }
  #  if (is.null(input$GeneSelection)) {
  #    metaColgene <- rownames(expr)[1]
  #  }
  #  LogChoice <- input$LogGeneSelection
  #  expr2 <- expr
  #  if (LogChoice == TRUE) {
  #    expr2 <- log2(as.matrix(expr) + 0.00001)
  #  }
  #  exprRange <- round(range(expr2[metaColgene,]),2)
  #  metaColkmean <- "Cluster"
  #  
  #  if (!is.null(input$GeneExprRange)) {
  #    exprMin <- as.numeric(gsub(" ","",strsplit(as.character(input$GeneExprRange),",")[[1]][1]))
  #    exprMax <- as.numeric(gsub(" ","",strsplit(as.character(input$GeneExprRange),",")[[1]][2]))
  #  }
  #  else if (is.null(input$GeneExprRange)) {
  #    exprMin <- as.numeric(gsub(" ","",strsplit(as.character(exprRange),",")[[1]][1]))
  #    exprMax <- as.numeric(gsub(" ","",strsplit(as.character(exprRange),",")[[1]][2]))
  #  }
  #  myPalette <- colorRampPalette(rev(brewer.pal(9, input$UMAPcolors)))
  #  
  #  ## Meta Annotation Plot
  #  if (!is.null(metaColanno)) {
  #    colnames(plot_df)[c(4,5)] <- c("AnnoName","GeneName")
  #    m <- plot_df %>%
  #      ggplot(aes(UMAP1, UMAP2, colour=GeneName,
  #                 text = paste("</br> <b>Sample Name:</b> ", SampleName,
  #                              "</br> <b>",metaColanno,":</b> ", AnnoName,
  #                              "</br> <b>",metaColgene," Gene Expression:</b> ", round(GeneName,4),
  #                              "</br> <b>Cluster:</b> ", Cluster,
  #                              sep = "")))
  #  }
  #  if (is.null(metaColanno)) {
  #    colnames(plot_df)[4] <- "GeneName"
  #    m <- plot_df %>%
  #      ggplot(aes(UMAP1, UMAP2, colour=GeneName,
  #                 text = paste("</br> <b>Sample Name:</b> ", SampleName,
  #                              "</br> <b>",metaColgene," Gene Expression:</b> ", round(GeneName,4),
  #                              "</br> <b>Cluster:</b> ", Cluster,
  #                              sep = "")))
  #  }
  #  m <- m + 
  #    geom_point(shape = 19,
  #               size = UMAPdotSize) +
  #    
  #    theme_minimal()
  #  
  #  m <- m + labs(x = "UMAP1",
  #                y = "UMAP2",
  #                color = metaColgene)
  #  
  #  m <-  m + scale_colour_gradientn(colours = rev(myPalette(100)), limits=c(exprMin, exprMax), oob = scales::squish)
  #  
  #  m <- m + theme(axis.text = element_text(size = UMAPaxisText),
  #                 axis.title = element_text(size = UMAPaxisText),
  #                 plot.title = element_text(size = UMAPtitleText),
  #                 legend.text=element_text(size=UMAPlegendText))
  #  
  #  if (length(sampSelected) > 0) {
  #    plotdata <- plot_df[sampSelected,]
  #    plotlabel <- rownames(plot_df[sampSelected,])
  #    if (!is.null(metaColanno)) {
  #      colnames(plotdata)[c(4,5)] <- c("AnnoName","GeneName")
  #    }
  #    if (is.null(metaColanno)) {
  #      colnames(plotdata)[4] <- "GeneName"
  #    }
  #    m <- m + geom_point(data = plotdata,
  #                        aes(x = UMAP1, y = UMAP2),
  #                        pch = 1,
  #                        color = "black",
  #                        size = UMAPdotSize,
  #                        stroke = .3)
  #  }
  #  m
  #  
  #})
  #
  ### UMAP Plot  - Base - expr
  #umap_plot_expr_react <- reactive({
  #  
  #  
  #  ## Variables}## Variables
  #  sampSelected <- input$UMAPsampSelect         # Sample Names to annotate
  #  metaColanno <- NULL
  #  
  #  plot_df <- UMAP_CoordTable_react()
  #  rownames(plot_df) <- plot_df[,1]
  #  if (is.null(input$UMAPmetaFile) == F) {
  #    if (input$UMAPannotateSamps != " ") {
  #      metaColanno <- input$UMAPannotateSamps
  #      if (input$UMAPannoContCheck != T) {
  #        plot_df[,metaColanno] <- as.factor(plot_df[,metaColanno])
  #      }
  #    }
  #  }
  #  expr <- umap_matrix_loaded_og()
  #  if (!is.null(input$GeneSelection)) {
  #    metaColgene <- input$GeneSelection
  #  }
  #  if (is.null(input$GeneSelection)) {
  #    metaColgene <- rownames(expr)[1]
  #  }
  #  LogChoice <- input$LogGeneSelection
  #  expr2 <- expr
  #  if (LogChoice == TRUE) {
  #    expr2 <- log2(as.matrix(expr) + 0.00001)
  #  }
  #  exprRange <- round(range(expr2[metaColgene,]),2)
  #  metaColkmean <- "Cluster"
  #  
  #  if (!is.null(input$GeneExprRange)) {
  #    exprMin <- as.numeric(gsub(" ","",strsplit(as.character(input$GeneExprRange),",")[[1]][1]))
  #    exprMax <- as.numeric(gsub(" ","",strsplit(as.character(input$GeneExprRange),",")[[1]][2]))
  #  }
  #  else if (is.null(input$GeneExprRange)) {
  #    exprMin <- as.numeric(gsub(" ","",strsplit(as.character(exprRange),",")[[1]][1]))
  #    exprMax <- as.numeric(gsub(" ","",strsplit(as.character(exprRange),",")[[1]][2]))
  #  }
  #  myPalette <- colorRampPalette(rev(brewer.pal(9, input$UMAPcolors)))
  #  
  #  m <- umap_plot_expr_react_base()
  #  
  #  if (length(sampSelected) > 0) {
  #    plotdata <- plot_df[sampSelected,]
  #    plotlabel <- rownames(plot_df[sampSelected,])
  #    if (!is.null(metaColanno)) {
  #      colnames(plotdata)[c(4,5)] <- c("AnnoName","GeneName")
  #    }
  #    if (is.null(metaColanno)) {
  #      colnames(plotdata)[4] <- "GeneName"
  #    }
  #  }
  #  
  #  ## Make it plotly
  #  m2 <- ggplotly(m,
  #                 tooltip = "text") %>% 
  #    
  #    config(displayModeBar = F)  %>% 
  #    
  #    layout(font=list(color="#black"),
  #           xaxis=list(title="UMAP1",zeroline=F),
  #           yaxis=list(title="UMAP2",zeroline=F)) %>% 
  #    
  #    hide_legend() %>%
  #    hide_colorbar()
  #  
  #  if (length(sampSelected) > 0) {
  #    m2 <- m2 %>%
  #      add_annotations(x = plotdata$UMAP1,
  #                      y = plotdata$UMAP2,
  #                      text = rownames(plotdata),
  #                      showarrow = TRUE,
  #                      arrowhead = 4,
  #                      arrowsize = .5)
  #  }
  #  
  #  
  #  m2
  #  
  #})
  #
  #umap_plot_kmean_react_base <- reactive({
  #  
  #  ## Variables
  #  sampSelected <- input$UMAPsampSelect         # Sample Names to annotate
  #  UMAPtitleText <- input$UMAPtitleTextSize     # Title text size
  #  UMAPaxisText <- input$UMAPaxisTextSize       # Axis text size
  #  UMAPlegendText <- input$UMAPlegendTextSize   # Legend text size
  #  UMAPdotSize <- input$UMAPdotSize             # Dot size
  #  metaColanno <- NULL
  #  
  #  plot_df <- UMAP_CoordTable_react()
  #  rownames(plot_df) <- plot_df[,1]
  #  if (is.null(input$UMAPmetaFile) == F) {
  #    if (input$UMAPannotateSamps != " ") {
  #      metaColanno <- input$UMAPannotateSamps
  #      if (input$UMAPannoContCheck != T) {
  #        plot_df[,metaColanno] <- as.factor(plot_df[,metaColanno])
  #      }
  #    }
  #  }
  #  expr <- umap_matrix_loaded_og()
  #  if (!is.null(input$GeneSelection)) {
  #    metaColgene <- input$GeneSelection
  #  }
  #  if (is.null(input$GeneSelection)) {
  #    metaColgene <- rownames(expr)[1]
  #  }
  #  metaColkmean <- "Cluster"
  #  
  #  ## Meta Annotation Plot
  #  if (!is.null(metaColanno)) {
  #    colnames(plot_df)[c(4,5)] <- c("AnnoName","GeneName")
  #    n <- plot_df %>%
  #      ggplot(aes(UMAP1, UMAP2, colour=Cluster,
  #                 text = paste("</br> <b>Sample Name:</b> ", SampleName,
  #                              "</br> <b>",metaColanno,":</b> ", AnnoName,
  #                              "</br> <b>",metaColgene," Gene Expression:</b> ", round(GeneName,4),
  #                              "</br> <b>Cluster:</b> ", Cluster,
  #                              sep = "")))
  #  }
  #  if (is.null(metaColanno)) {
  #    colnames(plot_df)[4] <- "GeneName"
  #    n <- plot_df %>%
  #      ggplot(aes(UMAP1, UMAP2, colour=Cluster,
  #                 text = paste("</br> <b>Sample Name:</b> ", SampleName,
  #                              "</br> <b>",metaColgene," Gene Expression:</b> ", round(GeneName,4),
  #                              "</br> <b>Cluster:</b> ", Cluster,
  #                              sep = "")))
  #  }
  #  n <- n + geom_point(shape = 19,
  #                      size = UMAPdotSize)+
  #    
  #    labs(x = "UMAP1",
  #         y = "UMAP2",
  #         color = "Cluster")  +
  #    
  #    theme_minimal()
  #  
  #  n <- n + theme(axis.text = element_text(size = UMAPaxisText),
  #                 axis.title = element_text(size = UMAPaxisText),
  #                 plot.title = element_text(size = UMAPtitleText),
  #                 legend.text=element_text(size=UMAPlegendText))
  #  
  #  if (length(sampSelected) > 0) {
  #    plotdata <- plot_df[sampSelected,]
  #    plotlabel <- rownames(plot_df[sampSelected,])
  #    if (!is.null(metaColanno)) {
  #      colnames(plotdata)[c(4,5)] <- c("AnnoName","GeneName")
  #    }
  #    if (is.null(metaColanno)) {
  #      colnames(plotdata)[4] <- "GeneName"
  #    }
  #    n <- n + geom_point(data = plotdata,
  #                        aes(x = UMAP1, y = UMAP2),
  #                        pch = 1,
  #                        color = "black",
  #                        size = UMAPdotSize,
  #                        stroke = .3)
  #  }
  #  n
  #  
  #})
  #
  ### UMAP Plot  - Base - kmean
  #umap_plot_kmean_react <- reactive({
  #  
  #  ## Variables
  #  sampSelected <- input$UMAPsampSelect         # Sample Names to annotate
  #  metaColanno <- NULL
  #  
  #  plot_df <- UMAP_CoordTable_react()
  #  rownames(plot_df) <- plot_df[,1]
  #  if (is.null(input$UMAPmetaFile) == F) {
  #    if (input$UMAPannotateSamps != " ") {
  #      metaColanno <- input$UMAPannotateSamps
  #      if (input$UMAPannoContCheck != T) {
  #        plot_df[,metaColanno] <- as.factor(plot_df[,metaColanno])
  #      }
  #    }
  #  }
  #  expr <- umap_matrix_loaded_og()
  #  if (!is.null(input$GeneSelection)) {
  #    metaColgene <- input$GeneSelection
  #  }
  #  if (is.null(input$GeneSelection)) {
  #    metaColgene <- rownames(expr)[1]
  #  }
  #  metaColkmean <- "Cluster"
  #  
  #  n <- umap_plot_kmean_react_base()
  #  
  #  if (length(sampSelected) > 0) {
  #    plotdata <- plot_df[sampSelected,]
  #    plotlabel <- rownames(plot_df[sampSelected,])
  #    if (!is.null(metaColanno)) {
  #      colnames(plotdata)[c(4,5)] <- c("AnnoName","GeneName")
  #    }
  #    if (is.null(metaColanno)) {
  #      colnames(plotdata)[4] <- "GeneName"
  #    }
  #  }
  #  
  #  ## Make it plotly
  #  n2 <- ggplotly(n,
  #                 tooltip = "text") %>% 
  #    
  #    config(displayModeBar = F)  %>% 
  #    
  #    layout(font=list(color="#black"),
  #           xaxis=list(title="UMAP1",zeroline=F),
  #           yaxis=list(title="UMAP2",zeroline=F)) %>% 
  #    
  #    hide_legend() %>%
  #    hide_colorbar()
  #  if (length(sampSelected) > 0) {
  #    n2 <- n2 %>%
  #      add_annotations(x = plotdata$UMAP1,
  #                      y = plotdata$UMAP2,
  #                      text = rownames(plotdata),
  #                      showarrow = TRUE,
  #                      arrowhead = 4,
  #                      arrowsize = .5)
  #  }
  #  
  #  
  #  n2
  #  
  #})
  #
  #output$UMAPallLegend <- renderPlot({
  #  
  #  k <- umap_plot_clin_react_base()
  #  m <- umap_plot_expr_react_base()
  #  n <- umap_plot_kmean_react_base()
  #  
  #  if (is.null(input$UMAPmetaFile) == F) {
  #    if (input$UMAPannotateSamps != " ") {
  #      legendk <- g_legend(k)
  #      legendm <- g_legend(m)
  #      legendn <- g_legend(n)
  #      legend_grid <- grid.arrange(legendk,legendm,legendn,nrow = 1)
  #    }
  #    else if (input$UMAPannotateSamps == " ") {
  #      legendm <- g_legend(m)
  #      legendn <- g_legend(n)
  #      legend_grid <- grid.arrange(legendm,legendn,nrow = 1)
  #    }
  #  }
  #  if (is.null(input$UMAPmetaFile) == T) {
  #    legendm <- g_legend(m)
  #    legendn <- g_legend(n)
  #    legend_grid <- grid.arrange(legendm,legendn,nrow = 1)
  #  }
  #  
  #  legend_grid
  #  
  #})
  #
  #UMAPplot_ALL_react <- reactive({
  #  
  #  k2 <- umap_plot_clin_react()
  #  m2 <- umap_plot_expr_react()
  #  n2 <- umap_plot_kmean_react()
  #  
  #  if (is.null(input$UMAPmetaFile) == F) {
  #    if (input$UMAPannotateSamps != " ") {
  #      UMAPtitlek <- paste("Annotated by",input$UMAPannotateSamps)
  #    }
  #    else if (input$UMAPannotateSamps == " ") {
  #      UMAPtitlek <- ""
  #    }
  #  }
  #  if (is.null(input$UMAPmetaFile) == T) {
  #    UMAPtitlek <- ""
  #  }
  #  expr <- umap_matrix_loaded_og()
  #  if (!is.null(input$GeneSelection)) {
  #    GeneSelec <- input$GeneSelection
  #  }
  #  if (is.null(input$GeneSelection)) {
  #    GeneSelec <- rownames(expr)[1]
  #  }
  #  
  #  if (input$LogGeneSelection == T) {
  #    GeneSelecExpr <- paste(GeneSelec,"Expression (Log2)")
  #  }
  #  if (input$LogGeneSelection == F) {
  #    GeneSelecExpr <- paste(GeneSelec,"Expression")
  #  }
  #  
  #  UMAPtitlem <- paste("Annotated by",GeneSelecExpr)
  #  UMAPtitlen <- paste("Annotated by",input$ClusterMethod,input$ClusterNumber,"Clusters")
  #  
  #  
  #  if (input$UMAPorientation == "Side-by-Side") {
  #    subplot_all_sbs <- subplot(k2, m2, n2, nrows = 1, shareY = TRUE, shareX = TRUE)
  #    plot_titles = list( 
  #      list( 
  #        x = 0,
  #        y = 1,
  #        text = UMAPtitlek,  
  #        xref = "x",  
  #        yref = "paper",  
  #        xanchor = "center",  
  #        yanchor = "bottom",  
  #        showarrow = FALSE 
  #      ),  
  #      list( 
  #        x = 0,
  #        y = 1,
  #        text = UMAPtitlem,  
  #        xref = "x2",
  #        yref = "paper",
  #        xanchor = "center",  
  #        yanchor = "bottom",    
  #        showarrow = FALSE 
  #      ),  
  #      list( 
  #        x = 0,
  #        y = 1,
  #        text = UMAPtitlen,  
  #        xref = "x3",
  #        yref = "paper",
  #        xanchor = "center",  
  #        yanchor = "bottom",  
  #        showarrow = FALSE 
  #      ))
  #  }
  #  else if (input$UMAPorientation == "Stacked") {
  #    subplot_all_sbs <- subplot(k2, m2, n2, nrows = 3, shareY = TRUE, shareX = TRUE, margin = 0.04)
  #    plot_titles = list( 
  #      list( 
  #        x = 0,
  #        y = 1,
  #        text = UMAPtitlek,  
  #        xref = "x",  
  #        yref = "paper",  
  #        xanchor = "center",  
  #        yanchor = "bottom",  
  #        showarrow = FALSE 
  #      ),  
  #      list( 
  #        x = 0,
  #        y = 0.65,
  #        text = UMAPtitlem,  
  #        xref = "x",
  #        yref = "paper",
  #        xanchor = "center",  
  #        yanchor = "bottom",    
  #        showarrow = FALSE 
  #      ),  
  #      list( 
  #        x = 0,
  #        y = 0.3,
  #        text = UMAPtitlen,  
  #        xref = "x",
  #        yref = "paper",
  #        xanchor = "center",  
  #        yanchor = "bottom",  
  #        showarrow = FALSE 
  #      ))
  #  }
  #  
  #  subplot_all_sbs <- subplot_all_sbs %>% layout(annotations = plot_titles)
  #  
  #  subplot_all_sbs
  #})
  #
  #output$UMAPplot_ALL <- renderPlotly({
  #  p_all <- UMAPplot_ALL_react()
  #  p_all
  #})
  #
  #####----PreC UMAP----####
  #
  #####----PreC UMAP----####
  #
  ### UMAP Plot  - PATH - Clin
  #umap_plot_PreC_clin_react_base <- reactive({
  #  
  #  req(input$UMAPmetaFile)
  #  ## Variables
  #  sampSelected <- input$UMAPsampSelect         # Sample Names to annotate
  #  UMAPtitleText <- input$UMAPtitleTextSize     # Title text size
  #  UMAPaxisText <- input$UMAPaxisTextSize       # Axis text size
  #  UMAPlegendText <- input$UMAPlegendTextSize   # Legend text size
  #  UMAPdotSize <- input$UMAPdotSize             # Dot size
  #  metaColanno <- NULL
  #  
  #  umap1 <- input$SelectPreCalc1
  #  umap2 <- input$SelectPreCalc2
  #  
  #  plot_df <- UMAP_PreC_CoordTable_react()
  #  rownames(plot_df) <- plot_df[,1]
  #  if (input$UMAPannotateSamps != " ") {
  #    metaColanno <- input$UMAPannotateSamps
  #    if (input$UMAPannoContCheck != T) {
  #      plot_df[,metaColanno] <- as.factor(plot_df[,metaColanno])
  #    }
  #  }
  #  if (!is.null(input$UMAPmatrixFile)) {
  #    expr <- umap_matrix_loaded_og()
  #    if (!is.null(input$GeneSelection)) {
  #      metaColgene <- input$GeneSelection
  #    }
  #    if (is.null(input$GeneSelection)) {
  #      metaColgene <- rownames(expr)[1]
  #    }
  #    metaColkmean <- "Cluster"
  #  }
  #  
  #  print(head(as_tibble(plot_df)))
  #  
  #  ## Meta Annotation Plot
  #  if (!is.null(input$UMAPmatrixFile)) {
  #    if (input$UMAPannotateSamps != " ") {
  #      colnames(plot_df)[c(4,5)] <- c("AnnoName","GeneName")
  #      print(head(as_tibble(plot_df)))
  #      k <- plot_df %>%
  #        ggplot(aes(UMAP1, UMAP2, colour=AnnoName,
  #                   text = paste("</br> <b>Sample Name:</b> ", SampleName,
  #                                "</br> <b>",metaColanno,":</b> ", AnnoName,
  #                                "</br> <b>",metaColgene," Gene Expression:</b> ", round(GeneName,4),
  #                                "</br> <b>Cluster:</b> ", Cluster,
  #                                sep = "")))
  #    }
  #    else if (input$UMAPannotateSamps == " ") {
  #      colnames(plot_df)[4] <- "GeneName"
  #      print(head(as_tibble(plot_df)))
  #      k <- plot_df %>%
  #        ggplot(aes(UMAP1, UMAP2,
  #                   text = paste("</br> <b>Sample Name:</b> ", SampleName,
  #                                "</br> <b>",metaColgene," Gene Expression:</b> ", round(GeneName,4),
  #                                "</br> <b>Cluster:</b> ", Cluster,
  #                                sep = "")))
  #    }
  #    
  #  }
  #  else if (is.null(input$UMAPmatrixFile)) {
  #    if (input$UMAPannotateSamps != " ") {
  #      colnames(plot_df)[4] <- c("AnnoName")
  #      print(head(as_tibble(plot_df)))
  #      k <- plot_df %>%
  #        ggplot(aes(UMAP1, UMAP2, colour=AnnoName,
  #                   text = paste("</br> <b>Sample Name:</b> ", SampleName,
  #                                "</br> <b>",metaColanno,":</b> ", AnnoName,
  #                                sep = "")))
  #    }
  #    else if (input$UMAPannotateSamps == " ") {
  #      k <- plot_df %>%
  #        ggplot(aes(UMAP1, UMAP2,
  #                   text = paste("</br> <b>Sample Name:</b> ", SampleName,
  #                                sep = "")))
  #    }
  #    
  #  }
  #  #colnames(plot_df)[c(4,5)] <- c("AnnoName","GeneName")
  #  #k <- plot_df %>%
  #  #  ggplot(aes(UMAP1, UMAP2, colour=AnnoName,
  #  #             text = paste("</br> <b>Sample Name:</b> ", SampleName,
  #  #                          "</br> <b>",metaColanno,":</b> ", AnnoName,
  #  #                          "</br> <b>",metaColgene," Gene Expression:</b> ", round(GeneName,4),
  #  #                          "</br> <b>Cluster:</b> ", Cluster,
  #  #                          sep = "")))
  #  
  #  k <- k + geom_point(shape = 19,
  #                      size = UMAPdotSize) +
  #    
  #    theme_minimal()
  #  
  #  if (input$UMAPannotateSamps != " ") {
  #    metaColanno <- input$UMAPannotateSamps
  #    k <- k + labs(x = umap1,
  #                  y = umap2,
  #                  color = metaColanno)
  #  }
  #  if (input$UMAPannotateSamps == " ") {
  #    k <- k + labs(x = umap1,
  #                  y = umap2,)
  #  }
  #  
  #  if (input$UMAPannotateSamps != " ") {
  #    if (input$UMAPannoContCheck == T) {
  #      myPalette <- colorRampPalette(rev(brewer.pal(9, input$UMAPcolors)))
  #      k <- k + scale_colour_gradientn(colours = rev(myPalette(100)))
  #    }
  #  }
  #  
  #  k <- k + theme(axis.text = element_text(size = UMAPaxisText),
  #                 axis.title = element_text(size = UMAPaxisText),
  #                 plot.title = element_text(size = UMAPtitleText),
  #                 legend.text=element_text(size=UMAPlegendText))
  #  
  #  if (length(sampSelected) > 0) {
  #    plotdata <- plot_df[sampSelected,]
  #    plotlabel <- rownames(plot_df[sampSelected,])
  #    #if (!is.null(metaColanno)) {
  #    #  if (input$UMAPmatrixFile == T) {
  #    #    colnames(plot_df)[c(4,5)] <- c("AnnoName","GeneName")
  #    #  }
  #    #  if (input$UMAPmatrixFile == F) {
  #    #    colnames(plot_df)[4] <- c("AnnoName","GeneName")
  #    #  }
  #    #}
  #    #if (!is.null(metaColanno)) {
  #    #  colnames(plotdata)[c(4,5)] <- c("AnnoName","GeneName")
  #    #}
  #    #if (is.null(metaColanno)) {
  #    #  colnames(plotdata)[4] <- "GeneName"
  #    #}
  #    k <- k + geom_point(data = plotdata,
  #                        aes(x = UMAP1, y = UMAP2),
  #                        pch = 1,
  #                        color = "black",
  #                        size = UMAPdotSize,
  #                        stroke = .3)
  #  }
  #  k
  #  
  #})
  #
  #umap_plot_PreC_clin_react <- reactive({
  #  
  #  req(input$UMAPmetaFile)
  #  ## Variables
  #  sampSelected <- input$UMAPsampSelect         # Sample Names to annotate
  #  metaColanno <- NULL
  #  
  #  plot_df <- UMAP_PreC_CoordTable_react()
  #  rownames(plot_df) <- plot_df[,1]
  #  
  #  umap1 <- input$SelectPreCalc1
  #  umap2 <- input$SelectPreCalc2
  #  
  #  if (input$UMAPannotateSamps != " ") {
  #    metaColanno <- input$UMAPannotateSamps
  #    if (input$UMAPannoContCheck != T) {
  #      plot_df[,metaColanno] <- as.factor(plot_df[,metaColanno])
  #    }
  #  }
  #  if (!is.null(input$UMAPmatrixFile)) {
  #    expr <- umap_matrix_loaded_og()
  #    if (!is.null(input$GeneSelection)) {
  #      metaColgene <- input$GeneSelection
  #    }
  #    if (is.null(input$GeneSelection)) {
  #      metaColgene <- rownames(expr)[1]
  #    }
  #    metaColkmean <- "Cluster"
  #  }
  #  
  #  
  #  k <- umap_plot_PreC_clin_react_base()
  #  
  #  if (length(sampSelected) > 0) {
  #    plotdata <- plot_df[sampSelected,]
  #    plotlabel <- rownames(plot_df[sampSelected,])
  #    if (!is.null(metaColanno)) {
  #      colnames(plotdata)[c(4,5)] <- c("AnnoName","GeneName")
  #    }
  #    if (is.null(metaColanno)) {
  #      colnames(plotdata)[4] <- "GeneName"
  #    }
  #  }
  #  print("anno df")
  #  print(head(as_tibble(k[["data"]])))
  #  
  #  ## Make it plotly
  #  k2 <- ggplotly(k,
  #                 tooltip = "text") %>% 
  #    
  #    config(displayModeBar = F)  %>% 
  #    
  #    layout(font=list(color="#black"),
  #           xaxis=list(title=umap1,zeroline=F),
  #           yaxis=list(title=umap2,zeroline=F)) 
  #  
  #  k2 <- k2 %>%
  #    hide_legend() %>%
  #    hide_colorbar()
  #  
  #  if (length(sampSelected) > 0) {
  #    k2 <- k2 %>%
  #      add_annotations(x = plotdata$UMAP1,
  #                      y = plotdata$UMAP2,
  #                      text = rownames(plotdata),
  #                      showarrow = TRUE,
  #                      arrowhead = 4,
  #                      arrowsize = .5)
  #  }
  #  
  #  k2
  #  
  #})
  #
  ### UMAP Plot  - PATH - Clin
  #umap_plot_PreC_expr_react_base <- reactive({
  #  
  #  req(input$UMAPmetaFile)
  #  req(input$UMAPmatrixFile)
  #  ## Variables
  #  sampSelected <- input$UMAPsampSelect         # Sample Names to annotate
  #  UMAPtitleText <- input$UMAPtitleTextSize     # Title text size
  #  UMAPaxisText <- input$UMAPaxisTextSize       # Axis text size
  #  UMAPlegendText <- input$UMAPlegendTextSize   # Legend text size
  #  UMAPdotSize <- input$UMAPdotSize             # Dot size
  #  metaColanno <- NULL
  #  
  #  umap1 <- input$SelectPreCalc1
  #  umap2 <- input$SelectPreCalc2
  #  
  #  plot_df <- UMAP_PreC_CoordTable_react()
  #  rownames(plot_df) <- plot_df[,1]
  #  if (input$UMAPannotateSamps != " ") {
  #    metaColanno <- input$UMAPannotateSamps
  #    if (input$UMAPannoContCheck != T) {
  #      plot_df[,metaColanno] <- as.factor(plot_df[,metaColanno])
  #    }
  #  }
  #  expr <- umap_matrix_loaded_og()
  #  if (!is.null(input$GeneSelection)) {
  #    metaColgene <- input$GeneSelection
  #  }
  #  if (is.null(input$GeneSelection)) {
  #    metaColgene <- rownames(expr)[1]
  #  }
  #  metaColkmean <- "Cluster"
  #  
  #  expr <- umap_matrix_loaded_og()
  #  if (!is.null(input$GeneSelection)) {
  #    metaColgene <- input$GeneSelection
  #  }
  #  if (is.null(input$GeneSelection)) {
  #    metaColgene <- rownames(expr)[1]
  #  }
  #  LogChoice <- input$LogGeneSelection
  #  expr2 <- expr
  #  if (LogChoice == TRUE) {
  #    expr2 <- log2(as.matrix(expr) + 0.00001)
  #  }
  #  exprRange <- round(range(expr2[metaColgene,]),2)
  #  
  #  if (!is.null(input$GeneExprRange)) {
  #    exprMin <- as.numeric(gsub(" ","",strsplit(as.character(input$GeneExprRange),",")[[1]][1]))
  #    exprMax <- as.numeric(gsub(" ","",strsplit(as.character(input$GeneExprRange),",")[[1]][2]))
  #  }
  #  else if (is.null(input$GeneExprRange)) {
  #    exprMin <- as.numeric(gsub(" ","",strsplit(as.character(exprRange),",")[[1]][1]))
  #    exprMax <- as.numeric(gsub(" ","",strsplit(as.character(exprRange),",")[[1]][2]))
  #  }
  #  myPalette <- colorRampPalette(rev(brewer.pal(9, input$UMAPcolors)))
  #  
  #  ## Expr Annotation Plot
  #  if (!is.null(input$UMAPmatrixFile)) {
  #    if (input$UMAPannotateSamps != " ") {
  #      colnames(plot_df)[c(4,5)] <- c("AnnoName","GeneName")
  #      print(head(as_tibble(plot_df)))
  #      k <- plot_df %>%
  #        ggplot(aes(UMAP1, UMAP2, colour=GeneName,
  #                   text = paste("</br> <b>Sample Name:</b> ", SampleName,
  #                                "</br> <b>",metaColanno,":</b> ", AnnoName,
  #                                "</br> <b>",metaColgene," Gene Expression:</b> ", round(GeneName,4),
  #                                "</br> <b>Cluster:</b> ", Cluster,
  #                                sep = "")))
  #    }
  #    else if (input$UMAPannotateSamps == " ") {
  #      colnames(plot_df)[4] <- "GeneName"
  #      print(head(as_tibble(plot_df)))
  #      k <- plot_df %>%
  #        ggplot(aes(UMAP1, UMAP2, colour=GeneName,
  #                   text = paste("</br> <b>Sample Name:</b> ", SampleName,
  #                                "</br> <b>",metaColgene," Gene Expression:</b> ", round(GeneName,4),
  #                                "</br> <b>Cluster:</b> ", Cluster,
  #                                sep = "")))
  #    }
  #    
  #  }
  #  
  #  
  #  
  #  #colnames(plot_df)[c(4,5)] <- c("AnnoName","GeneName")
  #  #k <- plot_df %>%
  #  #  ggplot(aes(UMAP1, UMAP2, colour=GeneName,
  #  #             text = paste("</br> <b>Sample Name:</b> ", SampleName,
  #  #                          "</br> <b>",metaColanno,":</b> ", AnnoName,
  #  #                          "</br> <b>",metaColgene," Gene Expression:</b> ", round(GeneName,4),
  #  #                          "</br> <b>Cluster:</b> ", Cluster,
  #  #                          sep = "")))
  #  #colnames(plot_df)[c(4,5)] <- c("AnnoName","GeneName")
  #  #k <- plot_df %>%
  #  #  ggplot(aes(UMAP1, UMAP2, colour=AnnoName,
  #  #             text = paste("</br> <b>Sample Name:</b> ", SampleName,
  #  #                          "</br> <b>",metaColanno,":</b> ", AnnoName,
  #  #                          "</br> <b>",metaColgene," Gene Expression:</b> ", round(GeneName,4),
  #  #                          "</br> <b>Cluster:</b> ", Cluster,
  #  #                          sep = "")))
  #  
  #  k <- k + geom_point(shape = 19,
  #                      size = UMAPdotSize) +
  #    
  #    theme_minimal()
  #  
  #  if (input$UMAPannotateSamps != " ") {
  #    metaColanno <- input$UMAPannotateSamps
  #    k <- k + labs(x = umap1,
  #                  y = umap2,
  #                  color = metaColgene)
  #  }
  #  if (input$UMAPannotateSamps == " ") {
  #    k <- k + labs(x = umap1,
  #                  y = umap2,)
  #  }
  #  k <- k + scale_colour_gradientn(colours = rev(myPalette(100)), limits=c(exprMin, exprMax), oob = scales::squish)
  #  
  #  
  #  k <- k + theme(axis.text = element_text(size = UMAPaxisText),
  #                 axis.title = element_text(size = UMAPaxisText),
  #                 plot.title = element_text(size = UMAPtitleText),
  #                 legend.text=element_text(size=UMAPlegendText),
  #                 legend.spacing.y = unit(0,"cm"))
  #  
  #  if (length(sampSelected) > 0) {
  #    plotdata <- plot_df[sampSelected,]
  #    plotlabel <- rownames(plot_df[sampSelected,])
  #    #if (!is.null(metaColanno)) {
  #    #  if (input$UMAPmatrixFile == T) {
  #    #    colnames(plot_df)[c(4,5)] <- c("AnnoName","GeneName")
  #    #  }
  #    #  if (input$UMAPmatrixFile == F) {
  #    #    colnames(plot_df)[4] <- c("AnnoName","GeneName")
  #    #  }
  #    #}
  #    #if (!is.null(metaColanno)) {
  #    #  colnames(plotdata)[c(4,5)] <- c("AnnoName","GeneName")
  #    #}
  #    #if (is.null(metaColanno)) {
  #    #  colnames(plotdata)[4] <- "GeneName"
  #    #}
  #    k <- k + geom_point(data = plotdata,
  #                        aes(x = UMAP1, y = UMAP2),
  #                        pch = 1,
  #                        color = "black",
  #                        size = UMAPdotSize,
  #                        stroke = .3)
  #  }
  #  k
  #  
  #})
  #
  #umap_plot_PreC_expr_react <- reactive({
  #  
  #  req(input$UMAPmetaFile)
  #  req(input$UMAPmatrixFile)
  #  ## Variables
  #  sampSelected <- input$UMAPsampSelect         # Sample Names to annotate
  #  metaColanno <- NULL
  #  
  #  plot_df <- UMAP_PreC_CoordTable_react()
  #  rownames(plot_df) <- plot_df[,1]
  #  
  #  umap1 <- input$SelectPreCalc1
  #  umap2 <- input$SelectPreCalc2
  #  
  #  if (input$UMAPannotateSamps != " ") {
  #    metaColanno <- input$UMAPannotateSamps
  #    if (input$UMAPannoContCheck != T) {
  #      plot_df[,metaColanno] <- as.factor(plot_df[,metaColanno])
  #    }
  #  }
  #  expr <- umap_matrix_loaded_og()
  #  if (!is.null(input$GeneSelection)) {
  #    metaColgene <- input$GeneSelection
  #  }
  #  if (is.null(input$GeneSelection)) {
  #    metaColgene <- rownames(expr)[1]
  #  }
  #  metaColkmean <- "Cluster"
  #  
  #  
  #  k <- umap_plot_PreC_expr_react_base()
  #  
  #  if (length(sampSelected) > 0) {
  #    plotdata <- plot_df[sampSelected,]
  #    plotlabel <- rownames(plot_df[sampSelected,])
  #    if (!is.null(metaColanno)) {
  #      colnames(plotdata)[c(4,5)] <- c("AnnoName","GeneName")
  #    }
  #    if (is.null(metaColanno)) {
  #      colnames(plotdata)[4] <- "GeneName"
  #    }
  #  }
  #  
  #  print("expr df")
  #  print(head(as_tibble(k[["data"]])))
  #  
  #  ## Make it plotly
  #  k2 <- ggplotly(k,
  #                 tooltip = "text") %>% 
  #    
  #    config(displayModeBar = F)  %>% 
  #    
  #    layout(font=list(color="#black"),
  #           xaxis=list(title=umap1,zeroline=F),
  #           yaxis=list(title=umap2,zeroline=F)) 
  #  
  #  k2 <- k2 %>%
  #    hide_legend() %>%
  #    hide_colorbar()
  #  
  #  if (length(sampSelected) > 0) {
  #    k2 <- k2 %>%
  #      add_annotations(x = plotdata$UMAP1,
  #                      y = plotdata$UMAP2,
  #                      text = rownames(plotdata),
  #                      showarrow = TRUE,
  #                      arrowhead = 4,
  #                      arrowsize = .5)
  #  }
  #  
  #  k2
  #  
  #})
  #
  ### UMAP Plot  - PATH - Clin
  #umap_plot_PreC_kmean_react_base <- reactive({
  #  
  #  req(input$UMAPmetaFile)
  #  req(input$UMAPmatrixFile)
  #  ## Variables
  #  sampSelected <- input$UMAPsampSelect         # Sample Names to annotate
  #  UMAPtitleText <- input$UMAPtitleTextSize     # Title text size
  #  UMAPaxisText <- input$UMAPaxisTextSize       # Axis text size
  #  UMAPlegendText <- input$UMAPlegendTextSize   # Legend text size
  #  UMAPdotSize <- input$UMAPdotSize             # Dot size
  #  metaColanno <- NULL
  #  
  #  umap1 <- input$SelectPreCalc1
  #  umap2 <- input$SelectPreCalc2
  #  
  #  plot_df <- UMAP_PreC_CoordTable_react()
  #  rownames(plot_df) <- plot_df[,1]
  #  if (input$UMAPannotateSamps != " ") {
  #    metaColanno <- input$UMAPannotateSamps
  #    if (input$UMAPannoContCheck != T) {
  #      plot_df[,metaColanno] <- as.factor(plot_df[,metaColanno])
  #    }
  #  }
  #  expr <- umap_matrix_loaded_og()
  #  if (!is.null(input$GeneSelection)) {
  #    metaColgene <- input$GeneSelection
  #  }
  #  if (is.null(input$GeneSelection)) {
  #    metaColgene <- rownames(expr)[1]
  #  }
  #  metaColkmean <- "Cluster"
  #  
  #  expr <- umap_matrix_loaded_og()
  #  if (!is.null(input$GeneSelection)) {
  #    metaColgene <- input$GeneSelection
  #  }
  #  if (is.null(input$GeneSelection)) {
  #    metaColgene <- rownames(expr)[1]
  #  }
  #  LogChoice <- input$LogGeneSelection
  #  expr2 <- expr
  #  if (LogChoice == TRUE) {
  #    expr2 <- log2(as.matrix(expr) + 0.00001)
  #  }
  #  exprRange <- round(range(expr2[metaColgene,]),2)
  #  
  #  if (!is.null(input$GeneExprRange)) {
  #    exprMin <- as.numeric(gsub(" ","",strsplit(as.character(input$GeneExprRange),",")[[1]][1]))
  #    exprMax <- as.numeric(gsub(" ","",strsplit(as.character(input$GeneExprRange),",")[[1]][2]))
  #  }
  #  else if (is.null(input$GeneExprRange)) {
  #    exprMin <- as.numeric(gsub(" ","",strsplit(as.character(exprRange),",")[[1]][1]))
  #    exprMax <- as.numeric(gsub(" ","",strsplit(as.character(exprRange),",")[[1]][2]))
  #  }
  #  myPalette <- colorRampPalette(rev(brewer.pal(9, input$UMAPcolors)))
  #  
  #  ## Cluster Annotation Plot
  #  if (!is.null(input$UMAPmatrixFile)) {
  #    if (input$UMAPannotateSamps != " ") {
  #      colnames(plot_df)[c(4,5)] <- c("AnnoName","GeneName")
  #      print(head(as_tibble(plot_df)))
  #      k <- plot_df %>%
  #        ggplot(aes(UMAP1, UMAP2, colour=Cluster,
  #                   text = paste("</br> <b>Sample Name:</b> ", SampleName,
  #                                "</br> <b>",metaColanno,":</b> ", AnnoName,
  #                                "</br> <b>",metaColgene," Gene Expression:</b> ", round(GeneName,4),
  #                                "</br> <b>Cluster:</b> ", Cluster,
  #                                sep = "")))
  #    }
  #    else if (input$UMAPannotateSamps == " ") {
  #      colnames(plot_df)[4] <- "GeneName"
  #      print(head(as_tibble(plot_df)))
  #      k <- plot_df %>%
  #        ggplot(aes(UMAP1, UMAP2, colour=Cluster,
  #                   text = paste("</br> <b>Sample Name:</b> ", SampleName,
  #                                "</br> <b>",metaColgene," Gene Expression:</b> ", round(GeneName,4),
  #                                "</br> <b>Cluster:</b> ", Cluster,
  #                                sep = "")))
  #    }
  #    
  #  }
  #  
  #  
  #  k <- k + geom_point(shape = 19,
  #                      size = UMAPdotSize) +
  #    
  #    theme_minimal()
  #  
  #  if (input$UMAPannotateSamps != " ") {
  #    metaColanno <- input$UMAPannotateSamps
  #    k <- k + labs(x = umap1,
  #                  y = umap2,
  #                  color = "Cluster")
  #  }
  #  if (input$UMAPannotateSamps == " ") {
  #    k <- k + labs(x = umap1,
  #                  y = umap2,)
  #  }
  #  #k <- k + scale_colour_gradientn(colours = rev(myPalette(100)), limits=c(exprMin, exprMax), oob = scales::squish)
  #  
  #  
  #  k <- k + theme(axis.text = element_text(size = UMAPaxisText),
  #                 axis.title = element_text(size = UMAPaxisText),
  #                 plot.title = element_text(size = UMAPtitleText),
  #                 legend.text=element_text(size=UMAPlegendText))
  #  
  #  if (length(sampSelected) > 0) {
  #    plotdata <- plot_df[sampSelected,]
  #    plotlabel <- rownames(plot_df[sampSelected,])
  #    #if (!is.null(metaColanno)) {
  #    #  if (input$UMAPmatrixFile == T) {
  #    #    colnames(plot_df)[c(4,5)] <- c("AnnoName","GeneName")
  #    #  }
  #    #  if (input$UMAPmatrixFile == F) {
  #    #    colnames(plot_df)[4] <- c("AnnoName","GeneName")
  #    #  }
  #    #}
  #    #if (!is.null(metaColanno)) {
  #    #  colnames(plotdata)[c(4,5)] <- c("AnnoName","GeneName")
  #    #}
  #    #if (is.null(metaColanno)) {
  #    #  colnames(plotdata)[4] <- "GeneName"
  #    #}
  #    k <- k + geom_point(data = plotdata,
  #                        aes(x = UMAP1, y = UMAP2),
  #                        pch = 1,
  #                        color = "black",
  #                        size = UMAPdotSize,
  #                        stroke = .3)
  #  }
  #  k
  #  
  #})
  #
  #umap_plot_PreC_kmean_react <- reactive({
  #  
  #  req(input$UMAPmetaFile)
  #  req(input$UMAPmatrixFile)
  #  ## Variables
  #  sampSelected <- input$UMAPsampSelect         # Sample Names to annotate
  #  metaColanno <- NULL
  #  
  #  plot_df <- UMAP_PreC_CoordTable_react()
  #  rownames(plot_df) <- plot_df[,1]
  #  
  #  umap1 <- input$SelectPreCalc1
  #  umap2 <- input$SelectPreCalc2
  #  
  #  if (input$UMAPannotateSamps != " ") {
  #    metaColanno <- input$UMAPannotateSamps
  #    if (input$UMAPannoContCheck != T) {
  #      plot_df[,metaColanno] <- as.factor(plot_df[,metaColanno])
  #    }
  #  }
  #  expr <- umap_matrix_loaded_og()
  #  if (!is.null(input$GeneSelection)) {
  #    metaColgene <- input$GeneSelection
  #  }
  #  if (is.null(input$GeneSelection)) {
  #    metaColgene <- rownames(expr)[1]
  #  }
  #  metaColkmean <- "Cluster"
  #  
  #  
  #  k <- umap_plot_PreC_kmean_react_base()
  #  
  #  if (length(sampSelected) > 0) {
  #    plotdata <- plot_df[sampSelected,]
  #    plotlabel <- rownames(plot_df[sampSelected,])
  #    if (!is.null(metaColanno)) {
  #      colnames(plotdata)[c(4,5)] <- c("AnnoName","GeneName")
  #    }
  #    if (is.null(metaColanno)) {
  #      colnames(plotdata)[4] <- "GeneName"
  #    }
  #  }
  #  print("Kmeans df")
  #  print(head(as_tibble(k[["data"]])))
  #  
  #  ## Make it plotly
  #  k2 <- ggplotly(k,
  #                 tooltip = "text") %>% 
  #    
  #    config(displayModeBar = F)  %>% 
  #    
  #    layout(font=list(color="#black"),
  #           xaxis=list(title=umap1,zeroline=F),
  #           yaxis=list(title=umap2,zeroline=F)) 
  #  
  #  k2 <- k2 %>%
  #    hide_legend()
  #  
  #  if (length(sampSelected) > 0) {
  #    k2 <- k2 %>%
  #      add_annotations(x = plotdata$UMAP1,
  #                      y = plotdata$UMAP2,
  #                      text = rownames(plotdata),
  #                      showarrow = TRUE,
  #                      arrowhead = 4,
  #                      arrowsize = .5)
  #  }
  #  
  #  k2
  #  
  #})
  #
  #output$UMAPprecLegend <- renderPlot({
  #  
  #  req(input$UMAPmetaFile)
  #  if (is.null(input$UMAPmatrixFile)) {
  #    if (input$UMAPannotateSamps != " ") {
  #      k <- umap_plot_PreC_clin_react_base()
  #      legendk <- g_legend(k)
  #      legend_grid <- legendk
  #    }
  #    else if (input$UMAPannotateSamps == " ") {
  #      legend_grid <- NULL
  #    }
  #  }
  #  if (!is.null(input$UMAPmatrixFile)) {
  #    if (input$UMAPannotateSamps != " ") {
  #      k <- umap_plot_PreC_clin_react_base()
  #      m <- umap_plot_PreC_expr_react_base()
  #      n <- umap_plot_PreC_kmean_react_base()
  #      legendk <- g_legend(k)
  #      legendm <- g_legend(m)
  #      legendn <- g_legend(n)
  #      legend_grid <- grid.arrange(legendk,legendm,legendn,nrow = 1)
  #    }
  #    else if (input$UMAPannotateSamps == " ") {
  #      m <- umap_plot_PreC_expr_react_base()
  #      n <- umap_plot_PreC_kmean_react_base()
  #      legendm <- g_legend(m)
  #      legendn <- g_legend(n)
  #      legend_grid <- grid.arrange(legendm,legendn,nrow = 1)
  #    }
  #  }
  #  
  #  legend_grid
  #  
  #})
  #
  #UMAPplot_PreC_ALL_react <- reactive({
  #  
  #  req(input$UMAPmetaFile)
  #  
  #  if (!is.null(input$UMAPmatrixFile)) {
  #    k2 <- umap_plot_PreC_clin_react()
  #    m2 <- umap_plot_PreC_expr_react()
  #    n2 <- umap_plot_PreC_kmean_react()
  #    
  #    if (input$UMAPannotateSamps != " ") {
  #      UMAPtitlek <- paste("Annotated by",input$UMAPannotateSamps)
  #    }
  #    else if (input$UMAPannotateSamps == " ") {
  #      UMAPtitlek <- ""
  #    }
  #    
  #    if (is.null(input$UMAPmetaFile) == T) {
  #      UMAPtitlek <- ""
  #    }
  #    expr <- umap_matrix_loaded_og()
  #    if (!is.null(input$GeneSelection)) {
  #      GeneSelec <- input$GeneSelection
  #    }
  #    if (is.null(input$GeneSelection)) {
  #      GeneSelec <- rownames(expr)[1]
  #    }
  #    
  #    if (input$LogGeneSelection == T) {
  #      GeneSelecExpr <- paste(GeneSelec,"Expression (Log2)")
  #    }
  #    if (input$LogGeneSelection == F) {
  #      GeneSelecExpr <- paste(GeneSelec,"Expression")
  #    }
  #    
  #    
  #    UMAPtitlem <- paste("Annotated by",GeneSelecExpr)
  #    UMAPtitlen <- paste("Annotated by",input$ClusterMethod,input$ClusterNumber,"Clusters")
  #    
  #    
  #    if (input$UMAPorientation == "Side-by-Side") {
  #      subplot_all_sbs <- subplot(k2, m2, n2, nrows = 1, shareY = TRUE, shareX = TRUE)
  #      plot_titles = list( 
  #        list( 
  #          x = 0,
  #          y = 1,
  #          text = UMAPtitlek,  
  #          xref = "x",  
  #          yref = "paper",  
  #          xanchor = "center",  
  #          yanchor = "bottom",  
  #          showarrow = FALSE 
  #        ),  
  #        list( 
  #          x = 0,
  #          y = 1,
  #          text = UMAPtitlem,  
  #          xref = "x2",
  #          yref = "paper",
  #          xanchor = "center",  
  #          yanchor = "bottom",    
  #          showarrow = FALSE 
  #        ),  
  #        list( 
  #          x = 0,
  #          y = 1,
  #          text = UMAPtitlen,  
  #          xref = "x3",
  #          yref = "paper",
  #          xanchor = "center",  
  #          yanchor = "bottom",  
  #          showarrow = FALSE 
  #        ))
  #    }
  #    else if (input$UMAPorientation == "Stacked") {
  #      subplot_all_sbs <- subplot(k2, m2, n2, nrows = 3, shareY = TRUE, shareX = TRUE, margin = 0.04)
  #      plot_titles = list( 
  #        list( 
  #          x = 0,
  #          y = 1,
  #          text = UMAPtitlek,  
  #          xref = "x",  
  #          yref = "paper",  
  #          xanchor = "center",  
  #          yanchor = "bottom",  
  #          showarrow = FALSE 
  #        ),  
  #        list( 
  #          x = 0,
  #          y = 0.65,
  #          text = UMAPtitlem,  
  #          xref = "x",
  #          yref = "paper",
  #          xanchor = "center",  
  #          yanchor = "bottom",    
  #          showarrow = FALSE 
  #        ),  
  #        list( 
  #          x = 0,
  #          y = 0.3,
  #          text = UMAPtitlen,  
  #          xref = "x",
  #          yref = "paper",
  #          xanchor = "center",  
  #          yanchor = "bottom",  
  #          showarrow = FALSE 
  #        ))
  #    }
  #    
  #    subplot_all_sbs <- subplot_all_sbs %>% layout(annotations = plot_titles)
  #  }
  #  
  #  if (is.null(input$UMAPmatrixFile)) {
  #    k2 <- umap_plot_PreC_clin_react()
  #    if (input$UMAPannotateSamps != " ") {
  #      UMAPtitlek <- paste("Annotated by",input$UMAPannotateSamps)
  #    }
  #    else if (input$UMAPannotateSamps == " ") {
  #      UMAPtitlek <- ""
  #    }
  #    
  #    if (is.null(input$UMAPmetaFile) == T) {
  #      UMAPtitlek <- ""
  #    }
  #    k2 <- k2 %>%
  #      layout(title = UMAPtitlek)
  #    subplot_all_sbs <- k2
  #  }
  #  
  #  subplot_all_sbs
  #})
  #
  #output$UMAPplot_PreC_ALL <- renderPlotly({
  #  req(input$UMAPmetaFile)
  #  p_all <- UMAPplot_PreC_ALL_react()
  #  p_all
  #})
  
  
  ####----Downloads----####
  
  ####----Jaccard----####
  
  output$dnldJaccHeatmapIntraPDF <- downloadHandler(
    filename = function() {
      gs.u <- input$UserPathwayFile
      TopFeat <- input$TopFeatureSelect
      file_name <- gs.u$name
      file_name <- tools::file_path_sans_ext(gs.u)
      paste(file_name,"_Top",TopFeat,"_JaccardHeatmap.pdf", sep = '')
    },
    content = function(file) {
      heat <- JaccHeatmapIntra_react()
      ggsave(file,heat, width = 10, height = 12)
    }
  )
  output$dnldJaccHeatmapIntraSVG <- downloadHandler(
    filename = function() {
      gs.u <- input$UserPathwayFile
      TopFeat <- input$TopFeatureSelect
      file_name <- gs.u$name
      file_name <- tools::file_path_sans_ext(gs.u)
      paste(file_name,"_Top",TopFeat,"_JaccardHeatmap.svg", sep = '')
    },
    content = function(file) {
      heat <- JaccHeatmapIntra_react()
      ggsave(file,heat, width = 10, height = 12)
    }
  )
  output$dnldJaccDendoPDF <- downloadHandler(
    filename = function() {
      gs.u <- input$UserPathwayFile
      TopFeat <- input$TopFeatureSelect
      file_name <- gs.u$name
      file_name <- tools::file_path_sans_ext(gs.u)
      paste(file_name,"_Top",TopFeat,"_JaccardCluster.pdf", sep = '')
    },
    content = function(file) {
      dendo <- JaccDendoIntra_react()
      ggsave(file, dendo, width = 10, height = 8)
    }
  )
  output$dnldJaccDendoSVG <- downloadHandler(
    filename = function() {
      gs.u <- input$UserPathwayFile
      TopFeat <- input$TopFeatureSelect
      file_name <- gs.u$name
      file_name <- tools::file_path_sans_ext(gs.u)
      paste(file_name,"_Top",TopFeat,"_JaccardCluster.svg", sep = '')
    },
    content = function(file) {
      dendo <- JaccDendoIntra_react()
      ggsave(file, dendo, width = 10, height = 8)
    }
  )
  ## Download Cluster Table - INTRA
  output$dnldClustTabIntra <- downloadHandler(
    filename = function() {
      
      gs.u <- input$UserPathwayFile
      TopFeat <- input$TopFeatureSelect
      file_name <- gs.u$name
      file_name <- tools::file_path_sans_ext(gs.u)
      paste(file_name,"_Top",TopFeat,"_Cluster_Result.txt", sep = '')
      
    },
    content = function(file) {
      
      ## Data matrix
      jacc_df <- jacc_react_Intra()
      jacc_mat <- as.matrix(jacc_df)
      clust_me <- input$ClustMethodIntra
      
      hm_res <- pheatmap::pheatmap(jacc_mat,
                                   clustering_method = clust_me,
                                   silent = T)
      
      jacc_df <- jacc_react_Intra()
      cut_k <- input$NumClusters
      clustTab <- data.frame(Pathways = rownames(jacc_df),
                             Cluster = cutree(hm_res$tree_row,k = cut_k))
      write_delim(clustTab,file,delim = '\t')
      
    }
  )
  output$dnldSIFTabIntra <- downloadHandler(
    filename = function() {
      
      gs.u <- input$UserPathwayFile
      TopFeat <- input$TopFeatureSelect
      file_name <- gs.u$name
      file_name <- tools::file_path_sans_ext(gs.u)
      paste(file_name,"_Top",TopFeat,"_sif_dend.sif", sep = '')
      
    },
    content = function(file) {
      
      jacc_df <- jacc_react_Intra()
      dist_cutoff <- input$JaccDistCutoff
      jacc_df <- as.data.frame(jacc_df)
      jacc_df_melt <- melt(jacc_df)
      jacc_cols <- colnames(jacc_df)
      jacc_df_melt$Pathway_B <- rep_len(jacc_cols,length.out = nrow(jacc_df_melt))
      colnames(jacc_df_melt)[c(1,2)] <- c("Pathway_A","Jaccard_Distance")
      jacc_df_melt <- jacc_df_melt[which(jacc_df_melt$Jaccard_Distance <= dist_cutoff),]
      jacc_df_melt <- jacc_df_melt[which(jacc_df_melt$Pathway_A != jacc_df_melt$Pathway_B),]
      jacc_df_melt <- jacc_df_melt[!duplicated(t(apply(jacc_df_melt,1,sort))),]
      write_delim(jacc_df_melt,file,delim = '\t')
      
      
    }
  )
  
  ## Download Jaccard output - Intra
  output$dnldJaccTabIntra <- downloadHandler(
    filename = function() {
      
      gs.u <- input$UserPathwayFile
      TopFeat <- input$TopFeatureSelect
      file_name <- gs.u$name
      file_name <- tools::file_path_sans_ext(gs.u)
      paste(file_name,"_Top",TopFeat,"_Jaccard_Connectivity.txt", sep = '')
      
    },
    content = function(file) {
      
      jacc_df <- jacc_react_Intra()
      jacc_df <- as.data.frame(jacc_df)
      jacc_df$Pathways <- rownames(jacc_df)
      jacc_df <- jacc_df %>%
        relocate(Pathways)
      write_delim(jacc_df,file,delim = '\t')
      
    }
  )
  
  ## Download Jaccard output - Intra
  output$dnldClusterTabAnno <- downloadHandler(
    filename = function() {
      
      gs.u <- input$UserPathwayFile
      TopFeat <- input$TopFeatureSelect
      file_name <- gs.u$name
      file_name <- tools::file_path_sans_ext(gs.u)
      paste(file_name,"_Top",TopFeat,"_Cluster_Annotation.txt", sep = '')
      
    },
    content = function(file) {
      
      df <- ClusterTabAnno_react()
      write_delim(df,file,delim = '\t')
      
    }
  )
  
  ####----Matrix Clustering----####
  
  output$dnldUserMatrixHeatmapPDF <- downloadHandler(
    filename = function() {
      gs.u <- input$UserMatrixFile
      file_name <- gs.u$name
      file_name <- tools::file_path_sans_ext(gs.u)
      paste(file_name,"_Heatmap.pdf", sep = '')
    },
    content = function(file) {
      heat <- UserMatrixHeatmap_react()
      ggsave(file, heat, width = 10, height = 12)
    }
  )
  output$dnldUserMatrixHeatmapSVG <- downloadHandler(
    filename = function() {
      gs.u <- input$UserMatrixFile
      file_name <- gs.u$name
      file_name <- tools::file_path_sans_ext(gs.u)
      paste(file_name,"_Heatmap.svg", sep = '')
    },
    content = function(file) {
      heat <- UserMatrixHeatmap_react()
      ggsave(file, heat, width = 10, height = 12)
    }
  )
  output$dnldUserMatixCLusteringPDF <- downloadHandler(
    filename = function() {
      gs.u <- input$UserMatrixFile
      file_name <- gs.u$name
      file_name <- tools::file_path_sans_ext(gs.u)
      paste(file_name,"_ClusterFigure.pdf", sep = '')
    },
    content = function(file) {
      dendo <- JaccDendoIntra_react()
      ggsave(file, dendo, width = 10, height = 8)
    }
  )
  output$dnldUserMatixCLusteringSVG <- downloadHandler(
    filename = function() {
      gs.u <- input$UserMatrixFile
      file_name <- gs.u$name
      file_name <- tools::file_path_sans_ext(gs.u)
      paste(file_name,"_ClusterFigure.svg", sep = '')
    },
    content = function(file) {
      dendo <- JaccDendoIntra_react()
      ggsave(file, dendo, width = 10, height = 8)
    }
  )
  ## Download Cluster Table - INTRA
  output$dnldClustTabMat <- downloadHandler(
    filename = function() {
      
      gs.u <- input$UserMatrixFile
      file_name <- gs.u$name
      file_name <- tools::file_path_sans_ext(gs.u)
      paste(file_name,"_Cluster_Result.txt", sep = '')
      
    },
    content = function(file) {
      
      ## Data matrix
      df <- user_matrix_loaded()
      rownames(df) <- df[,1]
      df <- df[,-1]
      df_mat <- as.matrix(df)
      clust_me <- input$ClustMethodMat
      
      hm_res <- pheatmap::pheatmap(df_mat,
                                   clustering_method = clust_me,
                                   silent = T)
      
      cut_k <- input$NumClustersMat
      clustTab <- data.frame(Pathways = rownames(df),
                             Cluster = cutree(hm_res$tree_row,k = cut_k))
      write_delim(clustTab,file,delim = '\t')
      
    }
  )
  
  output$dnldSIFTabMat <- downloadHandler(
    filename = function() {
      
      gs.u <- input$UserMatrixFile
      file_name <- gs.u$name
      file_name <- tools::file_path_sans_ext(gs.u)
      paste(file_name,"_sif_dend.sif", sep = '')
      
    },
    content = function(file) {
      
      df <- user_matrix_loaded()
      rownames(df) <- df[,1]
      df <- df[,-1]
      df <- as.data.frame(df)
      df_melt <- melt(df)
      df_cols <- colnames(df)
      df_melt$Sample_B <- rep_len(df_cols,length.out = nrow(df_melt))
      colnames(df_melt)[c(1,2)] <- c("Sample_A","Distance")
      df_melt <- df_melt[which(df_melt$Sample_A != df_melt$Sample_B),]
      df_melt <- df_melt[!duplicated(t(apply(df_melt,1,sort))),]
      write_delim(df_melt,file,delim = '\t')
      
      
    }
  )
  
  #####----UMAP----####
  #
  ### Cluster Table
  #
  #output$dnldUMAPclusterTab <- downloadHandler(
  #  filename = function() {
  #    gs.u <- input$UMAPmatrixFile
  #    n_clust <- input$ClusterNumber
  #    cluster_method <- input$ClusterMethod
  #    file_name <- gs.u$name
  #    file_name <- tools::file_path_sans_ext(gs.u)
  #    paste(file_name,"_",cluster_method,"_",n_clust,"_Clusters.txt", sep = '')
  #  },
  #  content = function(file) {
  #    tab <- umapClusterTable_react()
  #    write_delim(tab,file,delim = '\t')
  #  }
  #)
  #
  ### MVG Downloads
  #
  #output$dnldUMAPcoords_MVG <- downloadHandler(
  #  filename = function() {
  #    gs.u <- input$UMAPmatrixFile
  #    file_name <- gs.u$name
  #    file_name <- tools::file_path_sans_ext(gs.u)
  #    topNum <- input$TopNumMVG
  #    paste(file_name,"_",topNum,"MVG_UMAP_Coordinates.txt", sep = '')
  #  },
  #  content = function(file) {
  #    tdata_fit_df <- as.data.frame(UMAP_MVG_CoordTable_react())
  #    write_delim(tdata_fit_df,file,delim = '\t')
  #  }
  #)
  #
  #output$dnldUMAP_SVG_MVG_clin <- downloadHandler(
  #  filename = function() {
  #    gs.u <- input$UMAPmatrixFile
  #    file_name <- gs.u$name
  #    file_name <- tools::file_path_sans_ext(gs.u)
  #    topNum <- input$TopNumMVG
  #    metaCol <- input$UMAPannotateSamps
  #    metaCol <- gsub(" ","",input$UMAPannotateSamps)
  #    metaCol <- gsub("[[:punct:]]","_",metaCol)
  #    paste(file_name,"_",topNum,"MVG_",metaCol,"Anno_UMAP.svg", sep = '')
  #  },
  #  content = function(file) {
  #    plot <- umap_plot_MVG_react_clin_base()
  #    ggsave(file,plot, width = 8, height = 8)
  #  }
  #)
  #output$dnldUMAP_PDF_MVG_clin <- downloadHandler(
  #  filename = function() {
  #    gs.u <- input$UMAPmatrixFile
  #    file_name <- gs.u$name
  #    file_name <- tools::file_path_sans_ext(gs.u)
  #    topNum <- input$TopNumMVG
  #    metaCol <- input$UMAPannotateSamps
  #    metaCol <- gsub(" ","",input$UMAPannotateSamps)
  #    metaCol <- gsub("[[:punct:]]","_",metaCol)
  #    paste(file_name,"_",topNum,"MVG_",metaCol,"Anno_UMAP.pdf", sep = '')
  #  },
  #  content = function(file) {
  #    plot <- umap_plot_MVG_react_clin_base()
  #    ggsave(file,plot, width = 8, height = 8)
  #  }
  #)
  #
  #output$dnldUMAP_SVG_MVG_expr <- downloadHandler(
  #  filename = function() {
  #    gs.u <- input$UMAPmatrixFile
  #    file_name <- gs.u$name
  #    file_name <- tools::file_path_sans_ext(gs.u)
  #    topNum <- input$TopNumMVG
  #    metaCol <- paste(input$GeneSelection,"Expr",sep = "")
  #    paste(file_name,"_",topNum,"MVG_",metaCol,"Anno_UMAP.svg", sep = '')
  #  },
  #  content = function(file) {
  #    plot <- umap_plot_MVG_react_expr_base()
  #    ggsave(file,plot, width = 8, height = 8)
  #  }
  #)
  #output$dnldUMAP_PDF_MVG_expr <- downloadHandler(
  #  filename = function() {
  #    gs.u <- input$UMAPmatrixFile
  #    file_name <- gs.u$name
  #    file_name <- tools::file_path_sans_ext(gs.u)
  #    topNum <- input$TopNumMVG
  #    metaCol <- paste(input$GeneSelection,"Expr",sep = "")
  #    paste(file_name,"_",topNum,"MVG_",metaCol,"Anno_UMAP.pdf", sep = '')
  #  },
  #  content = function(file) {
  #    plot <- umap_plot_MVG_react_expr_base()
  #    ggsave(file,plot, width = 8, height = 8)
  #  }
  #)
  #
  #output$dnldUMAP_SVG_MVG_kmean <- downloadHandler(
  #  filename = function() {
  #    gs.u <- input$UMAPmatrixFile
  #    file_name <- gs.u$name
  #    file_name <- tools::file_path_sans_ext(gs.u)
  #    topNum <- input$TopNumMVG
  #    clusterMethod <- input$ClusterMethod
  #    ClusterNum <- input$ClusterNumber
  #    metaCol <- paste(clusterMethod,ClusterNum,"Clusters",sep = "")
  #    paste(file_name,"_",topNum,"MVG_",metaCol,"Anno_UMAP.svg", sep = '')
  #  },
  #  content = function(file) {
  #    plot <- umap_plot_MVG_react_kmean_base()
  #    ggsave(file,plot, width = 8, height = 8)
  #  }
  #)
  #output$dnldUMAP_PDF_MVG_kmean <- downloadHandler(
  #  filename = function() {
  #    gs.u <- input$UMAPmatrixFile
  #    file_name <- gs.u$name
  #    file_name <- tools::file_path_sans_ext(gs.u)
  #    topNum <- input$TopNumMVG
  #    clusterMethod <- input$ClusterMethod
  #    ClusterNum <- input$ClusterNumber
  #    metaCol <- paste(clusterMethod,ClusterNum,"Clusters",sep = "")
  #    paste(file_name,"_",topNum,"MVG_",metaCol,"Anno_UMAP.pdf", sep = '')
  #  },
  #  content = function(file) {
  #    plot <- umap_plot_MVG_react_kmean_base()
  #    ggsave(file,plot, width = 8, height = 8)
  #  }
  #)
  #
  ### Pathway downloads
  #
  #output$dnldUMAPcoords_PATH <- downloadHandler(
  #  filename = function() {
  #    gs.u <- input$UMAPmatrixFile
  #    file_name <- gs.u$name
  #    file_name <- tools::file_path_sans_ext(gs.u)
  #    rowsSelected <- input$PathSelectTable_rows_selected
  #    GeneSet <- gs_cat2[rowsSelected,3]
  #    paste(file_name,"_",GeneSet,"_UMAP_Coordinates.txt", sep = '')
  #  },
  #  content = function(file) {
  #    tdata_fit_df <- as.data.frame(UMAP_PATH_CoordTable_react())
  #    write_delim(tdata_fit_df,file,delim = '\t')
  #  }
  #)
  #
  #output$dnldUMAP_PATH_SVG_clin <- downloadHandler(
  #  filename = function() {
  #    gs.u <- input$UMAPmatrixFile
  #    file_name <- gs.u$name
  #    file_name <- tools::file_path_sans_ext(gs.u)
  #    rowsSelected <- input$PathSelectTable_rows_selected
  #    GeneSet <- gs_cat2[rowsSelected,3]
  #    metaCol <- input$UMAPannotateSamps
  #    metaCol <- gsub(" ","",input$UMAPannotateSamps)
  #    metaCol <- gsub("[[:punct:]]","_",metaCol)
  #    paste(file_name,"_",GeneSet,"_",metaCol,"Anno_UMAP.svg", sep = '')
  #  },
  #  content = function(file) {
  #    plot <- umap_plot_PATH_react_clin_base()
  #    ggsave(file,plot, width = 8, height = 8)
  #  }
  #)
  #output$dnldUMAP_PDF_PATH_clin <- downloadHandler(
  #  filename = function() {
  #    gs.u <- input$UMAPmatrixFile
  #    file_name <- gs.u$name
  #    file_name <- tools::file_path_sans_ext(gs.u)
  #    rowsSelected <- input$PathSelectTable_rows_selected
  #    GeneSet <- gs_cat2[rowsSelected,3]
  #    metaCol <- input$UMAPannotateSamps
  #    metaCol <- gsub(" ","",input$UMAPannotateSamps)
  #    metaCol <- gsub("[[:punct:]]","_",metaCol)
  #    paste(file_name,"_",GeneSet,"_",metaCol,"Anno_UMAP.pdf", sep = '')
  #  },
  #  content = function(file) {
  #    plot <- umap_plot_PATH_react_clin_base()
  #    ggsave(file,plot, width = 8, height = 8)
  #  }
  #)
  #
  #output$dnldUMAP_PATH_SVG_expr <- downloadHandler(
  #  filename = function() {
  #    gs.u <- input$UMAPmatrixFile
  #    file_name <- gs.u$name
  #    file_name <- tools::file_path_sans_ext(gs.u)
  #    rowsSelected <- input$PathSelectTable_rows_selected
  #    GeneSet <- gs_cat2[rowsSelected,3]
  #    metaCol <- paste(input$GeneSelection,"Expr",sep = "")
  #    paste(file_name,"_",GeneSet,"_",metaCol,"Anno_UMAP.svg", sep = '')
  #  },
  #  content = function(file) {
  #    plot <- umap_plot_PATH_react_expr_base()
  #    ggsave(file,plot, width = 8, height = 8)
  #  }
  #)
  #output$dnldUMAP_PDF_PATH_expr <- downloadHandler(
  #  filename = function() {
  #    gs.u <- input$UMAPmatrixFile
  #    file_name <- gs.u$name
  #    file_name <- tools::file_path_sans_ext(gs.u)
  #    rowsSelected <- input$PathSelectTable_rows_selected
  #    GeneSet <- gs_cat2[rowsSelected,3]
  #    metaCol <- paste(input$GeneSelection,"Expr",sep = "")
  #    paste(file_name,"_",GeneSet,"_",metaCol,"Anno_UMAP.pdf", sep = '')
  #  },
  #  content = function(file) {
  #    plot <- umap_plot_PATH_react_expr_base()
  #    ggsave(file,plot, width = 8, height = 8)
  #  }
  #)
  #
  #output$dnldUMAP_PATH_SVG_kmean <- downloadHandler(
  #  filename = function() {
  #    gs.u <- input$UMAPmatrixFile
  #    file_name <- gs.u$name
  #    file_name <- tools::file_path_sans_ext(gs.u)
  #    rowsSelected <- input$PathSelectTable_rows_selected
  #    GeneSet <- gs_cat2[rowsSelected,3]
  #    clusterMethod <- input$ClusterMethod
  #    ClusterNum <- input$ClusterNumber
  #    metaCol <- paste(clusterMethod,ClusterNum,"Clusters",sep = "")
  #    paste(file_name,"_",GeneSet,"_",metaCol,"Anno_UMAP.svg", sep = '')
  #  },
  #  content = function(file) {
  #    plot <- umap_plot_PATH_react_kmean_base()
  #    ggsave(file,plot, width = 8, height = 8)
  #  }
  #)
  #output$dnldUMAP_PDF_PATH_kmean <- downloadHandler(
  #  filename = function() {
  #    gs.u <- input$UMAPmatrixFile
  #    file_name <- gs.u$name
  #    file_name <- tools::file_path_sans_ext(gs.u)
  #    rowsSelected <- input$PathSelectTable_rows_selected
  #    GeneSet <- gs_cat2[rowsSelected,3]
  #    clusterMethod <- input$ClusterMethod
  #    ClusterNum <- input$ClusterNumber
  #    metaCol <- paste(clusterMethod,ClusterNum,"Clusters",sep = "")
  #    paste(file_name,"_",GeneSet,"_",metaCol,"Anno_UMAP.pdf", sep = '')
  #  },
  #  content = function(file) {
  #    plot <- umap_plot_PATH_react_kmean_base()
  #    ggsave(file,plot, width = 8, height = 8)
  #  }
  #)
  #
  ### Base downloads
  #
  #output$dnldUMAPcoords <- downloadHandler(
  #  filename = function() {
  #    gs.u <- input$UMAPmatrixFile
  #    file_name <- gs.u$name
  #    file_name <- tools::file_path_sans_ext(gs.u)
  #    paste(file_name,"_AllGenes_UMAP_Coordinates.txt", sep = '')
  #  },
  #  content = function(file) {
  #    tdata_fit_df <- as.data.frame(UMAP_CoordTable_react())
  #    write_delim(tdata_fit_df,file,delim = '\t')
  #  }
  #)
  #
  #output$dnldUMAP_SVG_clin <- downloadHandler(
  #  filename = function() {
  #    gs.u <- input$UMAPmatrixFile
  #    file_name <- gs.u$name
  #    file_name <- tools::file_path_sans_ext(gs.u)
  #    rowsSelected <- input$PathSelectTable_rows_selected
  #    GeneSet <- gs_cat2[rowsSelected,3]
  #    metaCol <- input$UMAPannotateSamps
  #    metaCol <- gsub(" ","",input$UMAPannotateSamps)
  #    metaCol <- gsub("[[:punct:]]","_",metaCol)
  #    paste(file_name,"_AllGenes_",metaCol,"Anno_UMAP.svg", sep = '')
  #  },
  #  content = function(file) {
  #    plot <- umap_plot_clin_react_base()
  #    ggsave(file,plot, width = 8, height = 8)
  #  }
  #)
  #output$dnldUMAP_PDF_clin <- downloadHandler(
  #  filename = function() {
  #    gs.u <- input$UMAPmatrixFile
  #    file_name <- gs.u$name
  #    file_name <- tools::file_path_sans_ext(gs.u)
  #    rowsSelected <- input$PathSelectTable_rows_selected
  #    GeneSet <- gs_cat2[rowsSelected,3]
  #    metaCol <- input$UMAPannotateSamps
  #    metaCol <- gsub(" ","",input$UMAPannotateSamps)
  #    metaCol <- gsub("[[:punct:]]","_",metaCol)
  #    paste(file_name,"_AllGenes_",metaCol,"Anno_UMAP.pdf", sep = '')
  #  },
  #  content = function(file) {
  #    plot <- umap_plot_clin_react_base()
  #    ggsave(file,plot, width = 8, height = 8)
  #  }
  #)
  #
  #output$dnldUMAP_SVG_expr <- downloadHandler(
  #  filename = function() {
  #    gs.u <- input$UMAPmatrixFile
  #    file_name <- gs.u$name
  #    file_name <- tools::file_path_sans_ext(gs.u)
  #    rowsSelected <- input$PathSelectTable_rows_selected
  #    GeneSet <- gs_cat2[rowsSelected,3]
  #    metaCol <- paste(input$GeneSelection,"Expr",sep = "")
  #    paste(file_name,"_AllGenes_",metaCol,"Anno_UMAP.svg", sep = '')
  #  },
  #  content = function(file) {
  #    plot <- umap_plot_expr_react_base()
  #    ggsave(file,plot, width = 8, height = 8)
  #  }
  #)
  #output$dnldUMAP_PDF_expr <- downloadHandler(
  #  filename = function() {
  #    gs.u <- input$UMAPmatrixFile
  #    file_name <- gs.u$name
  #    file_name <- tools::file_path_sans_ext(gs.u)
  #    rowsSelected <- input$PathSelectTable_rows_selected
  #    GeneSet <- gs_cat2[rowsSelected,3]
  #    metaCol <- paste(input$GeneSelection,"Expr",sep = "")
  #    paste(file_name,"_AllGenes_",metaCol,"Anno_UMAP.pdf", sep = '')
  #  },
  #  content = function(file) {
  #    plot <- umap_plot_expr_react_base()
  #    ggsave(file,plot, width = 8, height = 8)
  #  }
  #)
  #
  #output$dnldUMAP_SVG_kmean <- downloadHandler(
  #  filename = function() {
  #    gs.u <- input$UMAPmatrixFile
  #    file_name <- gs.u$name
  #    file_name <- tools::file_path_sans_ext(gs.u)
  #    rowsSelected <- input$PathSelectTable_rows_selected
  #    GeneSet <- gs_cat2[rowsSelected,3]
  #    clusterMethod <- input$ClusterMethod
  #    ClusterNum <- input$ClusterNumber
  #    metaCol <- paste(clusterMethod,ClusterNum,"Clusters",sep = "")
  #    paste(file_name,"_AllGenes_",metaCol,"Anno_UMAP.svg", sep = '')
  #  },
  #  content = function(file) {
  #    plot <- umap_plot_kmean_react_base()
  #    ggsave(file,plot, width = 8, height = 8)
  #  }
  #)
  #output$dnldUMAP_PDF_kmean <- downloadHandler(
  #  filename = function() {
  #    gs.u <- input$UMAPmatrixFile
  #    file_name <- gs.u$name
  #    file_name <- tools::file_path_sans_ext(gs.u)
  #    rowsSelected <- input$PathSelectTable_rows_selected
  #    GeneSet <- gs_cat2[rowsSelected,3]
  #    clusterMethod <- input$ClusterMethod
  #    ClusterNum <- input$ClusterNumber
  #    metaCol <- paste(clusterMethod,ClusterNum,"Clusters",sep = "")
  #    paste(file_name,"_AllGenes_",metaCol,"Anno_UMAP.pdf", sep = '')
  #  },
  #  content = function(file) {
  #    plot <- umap_plot_kmean_react_base()
  #    ggsave(file,plot, width = 8, height = 8)
  #  }
  #)
  #
  ### Pre-calculated downloads
  #
  #output$dnldUMAP_SVG_PreC_clin <- downloadHandler(
  #  filename = function() {
  #    gs.u <- input$UMAPmatrixFile
  #    file_name <- gs.u$name
  #    file_name <- tools::file_path_sans_ext(gs.u)
  #    umap1 <- input$SelectPreCalc1
  #    umap2 <- input$SelectPreCalc2
  #    rowsSelected <- input$PathSelectTable_rows_selected
  #    GeneSet <- gs_cat2[rowsSelected,3]
  #    metaCol <- input$UMAPannotateSamps
  #    metaCol <- gsub(" ","",input$UMAPannotateSamps)
  #    metaCol <- gsub("[[:punct:]]","_",metaCol)
  #    paste(file_name,"_",umap1,"_",umap2,"_",metaCol,"Anno_UMAP.svg", sep = '')
  #  },
  #  content = function(file) {
  #    plot <- umap_plot_PreC_clin_react_base()
  #    ggsave(file,plot, width = 8, height = 8)
  #  }
  #)
  #output$dnldUMAP_PDF_PreC_clin <- downloadHandler(
  #  filename = function() {
  #    gs.u <- input$UMAPmatrixFile
  #    file_name <- gs.u$name
  #    file_name <- tools::file_path_sans_ext(gs.u)
  #    umap1 <- input$SelectPreCalc1
  #    umap2 <- input$SelectPreCalc2
  #    rowsSelected <- input$PathSelectTable_rows_selected
  #    GeneSet <- gs_cat2[rowsSelected,3]
  #    metaCol <- input$UMAPannotateSamps
  #    metaCol <- gsub(" ","",input$UMAPannotateSamps)
  #    metaCol <- gsub("[[:punct:]]","_",metaCol)
  #    paste(file_name,"_",umap1,"_",umap2,"_",metaCol,"Anno_UMAP.pdf", sep = '')
  #  },
  #  content = function(file) {
  #    plot <- umap_plot_PreC_clin_react_base()
  #    ggsave(file,plot, width = 8, height = 8)
  #  }
  #)
  #
  #output$dnldUMAP_SVG_PreC_expr <- downloadHandler(
  #  filename = function() {
  #    gs.u <- input$UMAPmatrixFile
  #    file_name <- gs.u$name
  #    file_name <- tools::file_path_sans_ext(gs.u)
  #    umap1 <- input$SelectPreCalc1
  #    umap2 <- input$SelectPreCalc2
  #    rowsSelected <- input$PathSelectTable_rows_selected
  #    GeneSet <- gs_cat2[rowsSelected,3]
  #    metaCol <- paste(input$GeneSelection,"Expr",sep = "")
  #    paste(file_name,"_",umap1,"_",umap2,"_",metaCol,"Anno_UMAP.svg", sep = '')
  #  },
  #  content = function(file) {
  #    plot <- umap_plot_PreC_expr_react_base()
  #    ggsave(file,plot, width = 8, height = 8)
  #  }
  #)
  #output$dnldUMAP_PDF_PreC_expr <- downloadHandler(
  #  filename = function() {
  #    gs.u <- input$UMAPmatrixFile
  #    file_name <- gs.u$name
  #    file_name <- tools::file_path_sans_ext(gs.u)
  #    umap1 <- input$SelectPreCalc1
  #    umap2 <- input$SelectPreCalc2
  #    rowsSelected <- input$PathSelectTable_rows_selected
  #    GeneSet <- gs_cat2[rowsSelected,3]
  #    metaCol <- paste(input$GeneSelection,"Expr",sep = "")
  #    paste(file_name,"_",umap1,"_",umap2,"_",metaCol,"Anno_UMAP.pdf", sep = '')
  #  },
  #  content = function(file) {
  #    plot <- umap_plot_PreC_expr_react_base()
  #    ggsave(file,plot, width = 8, height = 8)
  #  }
  #)
  #
  #output$dnldUMAP_SVG_PreC_kmean <- downloadHandler(
  #  filename = function() {
  #    gs.u <- input$UMAPmatrixFile
  #    file_name <- gs.u$name
  #    file_name <- tools::file_path_sans_ext(gs.u)
  #    umap1 <- input$SelectPreCalc1
  #    umap2 <- input$SelectPreCalc2
  #    rowsSelected <- input$PathSelectTable_rows_selected
  #    GeneSet <- gs_cat2[rowsSelected,3]
  #    clusterMethod <- input$ClusterMethod
  #    ClusterNum <- input$ClusterNumber
  #    metaCol <- paste(clusterMethod,ClusterNum,"Clusters",sep = "")
  #    paste(file_name,"_",umap1,"_",umap2,"_",metaCol,"Anno_UMAP.svg", sep = '')
  #  },
  #  content = function(file) {
  #    plot <- umap_plot_PreC_kmean_react_base()
  #    ggsave(file,plot, width = 8, height = 8)
  #  }
  #)
  #output$dnldUMAP_PDF_PreC_kmean <- downloadHandler(
  #  filename = function() {
  #    gs.u <- input$UMAPmatrixFile
  #    file_name <- gs.u$name
  #    file_name <- tools::file_path_sans_ext(gs.u)
  #    umap1 <- input$SelectPreCalc1
  #    umap2 <- input$SelectPreCalc2
  #    rowsSelected <- input$PathSelectTable_rows_selected
  #    GeneSet <- gs_cat2[rowsSelected,3]
  #    clusterMethod <- input$ClusterMethod
  #    ClusterNum <- input$ClusterNumber
  #    metaCol <- paste(clusterMethod,ClusterNum,"Clusters",sep = "")
  #    paste(file_name,"_",umap1,"_",umap2,"_",metaCol,"Anno_UMAP.pdf", sep = '')
  #  },
  #  content = function(file) {
  #    plot <- umap_plot_PreC_kmean_react_base()
  #    ggsave(file,plot, width = 8, height = 8)
  #  }
  #)
  #
  #output$dnldKmeansClusterSubsetExpr <- downloadHandler(
  #  filename = function() {
  #    
  #    gs.u <- input$UMAPmatrixFile
  #    file_name <- gs.u$name
  #    file_name <- tools::file_path_sans_ext(gs.u)
  #    cluster_choice <- input$KmeansClusterSubsetnum
  #    
  #    paste(file_name,"_UMAPCluster",cluster_choice,"Subset.txt",sep = "")
  #  },
  #  content = function(file) {
  #    gs.u <- input$UMAPmatrixFile
  #    ext <- tools::file_ext(gs.u$datapath)
  #    req(gs.u)
  #    validate(need(ext == c("tsv","txt","csv","zip"), "Please upload .tsv, .txt, or .csv file"))
  #    
  #    if (ext == "csv") {
  #      expr <- as.data.frame(read_delim(gs.u$datapath, delim = ',', col_names = T))
  #    }
  #    else {
  #      expr <- as.data.frame(read_delim(gs.u$datapath, delim = '\t', col_names = T))
  #    }
  #    
  #    clusterTab <- umapClusterTable_react()
  #    cluster_choice <- input$KmeansClusterSubsetnum
  #    clusterTab_sub <- clusterTab[which(clusterTab$Cluster == cluster_choice),]
  #    samples <- clusterTab_sub$SampleName
  #    expr_sub <- expr[,c(colnames(expr)[1],samples)]
  #    
  #    write_delim(expr_sub,file,delim = '\t')
  #  }
  #)
  #
  #output$dnldKmeansClusterSubsetMeta <- downloadHandler(
  #  filename = function() {
  #    
  #    gs.u <- input$UMAPmetaFile
  #    file_name <- gs.u$name
  #    file_name <- tools::file_path_sans_ext(gs.u)
  #    cluster_choice <- input$KmeansClusterSubsetnum
  #    
  #    paste(file_name,"_UMAPCluster",cluster_choice,"Subset.txt",sep = "")
  #  },
  #  content = function(file) {
  #    
  #    meta <- umap_meta_loaded()
  #    
  #    clusterTab <- umapClusterTable_react()
  #    cluster_choice <- input$KmeansClusterSubsetnum
  #    clusterTab_sub <- clusterTab[which(clusterTab$Cluster == cluster_choice),]
  #    samples <- clusterTab_sub$SampleName
  #    meta_sub <- meta[which(meta$SampleName %in% samples),]
  #    
  #    write_delim(meta_sub,file,delim = '\t')
  #  }
  #)
  
  
}




















