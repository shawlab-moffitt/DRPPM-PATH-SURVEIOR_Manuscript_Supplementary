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




#increase file upload size
options(shiny.maxRequestSize=5000*1024^2)

shinytheme("sandstone")



ui <- 
  navbarPage("{ Jaccard Pathway Connectivity Index }",
             
             ####----Intra-Pathway Connectivity----####
             
             tabPanel("Pathway Connectivity",
                      fluidPage(
                        title = "Pathway Connectivity",
                        sidebarLayout(
                          sidebarPanel(
                            tabsetPanel(
                              
                              ####----Pathway Input----####
                              
                              tabPanel("Pathway Parameters",
                                       p(),
                                       p("File input pre-loaded with PAN ICI melanoma pre-treatment OS CoxH Immune Signatures ranked. This is the file output from DRPPM-PATH-SURVEIOR Pipeline."),
                                       h4("Upload Pathways of Interest"),
                                       fluidRow(
                                         column(9,
                                                fileInput("UserPathwayFile","Upload File (.gmt/.tsv/.txt)",
                                                          accept = c(".gmt",".tsv",".txt"),
                                                          placeholder = "ICI_iAtlas_Skin_Pre_ImmuneSig_OS_coxh_ranked_Filtered.txt"),
                                                uiOutput("rendTopFeatureSelect"),
                                         ),
                                         column(3,
                                                checkboxInput("HeaderCheckIntra","Header",value = T)
                                         )
                                       ),
                                       hr(),
                                       h4("Clustering Parameters"),
                                       fluidRow(
                                         column(4,
                                                selectInput("ClustMethodIntra","Clustering Method:",
                                                            choices = c("ward.D", "complete", "ward.D2", "single", "average", "mcquitty", "median", "centroid"))
                                         ),
                                         column(8,
                                                numericInput("NumClusters", step = 1, label = "Number of Clusters (Cut Tree with ~k)", value = 13)
                                         )
                                       ),
                                       checkboxInput("ViewClustTabIntra","View Cluster Results Table"),
                                       uiOutput("rendClustTabIntra"),
                                       downloadButton("dnldClustTabIntra","Download Cluster Result"),
                                       hr(),
                                       h4("SIF Download"),
                                       numericInput("JaccDistCutoff","Jaccard Distance Cutoff",
                                                    min = 0,max = 1, step = 0.1, value = 0.9, width = "200px"),
                                       checkboxInput("PrevSIF","Preview SIF File",value = F),
                                       uiOutput("rendSIFPreview"),
                                       downloadButton("dnldSIFTabIntra","Download SIF File")
                              ),
                              
                              ####----Figure Parameters----####
                              
                              tabPanel("Figure Parameters",
                                       h4("Heatmap Parameters"),
                                       selectInput("ColorPaletteIntra", "Select Color Palette:",
                                                   choices = c("Red/Blue" = "original",
                                                               "OmniBlueRed" = "OmniBlueRed",
                                                               "LightBlue/BlackRed" = "LightBlueBlackRed",
                                                               "Green/Black/Red" = "GreenBlackRed",
                                                               "Yellow/Green/Blue" = "YlGnBu","Inferno" = "Inferno",
                                                               "Viridis" = "Viridis","Plasma" = "Plasma",
                                                               "Reds" = "OrRd","Blues" = "PuBu","Greens" = "Greens")),
                                       fluidRow(
                                         column(6,
                                                checkboxInput("HeatColNamesIntra","Show Heatmap Column Names", value = F),
                                                numericInput("HeatColFontIntra", "Heatmap Column Font Size:",
                                                             min = 5, max = 75,
                                                             value = 12, step = 1),
                                                numericInput("HeatColDendHeight","Column Dendrogram Height:",
                                                             value = 50, step = 1)
                                         ),
                                         column(6,
                                                checkboxInput("HeatRowNamesIntra","Show Heatmap Row Names", value = F),
                                                numericInput("HeatRowFontIntra", "Heatmap Row Font Size:",
                                                             min = 5, max = 75,
                                                             value = 9, step = 1),
                                                numericInput("HeatRowDendHeight","Row Dendrogram Height:",
                                                             value = 50, step = 1)
                                         )
                                       ),
                                       hr(),
                                       h4("Connectivity Visualization Parameters"),
                                       selectInput("ConnView","View Connectivity as:",
                                                   choices = c("Phylogeny" = "phylogenic","Dendrogram" = "rectangle","Circular" = "circular")),
                                       uiOutput("rendPhyloLayout"),
                                       fluidRow(
                                         column(6,
                                                checkboxInput("ShowConnLabels","Show Labels:",value = T)
                                         ),
                                         column(6,
                                                numericInput("ConnFontSize","Font Size:",
                                                             value = 0.6, step = 0.1)
                                         )
                                       )
                              )
                            )
                          ),
                          mainPanel(
                            tabsetPanel(
                              
                              ####----Jaccard Table----####
                              
                              tabPanel("Jaccard Pathway Connectivity Table",
                                       p(),
                                       div(DT::dataTableOutput("JaccTableIntra"), style = "font-size:12px"),
                                       uiOutput("renddnldJaccTabIntra")
                              ),
                              
                              ####----Jaccard Heatmap----####
                              
                              tabPanel("Heatmap",
                                       withSpinner(jqui_resizable(plotOutput('JaccHeatmapIntra', width = "100%", height = "800px")), type = 6),
                                       downloadButton("dnldJaccHeatmapIntraPDF","Download as PDF"),
                                       downloadButton("dnldJaccHeatmapIntraSVG","Download as SVG")
                              ),
                              
                              ####----Jaccard Clustering----####
                              
                              tabPanel("Clustering",
                                       uiOutput("rendJaccDendo"),
                                       downloadButton("dnldJaccDendoPDF","Download as PDF"),
                                       downloadButton("dnldJaccDendoSVG","Download as SVG")
                              ),
                              
                              ####----Jaccard Cluster Annotation----####
                              
                              tabPanel("Gene Clusters and Annoation",
                                       fileInput("UserAnnotationFile","Upload Annotation File",
                                                 accept = c(".tsv",".txt",".csv")),
                                       div(DT::dataTableOutput("ClusterTabAnno"), style = "font-size:12px"),
                                       uiOutput("renddnldClusterTabAnno")
                              )
                            )
                          )
                        )
                      )
                      
                      
             ),
             tabPanel("Matrix Clustering",
                      fluidPage(
                        title = "Matrix Clustering",
                        sidebarLayout(
                          sidebarPanel(
                            tabsetPanel(
                              ####----Cluster Parameters----####
                              
                              tabPanel("Cluster Parameters",
                                       p(),
                                       fileInput("UserMatrixFile","Upload Matrix", accept = c(".txt",".tsv",".csv")),
                                       hr(),
                                       h4("Clustering Parameters"),
                                       fluidRow(
                                         column(4,
                                                selectInput("ClustMethodMat","Clustering Method:",
                                                            choices = c("ward.D", "complete", "ward.D2", "single", "average", "mcquitty", "median", "centroid"))
                                         ),
                                         column(8,
                                                numericInput("NumClustersMat", step = 1, label = "Number of Clusters (Cut Tree with ~k)", value = 10)
                                         )
                                       ),
                                       checkboxInput("ViewClustTabMat","View Cluster Results Table"),
                                       uiOutput("rendClustTabMat"),
                                       downloadButton("dnldClustTabMat","Download Cluster Result"),
                                       hr(),
                                       h4("SIF Download"),
                                       checkboxInput("PrevSIFMat","Preview SIF File",value = F),
                                       uiOutput("rendSIFPreviewMat"),
                                       downloadButton("dnldSIFTabMat","Download SIF File")
                              ),
                              tabPanel("Figure Parameters",
                                       p(),
                                       h4("Heatmap Parameters"),
                                       selectInput("ColorPaletteMat", "Select Color Palette:",
                                                   choices = c("Red/Blue" = "original",
                                                               "OmniBlueRed" = "OmniBlueRed",
                                                               "LightBlue/BlackRed" = "LightBlueBlackRed",
                                                               "Green/Black/Red" = "GreenBlackRed",
                                                               "Yellow/Green/Blue" = "YlGnBu","Inferno" = "Inferno",
                                                               "Viridis" = "Viridis","Plasma" = "Plasma",
                                                               "Reds" = "OrRd","Blues" = "PuBu","Greens" = "Greens")),
                                       fluidRow(
                                         column(6,
                                                checkboxInput("HeatColNamesMat","Show Heatmap Column Names", value = T),
                                                numericInput("HeatColFontMat", "Heatmap Column Font Size:",
                                                             min = 5, max = 75,
                                                             value = 12, step = 1),
                                                numericInput("HeatColDendHeightMat","Column Dendrogram Height:",
                                                             value = 50, step = 1)
                                         ),
                                         column(6,
                                                checkboxInput("HeatRowNamesMat","Show Heatmap Row Names", value = T),
                                                numericInput("HeatRowFontMat", "Heatmap Row Font Size:",
                                                             min = 5, max = 75,
                                                             value = 9, step = 1),
                                                numericInput("HeatRowDendHeightMat","Row Dendrogram Height:",
                                                             value = 50, step = 1)
                                         )
                                       ),
                                       hr(),
                                       h4("Connectivity Visualization Parameters"),
                                       selectInput("ConnViewMat","View Connectivity as:",
                                                   choices = c("Phylogeny" = "phylogenic","Dendrogram" = "rectangle","Circular" = "circular")),
                                       uiOutput("rendPhyloLayoutMat"),
                                       fluidRow(
                                         column(6,
                                                checkboxInput("ShowConnLabelsMat","Show Labels:",value = T)
                                         ),
                                         column(6,
                                                numericInput("ConnFontSizeMat","Font Size:",
                                                             value = 0.4, step = 0.1)
                                         )
                                       )
                              )
                            )
                          ),
                          mainPanel(
                            tabsetPanel(
                              tabPanel("Table Preview",
                                       p(),
                                       div(DT::dataTableOutput("UserMatrixPreview"), style = "font-size:12px")
                              ),
                              tabPanel("Heatmap",
                                       p(),
                                       withSpinner(jqui_resizable(plotOutput('UserMatrixHeatmap', width = "100%", height = "800px")), type = 6),
                                       downloadButton("dnldUserMatrixHeatmapPDF","Download as PDF"),
                                       downloadButton("dnldUserMatrixHeatmapSVG","Download as SVG")
                              ),
                              tabPanel("Clustering",
                                       uiOutput("rendUserMatixCLustering"),
                                       downloadButton("dnldUserMatixCLusteringPDF","Download as PDF"),
                                       downloadButton("dnldUserMatixCLusteringSVG","Download as SVG")
                              )
                            )
                          )
                        )
                      )
             )
             #tabPanel("UMAP Connectivity",
             #         fluidPage(
             #           title = "UMAP Connectivity",
             #           sidebarLayout(
             #             sidebarPanel(width = 3,
             #               conditionalPanel(condition = "input.MainUMAPpan == '1'",
             #                                h4("Upload Data"),
             #                                fileInput("UMAPmatrixFile","Expression Matrix", accept = c(".txt",".tsv",".csv",".zip")),
             #                                fileInput("UMAPmetaFile","Meta Data", accept = c(".txt",".tsv",".csv",".zip")),
             #                                h4("Transform Data"),
             #                                fluidRow(
             #                                  column(3,
             #                                         numericInput("UserSetSeed","Set Seed:", value = 101)
             #                                  ),
             #                                  column(4, style = 'padding-right:0px;',
             #                                         checkboxInput("LogUMAPmatrix","Log Matrix", value = T)
             #                                  ),
             #                                  column(5, style = 'padding-left:0px;',
             #                                         checkboxInput("NormUMAPmatrix","Scale Normalize Matrix", value = T)
             #                                  )
             #                                  )
             #                                ),
             #               conditionalPanel(condition = "input.MainUMAPpan != '1'",
             #                                tabsetPanel(
             #                                  id = "mainSidebar",
             #                                  tabPanel("Data Parameters",
             #                                           h4("Annotate UMAP"),
             #                                           uiOutput("rendUMAPsampSelect"),
             #                                           tabsetPanel(
             #                                             id = "AnnoUMAP",
             #                                             tabPanel("Clinical",
             #                                                      p(),
             #                                                      fluidRow(
             #                                                        column(8,
             #                                                               uiOutput("rendUMAPannotateSamps")
             #                                                        ),
             #                                                        column(4,
             #                                                               uiOutput("rendUMAPannoContCheck")
             #                                                        )
             #                                                      ),
             #                                                      value = 1),
             #                                             tabPanel("Expression",
             #                                                      p(),
             #                                                      fluidRow(
             #                                                        column(4, style = 'padding-right:2px;',
             #                                                               uiOutput("rendGeneSelection")
             #                                                        ),
             #                                                        column(2, style = 'padding-right:0px;padding-left:2px;',
             #                                                               checkboxInput("LogGeneSelection","Log2",value = T)
             #                                                        ),
             #                                                        column(6, style = 'padding-left:0px;',
             #                                                               uiOutput("rendGeneExprRange")
             #                                                        )
             #                                                      ),
             #                                                      hr(),
             #                                                      value = 2),
             #                                             tabPanel("K-Means Clustering",
             #                                                      p(),
             #                                                      fluidRow(
             #                                                        column(4,
             #                                                               selectInput("ClusterMethod", "Cluster Method:",
             #                                                                           choices = c("ward.D", "complete", "ward.D2", "single",
             #                                                                                       "average", "mcquitty", "median", "centroid")),
             #                                                               checkboxInput("ViewUMAPclusterTab","View Cluster Table", value = F)
             #                                                        ),
             #                                                        column(7,
             #                                                               numericInput("ClusterNumber","Number of Clusters (Cut Tree with ~k)",value = 3),
             #                                                               downloadButton("dnldUMAPclusterTab","Download Cluster Table")
             #                                                        )
             #                                                      ),
             #                                                      div(DT::dataTableOutput("umapClusterTable"), style = "font-size:10px"),
             #                                                      h4("Download Cluster Data"),
             #                                                      fluidRow(
             #                                                        column(4,
             #                                                               numericInput("KmeansClusterSubsetnum","Cluster",value = 1,step = 1)
             #                                                        ),
             #                                                        column(7,
             #                                                               downloadButton("dnldKmeansClusterSubsetExpr","Download Cluster Expression"),
             #                                                               downloadButton("dnldKmeansClusterSubsetMeta","Download Cluster Meta")
             #                                                        )
             #                                                      ),
             #                                                      hr(),
             #                                                      value = 3)
             #                                             ),
             #                                           uiOutput("rendUMAPtypeHeader"),
             #                                           uiOutput("rendRadioPathSelect"),
             #                                           textOutput("test_text"),
             #                                           fluidRow(
             #                                             column(8,
             #                                                    uiOutput("rendUserGSupload")
             #                                             ),
             #                                             column(4,
             #                                                    uiOutput("rendUserGSheaderCheck")
             #                                             )
             #                                           ),
             #                                           uiOutput("rendGeneSetCatSelect"),
             #                                           uiOutput("rendPathSelectTable"),
             #                                           fluidRow(
             #                                             column(4,
             #                                                    uiOutput("rendTopNumMVG")
             #                                             ),
             #                                             column(4,
             #                                                    uiOutput("rendVarMethodMVG")
             #                                             ),
             #                                             column(4,
             #                                                    uiOutput("rendVeiwMVGTable")
             #                                             )
             #                                           ),
             #                                           uiOutput("rendMVGlist"),
             #                                           uiOutput("renddnldMVGtab"),
             #                                           value = 1),
             #                                  tabPanel("Figure Parameters",
             #                                           p(),
             #                                           selectInput("UMAPcolors","Continuous Variable Color Palette:",
             #                                                       choices = c("Red" = "Reds","Blue" = "Blues","Purple" = "Purples",
             #                                                                   "Green" = "Greens","Orange" = "Oranges","Grey" = "Greys",
             #                                                                   "Blue/Green" = "BuGn","Yellow/Green" = "YlGn","Red/Purple" = "RdPu")),
             #                                           selectInput("UMAPorientation","UMAP Plot Trio Orientation",choices = c("Side-by-Side","Stacked")),
             #                                           fluidRow(
             #                                             column(6,
             #                                                    textInput("UMAPplotHeight","UMAP Plot Height:",value = "500px",
             #                                                              placeholder = "e.g. '500px','auto','100%'")
             #                                                    ),
             #                                             column(6,
             #                                                    textInput("UMAPplotWidth","UMAP Plot Width:",value = "100%",
             #                                                              placeholder = "e.g. '450px','auto','100%'")
             #                                                    )
             #                                           ),
             #                                           hr(),
             #                                           h4("Font Sizes"),
             #                                           fluidRow(
             #                                             column(3,
             #                                                    numericInput("UMAPtitleTextSize","Title:", value = 16)
             #                                             ),
             #                                             column(3,
             #                                                    numericInput("UMAPaxisTextSize","Axis:", value = 12)
             #                                             ),
             #                                             column(3,
             #                                                    numericInput("UMAPlegendTextSize","Legend:", value = 11)
             #                                             ),
             #                                             column(3,
             #                                                    numericInput("UMAPdotSize","Dot size:", value = 2)
             #                                             )
             #                                             ),
             #                                           hr(),
             #                                           fluidRow(
             #                                             column(6, style = 'padding-right:2px;',
             #                                                    textInput("UMAPxAxisLim","X-Axis Limits: min,max",value = "")
             #                                             ),
             #                                             column(6, style = 'padding-right:2px;padding-left:2px;',
             #                                                    textInput("UMAPyAxisLim","Y-Axis Limits: min,max",value = "")
             #                                             )
             #                                           ),
             #                                           value = 2)
             #                                  )
             #                                )
             #             ),
             #           mainPanel(
             #             tabsetPanel(
             #               id = "MainUMAPpan",
             #               tabPanel("Data Input Preview",
             #                        tabsetPanel(
             #                          tabPanel("Meta Data Preview",
             #                                   p(),
             #                                   div(DT::dataTableOutput("UMAPMetaPreview"), style = "font-size:12px")
             #                          ),
             #                          tabPanel("Expression Matrix Preview",
             #                                   p(),
             #                                   div(DT::dataTableOutput("UMAPMatrixPreview"), style = "font-size:12px")
             #                          )
             #                        ),
             #                        value = 1),
             #               tabPanel("Most Variable Genes UMAP",
             #                        p(),
             #                        fluidRow(
             #                          column(8,
             #                                 uiOutput("rendUMAPmvgMainTitle")
             #                                 ),
             #                          column(4,
             #                                 jqui_draggable(jqui_resizable(plotOutput('UMAPmvgLegend', width = "100%", height = "110px")))
             #                                 )
             #                        ),
             #                        p(),
             #                        uiOutput("rendUMAPplot_MVG_ALL"),
             #                        p(),
             #                        downloadButton("dnldUMAP_SVG_MVG_clin","Download Clinical Annotation UMAP"),
             #                        downloadButton("dnldUMAP_SVG_MVG_expr","Download Gene Expression UMAP"),
             #                        downloadButton("dnldUMAP_SVG_MVG_kmean","Download K-means Clustering UMAP"),
             #                        p(),
             #                        div(DT::dataTableOutput("UMAP_MVG_CoordTable"), style = "font-size:12px"),
             #                        downloadButton("dnldUMAPcoords_MVG","Download UMAP Coordinates"),
             #                        value = 4),
             #               tabPanel("Pathway UMAP",
             #                        p(),
             #                        fluidRow(
             #                          column(8,
             #                                 uiOutput("rendUMAPpathMainTitle")
             #                          ),
             #                          column(4,
             #                                 jqui_draggable(jqui_resizable(plotOutput('UMAPpathLegend', width = "100%", height = "110px")))
             #                          )
             #                        ),
             #                        p(),
             #                        uiOutput("rendUMAPplot_PATH_ALL"),
             #                        p(),
             #                        downloadButton("dnldUMAP_PATH_SVG_clin","Download Clinical Annotation UMAP"),
             #                        downloadButton("dnldUMAP_PATH_SVG_expr","Download Gene Expression UMAP"),
             #                        downloadButton("dnldUMAP_PATH_SVG_kmean","Download K-means Clustering UMAP"),
             #                        p(),
             #                        div(DT::dataTableOutput("UMAP_PATH_CoordTable"), style = "font-size:12px"),
             #                        downloadButton("dnldUMAPcoords_PATH","Download UMAP Coordinates"),
             #                        value = 3),
             #               tabPanel("All Genes UMAP",
             #                        p(),
             #                        fluidRow(
             #                          column(8,
             #                                 h3("UMAP of All Genes from Expression Matrix")
             #                          ),
             #                          column(4,
             #                                 jqui_draggable(jqui_resizable(plotOutput('UMAPallLegend', width = "100%", height = "110px")))
             #                          )
             #                        ),
             #                        p(),
             #                        uiOutput("rendUMAPplot_ALL"),
             #                        p(),
             #                        downloadButton("dnldUMAP_SVG_clin","Download Clinical Annotation UMAP"),
             #                        downloadButton("dnldUMAP_SVG_expr","Download Gene Expression UMAP"),
             #                        downloadButton("dnldUMAP_SVG_kmean","Download K-means Clustering UMAP"),
             #                        p(),
             #                        div(DT::dataTableOutput("UMAP_CoordTable"), style = "font-size:12px"),
             #                        downloadButton("dnldUMAPcoords","Download UMAP Coordinates"),
             #                        value = 2),
             #               tabPanel("Pre-Calculated UMAP",
             #                        p(),
             #                        fluidRow(
             #                          column(3,
             #                                 uiOutput("rendSelectPreCalc1")
             #                          ),
             #                          column(3,
             #                                 uiOutput("rendSelectPreCalc2")
             #                          ),
             #                          column(6,
             #                                 jqui_draggable(jqui_resizable(plotOutput('UMAPprecLegend', width = "100%", height = "110px")))
             #                          )
             #                        ),
             #                        uiOutput("rendUMAPplot_PreC_ALL"),
             #                        downloadButton("dnldUMAP_SVG_PreC_clin","Download Clinical Annotation UMAP"),
             #                        downloadButton("dnldUMAP_SVG_PreC_expr","Download Gene Expression UMAP"),
             #                        downloadButton("dnldUMAP_SVG_PreC_kmean","Download K-means Clustering UMAP"),
             #                        p(),
             #                        div(DT::dataTableOutput("UMAP_PreC_CoordTable"), style = "font-size:12px"),
             #                        value = 6)
             #               )
             #             )
             #           )
             #           )
             #         )
             )
             
             
             
             
             
             
             
             
             
             
             
             
             
             
             