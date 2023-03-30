


####----Input----####

ProjectName <- "PAN ICI Melanoma - Pre Treatment"

##--Advanced Setup--##
SurvPlot_Height <- "550px"
SurvPlot_Width <- "850px"




####----Install and load packages----####

## Check if Immune deconvolution package is installed
immudecon <- "immunedeconv"
immudecon_check <- immudecon %in% rownames(installed.packages())
if (immudecon_check == TRUE) {
  library(immunedeconv)
}
## Check for and load other packages
packages <- c("shiny","shinythemes","shinyjqui","gtsummary","tidyr","RColorBrewer",
              "dplyr","DT","ggplot2","ggpubr","tibble","survival","pheatmap",
              "readr","shinycssloaders","survminer","gridExtra","viridis","plotly","ggrepel")
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
invisible(lapply(packages, library, character.only = TRUE))
#bioconductor packages
bioCpacks <- c("GSVA","clusterProfiler")
installed_packages_BIOC <- bioCpacks %in% rownames(installed.packages())
if (any(installed_packages_BIOC == FALSE)) {
  BiocManager::install(bioCpacks[!installed_packages_BIOC], ask = F)
}
invisible(lapply(bioCpacks, library, character.only = TRUE))





####----UI----####

ui <-
  navbarPage(paste("{ ",ProjectName," Survival Analysis }",sep = ""),
             
             ####----Overall Survival Tab----####
             
             tabPanel("Survival Analysis",
                      fluidPage(
                        sidebarLayout(
                          
                          ####----Sidebar Panel----####
                          
                          sidebarPanel(
                            conditionalPanel(condition = "input.SurvPanels == '1' | input.SurvPanels == '2' | input.SurvPanels == '3' | input.SurvPanels == '4'",
                                             tabsetPanel(
                                               id = "survside",
                                               
                                               ##--Sample Parameters--##
                                               
                                               tabPanel("Sample Parameters",
                                                        p(),
                                                        uiOutput("rendSampleTypeSelection"),
                                                        uiOutput("rendFeatureSelection"),
                                                        uiOutput("rendSubFeatureSelection"),
                                                        fluidRow(
                                                          column(6,
                                                                 uiOutput("rendSurvivalType_time")
                                                          ),
                                                          column(6,
                                                                 uiOutput("rendSurvivalType_id")
                                                          )
                                                        ),
                                                        uiOutput("rendScoreMethodBox"),
                                                        tabsetPanel(
                                                          id = "GeneSetTabs",
                                                          tabPanel("Gene Sets",
                                                                   p(),
                                                                   uiOutput("rendGeneSetCat_Select"),
                                                                   uiOutput("rendGeneSetTable"),
                                                                   value = 1
                                                          ),
                                                          tabPanel("Single Genes",
                                                                   #radioButtons("RawOrSS","Survival Analysis By:",
                                                                   #             choices = c("Raw Gene Expression","Rank Normalized"),
                                                                   #             selected = "Raw Gene Expression", inline = T),
                                                                   uiOutput("rendGeneGeneSetTable"),
                                                                   value = 2
                                                          ),
                                                          tabPanel("User Gene Set",
                                                                   p(),
                                                                   radioButtons("UserGSoption","",choices = c("Gene Set Upload","Text Box Input"), inline = T),
                                                                   uiOutput("renduserGeneSet"),
                                                                   uiOutput("renduserGeneSetTextName"),
                                                                   uiOutput("renduserGeneSetText"),
                                                                   uiOutput("rendUserGeneSetTable"),
                                                                   value = 3
                                                          )
                                                        ),
                                                        uiOutput("rendViewGeneSetGenes"),
                                                        uiOutput("rendGenesInGeneSetTab")
                                               ),
                                               
                                               ##--Survival Parameters--##
                                               
                                               tabPanel("Risk Strat Parameters",
                                                        p(),
                                                        h4("Risk Stratification Plot Parameters"),
                                                        fluidRow(
                                                          column(6,
                                                                 numericInput("cutoffTime1","High-Risk Survival Time Cutoff:", value = 364, min = 0, step = 1),
                                                                 selectInput("survStatus1","Survival Status Below Cutoff:", choices = c("1","0","0/1"), selected = "1")
                                                          ),
                                                          column(6,
                                                                 numericInput("cutoffTime0","Low-Risk Survival Time Cutoff:", value = 365, min = 0, step = 1),
                                                                 selectInput("survStatus0","Survival Status Above Cutoff:", choices = c("1","0","0/1"), selected = "0")
                                                          )
                                                        )
                                               ),
                                               
                                               ##--Figure Parameters--##
                                               
                                               tabPanel("Figure Parameters",
                                                        p(),
                                                        h4("Survival Plot Parameters"),
                                                        fluidRow(
                                                          column(4,
                                                                 uiOutput("rendSurvXaxis")
                                                          ),
                                                          column(8,
                                                                 uiOutput("rendSurvPlotTitle")
                                                          )
                                                        ),
                                                        fluidRow(
                                                          column(3,
                                                                 selectInput("SurvLegendPos","Legend Position",choices = c("right","left","top","bottom","none"))
                                                          ),
                                                          column(3,
                                                                 checkboxInput("ShowPval","Show P.Value",value = T)
                                                          ),
                                                          column(3,
                                                                 checkboxInput("ShowConfInt","Show Confidence Interval",value = F)
                                                          ),
                                                          column(3,
                                                                 checkboxInput("ShowMedSurvLine","Show Median Survival Line",value = F)
                                                          )
                                                        ),
                                                        hr(),
                                                        h4("Boxplot Parameters"),
                                                        fluidRow(
                                                          column(6,
                                                                 numericInput("boxplotFont","Boxplot Font Size:", value = 15, step = 1),
                                                                 selectInput("boxoptselec","Boxplot Stat Compare Method:", choices = c("none","wilcox.test","t.test","kruskal.test","anova")) 
                                                          ),
                                                          column(6,
                                                                 numericInput("boxplotDot", "Boxplot Dot Size:", value = 0.75, step = 0.25),
                                                                 selectInput("boxplotTextAngle","X-Axis Text Orientation",
                                                                             choices = c("Horizontal (0 degrees)" = "0","Angled (45 degrees)" = "45","Vertical (90 degrees)" = "90","Stagger"))
                                                          )
                                                        ),
                                                        hr(),
                                                        h4("Heatmap Parameters"),
                                                        selectInput("ClusterMethod", "Select Cluster Method",
                                                                    choices = c("complete", "ward.D", "ward.D2", "single", "average", "mcquitty", "median", "centroid")),
                                                        fluidRow(
                                                          column(6,
                                                                 numericInput("heatmapFontR", "Heatmap Row Font Size:", value = 9, step = 1)
                                                          ),
                                                          column(6,
                                                                 numericInput("heatmapFontC", "Heatmap Column Font Size:", value = 10, step = 1)
                                                          )
                                                        ),
                                                        selectInput("ColorPaletteHeat", "Select Color Palette:",
                                                                    choices = c("Red/Blue" = "original",
                                                                                "OmniBlueRed" = "OmniBlueRed",
                                                                                "LightBlue/BlackRed" = "LightBlueBlackRed",
                                                                                "Green/Black/Red" = "GreenBlackRed",
                                                                                "Yellow/Green/Blue" = "YlGnBu","Inferno" = "Inferno",
                                                                                "Viridis" = "Viridis","Plasma" = "Plasma",
                                                                                "Reds" = "OrRd","Blues" = "PuBu","Greens" = "Greens")
                                                        ),
                                                        hr(),
                                                        h4("Forest Plot Parameters"),
                                                        numericInput("ForestFontSize","Font Size",value = 1),
                                                        hr(),
                                                        h4("Linearity Plot Parameters"),
                                                        fluidRow(
                                                          column(4,
                                                                 numericInput("linAxisFont","X/Y Axis Font Size",
                                                                              value = 14, step = 1)
                                                          ),
                                                          column(4,
                                                                 numericInput("linTickFont","Axis Tick Font Size",
                                                                              value = 10, step = 1)
                                                          ),
                                                          column(4,
                                                                 numericInput("linMainFont","Title Font Size",
                                                                              value = 16, step = 1)
                                                          )
                                                        )
                                               )
                                             )
                                             ),
                            conditionalPanel(condition = "input.SurvPanels == '5'",
                                             tabsetPanel(
                                               tabPanel("Input Parameters",
                                                        p(),
                                                        h4("Sample Selection"),
                                                        uiOutput("rendSampleTypeSelection_lasso"),
                                                        uiOutput("rendFeatureSelection_lasso"),
                                                        uiOutput("rendSubFeatureSelection_lasso"),
                                                        h4("Feature Selection"),
                                                        textInput("LassoModelName","Lasso Model Name:",value = "Custom_Lasso_Model"),
                                                        uiOutput("rendLassoFeatureSelection_lasso"),
                                                        h4("Lasso Parameters"),
                                                        fluidRow(
                                                          column(6,
                                                                 uiOutput("rendSurvTimeSelec_lasso")
                                                          ),
                                                          column(6,
                                                                 uiOutput("rendSurvIDSelect_lasso")
                                                          )
                                                        ),
                                                        fluidRow(
                                                          column(5,
                                                                 numericInput("LassoTrainProp","Training Sample Proportion:",value = 50, step = 1, max = 100,min = 1)
                                                          ),
                                                          column(4,
                                                                 numericInput("LassoAlpha","Set Lasso Alpha:", value = 1, min = 0, step = 0.1)
                                                          ),
                                                          column(3,
                                                                 numericInput("LassoSeedSelection","Set Seed:", value = 100, min = 1, step = 1)
                                                          )
                                                        ),
                                                        p(),
                                                        fluidRow(
                                                          column(8,
                                                                 radioButtons("viewLassoMinOrSE","Select Lambda to Generate Risk Score:", inline = T,choices = c("Lambda Min","Lambda SE","Custom")),
                                                                 ),
                                                          column(4,
                                                                 uiOutput("rendCustomLambda")
                                                                 )
                                                        ),
                                                        uiOutput("rednLassoCoefTable"),
                                                        #,
                                                        fluidRow(
                                                          column(6,
                                                                 downloadButton("dnldLassoModel","Download Lasso Model")
                                                                 ),
                                                          column(6,
                                                                 downloadButton("dnldLassoRunData","Download Lasso Run Data")
                                                                 )
                                                          )
                                                        ),
                                               tabPanel("Figure Parameters",
                                                        p(),
                                                        h4("Survival Plot Parameters"),
                                                        radioButtons("LassoPlotCutP","Plot Cut-Point", choices = c("Median","Quartile","Optimal","Quantile","User Specified"),
                                                                     inline = T),
                                                        uiOutput("rendCutPinput"),
                                                        fluidRow(
                                                          column(4,
                                                                 uiOutput("rendSurvXaxis_lasso")
                                                          ),
                                                          column(8,
                                                                 textInput("SurvPlotTitle_lasso","Survival Plot Title:",value = "")
                                                          )
                                                        ),
                                                        fluidRow(
                                                          column(3,
                                                                 selectInput("SurvLegendPos_lasso","Legend Position",choices = c("top","right","left","bottom","none"))
                                                          ),
                                                          column(3,
                                                                 checkboxInput("ShowPval_lasso","Show P.Value",value = T)
                                                          ),
                                                          column(3,
                                                                 checkboxInput("ShowConfInt_lasso","Show Confidence Interval",value = F)
                                                          ),
                                                          column(3,
                                                                 checkboxInput("ShowMedSurvLine_lasso","Show Median Survival Line",value = F)
                                                          )
                                                        )
                                                        )
                                               )
                                             )
                          ),
                          
                          ####----Main Panel----####
                          
                          mainPanel(
                            tabsetPanel(
                              id = "SurvPanels",
                              
                              ####----Survival Analysis Tab----####
                              
                              tabPanel("Pathway Level Survival Analysis",
                                       tabsetPanel(
                                         id = "SurvPanelsMain",
                                         
                                         ##--Median Cut Point--##
                                         
                                         tabPanel("Median Cut-Point Survival",
                                                  p(),
                                                  htmlOutput("BINSurvDescrip", style = "font-size:14px;"),
                                                  hr(),
                                                  withSpinner(jqui_resizable(plotOutput("SplotBIN", height = SurvPlot_Height, width = SurvPlot_Width)), type = 6),
                                                  fluidRow(
                                                    downloadButton("dnldSplotBIN_SVG","Download as SVG"),
                                                    downloadButton("dnldSplotBIN_PDF","Download as PDF")
                                                  ),
                                                  hr(),
                                                  withSpinner(jqui_resizable(plotOutput("ssgseaBINDensity", width = "650px", height = "400px")), type = 6),
                                                  fluidRow(
                                                    downloadButton("dnldssgseaBINDensity_SVG","Download as SVG"),
                                                    downloadButton("dnldssgseaBINDensity_PDF","Download as PDF")
                                                  ),
                                                  hr(),
                                                  h4("Cox Hazard Regression Analysis Summary"),
                                                  fluidRow(
                                                    column(6,
                                                           uiOutput("rendBINHRtab"),
                                                           style = 'border-right: 0.5px solid lightgray',
                                                    ),
                                                    column(6,
                                                           verbatimTextOutput("MedianCutPSumm")
                                                    )
                                                  ),
                                                  value = 1),
                                         
                                         ##--Quaritle Cutpoints--##
                                         
                                         tabPanel("Quartile Survival",
                                                  p(),
                                                  htmlOutput("QuartSurvDescrip", style = "font-size:14px;"),
                                                  hr(),
                                                  withSpinner(jqui_resizable(plotOutput("Splot", height = SurvPlot_Height, width = SurvPlot_Width)), type = 6),
                                                  fluidRow(
                                                    downloadButton("dnldSplot_SVG","Download as SVG"),
                                                    downloadButton("dnldSplot_PDF","Download as PDF")
                                                  ),
                                                  hr(),
                                                  withSpinner(jqui_resizable(plotOutput("ssgseaQuartDensity", width = "650px", height = "400px")), type = 6),
                                                  fluidRow(
                                                    downloadButton("dnldssgseaQuartDensity_SVG","Download as SVG"),
                                                    downloadButton("dnldssgseaQuartDensity_PDF","Download as PDF")
                                                  ),
                                                  hr(),
                                                  h4("Cox Hazard Regression Analysis Summary"),
                                                  fluidRow(
                                                    column(6,
                                                           uiOutput("rendQuartHRtab"),
                                                           style = 'border-right: 0.5px solid lightgray',
                                                    ),
                                                    column(6,
                                                           verbatimTextOutput("QuartileCutPSumm")
                                                    )
                                                  ),
                                                  value = 2),
                                         
                                         ##--Optimal Cut Point--##
                                         
                                         tabPanel("Optimal Cut-Point Survival",
                                                  p(),
                                                  htmlOutput("CutPSurvDescrip", style = "font-size:14px;"),
                                                  hr(),
                                                  withSpinner(jqui_resizable(plotOutput("ScutPointPlot", height = SurvPlot_Height, width = SurvPlot_Width)), type = 6),
                                                  fluidRow(
                                                    downloadButton("dnldScutPointPlot_SVG","Download as SVG"),
                                                    downloadButton("dnldScutPointPlot_PDF","Download as PDF")
                                                  ),
                                                  hr(),
                                                  withSpinner(jqui_resizable(plotOutput("ssgseaCutPDensity", width = "650px", height = "400px")), type = 6),
                                                  fluidRow(
                                                    downloadButton("dnldssgseaCutPDensity_SVG","Download as SVG"),
                                                    downloadButton("dnldssgseaCutPDensity_PDF","Download as PDF")
                                                  ),
                                                  hr(),
                                                  h4("Cox Hazard Regression Analysis Summary"),
                                                  fluidRow(
                                                    column(6,
                                                           uiOutput("rendCutPointHRtab"),
                                                           style = 'border-right: 0.5px solid lightgray',
                                                    ),
                                                    column(6,
                                                           verbatimTextOutput("OptimalCutPSumm")
                                                    )
                                                  ),
                                                  value = 3),
                                         
                                         ##--Quantile Cutpoints--##
                                         
                                         tabPanel("Top/Bottom Cut-Point Survival",
                                                  p(),
                                                  htmlOutput("QuantSurvDescrip", style = "font-size:14px;"),
                                                  hr(),
                                                  withSpinner(jqui_resizable(plotOutput("SquantPlot", height = SurvPlot_Height, width = SurvPlot_Width)), type = 6),
                                                  numericInput("QuantPercent","Top/Bottom Cut-Point Quantile Cutoff (%)", value = 25, min = 0, max = 100),
                                                  fluidRow(
                                                    downloadButton("dnldSquantPlot_SVG","Download as SVG"),
                                                    downloadButton("dnldSquantPlot_PDF","Download as PDF")
                                                  ),
                                                  hr(),
                                                  withSpinner(jqui_resizable(plotOutput("ssgseaQuantDensity", width = "650px", height = "400px")), type = 6),
                                                  fluidRow(
                                                    downloadButton("dnldssgseaQuantDensity_SVG","Download as SVG"),
                                                    downloadButton("dnldssgseaQuantDensity_PDF","Download as PDF")
                                                  ),
                                                  hr(),
                                                  h4("Cox Hazard Regression Analysis Summary"),
                                                  fluidRow(
                                                    column(6,
                                                           uiOutput("rendQuantHRtab"),
                                                           style = 'border-right: 0.5px solid lightgray',
                                                    ),
                                                    column(6,
                                                           verbatimTextOutput("QuantileCutPSumm")
                                                    )
                                                  ),
                                                  value = 4),
                                         
                                         ##--User Cut Point--##
                                         
                                         tabPanel("User Cut-Point Survival",
                                                  p(),
                                                  htmlOutput("Quant2SurvDescrip", style = "font-size:14px;"),
                                                  hr(),
                                                  withSpinner(jqui_resizable(plotOutput("SquantPlot2", height = SurvPlot_Height, width = SurvPlot_Width)), type = 6),
                                                  numericInput("QuantPercent2","Above/Below User Quantile Cut-Point (%)", value = 25, min = 0, max = 100),
                                                  fluidRow(
                                                    downloadButton("dnldSquantPlot2_SVG","Download as SVG"),
                                                    downloadButton("dnldSquantPlot2_PDF","Download as PDF")
                                                  ),
                                                  hr(),
                                                  withSpinner(jqui_resizable(plotOutput("ssgseaQuant2Density", width = "650px", height = "400px")), type = 6),
                                                  fluidRow(
                                                    downloadButton("dnldssgseaQuant2Density_SVG","Download as SVG"),
                                                    downloadButton("dnldssgseaQuant2Density_PDF","Download as PDF")
                                                  ),
                                                  hr(),
                                                  h4("Cox Hazard Regression Analysis Summary"),
                                                  fluidRow(
                                                    column(6,
                                                           uiOutput("rendQuantHRtab2"),
                                                           style = 'border-right: 0.5px solid lightgray',
                                                    ),
                                                    column(6,
                                                           verbatimTextOutput("UserCutPSumm")
                                                    )
                                                  ),
                                                  value = 5)
                                       ),
                                       value = 1),
                              
                              ####----Univariate Survival----####
                              
                              tabPanel("Univariate Survival Analysis",
                                       p(),
                                       fluidRow(
                                         column(5,
                                                uiOutput("rendSurvivalFeatureSingle"),
                                                ## Allows all select inputs to be wide enough to read the contents
                                                tags$head(
                                                  tags$style(HTML('
                                                                  .selectize-input {
                                                                      white-space: nowrap;
                                                                  }
                                                                  .selectize-dropdown {
                                                                      width: 500px !important;
                                                                  }'
                                                  )
                                                  )
                                                ),
                                                fluidRow(
                                                  column(3,
                                                         checkboxInput("UniVarNAcheck","Remove NA/Unknown/Inf",value = T)
                                                  ),
                                                  column(3,
                                                         checkboxInput("UniVarContCheck","Continuous Feature",value = F)
                                                  ),
                                                  column(3,
                                                         uiOutput("rendUniVarContHiLoCheck")
                                                  )
                                                ),
                                                uiOutput("rendSurvFeatVariableUni"),
                                         ),
                                         column(7,
                                                #htmlOutput("UnivarSummExpl", style = "font-size:14px;"),
                                         )
                                       ),
                                       tabsetPanel(
                                         id = "UniVarPlots",
                                         
                                         ##--Survival Plot--##
                                         
                                         tabPanel("Survival Plot",
                                                  p(),
                                                  withSpinner(jqui_resizable(plotOutput("featSplot", width = SurvPlot_Width, height = SurvPlot_Height)), type = 6),
                                                  fluidRow(
                                                    downloadButton("dnldfeatSplot_SVG","Download as SVG"),
                                                    downloadButton("dnldfeatSplot_PDF","Download as PDF")
                                                  )
                                                  
                                         ),
                                         
                                         ##--Coxh Tables--##
                                         
                                         tabPanel("Coxh Table",
                                                  p(),
                                                  fluidRow(
                                                    column(6,
                                                           div(withSpinner(tableOutput("SSingleFeatureHRtab"), type = 7, size = 0.5), style = "font-size:12px; width:500px; overflow-X: scroll")
                                                    ),
                                                    column(6,
                                                           verbatimTextOutput("UnivarSummary")
                                                    )
                                                  )
                                         ),
                                         
                                         ##--Forest Plot--##
                                         
                                         tabPanel("Forest Plot",
                                                  p(),
                                                  withSpinner(jqui_resizable(plotOutput("SinglevarForestPlot", width = "100%", height = "800px")), type = 6),
                                                  fluidRow(
                                                    downloadButton("dnldUniVarForestplot_SVG","Download as SVG"),
                                                    downloadButton("dnldUniVarForestplot_PDF","Download as PDF")
                                                  )
                                         ),
                                         
                                         ##--Linearity Check--##
                                         
                                         tabPanel("Linearity Check",
                                                  p(),
                                                  fluidRow(
                                                    column(3,
                                                           selectInput("ResidualTypeUni","Select Residual Type",
                                                                       choices = c("deviance", "martingale", "score", "schoenfeld", "dfbeta", "dfbetas", "scaledsch", "partial"))
                                                    ),
                                                    column(3,
                                                           selectInput("linPredict1", "X-axis Scale:",
                                                                       choices = c("linear.predictions","observation.id","time"))
                                                    )
                                                  ),
                                                  uiOutput("timewarnmessage1"),
                                                  withSpinner(jqui_resizable(plotOutput("UnivarLinearityPlot", width = "100%", height = "500px")), type = 6),
                                                  fluidRow(
                                                    downloadButton("dnldUniVarLinplot_SVG","Download as SVG"),
                                                    downloadButton("dnldUniVarLinplot_PDF`","Download as PDF")
                                                  )
                                         )
                                       ),
                                       value = 2),
                              
                              ####----Multivariate Survival----####
                              
                              tabPanel("Multivariate Coxh Analysis",
                                       tabsetPanel(
                                         id = "multivariate",
                                         
                                         ####----Bivariate Additive Survival----####
                                         
                                         tabPanel("Bivariate Additive Survival Analysis",
                                                  p(),
                                                  fluidRow(
                                                    column(4,
                                                           uiOutput("rendSurvivalFeatureBi1"),
                                                           fluidRow(
                                                             column(4,
                                                                    checkboxInput("BiVarAddNAcheck1","Remove NA/Unknown/Inf",value = T)
                                                             ),
                                                             column(4,
                                                                    checkboxInput("BiVarAddContCheck1","Continuous Feature",value = F)
                                                             ),
                                                             column(4,
                                                                    uiOutput("rendBiVarAddContHiLoCheck1")
                                                             )
                                                           ),
                                                           uiOutput("rendSurvFeatVariableBi1")
                                                    ),
                                                    column(4,
                                                           uiOutput("rendSurvivalFeatureBi2"),
                                                           fluidRow(
                                                             column(4,
                                                                    checkboxInput("BiVarAddNAcheck2","Remove NA/Unknown/Inf",value = T)
                                                             ),
                                                             column(4,
                                                                    checkboxInput("BiVarAddContCheck2","Continuous Feature",value = F)
                                                             ),
                                                             column(4,
                                                                    uiOutput("rendBiVarAddContHiLoCheck2")
                                                             )
                                                           ),
                                                           uiOutput("rendSurvFeatVariableBi2")
                                                    ),
                                                    column(4,
                                                           htmlOutput("BivarAddSummExpl", style = "font-size:14px;")
                                                    )
                                                  ),
                                                  tabsetPanel(
                                                    id = "BiVarPlots",
                                                    
                                                    ##--Coxh Tables--##
                                                    
                                                    tabPanel("Cox HR Table",
                                                             p(),
                                                             fluidRow(
                                                               column(6,
                                                                      div(withSpinner(tableOutput("BiFeatureHRtab"), type = 7, size = 0.5), style = "font-size:12px; width:500px; overflow-X: scroll")
                                                               ),
                                                               column(6,
                                                                      verbatimTextOutput("bivarSummary"),
                                                                      fluidRow(
                                                                        column(6,
                                                                               verbatimTextOutput("bivarAnova1")
                                                                        ),
                                                                        column(6,
                                                                               verbatimTextOutput("bivarAnova2")
                                                                        )
                                                                      )
                                                               )
                                                             )
                                                    ),
                                                    
                                                    ##--Forest Plot--##
                                                    
                                                    tabPanel("Forest Plot",
                                                             p(),
                                                             withSpinner(jqui_resizable(plotOutput("BivarForestPlot", width = "100%", height = "800px")), type = 6),
                                                             fluidRow(
                                                               downloadButton("dnldBiVarAddForest_SVG","Download as SVG"),
                                                               downloadButton("dnldBiVarAddForest_PDF","Download as PDF")
                                                             )
                                                    ),
                                                    
                                                    ##--Linearity Check--##
                                                    
                                                    tabPanel("Linearity Check",
                                                             p(),
                                                             fluidRow(
                                                               column(3,
                                                                      selectInput("ResidualTypeBi","Select Residual Type",
                                                                                  choices = c("deviance", "martingale", "score", "schoenfeld", "dfbeta", "dfbetas", "scaledsch", "partial"))
                                                               ),
                                                               column(3,
                                                                      selectInput("linPredict2", "X-axis Scale:",
                                                                                  choices = c("linear.predictions","observation.id","time"))
                                                               )
                                                             ),
                                                             uiOutput("timewarnmessage2"),
                                                             withSpinner(jqui_resizable(plotOutput("BivarLinearityPlot", width = "100%", height = "500px")), type = 6),
                                                             fluidRow(
                                                               downloadButton("dnldBiVarAddLinplot_SVG","Download as SVG"),
                                                               downloadButton("dnldBiVarAddLinplot_PDF","Download as PDF")
                                                             )
                                                    )
                                                  )
                                         ),
                                         
                                         ####----Bivariate Interaction Survival----####
                                         
                                         tabPanel("Bivariate Interaction Survival Analysis",
                                                  p(),
                                                  fluidRow(
                                                    column(4,
                                                           uiOutput("rendSurvivalFeatureBi1Inter"),
                                                           fluidRow(
                                                             column(4,
                                                                    checkboxInput("BiVarIntNAcheck1","Remove NA/Unknown/Inf",value = T)
                                                             ),
                                                             column(4,
                                                                    checkboxInput("BiVarIntContCheck1","Continuous Feature",value = F)
                                                             ),
                                                             column(4,
                                                                    uiOutput("rendBiVarIntContHiLoCheck1")
                                                             )
                                                           ),
                                                           uiOutput("rendSurvFeatVariableBi1Inter")
                                                           
                                                    ),
                                                    column(4,
                                                           uiOutput("rendSurvivalFeatureBi2Inter"),
                                                           fluidRow(
                                                             column(4,
                                                                    checkboxInput("BiVarIntNAcheck2","Remove NA/Unknown/Inf",value = T)
                                                             ),
                                                             column(4,
                                                                    checkboxInput("BiVarIntContCheck2","Continuous Feature",value = F)
                                                             ),
                                                             column(4,
                                                                    uiOutput("rendBiVarIntContHiLoCheck2")
                                                             )
                                                           ),
                                                           uiOutput("rendSurvFeatVariableBi2Inter")
                                                    ),
                                                    column(4,
                                                           verbatimTextOutput("BivarIntSummExpl")
                                                    )
                                                  ),
                                                  tabsetPanel(
                                                    id = "BiVarInterTabs",
                                                    
                                                    ##--Survival Plot--##
                                                    
                                                    tabPanel("Survival Plot",
                                                             p(),
                                                             withSpinner(jqui_resizable(plotOutput("featSplotBi", width = SurvPlot_Width, height = SurvPlot_Height)), type = 6),
                                                             fluidRow(
                                                               downloadButton("dnldfeatSplotBi_SVG","Download as SVG"),
                                                               downloadButton("dnldfeatSplotBi_PDF","Download as PDF")
                                                             )
                                                    ),
                                                    
                                                    ##--Coxh Tables--##
                                                    
                                                    tabPanel("Cox HR Table",
                                                             p(),
                                                             fluidRow(
                                                               column(6,
                                                                      div(withSpinner(tableOutput("BiFeatureHRtabInter"), type = 7, size = 0.5), style = "font-size:12px; width:500px; overflow-X: scroll")
                                                               ),
                                                               column(6,
                                                                      verbatimTextOutput("bivarSummaryInter"),
                                                                      verbatimTextOutput("bivarAnovaInter1")
                                                               )
                                                             )
                                                    ),
                                                    
                                                    ##--Forest Plot--##
                                                    
                                                    tabPanel("Forest Plot",
                                                             p(),
                                                             withSpinner(jqui_resizable(plotOutput("BivarForestPlotInter", width = "100%", height = "800px")), type = 6),
                                                             fluidRow(
                                                               downloadButton("dnldBiVarIntForest_SVG","Download as SVG"),
                                                               downloadButton("dnldBiVarIntForest_PDF","Download as PDF")
                                                             )
                                                    ),
                                                    
                                                    ##--Linearity Check--##
                                                    
                                                    tabPanel("Linearity Check",
                                                             p(),
                                                             fluidRow(
                                                               column(3,
                                                                      selectInput("ResidualTypeInter","Select Residual Type",
                                                                                  choices = c("deviance", "martingale", "score", "schoenfeld", "dfbeta", "dfbetas", "scaledsch", "partial"))
                                                               ),
                                                               column(3,
                                                                      selectInput("linPredict3", "X-axis Scale:",
                                                                                  choices = c("linear.predictions","observation.id","time"))
                                                               )
                                                             ),
                                                             uiOutput("timewarnmessage3"),
                                                             withSpinner(jqui_resizable(plotOutput("BivarLinearityPlotInter", width = "100%", height = "500px")), type = 6),
                                                             fluidRow(
                                                               downloadButton("dnldBiVarIntLinplot_SVG","Download as SVG"),
                                                               downloadButton("dnldBiVarIntLinplot_PDF","Download as PDF")
                                                             )
                                                    )
                                                    
                                                  )
                                         ),
                                         
                                         ####----Multivariate Survival----####
                                         
                                         tabPanel("Multivariate Coxh Analysis",
                                                  p(),
                                                  fluidRow(
                                                    column(4,
                                                           uiOutput("rendSurvivalFeature"),
                                                           checkboxInput("MultiVarNAcheck","Remove NA/Unknown",value = T)
                                                    )
                                                  ),
                                                  tabsetPanel(
                                                    id = "multivartabstwo",
                                                    
                                                    ##--Coxh Tables--##
                                                    
                                                    tabPanel("Coxh Tables",
                                                             fluidRow(
                                                               column(6,
                                                                      h4("Coxh Hazard Ratio (Categorical)"),
                                                                      verbatimTextOutput("multivarSummaryCat"),
                                                                      div(withSpinner(tableOutput("SFeatureHRtabCat"), type = 7, size = 0.5), style = "font-size:12px; width:500px; overflow-X: scroll")
                                                               ),
                                                               column(6,
                                                                      h4("Coxh Hazard Ratio (Continuous)"),
                                                                      verbatimTextOutput("multivarSummaryCont"),
                                                                      div(withSpinner(tableOutput("SFeatureHRtabCont"), type = 7, size = 0.5), style = "font-size:12px; width:500px; overflow-X: scroll")
                                                               )
                                                             )
                                                    ),
                                                    
                                                    ##--Forest Plot--##
                                                    
                                                    tabPanel("Forest Plot",
                                                             p(),
                                                             withSpinner(jqui_resizable(plotOutput("MultivarForestPlot", width = "100%", height = "800px")), type = 6),
                                                             fluidRow(
                                                               downloadButton("dnldMultiVarForest_SVG","Download as SVG"),
                                                               downloadButton("dnldMultiVarForest_PDF","Download as PDF")
                                                             )
                                                    )
                                                  )
                                         )
                                       ),
                                       value = 3),
                              
                              ####----Lasso Analysis----####
                              
                              #tabPanel("Lasso Cox Model",
                              #         p(),
                              #         fluidRow(
                              #           column(2, 
                              #                  br(),
                              #                  actionButton("RunLassoModelGen","Generate Lasso Model")
                              #                  ),
                              #           column(3,
                              #                  uiOutput("rendSurvTimeSelecView_lasso")
                              #                  ),
                              #           column(3,
                              #                  uiOutput("rendSurvIDSelecView_lasso")
                              #                  )
                              #         ),
                              #         fluidRow(
                              #           column(6,
                              #                  withSpinner(jqui_resizable(plotOutput("Lasso_Train_Splot", width = "100%", height = "500px")),type = 6),
                              #                  downloadButton("dnldSplotLassoTrain_SVG","Download as SVG"),
                              #                  hr(),
                              #                  h4("Path of Coefficients"),
                              #                  jqui_resizable(plotOutput("Lasso_CoeffPlot", width = "100%", height = "400px")),
                              #                  hr(),
                              #                  h4("Training Cox Hazard Regression Analysis Summary"),
                              #                  uiOutput("rendLassoTrainHRtab"),
                              #                  verbatimTextOutput("LassoTrainCoxSumm")
                              #           ),
                              #           column(6,
                              #                  withSpinner(jqui_resizable(plotOutput("Lasso_Test_Splot", width = "100%", height = "500px")),type = 6),
                              #                  downloadButton("dnldSplotLassoTest_SVG","Download as SVG"),
                              #                  hr(),
                              #                  h4("Lambda Cross-Validation"),
                              #                  jqui_resizable(plotOutput("Lasso_LambdaPlot", width = "100%", height = "400px")),
                              #                  hr(),
                              #                  h4("Testing Cox Hazard Regression Analysis Summary"),
                              #                  uiOutput("rendLassoTestHRtab"),
                              #                  verbatimTextOutput("LassoTestCoxSumm")
                              #           )
                              #           ),
                              #         value = 5),
                              
                              ####----Data Exploration----####
                              
                              tabPanel("Data Exploration",
                                       tabsetPanel(
                                         id = "DataExploration",
                                         tabPanel("Download Survival Data",
                                                  p(),
                                                  uiOutput("rendMetaTableCols"),
                                                  uiOutput("rendMetaTable"),
                                                  fluidRow(
                                                    column(6,
                                                           uiOutput("DnldMetaButon")
                                                    ),
                                                    column(6,
                                                           uiOutput("DnldExprButon")
                                                    )
                                                  ),
                                                  value = 5),
                                         tabPanel("Score Density",
                                                  p(),
                                                  fluidRow(
                                                    numericInput("densityPercent","User Defined Percentile (Red)",value = 15, width = "200px"),
                                                    checkboxInput("QuartileLinesCheck","Show Quartile Lines (Blue)",value = T)
                                                  ),
                                                  withSpinner(jqui_resizable(plotOutput("ssgseaDensity", width = "100%", height = "500px")), type = 6),
                                                  fluidRow(
                                                    downloadButton("dnldssgseaDensity_SVG","Download as SVG"),
                                                    downloadButton("dnldssgseaDensity_PDF","Download as PDF")
                                                  ),
                                                  div(DT::dataTableOutput("ssgseaDensityTable"), style = "font-size:12px"),
                                                  downloadButton("dnldssgseaDensityTable","Download Table"),
                                                  value = 6),
                                         tabPanel("Feature Comparison",
                                                  p(),
                                                  fluidRow(
                                                    column(3,
                                                           uiOutput("rendScatterFeature")
                                                    ),
                                                    column(2,
                                                           radioButtons("ColorScatterChoice","Color Plot by:",choices = c("Feature","Single Color"))
                                                    ),
                                                    column(3,
                                                           uiOutput("rendScatterColor")
                                                    ),
                                                    column(2,
                                                           checkboxGroupInput("ScatterLog","", choices = c("Log x-axis","Log y-axis"))
                                                    ),
                                                  ),
                                                  withSpinner(jqui_resizable(plotlyOutput("FeatCompScatterPlot", width = "100%", height = "500px")), type = 6),
                                                  fluidRow(
                                                    downloadButton("dnldFeatCompScatter_SVG","Download as SVG"),
                                                    downloadButton("dnldFeatCompScatter_PDF","Download as PDF")
                                                  ),
                                                  p(),
                                                  div(DT::dataTableOutput("FeatCompScatterTable"), style = "font-size:12px"),
                                                  downloadButton("dnldFeatCompScatterTable","Download Table"),
                                                  value = 7),
                                         tabPanel("Risk Stratification",
                                                  tabsetPanel(
                                                    tabPanel("Risk Straification Boxplot",
                                                             p("Users may adjust risk-stratification cutoff parameters under the 'Risk Strat Parameters' tab on the sidebar panel."),
                                                             checkboxInput("SBoxLog", "Log Transform Score", value = T),
                                                             withSpinner(jqui_resizable(plotOutput("Sboxplot", width = "100%", height = "500px")), type = 6),
                                                             fluidRow(
                                                               downloadButton("dnldSboxplot_SVG","Download as SVG"),
                                                               downloadButton("dnldSboxplot_PDF","Download as PDF")
                                                             ),
                                                             div(DT::dataTableOutput("SboxplotTable"), style = "font-size:12px; height:450px; overflow-Y: scroll"),
                                                             p(),
                                                             downloadButton("dnldSBoxplotTab","Download Table")
                                                    ),
                                                    tabPanel("Risk Straification Heatmap",
                                                             p("Users may adjust risk-stratification cutoff parameters under the 'Risk Strat Parameters' tab on the sidebar panel."),
                                                             uiOutput("heatmap_error_message"),
                                                             withSpinner(jqui_resizable(plotOutput("Sheatmap", width = "100%", height = "2000px")), type = 6),
                                                             fluidRow(
                                                               downloadButton("dnldSheatmap_SVG","Download as SVG"),
                                                               downloadButton("dnldSheatmap_PDF","Download as PDF"),
                                                               downloadButton("dnldSheatmapexpr","Download Expression Matrix From Heatmap")
                                                             )
                                                    )
                                                  )
                                         ),
                                         tabPanel("Feature Stratification",
                                                  tabsetPanel(
                                                    tabPanel("Feature Boxplot",
                                                             p(),
                                                             fluidRow(
                                                               column(4,
                                                                      uiOutput("rendBoxplotFeature")
                                                               ),
                                                               column(3,
                                                                      checkboxInput("BoxPRemoveNA","Remove NA/Unknowns",value = T)
                                                               ),
                                                               column(3,
                                                                      checkboxInput("FBoxLog", "Log Transform Score", value = T)
                                                               )
                                                             ),
                                                             withSpinner(jqui_resizable(plotOutput("Featureboxplot", width = "100%", height = "500px")), type = 6),
                                                             fluidRow(
                                                               downloadButton("dnldFboxplot_SVG","Download as SVG"),
                                                               downloadButton("dnldFboxplot_PDF","Download as PDF")
                                                             ),
                                                             div(DT::dataTableOutput("FeatureboxplotTable"), style = "font-size:12px; height:450px; overflow-Y: scroll"),
                                                             p(),
                                                             downloadButton("dnldFeatureboxplotTab","Download Table")
                                                    ),
                                                    tabPanel("Feature Heatmap",
                                                             p(),
                                                             uiOutput("heatmap_error_message2"),
                                                             fluidRow(
                                                               column(4,
                                                                      uiOutput("rendHeatmapFeature")
                                                               ),
                                                               column(3,
                                                                      checkboxInput("HeatRemoveNA","Remove NA/Unknowns",value = T)
                                                               )
                                                             ),
                                                             withSpinner(jqui_resizable(plotOutput("FeatureHeatmap", width = "100%", height = "2000px")), type = 6),
                                                             fluidRow(
                                                               downloadButton("dnldFheatmap_SVG","Download as SVG"),
                                                               downloadButton("dnldFheatmap_PDF","Download as PDF"),
                                                               downloadButton("dnldFheatmapexpr","Download Expression Matrix From Heatmap")
                                                             )
                                                    )
                                                  )
                                         )
                                       ),
                                       value = 4)
                            )
                          )
                        )
                      )
             ),
             tabPanel("About",
                      fluidPage(
                        mainPanel(
                          tabPanel("Purpose and Methods",
                                   uiOutput("rendPurposeAndMethodsMD")
                                   #tabsetPanel(
                                   #  tabPanel("Purpose and Methods",
                                   #           uiOutput("rendPurposeAndMethodsMD")),
                                   #  tabPanel("Tutorial",
                                   #           uiOutput("rendTutorialMD"))
                                   #)
                          )
                        )
                      )
             )
  )
