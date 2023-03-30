####----Install and load packages----####

packages <- c("shiny","shinythemes","shinyjqui","shinycssloaders","dplyr","DT","readr","ggplot2","shinyjs")

installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
invisible(lapply(packages, library, character.only = TRUE))
#bioconductor packages
bioCpacks <- c("clusterProfiler","enrichplot")
installed_packages_BIOC <- bioCpacks %in% rownames(installed.packages())
if (any(installed_packages_BIOC == FALSE)) {
  BiocManager::install(bioCpacks[!installed_packages_BIOC], ask = F)
}
invisible(lapply(bioCpacks, library, character.only = TRUE))

ui <- 
  
  tagList(
    
    useShinyjs(),
    div(
      id = "loading_page",
      h1("Loading...")
    ),
    hidden(
      div(
        id = "main_content",
        
        navbarPage("{ Hazard Ratio Ranked GSEA }",
                   
                   ####----Enrichment Table----####
                   
                   tabPanel("Enrichment Table",
                            fluidPage(
                              title = "Enrichment Table",
                              sidebarPanel(
                                width = 3,
                                p("File input pre-loaded with PAN ICI melanoma pre-treatment OS CoxH genes ranked inverse. This is the file output from DRPPM-PATH-SURVEIOR Pipeline and the inverse of the hazard ration is used."),
                                fluidRow(
                                  column(9,
                                         fileInput("UserGeneList","Upload Ranked Gene List", placeholder = "ICI_iAtlas_Skin_Pre_OS_Genes_coxh_ranked_Inverse.txt")
                                  ),
                                  column(3,
                                         checkboxInput("UserGeneListheaderCheck","Header",value = T)
                                  )
                                ),
                                fluidRow(
                                  column(9,
                                         numericInput("gseaPval","GSEA Adj Pvalue Cutoff",value = 1)
                                  ),
                                  column(3,
                                         numericInput("UserSeedSet","Set Seed:",value = 101)
                                  )
                                ),
                                uiOutput("rendSpecimenType"),
                                uiOutput("rendGeneSetChosen"),
                                #uiOutput("rendUserUploadCheck"),
                                checkboxInput("UserUploadCheck","Upload Gene Set"),
                                fluidRow(
                                  column(9,
                                         uiOutput("rendUserGSupload")
                                  ),
                                  column(3,
                                         uiOutput("rendUserGSheaderCheck")
                                  )
                                ),
                                checkboxInput("PreviewGeneSet","Preview Selected Gene Set"),
                                uiOutput("rendGeneSetPreview"),
                              ),
                              mainPanel(
                                p("Please note this may take excess time depending on the number of gene sets being analyzed."),
                                uiOutput("rendEnrichTable"),
                                uiOutput("renddnldEnrichTable"),
                              )
                            )
                   ),
                   tabPanel("Enrichment Plot",
                            fluidPage(
                              title = "Enrichment Plot",
                              sidebarPanel(
                                h4("Select Gene Set:"),
                                div(DT::dataTableOutput("GeneSetTable"), style = "font-size:10px")
                              ),
                              mainPanel(
                                h3("GSEA Enrichment Plot"),
                                verbatimTextOutput("NESandPval"),
                                withSpinner(plotOutput("enrichplot", width = "500px", height = "450px"), type = 6),
                                fluidRow(
                                  downloadButton("dnldPlotSVG_gsea","Download as SVG"),
                                  downloadButton("dnldPlotPDF_gsea","Download as PDF")
                                ),
                                h3("Leading Edge Genes"),
                                downloadButton("LEGdownload", "Download Leading Edge Gene Table"),
                                div(DT::dataTableOutput("LeadingEdgeGenes"), style = "font-size:12px; height:500px; width:400px")
                                )
                              )
                            )
                   )
        )
      )
    )
  
  
  