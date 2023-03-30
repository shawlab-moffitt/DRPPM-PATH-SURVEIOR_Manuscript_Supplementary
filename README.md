# DRPPM-PATH-SURVEIOR Manuscript Supplementary

To ensure reproducibility, transparency, and accessibility for the work we have done to develop the DRPPM-PATH-SURVEIOR suite of tools, we have arranged this github in tandem with our manuscript that we are working towards publication. Users may explore the code and figures from the use cases we worked on, including the TARGET AML drug screening and identification of immune pathways in melanoma patients. To access the input and output from our analysis, users can use the URL links below. Please note that some files in the Github are zipped, compared to the files found hosted on the external server due to compliance with Github's file size requirements.

The source and url links to the data below were also added to our existing GitHub page, https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR-Suite. Here you can find up-to-date source code to performe this analysis.

* For the parent directory of entire analysis please go here:
  * http://shawlab.science/shiny/DRPPM_PATH_SURVEIOR_ExampleUseCases/

# PAN ICI Melanoma Immune Pathway Identification Use Case Example

* Parent directory for this analysis:
  * http://shawlab.science/shiny/DRPPM_PATH_SURVEIOR_ExampleUseCases/PANICI_Melanoma_UseCase/
  * https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR_Manuscript_Supplementary/tree/main/PANICI_Melanoma_UseCase

## Step 1 - Pipeline mode

* This was used to derive a ranked gene and immune signatures list associated with patient survival
* Figures 1 and 3A from the manscript provide a visual representation of using the pipeline as well as incorporating it with hazard ratio ranked GSEA
* The work for this analysis can be found here:
  * http://shawlab.science/shiny/DRPPM_PATH_SURVEIOR_ExampleUseCases/PANICI_Melanoma_UseCase/Fig3A_Step1_Pipeline_Mode/
  * https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR_Manuscript_Supplementary/tree/main/PANICI_Melanoma_UseCase/Step1_Pipeline_Mode

## Step 2 - Immune Deconvolution

* To examine immune features we derived immune deconvolution scores through Xcell and Cibersort and performed a bivariate survival analysis.
* Visualization of this analysis can be seen in Figure 4
* The code used to derive and append the immune deconvolution scores to the patient meta data can be found here:
  * http://shawlab.science/shiny/DRPPM_PATH_SURVEIOR_ExampleUseCases/PANICI_Melanoma_UseCase/Fig4_Step2_Immune_Deconvolution/
  * https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR_Manuscript_Supplementary/tree/main/PANICI_Melanoma_UseCase/Step2_Immune_Deconvolution

## Step 3 - On The Fly Analysis Mode

* This is the app used to generate the figures and perform much of the visualization and real time analysis for the manuscript.
* Here we host a live version of the R Shiny application for users to interact with.
* A schematic of the application can be seen in figure 1
* There are many figures from the manuscript that provide a small visualization of the UI from the application, such as Figures 2C, 4B, 5, S1, and S2B & C.
* The application is available here:
  * http://shawlab.science/shiny/DRPPM_PATH_SURVEIOR_ExampleUseCases/PANICI_Melanoma_UseCase/Fig2_FigS1_Step3_On_The_Fly_Analysis_Mode/
  * https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR_Manuscript_Supplementary/tree/main/PANICI_Melanoma_UseCase/Step3_On_The_Fly_Analysis_Mode/PAN_ICI_Melanoma_Pre_SurvivalApp

## Step 4 - Hazard Ratio Ranked Gene Set Enrichment Analysis (GSEA)

* The ranked gene list derived from the step 1 pipeline mode was used as input to the hazard ratio ranked GSEA Shiny application.
* The ranked gene list associated with overall survival data was used.
* This analysis can be seen in Figures 3 and S2.
* The code used for this analysis and a ready to use hazard ration ranked GSEA App is available here:
  * http://shawlab.science/shiny/DRPPM_PATH_SURVEIOR_ExampleUseCases/PANICI_Melanoma_UseCase/Fig3_FigS2_Step4_HazardRatio_Ranked_GSEA/
  * https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR_Manuscript_Supplementary/tree/main/PANICI_Melanoma_UseCase/Step4_HR_Ranked_GSEA

## Step 5 - Pathway Connectivivty

* The ranked immune signatures output from the step 1 pipeline mode was used as input to the pathway connectivit Shiny application.
* The connectivity was calculated through Jaccard distance based on overlapping genes between signatures.
* This analysis can be seen in Figures 5 and S4.
* The code used for this analysis and a ready to use pathway connectivity app is available here:
  * http://shawlab.science/shiny/DRPPM_PATH_SURVEIOR_ExampleUseCases/PANICI_Melanoma_UseCase/Fig5_FigS4_Step5_Pathway_Connectivity/
  * https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR_Manuscript_Supplementary/tree/main/PANICI_Melanoma_UseCase/Step5_Pathway_Connectivity

# TARGET AML 1031 Drug Discovery Use Case Example

* Parent directory for this analysis:
  * http://shawlab.science/shiny/DRPPM_PATH_SURVEIOR_ExampleUseCases/TARGET_AML_UseCase/
  * https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR_Manuscript_Supplementary/tree/main/TARGET_AML_UseCase

## Interactive DRPPM-PATH-SURVEIOR R Shiny App featuring the TARGET AML 1031 data

* An interactive Shiny app that was used to generate a majority of the figures for this use case can be found here:
  * http://shawlab.science/shiny/DRPPM_PATH_SURVEIOR_ExampleUseCases/TARGET_AML_UseCase/TARGET_AML_KMT2A_SurvivalApp/

## Step 1 - DRPPM PATH SURVEIOR Pipeline

* A systematic drug screening was performed using LINCS L1000 drug induced downregulated genes based on the gene reversion strategy.
* This was performed through the pipeline visualized the schematic seen in Figure 1 and S5.
* Two of the top hits from this pipeline as well as a summary of enriched Cmap names and mechanisms of action (MOAs) can be seen in Figure 6.
* The code for the analysis can be found here:
  * http://shawlab.science/shiny/DRPPM_PATH_SURVEIOR_ExampleUseCases/TARGET_AML_UseCase/Fig6_Step1_DRPPM_PATH_SURVEIOR_Pipeline/
  * https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR_Manuscript_Supplementary/tree/main/TARGET_AML_UseCase/Step1_DRPPM_PATH_SURVEIOR_Pipeline

## Step 2 - Subset data for bootstrap analysis

* To examine the robustness of our results we performed a resampling bootstrap of 75% of the patient cohort and reperformed the analysis.
* The figure depicting the results we derived from this analysis can be seen in Figure S6.
* The code for the resampling can be found here:
  * http://shawlab.science/shiny/DRPPM_PATH_SURVEIOR_ExampleUseCases/TARGET_AML_UseCase/FigS6_Step2_Subset_Data_For_Bootstrap/
  * https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR_Manuscript_Supplementary/tree/main/TARGET_AML_UseCase/Step2_Subset_Data_For_Bootstrap

## Step 3 - DRPPM PATH SURVEIOR Pipeline of 75% Bootstrap patient cohort

* A systematic drug screening using the resampled 75% AML KMT2A-fusion positive patients was performed using the same pipeline that was used in step 1 for this cohort.
* The figure depicting the results we derived from this analysis can be seen in Figure S6.
* The code for the resampling can be found here:
  * http://shawlab.science/shiny/DRPPM_PATH_SURVEIOR_ExampleUseCases/TARGET_AML_UseCase/FigS6_Step3_DRPPM_PATH_SURVEIOR_Pipeline_Bootstrap/
  * https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR_Manuscript_Supplementary/tree/main/TARGET_AML_UseCase/Step3_DRPPM_PATH_SURVEIOR_Pipeline_Bootstrap

## Step 4 - LINCS L1000 Drug Prioritization

* The results of the Step 1 and Step 3 pipelines were ran through a screening and prioritization script to derive plots of some of the top enriched Cmap names and MOAs
* The prioritization of the bootstrap analysis also allowed us to validat and visualiza the robustness of our results.
* This step was used in producing Figures 6 and S6.
* Code for this analysis can be found here:
  * http://shawlab.science/shiny/DRPPM_PATH_SURVEIOR_ExampleUseCases/TARGET_AML_UseCase/Fig6_FigS6_Step4_LINCSL1000_Drug_Prioritization/
  * https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR_Manuscript_Supplementary/tree/main/TARGET_AML_UseCase/Step4_LINCSL1000_Drug_Prioritization



