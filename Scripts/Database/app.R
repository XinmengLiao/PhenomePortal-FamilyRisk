# Install and load necessary packages
if (!require("shiny")) install.packages("shiny")
if (!require("shinythemes")) install.packages("shinythemes")
if (!require("shinyBS")) install.packages("shinyBS")
if (!require("DT")) install.packages("DT")
if (!require("readxl")) install.packages("readxl")
if (!require("dplyr")) install.packages("dplyr")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("jsonlite")) install.packages("jsonlite")
if (!require("stringr")) install.packages("stringr")
if (!require("tidyr")) install.packages("tidyr")
if (!require("UpSetR")) install.packages("UpSetR")
if (!require("plotly")) install.packages("plotly")
if (!require("shinyjs")) install.packages("shinyjs")
if (!require("VennDiagram")) install.packages("VennDiagram")
if (!require("grid")) install.packages("grid")

library(shiny)
library(shinythemes)
library(shinyBS)
library(DT)
library(readxl)
library(dplyr)
library(ggplot2)
library(jsonlite)
library(stringr)
library(tidyr)
library(UpSetR)
library(plotly)
library(tidyverse)
library(shinyjs)
library(VennDiagram)
library(grid)

# Define CSS for homepage styling
css_code <- "
.project-grid {
  display: grid;
  grid-template-columns: repeat(4, 1fr);
  gap: 20px;
  margin: 30px auto;
  max-width: 1400px;
  padding: 0 40px;
  grid-auto-rows: minmax(300px, auto);
  width: 100%;
  box-sizing: border-box;
}

.project-card {
  border: 1px solid #ddd;
  border-radius: 8px;
  padding: 20px;
  text-align: center;
  background-color: #f8f9fa;
  cursor: pointer;
  transition: all 0.3s ease;
  box-shadow: 0 2px 4px rgba(0,0,0,0.1);
  min-height: 280px;
  display: flex;
  flex-direction: column;
  justify-content: center;
  align-items: center;
}

.col-md-3 {
  margin-bottom: 20px;
}

.home-container {
  width: 80%;
  #margin: 0 40px;
  margin: 0 auto;
  padding: 20px 40px;
}

.project-card:hover {
  box-shadow: 0 4px 12px rgba(0,0,0,0.2);
  transform: translateY(-5px);
  background-color: #e8f4f8;
}

.project-logo {
  width: 150px;
  height: 120px;
  margin: 0 auto 15px;
  display: flex;
  align-items: center;
  justify-content: center;
  background-color: #f8f9fa;
  border-radius: 8px;
}

.project-logo img {
  max-width: 100%;
  max-height: 100%;
  object-fit: contain;
}

.project-name {
  font-size: 16px;
  font-weight: bold;
  color: #333;
  margin: 8px 0;
}

.project-description {
  font-size: 14px;
  color: #666;
  margin-top: 10px;
}

.screening-results-container {
  width: 80%;
  margin: 0 auto;
  padding: 20px;
}

.generator-container {
  width: 80%;
  margin: 0 auto;
  padding: 20px;
}

.detail-panel-container {
  width: 80%;
  margin: 0 auto;
  text-align: center;
}

.section-title {
  display: inline-block;
  background: linear-gradient(135deg, #2c3e50 0%, #34495e 100%);
  color: white;
  padding: 12px 20px;
  border-radius: 8px;
  font-size: 24px;
  font-weight: bold;
  margin: 15px 0;
}

/* DT Filter Styling */
.dataTables_filter {
  margin-bottom: 15px;
}

.dataTables_filter input {
  width: 300px;
  padding: 8px 12px;
  border: 1px solid #ddd;
  border-radius: 4px;
  font-size: 14px;
}

thead input, thead select {
  width: 100%;
  padding: 5px;
  box-sizing: border-box;
  font-size: 13px;
  border: 1px solid #ddd;
  border-radius: 3px;
}

thead input:focus, thead select:focus {
  outline: none;
  border-color: #4CAF50;
  box-shadow: 0 0 5px rgba(76, 175, 80, 0.3);
}

.dt-filter-select {
  max-height: 200px;
  overflow-y: auto;
}

/* 改进表格外观 */
.dataTables_wrapper .dataTables_length,
.dataTables_wrapper .dataTables_filter,
.dataTables_wrapper .dataTables_info,
.dataTables_wrapper .dataTables_paginate {
  margin: 10px 0;
}

table.dataTable thead th {
  background-color: #f5f5f5;
  border-bottom: 2px solid #ddd;
  font-weight: 600;
  color: #333;
}

/* Collapsible panel styling */
.panel-group {
  margin-bottom: 15px;
  border-radius: 4px;
  box-shadow: 0 1px 1px rgba(0,0,0,.05);
}

.panel {
  background-color: #fff;
  border: 1px solid #ddd;
  border-radius: 4px;
  box-shadow: 0 1px 1px rgba(0,0,0,.05);
  margin-bottom: 0;
}

.panel-heading {
  background-color: #f5f5f5;
  border: 1px solid #ddd;
  border-radius: 4px 4px 0 0;
  padding: 0;
  transition: all 0.3s ease;
}

.panel-heading:hover {
  background-color: #efefef;
}

.panel-title {
  font-size: 15px;
  margin: 0;
}

.panel-title a {
  cursor: pointer;
}

.panel-body {
  padding: 15px;
  border-top: 1px solid #ddd;
}

.panel-collapse.collapse {
  display: none;
}

.panel-collapse.collapse.in {
  display: block;
}
"

# Define UI
ui <- navbarPage(
  "FamilyRisk Database",
  theme = shinytheme("flatly"),
  tags$head(
    useShinyjs(),
    tags$style(HTML(css_code)),
    tags$script(HTML("
      $(function() {
        // Handle project card clicks
        $(document).on('click', '.project-card', function() {
          var projectId = $(this).data('project-id');
          Shiny.setInputValue('selected_project_id', projectId, {priority: 'event'});
          setTimeout(function() {
            $('a[data-value=\"NBSeq Details\"]').tab('show');
          }, 100);
        });
        
        // Handle NBSeq Details tab click to reset project selection
        $(document).on('shown.bs.tab', 'a[data-value=\"NBSeq Details\"]', function() {
          // When NBSeq Details tab is shown, check if we should reset
          var currentSelected = $('#selected_project_id').val();
          // If a project is selected and user clicked the tab, reset it
          Shiny.setInputValue('check_tab_state', Date.now(), {priority: 'event'});
        });
        
        // Handle collapse toggle for Advanced Filters with Bootstrap
        $(document).on('click', '.panel-title a[data-toggle=\"collapse\"]', function(e) {
          e.preventDefault();
          var targetPanel = $(this).attr('href');
          $(targetPanel).collapse('toggle');
          // Update arrow icon
          var icon = $(this).find('.collapse-icon');
          if ($(targetPanel).hasClass('in')) {
            $(this).html('▼ Advanced Filters');
          } else {
            $(this).html('▶ Advanced Filters');
          }
        });
      });
    "))
  ),
  tabPanel(
    "Home",
    value = "Home",
    div(
      class = "home-container",
      h1("Welcome to the FamilyRisk Database!", style = "text-align: left; margin-top: 30px;"),
      br(),
      p("FamilyRisk Database is a unified knowledge base for newborn and carrier screening. It brings together curated gene–disease screening lists from six leading NBSeq projects and ACMG recommendations. For each gene-disease associaiton, the corresponding validity from GenCC is provided. More details about the projects and screening lists could be checked in “Screening List Details”. This database also aggretateds validated screening results from relevant publications, which could be checked in “NBSeq Screening Results”. ", 
        style = "text-align: justify; font-size: 20px; color: #333; margin-bottom: 20px; line-height: 1.4;"),
      
      div(
        style = "text-align: center; margin-bottom: 30px;",
        img(src = "logo/DB_Overview.jpg", alt = "Database Overview", style = "max-width: 90%; height: auto;")
      ),

      br(),
      
      
      fluidRow(
        column(12, 
               h3("Disease Ontology", class = "section-title"),
               p("The Database collects 3115 gene-disease associations from 8 different screening lists, covering 13 disease ontologies. Disease with more than one ontology is classified as \"Complex genetic syndromes\". See the plot below of ontology distribution across different projects/screening lists.",
                 style = "text-align: justify; color: #333; font-size: 20px; line-height: 1.4; margin-bottom: 15px;"),
               div(
                 plotlyOutput("ontology_chart", height = "600px", width = "80%"),
                 style = "display: flex; justify-content: center; width: 100%;"
               )
        ),
        column(12, 
               br(),
               h3("Screening List Intersection", class = "section-title"),
               p("Some gene-disease associations for screening is overlapped between different lists. See the plot below of the statistics of the intersections across projects/screening lists.",
                 style = "text-align: justify; color: #333; font-size: 20px; line-height: 1.4; margin-bottom: 15px;"),
               div(
                 plotlyOutput("upset_chart", height = "600px", width = "85%"),
                 style = "display: flex; justify-content: center; width: 100%;"
               )
        )
      ),
      
      br(),br(),br(),
      
      h3("Screening List Generator", style = "text-align: left; margin-top: 30px;", class = "section-title"),
      p("This tools allows users to create and download a tailored screening list, which could be applied directly in FamilyRisk. Filter and select the screening gene-disease association of interest by project/list, gene, disease, inheritance, ontology, and GenCC evidence, or add additional gene-disease association. Export your screening panel for downstream genomics workflows.", 
        style = "text-align: justify; font-size: 20px; color: #333; margin-bottom: 20px; line-height: 1.4;"),

      br()
      
      # h2("NBSeq Projects and Expanded Carrier Screening List", style = "text-align: left; margin-top: 30px;"),
      # p("Click on a project to view the project information and detailed screening list.", 
      #   style = "text-align: left; color: #666; margin-bottom: 20px;"),
      # # Projects section - 2 rows x 4 columns
      # uiOutput("project_cards"),
      # br()
    )
  ),
  tabPanel(
    "Screening List Details",
    value = "NBSeq Details",
    uiOutput("project_detail_panel")
  ),
  tabPanel(
    "NBSeq Screening Results",
    value = "NBSeq Screening Results",
    div(
      class = "screening-results-container",
      h2("Overview of the screening results"),
      p("FamilyRisk database collects published screening results of positive monogenic disease from four NBSeq projects: BabySeq, BabyDetect, EarlyCheck, and Guardian. Only BabySeq and EarlyCheck disclosed positive disease carriers in their publications.",
        tags$br(),
        "The screening rate of each NBSeq projects are:", tags$br(), tags$br(),
        "GUARDIAN: 3.7%", tags$br(),
        "BabyDetect: 1.8%", tags$br(),
        "BabyScreen+: 1.6%", tags$br(),
        "BabySeq: 9.4%", tags$br(),
        "EarlyCheck: 3.8%." , tags$br(), tags$br(),
        "See the plots below for statistics of overlapping mutated genes associated with positive screening results.",
        style = "text-align: justified; color: #666; font-size: 18px"),
      br(),br(),
      fluidRow(
        column(4,
               h4("Number of mutated genes of screen-positive rare dieases", style = "text-align: justified; color: #333; font-size: 18px"),
               p(HTML("* G6PD is shared among GUARDIAN, BabyDetect, BabyScreen+. <br/>* G6PD and ACADS are shared between GUARDIAN and BabyDetect.<br/>* BRCA2 is shared between BabySeq and EarlyCheck.")),
               div(
                 plotlyOutput("positive_rate_plot", height = "500px"),
                 style = "width: 80%"
               )
        ),
        column(4,
          h4("Number of the intersected mutated genes of disease carrier status", style = "text-align: justified; color: #333; font-size: 18px"),
          div(
            plotOutput("carrier_venn_plot", height = "400px", width = "400px"),
            style = "width: 60%;"
          )
        ),
        column(4,
          h4("Number of the intersected mutated genes of screen-positive rare diseases", style = "text-align: justified; color: #333; font-size: 18px"),
          div(
            plotOutput("positive_venn_plot", height = "500px", width = "500px"),
            style = "width: 60%; "
          )
        )
          ),

      br(),br(),
      
      h2("Details of the Screening results. "),
      p("'Positives' refers to the screen-positive rare diseases. 'Carriers' refers to the disease carrier-status.",
        style = "font-size: 18px"),
      br(),
      uiOutput("screening_results_tabs")
    )
  ),
  tabPanel(
    "Screening List Generator",
    value = "Screening List Generator",
    div(
      class = "generator-container",
      uiOutput("generator_ui")
    )
  )
)

# Define server logic
server <- function(input, output, session) {

  # Add resource path for logo folder
  shiny::addResourcePath("logo", file.path(getwd(), "logo"))

  # Define featured projects (8 main projects) with mapping to data IDs
  featured_projects <- list(
    list(id = "BabySeq", name = "BabySeq", short_name = "BabySeq", logo = "BabySeq.jpg",
         data_ids = c("BabySeq_GroupA", "BabySeq_GroupB"),
         description = "The BabySeq Project is a pilot randomized clinical trial that explores the medical, behavioral, and economic impacts of nGS in well newborns and those admitted to a neonatal intensive care unit (NICU).",
         website = "https://www.babyseqproject.org/", publication = "https://doi.org/10.1016/j.ajhg.2018.11.016",
         criteria = "1. Variants classified as disease-causing in public databases with minor allele frequency < 5.0% in gnomAD.\n2. Predicted loss-of-function variants (nonsense, frameshift, ±1/2 splice-site) in disease-associated genes with minor allele frequency ≤ 0.1% in gnomAD.",
         returned_results = "1.	Monogenic disease risk: (likely) pathogenic variants in genes associated with dominantly inherited diseases or as bi-allelic (likely) pathogenic variants in genes associated with recessively inherited diseases) that present or are manageable during childhood. \n2.	Monogenic disease carrier risk: the carrier status for any gene meeting the monogenic disease risk reporting criteria. \n3.	PGx-associated genes (DPYD, TPMT, G6PD) that are related to atypical reactions to medications used in the pediatric population."),
    list(id = "NBScreening", name = "NBScreening", short_name = "NBScreening", logo = "NBScreening.png",
         data_ids = c("NBScreening"),
         description = "NBScreening is one of the subprojects of the Anatolian Precision Medicine Initiative (APMI). It was initiated in 2023 and represents the first NBSeq project in Türkiye. ",
         website = "https://mphenome.org/APMI", publication = "Coming soon.",
         criteria = "1.	ClinVar (likely) pathogenic and ACMG (likely) pathogenic variants.\n2.	Common variants (population-level minor allele frequency > 5%) within Turkish population are excluded.\n3.	The minor allele frequency < 1%.",
         returned_results = "1.	Monogenic disease risk: One (likely) pathogenic variant for autosomal dominant diseases. Hemizygous (likely) pathogenic variants in male neonates and homozygous (likely) pathogenic variants in female neonates for X-linked diseases.\n2.	Monogenic disease carrier risk: Two or more (likely) pathogenic variants for autosomal recessive diseases.\n3.	PGx-associated genes: CACNA1S, CFTR, CYP2B6, CYP2C19, CYP2C19_GKB, CYP2C9, CYP2D6, DPYD, G6PD, HLA-A, HLA-B, IFNL3, MT-RNR1, NAT2, NUDT15, RYR1, SLCO1B1, TPMT, UGT1A1, VKORC1."),
    list(id = "BabyDetect", name = "BabyDetect", short_name = "BabyDetect", logo = "BabyDetect.png",
         data_ids = c("BabyDetect"),
         description = "The BabyDetect Project is designed to enroll babies born in the maternities of La Citadelle or Notre-Dame des Bruyèresuntil summer 2025.",
         website = "https://babydetect.com/en/", publication = "https://doi.org/10.1038/s41591-024-03465-x",
         criteria = "1.	ACMG class 4 (likely pathogenic) or class 5 (pathogenic) variants.\n2.	ClinVar variants and the ones included in their own curated list. ClinVar variants are are subsequently reviewed with particular caution using the Franklin and VarSome platforms, which have an advanced artificial intelligence-driven engine designed to prioritize and interpret variant data.\n3.	The minor allele frequency threshold for defining rare variants is not available.",
         returned_results = "Monogenic disease risk: One (likely) pathogenic variant for autosomal dominant diseases. Two (likely) pathogenic variants for autosomal recessive diseases. Hemizygous (likely) pathogenic variants in male neonates and homozygous or possible compound heterozygous (likely) pathogenic variants in female neonates for X-linked diseases."),
    list(id = "BabyScreen+", name = "BabyScreen+", short_name = "BabyScreen+", logo = "BabyScreen+.jpg",
         data_ids = c("BabyScreen+"),
         description = "BabyScreen+ is a research study where parents can choose to have their baby’s heel prick sample screened for over 500 additional treatable, childhood-onset conditions. This is called genomic newborn screening. This study is investigating the best way to deliver genomic newborn screening for babies across Victoria.",
         website = "https://babyscreen.mcri.edu.au/", publication = "https://doi.org/10.1038/s41591-025-03986-z",
         criteria = "1.	(Likely) pathogenic variants.\n2.	For known pathogenic variants, minor allele frequency (MAF) < 2% in gnomAD; for high impact variant, MAF < 0.01% for monoallelic (dominant) disease or MAF < 0.1% for biallelic (recessive) disease.",
         returned_results = "Only Monogenic disease risk will be reported. Carrier status and adult-onset conditions will not be reported."),
    list(id = "EarlyCheck", name = "EarlyCheck", short_name = "EarlyCheck", logo = "EarlyCheck.png",
         data_ids = c("EarlyCheck_Group1", "EarlyCheck_Group2", "EarlyCheck_Group3", "EarlyCheck_Group4"),
         description = "The NC NEXUS study seeks to assess the technical possibilities and limitations of NGS-NBS, devise and evaluate a framework to convey various types of genetic information, and develop best practices for incorporating NGS-NBS into clinical care. Early Check is a current research study that builds on NC NEXUS (Newborn Exome Sequencing for Universal Screening), open to newborns in North Carolina to provide additional, voluntary genetic screening beyond the standard newborn screening program.",
         website = "https://earlycheck.org", publication = "https://doi.org/10.1016/j.ajhg.2020.08.001",
         criteria = "Criteria are not publicly available.",
         returned_results = "Monogenic disease risk: heterozygous (likely) pathogenic variants in a gene associated with a dominant condition; two or more (likely) pathogenic variants in a gene associated with a recessive condition."),
    list(id = "GUARDIAN", name = "GUARDIAN", short_name = "GUARDIAN", logo = "Guardian.png",
         data_ids = c("Guardian_Group1", "Guardian_Group2"),
         description = "The Genomic Uniform-screening Against Rare Disease in All Newborns (GUARDIAN) study was initiated with the aim of screening a diverse population by parent-reported race and ethnicity of newborns in New York City within the context of the New York State Department of Health Newborn Screening Program.",
         website = "https://guardian-study.org/", publication = "https://jamanetwork.com/journals/jama/fullarticle/2825327",
         criteria = "1.	CFTR p. D1270N and p.R74W variants are excluded.\n2.	The minor allele frequency threshold for defining rare variants is not available.",
         returned_results = "Monogenic disease risk: heterozygous (likely) pathogenic variants in a gene associated with a dominant condition. Two or more (likely) pathogenic variants, or variants of uncertain significance co-occurred with (likely) pathogenic variants in a gene associated with a recessive condition. Hemizygous (likely) pathogenic variants in male neonates and homozygous (likely) pathogenic variants in female neonates for X-linked diseases."),
    list(id = "ACMG_Carrier", name = "ACMG Carrier Screening", short_name = "ACMG Carrier", logo = "ACMG.jpg",
         data_ids = c("ACMG_Carrier_Tier_1", "ACMG_Carrier_Tier_2", "ACMG_Carrier_Tier_3", "ACMG_Carrier_Tier_4", "ACMG_Carrier_Outside_gnomAD"),
         description = "ACMG expanded carrier screening is a genetic test that checks for a wide range of inherited recessive disorders to help people make informed reproductive decisions, especially those who are pregnant or planning a pregnancy.",
         website = "https://www.acmg.net", publication = "https://doi.org/10.1038/s41436-021-01203-z",
         criteria = "",
         returned_results = ""),
    list(id = "ACMGv3.3", name = "ACMG v3.3", short_name = "ACMG v3.3", logo = "ACMG.jpg",
         data_ids = c("ACMGv3.3"),
         description = "The ACMG secondary findings list contains genes associated with serious, actionable genetic diseases. These genes are recommended for opportunistic screening during clinical exome or genome sequencing, even if the patient's primary reason for testing was different. ",
         website = "https://www.acmg.net/", publication = "https://doi.org/10.1016/j.gim.2025.101454",
         criteria = "",
         returned_results = ""),
    list(id = "Genomics101", name = "Genomics101", short_name = "Genomics101", logo = "Genomics101.png",
         data_ids = c("Genomic101"),
         description = "Genomics101, also known as the Generation Study, is the research project carried out by Genomics England in partnership with the NHS, to conduct whole-genome sequencing on 100,0000 newborns. It aims to identify treatable rare conditions in babies earlier, support research on rare genetic conditions, and explore the benefits and risks of stroing a baby’s genome throughout their lifetime.",
         website = "https://www.genomicsengland.co.uk/blog/genomics-101-what-is-the-genome", publication = "Coming soon.")
  )

  # Project information mapping (can be replaced by JSON config)
  get_project_info <- function(project_id) {
    for (proj in featured_projects) {
      if (proj$id == project_id) {
        return(proj$description)
      }
    }
    return("Project information not available.")
  }

  # Load Preset_screening_list.txt
  preset_data <- reactive({
    req(file.exists("Preset_screening_list_GenCC_20251125_WebPage.txt"))
    df <- read.delim("Preset_screening_list_GenCC_20251125_WebPage.txt", stringsAsFactors = FALSE, header = T)
    # Clean up trailing/leading whitespace in all columns
    df <- df %>% 
      mutate(across(everything(), ~trimws(.x)))
    df
  })

  # Generate ontology bar plot
  output$ontology_chart <- renderPlotly({
    df <- preset_data()
    bar.df <- df %>% 
      mutate(Disorder_Group = if_else(grepl(";", Disorder_Group), 
                                      "Complex genetic syndromes", Disorder_Group)) %>% unique() %>% 
      filter(!is.na(Genes)) %>% 
      mutate(Project = case_when(
        grepl("ACMG_Carrier", Project) ~ "ACMG_Carrier",
        grepl("BabySeq", Project) ~ "BabySeq",
        grepl("EarlyCheck", Project) ~ "EarlyCheck",
        grepl("Guardian", Project) ~ "GUARDIAN",
        TRUE ~ Project)) %>%
      mutate(cor = paste(Genes, Disease, sep = "_")) %>% 
      select(cor, Disorder_Group, Project) %>% 
      mutate(Disorder_Group = trimws(Disorder_Group)) %>% unique() 
    
    bar.df.cor <- bar.df %>%
      group_by(Disorder_Group, Project) %>%
      summarise(n_cor = n_distinct(cor), .groups = 'drop') %>% ungroup() %>% arrange(desc(n_cor)) %>% 
      mutate(Disorder_Group = if_else(Disorder_Group == "Neuromuscular and musculoskeletal disorders",
                                      str_wrap(Disorder_Group, width = 30),
                                      str_wrap(Disorder_Group, width = 40)))
    
    sum_order = bar.df.cor %>% group_by(Disorder_Group) %>% summarise(total = sum(n_cor), .groups = 'drop') %>% arrange(desc(total))
    bar.df.cor$Disorder_Group <- factor(bar.df.cor$Disorder_Group, levels = rev(sum_order$Disorder_Group))
    
    scientific_colors <- c(
      "#A8D5BA", "#F7C59F", "#C9B6E4", "#F4A7B9", "#DCEFB7", "#F9E79F", "#D7CCC8", "#AED6F1", "#F5B7B1"
    )
    
    p <- ggplot(bar.df.cor, aes(y = Disorder_Group, fill = Project, x = n_cor)) +
      geom_bar(position = "stack", stat = "identity") +
      theme_classic(base_size = 14) +
      scale_fill_manual(values = scientific_colors) +
      theme(axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 11),
            axis.title.x = element_text(size = 13),
            axis.title.y = element_text(size = 13),
            legend.title = element_text(size = 11),
            legend.text = element_text(size = 10),
            plot.title = element_text(size = 15, face = "bold")) +
      labs(y = NULL, x = "Count", fill = "Project") +
      theme(legend.position = "bottom")
    
    ggplotly(p, tooltip = c("x", "fill")) %>%
      layout(hovermode = "y unified")
  })

  # Generate UpSet plot (interactive with plotly)
  output$upset_chart <- renderPlotly({
    df <- preset_data()
    upset.df.main <- df %>% 
      filter(!is.na(Genes)) %>% 
      mutate(Project = case_when(
        grepl("ACMG_Carrier", Project) ~ "ACMG_Carrier",
        grepl("BabySeq", Project) ~ "BabySeq",
        grepl("EarlyCheck", Project) ~ "EarlyCheck",
        grepl("Guardian", Project) ~ "GUARDIAN",
        TRUE ~ Project)) %>%
      mutate(cor = paste(Genes, Disease, sep = "_")) %>% 
      select(cor, Project) %>% unique() 
    
    # Convert to binary matrix for UpSet plot
    upset_matrix <- upset.df.main %>%
      mutate(value = 1) %>%
      pivot_wider(names_from = Project, values_from = value, values_fill = 0) %>%
      column_to_rownames(var = "cor") %>%
      as.data.frame()
    
    # Sort columns
    upset_matrix <- upset_matrix[, sort(colnames(upset_matrix), decreasing = TRUE)]
    
    # Calculate intersection sizes
    intersections <- list()
    
    # Add single sets
    for (i in seq_len(ncol(upset_matrix))) {
      set_name <- colnames(upset_matrix)[i]
      count <- sum(upset_matrix[[i]])
      intersections[[set_name]] <- count
    }
    
    # Add pairwise intersections
    if (ncol(upset_matrix) >= 2) {
      for (i in 1:(ncol(upset_matrix)-1)) {
        for (j in (i+1):ncol(upset_matrix)) {
          set1 <- colnames(upset_matrix)[i]
          set2 <- colnames(upset_matrix)[j]
          count <- sum(upset_matrix[[i]] & upset_matrix[[j]])
          if (count > 0) {
            intersections[[paste(set1, set2, sep = " ∩ ")]] <- count
          }
        }
      }
    }
    
    # Convert to dataframe and sort
    intersection_df <- data.frame(
      intersection = names(intersections),
      count = unlist(intersections),
      row.names = NULL,
      stringsAsFactors = FALSE
    )
    intersection_df <- intersection_df[order(intersection_df$count, decreasing = TRUE), ]
    
    # Create interactive plotly visualization
    p <- plot_ly(
      intersection_df,
      x = ~reorder(intersection, count),
      y = ~count,
      type = 'bar',
      marker = list(
        color = ~count,
        colorscale = list(
          c(0, '#89c2d9'),  # Light blue for low values
          c(1, '#2a6f97')   # Dark blue for high values
        ),
        showscale = TRUE,
        colorbar = list(title = "Count")
      ),
      hovertemplate = '<b>%{x}</b><br>Count: %{y}<extra></extra>'
    ) %>%
      layout(
        #title = "Gene-Disease Intersections Across Screening Lists",
        xaxis = list(
          title = "Project Intersections",
          tickangle = -45,
          tickfont = list(size = 10)
        ),
        yaxis = list(
          title = "Number of Shared Gene-Disease Pairs",
          tickfont = list(size = 10)
        ),
        margin = list(b = 120, l = 80, r = 60, t = 80),
        hovermode = "closest",
        font = list(size = 11)
      )
    p
  })


  # Render project cards on homepage
  output$project_cards <- renderUI({
    # Create 4 columns per row using Bootstrap grid
    cards <- lapply(seq_along(featured_projects), function(i) {
      project <- featured_projects[[i]]
      div(
        class = "col-md-3 col-sm-6",
        div(
          class = "project-card",
          `data-project-id` = project$id,
          div(class = "project-logo", 
              img(src = paste0("logo/", project$logo), alt = project$name, style = "max-width: 100%; max-height: 100%; object-fit: contain;")
          ),
          div(class = "project-name", project$name)
        )
      )
    })
    
    # Organize into rows of 4 cards
    row_list <- list()
    for (i in seq(1, length(cards), by = 4)) {
      row_cards <- cards[i:min(i+3, length(cards))]
      row_list[[length(row_list) + 1]] <- div(class = "row", row_cards)
    }
    do.call(tagList, row_list)
  })

  # Initialize reactive value to store selected project
  selected_project_val <- reactiveVal(NULL)

  # Handle project selection
  observeEvent(input$selected_project_id, {
    if (!is.null(input$selected_project_id)) {
      selected_project_val(input$selected_project_id)
    }
  })

  # Handle reset project selection when clicking NBSeq Details tab
  observeEvent(input$reset_project_selection, {
    selected_project_val(NULL)
  })
  
  # Handle back button in project detail panel
  observeEvent(input$back_to_projects, {
    selected_project_val(NULL)
  })

  # Store selected project
  selected_project <- reactive({
    selected_project_val()
  })

  # Output selected project title
  output$detail_project_title <- renderText({
    project_id <- selected_project()
    if (is.null(project_id)) return("")
    project <- featured_projects[sapply(featured_projects, function(p) p$id == project_id)][[1]]
    if (is.null(project)) return("")
    project$name
  })

  # Output project links
  output$project_links <- renderUI({
    project_id <- selected_project()
    if (is.null(project_id)) return(NULL)
    project <- featured_projects[sapply(featured_projects, function(p) p$id == project_id)][[1]]
    if (is.null(project)) return(NULL)
    
    tagList(
      a("🌐 Website", href = project$website, target = "_blank", class = "btn btn-sm btn-info", style = "display: block; margin: 5px 0;"),
      a("📄 Publication", href = "#", class = "btn btn-sm btn-primary", style = "display: block; margin: 5px 0;")
    )
  })

  # Output project info
  output$project_info_text <- renderText({
    project_id <- selected_project()
    if (is.null(project_id)) return("")
    get_project_info(project_id)
  })

  # Output filtering criteria
  output$filtering_criteria_text <- renderText({
    project_id <- selected_project()
    if (is.null(project_id)) return("")
    for (proj in featured_projects) {
      if (proj$id == project_id) {
        return(proj$criteria)
      }
    }
    return("")
  })

  # Output returned results
  output$returned_results_text <- renderText({
    project_id <- selected_project()
    if (is.null(project_id)) return("")
    for (proj in featured_projects) {
      if (proj$id == project_id) {
        return(proj$returned_results)
      }
    }
    return("")
  })

  # Render project detail panel
  output$project_detail_panel <- renderUI({
    project_id <- selected_project()
    
    # If no project is selected, show all project cards
    if (is.null(project_id) || project_id == "") {
      # Create project cards (same as Home page)
      cards <- lapply(seq_along(featured_projects), function(i) {
        project <- featured_projects[[i]]
        div(
          class = "col-md-3 col-sm-6",
          div(
            class = "project-card",
            `data-project-id` = project$id,
            div(class = "project-logo", 
                img(src = paste0("logo/", project$logo), alt = project$name, style = "max-width: 100%; max-height: 100%; object-fit: contain;")
            ),
            div(class = "project-name", project$name)
          )
        )
      })
      
      # Organize into rows of 4 cards
      row_list <- list()
      for (i in seq(1, length(cards), by = 4)) {
        row_cards <- cards[i:min(i+3, length(cards))]
        row_list[[length(row_list) + 1]] <- div(class = "row", row_cards)
      }
      
      return(
        div(
          class = "detail-panel-container",
          style = "margin: 40px;",
          h2("NBSeq Projects and Expanded Carrier Screening List", style = "text-align: left;"),
          p("Click on a project to view the project information and detailed screening list.", style = "text-align: left; color: #666;"),
          br(),
          do.call(tagList, row_list)
        )
      )
    }
    
    # Show screening list when project is selected
    div(
      class = "detail-panel-container",
      sidebarLayout(
        sidebarPanel(
          h3(textOutput("detail_project_title")),
          hr(),
          h4("Project Information:"),
          textOutput("project_info_text"),
          hr(),
          h4("Filtering Criteria:"),
          htmlOutput("filtering_criteria_text"),
          hr(),
          h4("Returned Results:"),
          htmlOutput("returned_results_text"),
          hr(),
          h4("Quick Links:"),
          uiOutput("project_links"),
          width = 3
        ),
        mainPanel(
          fluidRow(
            column(10,
              h3("Screening List", style = "margin-top: 0; margin-left: 20px;")
            ),
            column(2,
                   actionButton("back_to_projects", "Back to all projects", class = "btn btn-secondary btn-sm", style = "width: 100%; margin-bottom: 10px;")
            )
          ),
          br(),
          DTOutput("detail_gene_disease_table")
        )
      )
    )
  })

  # Project-specific data for the detail tab
  detail_project_data <- reactive({
    df <- preset_data()
    project_id <- selected_project()
    
    if (is.null(project_id)) return(df[0, ])
    
    # Find the project definition
    project_def <- NULL
    for (proj in featured_projects) {
      if (proj$id == project_id) {
        project_def <- proj
        break
      }
    }
    
    if (is.null(project_def)) return(df[0, ])
    
    # Filter by all data_ids for this project
    df[df$Project %in% project_def$data_ids, ]
  })

  # Filter data - now just returns the project data (filtering done in table)
  detail_filtered_data <- reactive({
    detail_project_data()
  })

  # Render gene-disease table with column filters
  output$detail_gene_disease_table <- renderDT({
    df <- detail_filtered_data()
    datatable(
      df,
      options = list(
        pageLength = 15,
        scrollX = TRUE,
        columnDefs = list(
          list(targets = "_all", searchable = TRUE)
        )
      ),
      filter = "top",
      extensions = "Buttons",
      selection = "multiple"
    )
  })

  # Save selected rows
  observeEvent(input$detail_save_button, {
    selected_rows <- input$detail_gene_disease_table_rows_selected
    if (length(selected_rows) > 0) {
      selected_data <- detail_filtered_data()[selected_rows, ]
      file_name <- paste0("selected_rows_", gsub(" ", "_", selected_project()), ".txt")
      write.table(selected_data, file_name, sep = "\t", row.names = FALSE, quote = FALSE)
      showNotification(paste("Saved selected rows to", file_name), type = "message")
    } else {
      showNotification("No rows selected.", type = "warning")
    }
  })

  # Load NBSeq_Results_WebPage.xlsx (if exists) - for Screening Results tab
  screening_data <- reactive({
    if (file.exists("NBSeq_Results_WebPage.xlsx")) {
      excel_sheets("NBSeq_Results_WebPage.xlsx")
    } else {
      character(0)
    }
  })

  # NBSeq Screening Results Tab
  output$screening_results_tabs <- renderUI({
    sheet_names <- screening_data()
    if (length(sheet_names) == 0) {
      return(tagList(
        h3("No screening results file found"),
        p("Please add NBSeq_Results_WebPage.xlsx to the application directory.")
      ))
    }
    
    # Filter out All-Positive and All-Carrier sheets
    sheet_names <- sheet_names[!sheet_names %in% c("All-Positive", "All-Carrier", "Positive-rates")]
    
    tabs <- lapply(sheet_names, function(sheet) {
      tabPanel(
        title = sheet,
        DTOutput(paste0("table_", sheet))
      )
    })
    do.call(tabsetPanel, tabs)
  })

  # Render Venn plots and positive rate plot
  output$carrier_venn_plot <- renderPlot({
    if (!file.exists("NBSeq_Results_WebPage.xlsx")) return(NULL)
    
    tryCatch({
      car.res <- read_excel("NBSeq_Results_WebPage.xlsx", sheet = "All-Carrier")
      
      # Create Venn diagram for carrier-status results
      grid.newpage()
      venn.plot <- venn.diagram(
        x = list(
          BabySeq = car.res %>% filter(Project == "BabySeq") %>% pull(Genes) %>% unique(),
          EarlyCheck = car.res %>% filter(Project == "EarlyCheck") %>% pull(Genes) %>% unique()
        ),
        category.names = c("BabySeq", "EarlyCheck"),
        filename = NULL,
        output = TRUE,
        lwd = 2,
        lty = 'solid',
        fill = c("#FF9999", "#99CCFF"),
        cex = 1.5,
        fontface = "bold",
        fontfamily = "sans",
        cat.cex = 1,
        cat.fontface = "bold",
        #cat.default.pos = "outer",
        cat.fontfamily = "sans",
        cat.pos = c(0, 65),      # 统一角度
        cat.dist = c(-0.05, -0.05)  # 负值 = 圆内
      )
      grid.draw(venn.plot)
    }, error = function(e) {
      plot.new()
      text(0.5, 0.5, "Unable to load carrier-status data", cex = 1.2)
    })
  })

  output$positive_venn_plot <- renderPlot({
    if (!file.exists("NBSeq_Results_WebPage.xlsx")) return(NULL)
    
    tryCatch({
      po.res <- read_excel("NBSeq_Results_WebPage.xlsx", sheet = "All-Positive") %>% 
        filter(Genes != "BRCA2")
      
      # Create Venn diagram for screen-positive results
      grid.newpage()
      venn.plot <- venn.diagram(
        x = list(
          BabySeq = po.res %>% filter(Project == "BabySeq") %>% pull(Genes) %>% unique(),
          BabyDetect = po.res %>% filter(Project == "BabyDetect") %>% pull(Genes) %>% unique(),
          EarlyCheck = po.res %>% filter(Project == "EarlyCheck") %>% pull(Genes) %>% unique(),
          Guardian = po.res %>% filter(Project == "Guardian") %>% pull(Genes) %>% unique(),
          BabyScreen = po.res %>% filter(Project == "BabyScreen+") %>% pull(Genes) %>% unique()
        ),
        category.names = c("BabySeq", "BabyDetect", "EarlyCheck", "GUARDIAN", "BabyScreen+"),
        filename = NULL,
        output = TRUE,
        lwd = 2,
        lty = 'solid',
        fill = c("#FF9999", "#99FF99", "#99CCFF", "#FFFF99","#cdb4db"),
        cex = 1.2,
        fontface = "bold",
        fontfamily = "sans",
        cat.cex = 0.9,
        cat.fontface = "bold",
        #cat.default.pos = "outer",
        cat.fontfamily = "sans",
        cat.pos = c(0, 0, 0, 0, 0),      # 统一角度
        cat.dist = c(0.05, 0.05, -0.05, -0.05, 0.05)  # 负值 = 圆内
      )
      grid.draw(venn.plot)
    }, error = function(e) {
      plot.new()
      text(0.5, 0.5, "Unable to load screen-positive data", cex = 1.2)
    })
  })

  output$positive_rate_plot <- renderPlotly({
    if (!file.exists("NBSeq_Results_WebPage.xlsx")) return(NULL)
    
    tryCatch({
      po.res.rate <- read_excel("NBSeq_Results_WebPage.xlsx", sheet = "Positive-rates") %>% 
        mutate(Total = c(17, 14,10, 12, 7)) %>%
        mutate(Project = factor(Project, levels = rev(Project)))
      
      p <- ggplot(po.res.rate, aes(y = Project, x = Total, fill = Project)) +
        geom_bar(stat = "identity", fill = "#cdb4db") +
        scale_x_continuous(breaks = seq(0, max(po.res.rate$Total), by = 2)) +
        theme_classic() +
        labs(y = NULL, x = "No. of mutated genes") +
        theme(
          legend.position = "none",
          axis.text.x = element_text(size = 12),
          axis.title.x = element_text(size = 12),
          axis.text.y = element_text(size = 12)
        )
      
      ggplotly(p, tooltip = c("x", "y"))
    }, error = function(e) {
      plot.new()
      text(0.5, 0.5, "Unable to load positive rate data", cex = 1.5)
    })
  })

  # Render tables for each screening results tab
  observe({
    sheet_names <- screening_data()
    if (length(sheet_names) > 0) {
      lapply(sheet_names, function(sheet) {
        output[[paste0("table_", sheet)]] <- renderDT({
          df <- read_excel("NBSeq_Results_WebPage.xlsx", sheet = sheet)
          datatable(df, options = list(pageLength = 10))
        })
      })
    }
  })

  # ============ Screening List Generator Tab ============
  
  # Store additional custom rows
  custom_rows <- reactiveVal(data.frame())
  
  # Load preset screening list for generator
  generator_preset_data <- reactive({
    file_path <- "Preset_screening_list_GenCC_20251125_WebPage.txt"
    if (file.exists(file_path)) {
      df <- read.delim(file_path, stringsAsFactors = FALSE)
      # Clean up trailing/leading whitespace in all columns
      df <- df %>% 
        mutate(across(everything(), ~trimws(.x)))
      df
    } else {
      data.frame()
    }
  })
  
  # Combined data (preset + custom rows)
  generator_combined_data <- reactive({
    preset <- generator_preset_data()
    custom <- custom_rows()
    
    if (nrow(custom) == 0) {
      preset
    } else {
      # Ensure custom rows have all columns
      for (col in names(preset)) {
        if (!col %in% names(custom)) {
          custom[[col]] <- ""
        }
      }
      rbind(preset, custom[names(preset)])
    }
  })
  
  # Generator UI
  output$generator_ui <- renderUI({
    data <- generator_preset_data()
    if (nrow(data) == 0) {
      return(tagList(
        h3("Preset Screening List not found"),
        p("Please ensure Preset_screening_list_GenCC_20251125_WebPage.txt exists.")
      ))
    }
    
    col_names <- names(data)
    
    tagList(
      div(
        style = "margin: 20px;",
        h2("Screening List Generator"),
        p("View the screening lists in FamilyRisk Database. Select gene-disease associations to download, or add customized associations.",
          tags$br(),
          "The customized association will be added to the last raw. Add first and download next.",
          style = "text-align: justified; color: #666; font-size: 18px"),
        
        hr(),
        
        # Add Custom Gene/Disease Entry above table
        div(
          class = "well",
          style = "margin-bottom: 15px; padding: 12px;",
          h4("Quick Add Custom Entry:", style = "margin-top: 0; margin-bottom: 10px;"),
          div(
            class = "row",
            lapply(head(names(generator_preset_data()), 6), function(col) {
              div(
                class = "form-group col-md-2",
                textInput(
                  paste0("custom_", col),
                  label = NULL,
                  placeholder = col,
                  width = "100%"
                )
              )
            }),
            div(
              class = "form-group col-md-2",
              actionButton("add_custom_row", "Add Row", class = "btn btn-success btn-sm", style = "width: 100%; margin-top: 25px;")
            )
          ),
          uiOutput("custom_rows_display")
        ),
        
        br(),br(),
        
        fluidRow(
          column(6, h3("Preset Screening List:", style = "margin-top: 0; margin-bottom: 0;")),
          column(3),
          column(3,
                 div(
                   downloadButton("download_data", "Download Selected", class = "btn btn-sm btn-primary", style = "width: 100%; margin-bottom: 5px;"),
                   br(),
                   downloadButton("download_all_data", "Download All", class = "btn btn-sm btn-info", style = "width: 100%;")
                 )
          ),
          style = "margin-bottom: 10px;"
        ),
        
        # Collapsible Advanced Filters panel
        div(
          class = "panel-group",
          div(
            class = "panel panel-default",
            div(
              class = "panel-heading",
              h4(
                class = "panel-title",
                style = "margin: 0;",
                a(
                  "data-toggle" = "collapse",
                  "href" = "#filterPanel",
                  "▶ Advanced Filters",
                  style = "text-decoration: none; color: #333; display: block; padding: 10px;"
                )
              )
            ),
            div(
              id = "filterPanel",
              class = "panel-collapse collapse",
              div(
                class = "panel-body",
                style = "padding: 15px;",
                fluidRow(
                  column(5,
                    h5("Project"),
                    div(
                      style = "margin-bottom: 10px;",
                      actionButton("select_all_project", "Select All", class = "btn btn-xs btn-default", style = "margin-right: 5px;"),
                      actionButton("deselect_all_project", "Deselect All", class = "btn btn-xs btn-default")
                    ),
                    uiOutput("filter_project_ui")
                  ),
                  column(5,
                    h5("Inheritance"),
                    div(
                      style = "margin-bottom: 10px;",
                      actionButton("select_all_inheritance", "Select All", class = "btn btn-xs btn-default", style = "margin-right: 5px;"),
                      actionButton("deselect_all_inheritance", "Deselect All", class = "btn btn-xs btn-default")
                    ),
                    uiOutput("filter_inheritance_ui")
                  ),
                  column(2,
                    div(
                      style = "margin-top: 25px;",
                      actionButton("apply_filters", "Apply", class = "btn btn-primary btn-sm", style = "width: 100%; margin-bottom: 5px;"),
                      actionButton("reset_filters", "Reset", class = "btn btn-secondary btn-sm", style = "width: 100%;")
                    )
                  )
                )
              )
            )
          )
        ),
        
        DTOutput("generator_table"),
        
        br(),br()
      )
    )
  })
  
  # Render the main data table with advanced filtering
  output$generator_table <- renderDT({
    # Always include custom rows, even if table is filtered
    custom <- custom_rows()
    is_filtered <- !is.null(filtered_generator_data())
    
    if (is_filtered) {
      data <- filtered_generator_data()
      # Add custom rows at the bottom
      if (nrow(custom) > 0) {
        # Ensure custom rows have all columns in correct order
        for (col in names(data)) {
          if (!col %in% names(custom)) {
            custom[[col]] <- ""
          }
        }
        data <- rbind(data, custom[names(data)])
      }
    } else {
      data <- generator_combined_data()
    }
    
    if (nrow(data) == 0) {
      return(datatable(
        data,
        options = list(pageLength = 5),
        filter = "top"
      ))
    }
    
    # Create the datatable with improved options
    # If filtered, auto-select all rows
    dt <- datatable(
      data,
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        selection = list(
          mode = 'multiple', 
          selected = if(is_filtered) seq_len(nrow(data)) else NULL,
          style = 'multi',
          bluelight = TRUE
        ),
        dom = 'lfrtip',
        autoWidth = FALSE,
        columnDefs = list(
          list(width = '150px', targets = "_all")
        ),
        language = list(
          search = "Quick search:",
          lengthMenu = "Display _MENU_ rows per page",
          info = "Showing _START_ to _END_ of _TOTAL_ entries"
        )
      ),
      filter = "top",
      extensions = c('Select', 'FixedColumns', 'FixedHeader'),
      callback = JS("
        table.on('draw', function() {
          Shiny.setInputValue('generator_table_initialized', Math.random(), {priority: 'event'});
        });
      ")
    )
    
    dt
  })
  
  # Reactive filtered data for generator table
  filtered_generator_data <- reactiveVal(NULL)
  
  # Create filter UI elements dynamically
  output$filter_project_ui <- renderUI({
    data <- generator_combined_data()
    projects <- unique(data$Project)
    projects <- projects[!is.na(projects) & projects != ""]
    
    checkboxGroupInput(
      "filter_project_selected",
      label = NULL,
      choices = sort(projects),
      selected = sort(projects)
    )
  })
  
  output$filter_inheritance_ui <- renderUI({
    data <- generator_combined_data()
    inheritance <- unique(data$Inheritance)
    inheritance <- inheritance[!is.na(inheritance) & inheritance != ""]
    
    checkboxGroupInput(
      "filter_inheritance_selected",
      label = NULL,
      choices = sort(inheritance),
      selected = sort(inheritance)
    )
  })
  
  # Select All / Deselect All for Project
  observeEvent(input$select_all_project, {
    data <- generator_combined_data()
    projects <- unique(data$Project)
    projects <- projects[!is.na(projects) & projects != ""]
    
    updateCheckboxGroupInput(session, "filter_project_selected", selected = sort(projects))
  })
  
  observeEvent(input$deselect_all_project, {
    updateCheckboxGroupInput(session, "filter_project_selected", selected = character(0))
  })
  
  # Select All / Deselect All for Inheritance
  observeEvent(input$select_all_inheritance, {
    data <- generator_combined_data()
    inheritance <- unique(data$Inheritance)
    inheritance <- inheritance[!is.na(inheritance) & inheritance != ""]
    
    updateCheckboxGroupInput(session, "filter_inheritance_selected", selected = sort(inheritance))
  })
  
  observeEvent(input$deselect_all_inheritance, {
    updateCheckboxGroupInput(session, "filter_inheritance_selected", selected = character(0))
  })
  
  # Apply filters
  observeEvent(input$apply_filters, {
    data <- generator_combined_data()
    
    # Get selected values, handle NULL case
    selected_projects <- input$filter_project_selected
    selected_inheritance <- input$filter_inheritance_selected
    
    if (is.null(selected_projects)) {
      selected_projects <- character(0)
    }
    if (is.null(selected_inheritance)) {
      selected_inheritance <- character(0)
    }
    
    # Apply filters
    filtered <- data
    
    # Determine if all projects/inheritances are selected
    all_projects <- unique(data$Project[!is.na(data$Project) & data$Project != ""])
    all_inheritance <- unique(data$Inheritance[!is.na(data$Inheritance) & data$Inheritance != ""])
    
    all_proj_selected <- length(selected_projects) == length(all_projects)
    all_inher_selected <- length(selected_inheritance) == length(all_inheritance)
    
    # If everything is selected, we still want to treat it as a filter
    # to ensure all rows are selected for download.
    if (all_proj_selected && all_inher_selected) {
      filtered <- data # Use the full dataset
    } else {
      if (length(selected_projects) > 0) {
        filtered <- filtered %>%
          filter(Project %in% selected_projects)
      }
      
      if (length(selected_inheritance) > 0) {
        filtered <- filtered %>%
          filter(Inheritance %in% selected_inheritance)
      }
    }
    
    filtered_generator_data(filtered)
    message_text <- paste0("Filters applied! Found ", nrow(filtered), " matching row(s).")
    showNotification(message_text, type = "message", duration = 3)
  })
  
  # Reset filters
  observeEvent(input$reset_filters, {
    filtered_generator_data(NULL)
    showNotification("Filters reset!", type = "message")
  })
  
  # Display custom rows info
  output$custom_rows_display <- renderUI({
    custom <- custom_rows()
    if (nrow(custom) > 0) {
      tagList(
        div(
          style = "margin-top: 10px; padding: 10px; background-color: #d4edda; border: 1px solid #c3e6cb; border-radius: 4px;",
          h5(paste0("✓ ", nrow(custom), " Custom Entr", if(nrow(custom) > 1) "ies" else "y", " Added:"), style = "color: #155724; margin: 0 0 8px 0;"),
          div(
            style = "max-height: 120px; overflow-y: auto;",
            lapply(seq_len(nrow(custom)), function(i) {
              genes <- if(!is.na(custom[i, "Genes"])) custom[i, "Genes"] else ""
              disease <- if(!is.na(custom[i, "Disease"])) custom[i, "Disease"] else ""
              
              div(
                style = "padding: 6px; font-size: 12px; color: #155724;",
                if(genes != "") tags$b(genes),
                if(disease != "" && genes != "") " - ",
                if(disease != "") disease
              )
            })
          ),
          br(),
          actionButton("clear_custom", "Clear Custom Rows", class = "btn btn-warning btn-sm")
        )
      )
    }
  })
  
  # Add custom row
  observeEvent(input$add_custom_row, {
    data <- generator_preset_data()
    col_names <- names(data)
    
    # Collect values from inputs
    new_row <- data.frame(t(sapply(col_names, function(col) {
      val <- input[[paste0("custom_", col)]]
      if (is.null(val)) "" else val
    })), stringsAsFactors = FALSE)
    names(new_row) <- col_names
    
    # Check if any field has content
    if (any(new_row != "")) {
      # Add to custom rows
      current <- custom_rows()
      new_data <- if (nrow(current) == 0) new_row else rbind(current, new_row)
      custom_rows(new_data)
      
      # Clear input fields
      for (col in col_names) {
        shinyjs::runjs(paste0("$('#custom_", col, "').val('');"))
      }
      
      showNotification("Custom row added!", type = "message")
    } else {
      showNotification("Please fill in at least one field!", type = "warning")
    }
  })
  
  # Clear custom rows
  observeEvent(input$clear_custom, {
    custom_rows(data.frame())
    showNotification("Custom rows cleared!", type = "message")
  })
  
  # Download handler for selected data
  output$download_data <- downloadHandler(
    filename = function() {
      paste0("screening_list_selected_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")
    },
    content = function(file) {
      # Check if rows are selected
      if (is.null(input$generator_table_rows_selected) || length(input$generator_table_rows_selected) == 0) {
        # If no rows selected, write empty file with headers
        empty_data <- generator_combined_data()[0, ]
        write.table(empty_data, file = file, sep = "\t", row.names = FALSE, quote = FALSE)
        return()
      }
      
      # Build the data table as it appears in the UI (with filtered data + custom rows)
      custom <- custom_rows()
      is_filtered <- !is.null(filtered_generator_data())
      
      if (is_filtered) {
        data <- filtered_generator_data()
        # Add custom rows at the bottom
        if (nrow(custom) > 0) {
          for (col in names(data)) {
            if (!col %in% names(custom)) {
              custom[[col]] <- ""
            }
          }
          data <- rbind(data, custom[names(data)])
        }
      } else {
        data <- generator_combined_data()
      }
      
      # Get selected rows indices - must be numeric and valid
      selected_indices <- input$generator_table_rows_selected
      
      if (!is.null(selected_indices) && length(selected_indices) > 0) {
        # Filter to valid indices
        valid_indices <- selected_indices[selected_indices > 0 & selected_indices <= nrow(data)]
        
        if (length(valid_indices) > 0) {
          selected_data <- data[valid_indices, ]
          write.table(selected_data, file = file, sep = "\t", row.names = FALSE, quote = FALSE)
        } else {
          # No valid indices
          empty_data <- data[0, ]
          write.table(empty_data, file = file, sep = "\t", row.names = FALSE, quote = FALSE)
        }
      } else {
        # No rows selected
        empty_data <- data[0, ]
        write.table(empty_data, file = file, sep = "\t", row.names = FALSE, quote = FALSE)
      }
    }
  )
  
  # Download all rows
  output$download_all_data <- downloadHandler(
    filename = function() {
      paste0("screening_list_all_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")
    },
    content = function(file) {
      # Build the data table as it appears in the UI (with filtered data + custom rows)
      custom <- custom_rows()
      is_filtered <- !is.null(filtered_generator_data())
      
      if (is_filtered) {
        data <- filtered_generator_data()
        # Add custom rows at the bottom
        if (nrow(custom) > 0) {
          for (col in names(data)) {
            if (!col %in% names(custom)) {
              custom[[col]] <- ""
            }
          }
          data <- rbind(data, custom[names(data)])
        }
      } else {
        data <- generator_combined_data()
      }
      
      write.table(data, file = file, sep = "\t", row.names = FALSE, quote = FALSE)
    }
  )

}

# Run the application
# Automatically open in browser
shinyApp(ui = ui, server = server, options = list(launch.browser = TRUE))
