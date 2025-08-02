library(shiny)
library(bslib)
library(DT)
library(shinyWidgets)
library(dplyr)
library(shinythemes)
library(kinship2)
library(ggplot2)
library(plotly)
library(ggthemes)
library(viridis)
library(hrbrthemes)
library(patchwork)
library(tidyverse)
library(visNetwork)
library(shinyjs)

for (i in c("select","filter","mutate","rename","left_join","slice")){
  conflicted::conflict_prefer(i,"dplyr")
}
conflicted::conflicts_prefer(stats::sd)
rm(i)

## ---- UI ----
ui <- fluidPage(
  # CSS
  tags$head(
    tags$style(HTML("
      body {
        background-color: #f5f5f5;
      }
      
      .main-container {
        padding: 20px;
        font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
        background-color: #f5f5f5;
      }
      
      .section-header {
        background: white;
        color: #000000;
        padding: 15px 25px;
        border-radius: 10px 10px 0 0;
        margin: 0;
        font-size: 18px;
        font-weight: 600;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        border: 1px solid #e0e0e0;
        border-bottom: none;
      }
      
      .section-content {
        background: white;
        border: 1px solid #e0e0e0;
        border-top: none;
        border-radius: 0 0 10px 10px;
        padding: 20px;
        box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        margin-bottom: 25px;
      }
      
      .tab-navigation {
        background: #f5f5f5;
        border-radius: 0;
        padding: 0;
        margin: 0;
        box-shadow: none;
      }
      
      .nav-tabs .nav-link {
        border-radius: 8px 8px 0 0;
        margin: 0 5px;
        padding: 15px 25px;
        font-weight: 500;
        transition: all 0.3s ease;
        border: 1px solid #e0e0e0;
        background: white;
        color: #333333;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
      }
      
      .nav-tabs .nav-link.active {
        background: #007bff;
        color: white;
        border-color: #007bff;
        box-shadow: 0 4px 12px rgba(0, 123, 255, 0.3);
      }
      
      .nav-tabs .nav-link:hover:not(.active) {
        background: #e9ecef;
        color: #495057;
        border-color: #e0e0e0;
      }
      
      .nav-tabs {
        border-bottom: none;
        margin-bottom: 20px;
      }
      
      .tab-content {
        background-color: #f5f5f5;
        padding: 0;
      }
      
      .pedigree-carousel {
        background: white;
        border-radius: 10px;
        padding: 20px;
        margin: 15px 0;
        border: 1px solid #e0e0e0;
        overflow-x: auto;
        white-space: nowrap;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
      }
      
      .pedigree-item {
        display: inline-block;
        margin: 0 15px;
        text-align: center;
        vertical-align: top;
        background: white;
        border-radius: 8px;
        padding: 15px;
        box-shadow: 0 2px 6px rgba(0,0,0,0.1);
        min-width: 300px;
        border: 1px solid #e0e0e0;
      }
      
      .pedigree-item img {
        border-radius: 8px;
        border: 2px solid #e9ecef;
        transition: transform 0.3s ease;
        max-width: 100%;
        height: auto;
        cursor: pointer;
      }
      
      .pedigree-item img:hover {
        transform: scale(1.05);
        border-color: #007bff;
      }
      
      .modal-overlay {
        position: fixed;
        top: 0;
        left: 0;
        width: 100%;
        height: 100%;
        background-color: rgba(0, 0, 0, 0.8);
        z-index: 10000;
        display: none;
        justify-content: center;
        align-items: center;
      }
      
      .modal-content {
        position: relative;
        max-width: 90%;
        max-height: 90%;
        background: white;
        padding: 20px;
        border-radius: 10px;
        box-shadow: 0 4px 20px rgba(0, 0, 0, 0.3);
      }
      
      .modal-content img {
        max-width: 100%;
        max-height: 70vh;
        object-fit: contain;
      }
      
      .modal-close {
        position: absolute;
        top: 10px;
        right: 15px;
        font-size: 24px;
        font-weight: bold;
        color: #666;
        cursor: pointer;
        background: none;
        border: none;
        padding: 5px;
      }
      
      .modal-close:hover {
        color: #000;
      }
      
      .modal-download {
        margin-top: 15px;
        text-align: center;
      }
      
      .download-btn {
        background-color: #007bff;
        color: white;
        border: none;
        padding: 10px 20px;
        border-radius: 5px;
        cursor: pointer;
        font-size: 14px;
        text-decoration: none;
        display: inline-block;
      }
      
      .download-btn:hover {
        background-color: #0056b3;
        color: white;
        text-decoration: none;
      }
      
      .pedigree-title {
        font-weight: 600;
        color: #000000;
        margin-bottom: 10px;
        font-size: 14px;
      }
      
      /* 加载overlay */
      .loading-overlay {
        position: absolute;
        top: 0;
        left: 0;
        right: 0;
        bottom: 0;
        background: rgba(255,255,255,0.9);
        display: flex;
        align-items: center;
        justify-content: center;
        z-index: 1000;
        border-radius: 10px;
      }
      
      .spinner {
        border: 4px solid #f3f3f3;
        border-top: 4px solid #007bff;
        border-radius: 50%;
        width: 40px;
        height: 40px;
        animation: spin 1s linear infinite;
      }
      
      @keyframes spin {
        0% { transform: rotate(0deg); }
        100% { transform: rotate(360deg); }
      }
      
      .dataTables_wrapper {
        margin-top: 15px;
        background: white;
        border-radius: 8px;
      }
      
      .alert-info {
        background-color: white;
        border: 1px solid #bee5eb;
        color: #0c5460;
        border-radius: 8px;
        padding: 15px;
        margin: 15px 0;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
      }
      
      h1, h3, h4, h5 {
        color: #000000 !important;
      }
      
      .nav-tabs-custom {
        background-color: #f5f5f5;
        border: none;
        box-shadow: none;
      }
    "))
  ),
  
  # Main content
  div(class = "main-container",
      # Top navigation tabs for Summary and Variant Details
      tabsetPanel(
        tabPanel("Summary", value = "risk_summary",
                 div(class = "section-header", "Family Information"),
                 div(class = "section-content",
                     div(id = "family_loading", class = "loading-overlay", style = "display: none;",
                         div(class = "spinner")
                     ),
                     uiOutput("variant_summary_title"),
                     DT::dataTableOutput("family_table")
                 ),
                 
                 div(class = "section-header",
                     "Screen-positive Neonatal Diseases and Disease carrier status",
                     div(h5("Caused by ClinVar and ACMG (likely) pathogenic variants"), class = "sub-header")
                 ),
                 div(class = "section-content",
                     div(id = "brief_loading", class = "loading-overlay", style = "display: none;",
                         div(class = "spinner")
                     ),
                     DT::dataTableOutput("brief_table")
                 ),
                 
                 div(class = "section-header", "Variant Pedigree Analysis"),
                 div(class = "section-content",
                     div(id = "pedigree_loading", class = "loading-overlay", style = "display: none;",
                         div(class = "spinner")
                     ),
                     conditionalPanel(
                       condition = "output.has_variants",
                       div(div(class = "alert-info",
                               icon("info-circle"), " Click on any pedigree plot to view in full size and download"
                       ),
                       div(class = "pedigree-carousel",
                           uiOutput("pedigree_plots_ui")
                       )
                       )
                     ),
                     conditionalPanel(
                       condition = "!output.has_variants",
                       div(class = "alert-info",
                           h5("No variants available for pedigree analysis"),
                           p("Please ensure that variant data has been loaded and processed.")
                       )
                     )
                 )
        ),
        
        tabPanel("Variant Details", value = "variant_details",
                 div(class = "section-header", "Detailed Variant Information"),
                 div(class = "section-content",
                     div(id = "detail_loading", class = "loading-overlay", style = "display: none;",
                         div(class = "spinner")
                     ),
                     DT::dataTableOutput("detail_table")
                 )
        )
      )
  ),
  
  # JavaScript for loading states and modal functionality
  tags$script(HTML("
    Shiny.addCustomMessageHandler('showSectionLoading', function(section) {
      $('#' + section + '_loading').show();
    });
    Shiny.addCustomMessageHandler('hideSectionLoading', function(section) {
      $('#' + section + '_loading').hide();
    });
    
    // Zoom in plot
    $(document).on('click', '.pedigree-item img', function() {
      var src = $(this).attr('src');
      var alt = $(this).attr('alt');
      var filename = src.split('/').pop();
      
      var modal = $('<div class=\"modal-overlay\" style=\"display: flex;\">').html(
        '<div class=\"modal-content\">' +
          '<button class=\"modal-close\">&times;</button>' +
          '<img src=\"' + src + '\" alt=\"' + alt + '\">' +
          '<div class=\"modal-download\">' +
            '<a href=\"' + src + '\" download=\"' + filename + '\" class=\"download-btn\">' +
              '<i class=\"fa fa-download\"></i> Download Image' +
            '</a>' +
          '</div>' +
        '</div>'
      );
      
      $('body').append(modal);
      
      modal.on('click', function(e) {
        if (e.target === this || $(e.target).hasClass('modal-close')) {
          modal.remove();
        }
      });
      
      $(document).on('keydown.modal', function(e) {
        if (e.keyCode === 27) {
          modal.remove();
          $(document).off('keydown.modal');
        }
      });
    });
  "))
)

# Server
server <- function(input, output, session) {
  
  #### --- Load files ----
  vep_result_file <- "/Users/xinmengliao/Documents/Project/20250710_NewbornRisk/Family/F4/F4.txt"
  ped <- "/Users/xinmengliao/Documents/Project/20250710_NewbornRisk/Family/F4/F4.ped"
  
  # Make www directory for plots
  if (!dir.exists("www")) {
    dir.create("www")
  }
  
  # Loading status
  loading_states <- reactiveValues(
    family = FALSE,
    brief = FALSE,
    pedigree = FALSE,
    detail = FALSE
  )
  
  # Loading functions
  showSectionLoading <- function(section) {
    loading_states[[section]] <- TRUE
    session$sendCustomMessage("showSectionLoading", section)
  }
  
  hideSectionLoading <- function(section) {
    loading_states[[section]] <- FALSE
    session$sendCustomMessage("hideSectionLoading", section)
  }
  
  # choose the file
  selected_vep_file <- reactive({
    return(vep_result_file)
  })
  
  # original file
  original_data <- reactive({
    if (!file.exists(selected_vep_file())) {
      return(data.frame(Message = "No result file found. Please upload and analyze a variant or gene first."))
    }
    
    df <- read.csv(selected_vep_file(), header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
    return(df)
  })
  
  #### ---- Family Table ----
  family_data <- reactive({
    showSectionLoading("family")
    
    if (!file.exists(ped)) {
      hideSectionLoading("family")
      return(data.frame(Message = "No pedigree file found"))
    }
    
    ped_data <- read.csv(ped, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    colnames(ped_data) <- c("FID", "Individual", "Father", "Mother", "Sex", "Genotype")
    
    display_data <- data.frame()
    fatherid <- ped_data %>% filter(Father == 0 & Sex == 1) %>% select(Individual) %>% unique() %>% unlist()
    motherid <- ped_data %>% filter(Mother == 0 & Sex == 2) %>% select(Individual) %>% unique() %>% unlist()
    kidsid <- ped_data %>% filter(Mother != 0 & Father != 0) %>% select(Individual) %>% unique() %>% unlist()
    
    display_data <- data.frame(
      Relationship = c(
        paste(as.character(icon("mars")), "Father"),
        paste(as.character(icon("venus")), "Mother")
      ),
      SampleID = c(fatherid, motherid),
      Sex = c("Male", "Female")
    )
    
    if (length(kidsid) > 0) {
      for (i in 1:length(kidsid)) {
        kid_sex <- ifelse(ped_data$Sex[which(ped_data$Individual == kidsid[i])] == 2, "Female", "Male")
        kid_icon <- ifelse(kid_sex == "Female", as.character(icon("baby")), as.character(icon("baby")))
        
        display_data <- rbind(display_data, data.frame(
          Relationship = paste(kid_icon, paste0("Kid", i)),
          SampleID = kidsid[i],
          Sex = kid_sex
        ))
      }
    }
    
    hideSectionLoading("family")
    return(display_data)
  })
  
  #### ---- Brief table ----
  brief_data <- reactive({
    showSectionLoading("brief")
    
    req(original_data())
    if ("Message" %in% colnames(original_data())) {
      hideSectionLoading("brief")
      return(data.frame(Message = "No data available"))
    }
    
    df <- original_data()
    
    index <- grep("Kid", colnames(df))
    kid_pattern_col <- index[grep("Pattern", colnames(df)[index])]
    kid_genotype_col <- index[grep("Genotype", colnames(df)[index])]
    kid_cols_to_keep <- colnames(df)[c(kid_genotype_col, kid_pattern_col)]
    
    brief <- df %>%
      filter(Inheritance != "Unknown") %>%
      filter(rowSums(!is.na(select(., all_of(kid_pattern_col)))) > 0) %>%
      filter(MAX_AF < 0.05 | MAX_AF == "" | is.na(MAX_AF)) %>%
      filter((grepl("Pathogenic|Likely_pathogenic", ClinVar_CLNSIG) &
                grepl("Pathogenic|Likely pathogenic", acmg_classification))) %>%
      mutate(
        ref = sapply(strsplit(variant_info, split = "_"), `[`, 3),
        alt = sapply(strsplit(variant_info, split = "_"), `[`, 4),
        FatherGenotype = mapply(function(fg, r, a) gsub("1", a, gsub("0", r, fg)), FatherGenotype, ref, alt),
        MotherGenotype = mapply(function(mg, r, a) gsub("1", a, gsub("0", r, mg)), MotherGenotype, ref, alt),
        across(all_of(kid_genotype_col),
               ~ mapply(function(kid, r, a) gsub("1", a, gsub("0", r, kid)), ., ref, alt),
               .names = "{.col}"),
        ClinVar_CLNSIG = gsub("_", " ", ClinVar_CLNSIG),
        acmg_classification = gsub("_", " ", acmg_classification),
        ClinVar_CLNSIG = sapply(strsplit(split = "&", ClinVar_CLNSIG), `[`, 1),
        ClinVar_CLNSIG = if_else(ClinVar_CLNSIG == "Pathogenic/Likely pathogenic/Pathogenic",
                                 "Pathogenic/Likely pathogenic",
                                 ClinVar_CLNSIG),
        Existing_variation = sapply(strsplit(split = "&", Existing_variation), `[`, 1)
      ) %>%
      select(Disease, Genes, Inheritance, variant_info, Existing_variation, ClinVar_CLNSIG, acmg_classification, FatherGenotype, MotherGenotype, all_of(kid_cols_to_keep)) %>%
      rename(`Existing variant` = Existing_variation, ClinVar = ClinVar_CLNSIG, ACMG = acmg_classification) %>%
      distinct()
    
    hideSectionLoading("brief")
    return(brief)
  })
  
  #### ---- Detailed table ----
  detail_data <- reactive({
    showSectionLoading("detail")
    
    req(original_data())
    if ("Message" %in% colnames(original_data())) {
      hideSectionLoading("detail")
      return(data.frame(Message = "No data available"))
    }
    
    detail_df <- original_data() %>% select(-variant_info,-INFO)
    hideSectionLoading("detail")
    return(detail_df)
  })
  
  #### ---- Pedigree plots ----
  generate_pedigree_plots <- reactive({
    showSectionLoading("pedigree")
    
    brief <- brief_data()
    df <- original_data()
    
    if (!file.exists(ped)) {
      hideSectionLoading("pedigree")
      return(list())
    }
    
    ped_data <- read.csv(ped, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    colnames(ped_data) <- c("FID", "Individual", "Father", "Mother", "Sex", "Genotype")
    familyid <- ped_data$FID %>% unique()
    
    plot_info <- list()
    
    for (i in 1:nrow(brief)) {
      variant <- brief$variant_info[i]
      row <- df %>% filter(variant_info == variant)
      
      if (nrow(row) == 0) next
      
      transcript <- tryCatch({
        unlist(strsplit(row$HGVSc, split = ":"))[2]
      }, error = function(e) "Unknown")
      
      protein <- tryCatch({
        unlist(strsplit(row$HGVSp, split = ":"))[2]
      }, error = function(e) "Unknown")
      
      gene <- row$SYMBOL %>% unique()
      if (length(gene) == 0) gene <- "Unknown"
      
      tryCatch({
        kid_gender_col <- colnames(row)[grep("Kid.*Gender",colnames(row))]
        kid_genotype_col <- colnames(row)[grep("Kid.*Genotype",colnames(row))]
        
        ped_plot <- ped_data
        for (k in 1:length(kid_genotype_col)) {
          kidid <- unlist(strsplit(kid_gender_col[k], split = "Gender_"))[2]
          genotype = row[1,kid_genotype_col[k]]
          ped_plot$Genotype[which(ped_plot$Individual == kidid)] = genotype
        }
        
        f_genotype = row[1,"FatherGenotype"]
        ped_plot$Genotype[which(ped_plot$Father == 0 & ped_plot$Mother == 0 & ped_plot$Sex == 1)] = f_genotype
        m_genotype = row[1,"MotherGenotype"]
        ped_plot$Genotype[which(ped_plot$Father == 0 & ped_plot$Mother == 0 & ped_plot$Sex == 2)] = m_genotype
        ped_plot$Father[ped_plot$Father == "0"] <- NA
        ped_plot$Mother[ped_plot$Mother == "0"] <- NA
        ped_plot$Sex <- as.integer(ped_plot$Sex)
        ped_plot <- ped_plot %>% separate(Genotype, into = c("affected","avail"),sep = "[/|]")
        ped_plot$Mother[ped_plot$Mother == "0"] <- NA
        ped_plot$affected <- as.numeric(ped_plot$affected)
        ped_plot$avail <- as.numeric(ped_plot$avail)
        
        genotype <- brief %>% filter(variant_info == variant) %>%
          pivot_longer(cols = c(FatherGenotype, MotherGenotype, kid_genotype_col), values_to = "Genotype") %>%
          select(variant_info, name, Genotype) %>%
          mutate(Individual = sapply(strsplit(split = "Genotype_",name,), `[`, 2))
        
        genotype$Individual[1:2] <- c(ped_data %>% filter(Sex == 1 & Father == 0) %>% select(Individual),
                                      ped_data %>% filter(Sex == 2 & Mother == 0) %>% select(Individual))
        
        ped_plot1 <- ped_plot %>% left_join(., genotype %>% mutate(Individual = as.character(Individual)),
                                            by = c("Individual" = "Individual"))
        
        ped <- pedigree(id = ped_plot$Individual,
                        dadid = ped_plot$Father,
                        momid = ped_plot$Mother,
                        sex = ped_plot$Sex,
                        affected = cbind(ped_plot$affected, ped_plot$avail),
                        famid = ped_plot$FID)
        
        plot_filename <- paste0("pedigree_plot_", i, ".png")
        plot_path <- file.path("www", plot_filename)
        
        png(plot_path, width = 1200, height = 1350, res = 300)
        
        ped1 <- ped[familyid]
        note <- ped_plot1 %>% mutate(combine = paste0(Individual, " (", Genotype, ")")) %>%
          summarise(note = paste(combine, collapse = "\n")) %>% pull(note)
        
        plot(ped1, col=ifelse(ped_plot$avail, 2, 1), cex=0.6)
        title(main="Pedigree analysis",cex.main = 0.7)
        mtext(paste(gene, transcript, paste0("(", protein, ")")),
              side=3, line=0.5, cex=0.5)
        mtext(note, side=1, line=0.5, cex=0.5, adj = 1, col = "black")
        dev.off()
        
        plot_info[[i]] <- list(
          filename = plot_filename,
          path = plot_path, 
          gene = gene,
          transcript = transcript,
          protein = protein,
          variant = variant
        )
      }, error = function(e) {
        print(paste("Error generating plot for variant", variant, ":", e$message))
      })
    }
    
    hideSectionLoading("pedigree")
    return(plot_info)
  })
  
  ### ---- Show Family Table ----
  output$family_table <- DT::renderDataTable({
    data <- family_data()
    if ("Message" %in% colnames(data)) {
      return(datatable(data, options = list(dom = 't'), rownames = FALSE))
    }
    
    datatable(data,
              options = list(
                pageLength = 20,
                dom = 't',
                ordering = FALSE,
                columnDefs = list(
                  list(className = 'dt-center', targets = '_all')
                )
              ),
              rownames = FALSE,
              escape = FALSE) %>%  
      formatStyle(columns = "Relationship",
                  backgroundColor = styleEqual(
                    paste(as.character(icon("male")), "Father"),
                    "#e3f2fd")) %>%
      formatStyle(columns = "Relationship",
                  backgroundColor = styleEqual(
                    paste(as.character(icon("female")), "Mother"),
                    "#fce4ec")) %>%
      formatStyle(columns = "Relationship",
                  backgroundColor = styleEqual(
                    c(paste(as.character(icon("child")), "Kid1"),
                      paste(as.character(icon("child")), "Kid2"),
                      paste(as.character(icon("child")), "Kid3"),
                      paste(as.character(icon("child")), "Kid4"),
                      paste(as.character(icon("child")), "Kid5")),
                    rep("#f3e5f5", 5)))
  }, server = FALSE)
  
  ### ---- Show Brief Table ----
  output$brief_table <- DT::renderDataTable({
    data <- brief_data()
    if ("Message" %in% colnames(data)) {
      return(datatable(data, options = list(dom = 't'), rownames = FALSE))
    }
    
    datatable(data,
              options = list(
                scrollX = TRUE,
                pageLength = 10,
                dom = 'Bfrtip',
                buttons = c('selectAll', 'selectNone'),
                lengthMenu = c(5, 10, 25, 50),
                autoWidth = TRUE
              ),
              extensions = c('Buttons', 'Select'),
              class = "display compact stripe hover",
              filter = "top",
              rownames = FALSE,
              selection = list(mode = 'multiple', target = 'row')) %>%
      formatStyle(columns = c("variant_info", "Genes"),
                  fontWeight = "bold") %>%
      formatStyle(
        columns = "ACMG",
        target = "cell",
        backgroundColor = styleEqual(
          c("Pathogenic", "Pathogenic/Likely pathogenic", "Likely pathogenic"),
          c("#dc3545", "#e74c3c", "#f39c12")
        ),
        color = styleEqual(
          c("Pathogenic", "Pathogenic/Likely pathogenic", "Likely pathogenic"),
          c("white", "white", "white")
        ),
        `white-space` = "nowrap",
        fontWeight = "bold",
        textAlign = "center",
        padding = "4px 8px",
        borderRadius = "10px",
        `font-size` = "12px",
        display = "inline-block"
      ) %>% 
      formatStyle(
        columns = "ACMG",
        target = "cell",
        backgroundColor = styleEqual(
          c("Pathogenic", "Pathogenic/Likely pathogenic", "Likely pathogenic"),
          c("#dc3545", "#e74c3c", "#f39c12")
        ),
        color = styleEqual(
          c("Pathogenic", "Pathogenic/Likely pathogenic", "Likely pathogenic"),
          c("white", "white", "white")
        ),
        `white-space` = "nowrap",
        fontWeight = "bold",
        textAlign = "center",
        padding = "4px 8px",
        borderRadius = "10px",
        `font-size` = "12px"
      ) %>%
      formatStyle(
        columns = "ClinVar",
        target = "cell",
        backgroundColor = styleEqual(
          c("Pathogenic", "Pathogenic/Likely pathogenic", "Likely pathogenic"),
          c("#dc3545", "#e74c3c", "#f39c12")
        ),
        color = styleEqual(
          c("Pathogenic", "Pathogenic/Likely pathogenic", "Likely pathogenic"),
          c("white", "white", "white")
        ),
        `white-space` = "nowrap",
        fontWeight = "bold",
        textAlign = "center",
        padding = "20px 8px",
        borderRadius = "10px",
        `font-size` = "12px"
      )
  })
  
  #### ---- Show Detailed Table ----
  output$detail_table <- DT::renderDataTable({
    data <- detail_data()
    if ("Message" %in% colnames(data)) {
      return(datatable(data, options = list(dom = 't'), rownames = FALSE))
    }
    
    datatable(data,
              options = list(
                scrollX = TRUE,
                pageLength = 15,
                dom = 'Bfrtip',
                buttons = c('copy', 'csv', 'excel'),
                lengthMenu = c(10, 15, 25, 50),
                autoWidth = TRUE
              ),
              extensions = c('Buttons'),
              class = "display compact stripe hover",
              filter = "top",
              rownames = FALSE)
  })
  
  output$has_variants <- reactive({
    data <- brief_data()
    return(!"Message" %in% colnames(data) && nrow(data) > 0)
  })
  outputOptions(output, "has_variants", suspendWhenHidden = FALSE)
  
  output$pedigree_plots_ui <- renderUI({
    plots <- generate_pedigree_plots()
    
    if (length(plots) == 0) {
      return(div(class = "alert-info", "No pedigree plots available"))
    }
    
    plot_elements <- lapply(plots, function(plot_info) {
      div(class = "pedigree-item",
          div(class = "pedigree-title",
              paste(plot_info$gene, "-", plot_info$transcript)
          ),
          if (file.exists(plot_info$path)) {
            div(
              tags$img(src = plot_info$filename, 
                       style = "max-width: 280px; height: auto;",
                       alt = paste("Pedigree plot for", plot_info$variant),
                       class = "clickable-image"),
              div(style = "margin-top: 10px;",
                  tags$a(href = plot_info$filename, 
                         download = paste0(plot_info$gene, "_", gsub("[^A-Za-z0-9]", "_", plot_info$variant), "_pedigree.png"),
                         class = "btn btn-sm btn-outline-primary",
                         style = "font-size: 12px;",
                         icon("download"), " Download")
              )
            )
          } else {
            div("Plot not available", style = "color: #6c757d; font-style: italic;")
          },
          div(style = "margin-top: 10px; font-size: 12px; color: #6c757d;",
              paste("Variant:", plot_info$variant)
          )
      )
    })
    
    return(do.call(tagList, plot_elements))
  })
}

shinyApp(ui = ui, server = server)