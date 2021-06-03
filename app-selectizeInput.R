############################################################################################################# #
# Project: A Shiny app to interact with the ChEP time course through mitotic entry (Earnshaw lab)
# Authors: Georg Kustatscher and Natalia Kochanova
# Date: May 2021
############################################################################################################# #

#### Load necessary libraries ####
library(shiny)
library(shinydashboard)
library(ggplot2)
library(data.table)
library(shinycssloaders)
library(shinyWidgets)
library(colourpicker)
library(plotly)
library(DT)

#### Load and prep the data ####

# Load the app input files 
DT9 <- fread("DT9app.csv", colClasses = c("Ratio_log2_repA_G2" = "numeric",    # Make sure these are not read in as integer
                                          "Ratio_log2_repB_G2" = "numeric", 
                                          "Ratio_log2_avg_G2"  = "numeric"))
DT9_download_version <- fread("mitoChEP_data.csv")
cl_GO <- fread("cl_GO_app.csv")

# Melt the table and prepare replicate and timepoint columns
mDT9 <- melt(DT9, measure.vars = grep("Ratio_log2", names(DT9), value = TRUE))
mDT9[, replicate := gsub("Ratio_log2_(.+)_.+$", "\\1", variable) ]
mDT9[, timepoint := gsub("Ratio_log2_.+_", "", variable) ]

# Turn timepoints into integers
mDT9[ timepoint == "G2"   , timepoint := 1 ]
mDT9[ timepoint == "0min" , timepoint := 2 ]
mDT9[ timepoint == "5min" , timepoint := 3 ]
mDT9[ timepoint == "10min", timepoint := 4 ]
mDT9[ timepoint == "15min", timepoint := 5 ]
mDT9[ timepoint == "20min", timepoint := 6 ]
mDT9[ timepoint == "25min", timepoint := 7 ]
mDT9[, timepoint := as.integer(timepoint) ]

# Make sure table is ordered by time point
setkey(mDT9, timepoint)


#### The user interface ####

ui <- dashboardPage( title = "Chromatin changes in mitosis",
                     dashboardHeader( title = tags$h4("Chromatin transactions during mitotic entry", style = "text-align: left" ), titleWidth = "100%"),
                     dashboardSidebar(
                       sidebarMenu( id = "tabs",
                                    menuItem("Home", tabName = "LandingPage", icon = icon("home")),
                                    menuItem("Help", tabName = "HowToPage", icon = icon("info-circle")),
                                    menuItem("Proteins", tabName = "T1", icon = icon("chart-area")),
                                    menuItem("Clusters", tabName = "T2", icon = icon("chart-area")),
                                    menuItem("Download", tabName = "T3", icon = icon("download")))),
                     
                     dashboardBody(
                     tags$head(tags$style(HTML('
                     .sidebar-menu li { margin-bottom: 5%; }
                     /* logo */ .skin-blue .main-header .logo { background-color: #001f3f; }
                     /* logo when hovered */ .skin-blue .main-header .logo:hover { background-color: #001f3f; }
                     /* navbar (rest of the header) */ .skin-blue .main-header .navbar { background-color: #001f3f; }
                     /* main sidebar */ .skin-blue .main-sidebar { background-color: #001f3f; font-size: 125%; }
                     /* active selected tab in the sidebarmenu */ .skin-blue .main-sidebar .sidebar .sidebar-menu .active a{ background-color: #2C7DA5; }
                     /* other links in the sidebarmenu */ .skin-blue .main-sidebar .sidebar .sidebar-menu a{ background-color: #001f3f; color: white;}
                     /* other links in the sidebarmenu when hovered */ .skin-blue .main-sidebar .sidebar .sidebar-menu a:hover{ background-color: #2C7DA5;}
                                               /* toggle button when hovered  */ .skin-blue .main-header .navbar .sidebar-toggle:hover{ background-color: NOT-SET;} '))),
    
      tabItems(
        
        tabItem(tabName = "LandingPage",
                fluidRow( column(7),
                          column(3, align = "center",
                                 img(src = "wcb_logo.png", style="padding-bottom: 5%; width: 200px"),
                                 img(src = "UoE_logo.png", style="width: 140px")),
                          column(2)),
                fluidRow( column(1),
                          column(9, align = "left", 
                                 tags$h4(tags$b("ONLINE RESOURCE")),
                                 tags$h2(tags$b("Mapping the invisible chromatin transactions of prophase chromosome remodelling"))),
                          column(2)),
                fluidRow( column(1),
                          column(9,
                                 tags$h4(tags$i("Itaru Samejima, Christos Spanos, Kumiko Samejima, Juri Rappsilber, Georg Kustatscher and William C. Earnshaw (submitted, 2021)")),
                                 br(),
                                 tags$blockquote("This is an R Shiny app designed to interactively explore the data presented in our manuscript", style = "border-color: #1f5f80; background: #abcbdb;")),
                          column(2)),
                
                fluidRow( column(1),
                          column(5,
                                 tags$p(align = "justify", style = "font-size: 125%",
                                        tags$b("SUMMARY.")," We have used a combination of chemical genetics, chromatin proteomics and imaging to map the earliest chromatin transactions during vertebrate cell entry into mitosis. Chicken DT40 Cdk1as cells undergo synchronous mitotic entry within 15 minutes following release from a 1NM-PP1-induced arrest in late G2. In addition to changes in chromatin association with nuclear pores and the nuclear envelope, earliest prophase is dominated by changes in the association of ribonucleoproteins with chromatin, particularly in the nucleolus, where pre-rRNA processing factors leave chromatin significantly before RNA polymerase I. Nuclear envelope barrier function is lost early in prophase and cytoplasmic proteins begin to accumulate on the chromatin. As a result, outer kinetochore assembly appears complete by nuclear envelope breakdown (NEBD). Most interphase chromatin proteins remain associated with chromatin until NEBD, after which their levels drop sharply. This website presents an interactive proteomic map of these chromatin transactions during mitotic entry."
                                 )),
                          column(4, img(src = "Figure_1.jpg", style = "display: block; margin-left: auto; margin-right: auto; width: 100%;")),
                          column(1)),
                
                fluidRow( column(4), 
                          column(4, align = "center", style = "font-size: 125%",
                                 br(),
                                 actionButton( inputId = "bttn_switch_X", label = "Look up a protein", icon = icon("share"), style = "color: white; background-color: #2C7DA5; border-color: #1f5f80; font-size: 100%; padding: 8px" ),
                                 actionButton( inputId = "bttn_switch_Y", label = "See how to use this app", icon = icon("info-circle"), style = "color: white; background-color: #2C7DA5; border-color: #1f5f80; font-size: 100%; padding: 8px" ),
                                 ),
                          column(4))
                ),
        
        tabItem(tabName = "HowToPage",        
                fluidRow( column(1),
                          column(6, align = "center", 
                                 tags$h2(tags$b("A brief introduction: How to use this app")),
                                 br(),
                          column(5))),
                
                fluidRow( column(1),
                          column(6, tags$p(align = "justify", style = "font-size: 125%",
                          
                          "The purpose of this app is to enable users to analyse our data on chromatin transactions during mitotic entry in a straightforward manner. 
                          There are two main options: One can search for some proteins of interest and find out how their chromatin association changes as cell enter mitotis. 
                          Alternatively, it is possible to interactively analyse the results of our hierarchical clustering analysis. 
                          Watch the brief video below to find out more, and please have a look at our manuscript for additional information.")),
                          
                          column(5)),
                
                fluidRow(
                  column(1),
                  column(6, align = "center",
                         br(),
                         HTML('<iframe src="https://player.vimeo.com/video/556647752?title=0&amp;byline=0&amp;portrait=0&amp;speed=0&amp;badge=0&amp;autopause=0&amp;player_id=0&amp;app_id=58479" width="750" height="438" frameborder="0" allow="autoplay; fullscreen; picture-in-picture" allowfullscreen title="How to use this app"></iframe>'))),
                
                fluidRow(
                  column(1),
                  column(6, align = "center", style = "font-size: 125%",
                         br(),
                         actionButton( inputId = "bttn_switch_1", label = "Look up a protein", icon = icon("share"), style = "color: white; background-color: #2C7DA5; border-color: #1f5f80; font-size: 100%; padding: 8px" ),
                         actionButton( inputId = "bttn_switch_2", label = "Explore protein clusters", icon = icon("share"), style = "color: white; background-color: #2C7DA5; border-color: #1f5f80; font-size: 100%; padding: 8px" )
                         ))
                ),
        
        tabItem(tabName = "T1",
                fluidRow( 
                  box( solidHeader = TRUE, background = "navy", width = 3, 
                       br(),
                       selectizeInput(inputId = "highlight_1",
                                      label = "Select protein(s) to highlight",
                                      choices =  DT9$plot_label,
                                      options = list( placeholder = "Search by gene name, protein name or UniProt ID"),
                                      multiple = TRUE),
                     colourInput("col_highlight_1", "Select colour", "#FF00FF", allowTransparent = TRUE),
                     actionButton("add_trace", "   Show selected", icon = icon("share")),
                     actionButton("delete_trace", "  Clear selection", icon = icon("undo"))),
                
                box( solidHeader = TRUE, background = "navy", collapsible = FALSE, width = 4, 
                     withSpinner( plotlyOutput("tSNE_plot") )),
                
                box( solidHeader = TRUE, background = "navy", width = 5,
                     withSpinner( plotlyOutput("lineplot") ), 
                     tags$hr( style = "border-color: white;" ),
                     tags$h5( tags$b( "Additional plotting options:")),
                     awesomeCheckbox(inputId = "showHistones", label = "Show core histones", value = FALSE),
                     awesomeCheckbox(inputId = "showMedian", label = "Show overall median", value = FALSE),
                     awesomeCheckbox(inputId = "showICP", label = "Show interphase chromatin proteins", value = FALSE),
                     awesomeCheckbox(inputId = "show_gene_names", label = "Show gene names for highlighted proteins", value = TRUE),
                     tags$hr( style = "border-color: white;" ),
                     radioGroupButtons(inputId = "dataset", label = "Dataset to display:", justified = TRUE,
                                       choices = c("Average" = "avg", "Replicate A" = "repA", "Replicate B" = "repB"), 
                                       selected = "avg",
                                       checkIcon = list(yes = tags$i(class = "fa fa-circle",  style = "color: steelblue"),
                                                        no = tags$i(class = "fa fa-circle-o", style = "color: steelblue"))))),
                column(width = 12,
                       infoBox( "Info", width = NULL, value = "How to find your proteins of interest",
                                subtitle = tags$p( 
                                  "Use the search mask in the top left to find your proteins of interest and display their behaviour during the mitotic time course.
                                  You can search by gene name, protein name or UniProt ID. Choose a display colour and click `Show selected`. Alternatively, it is 
                                  also possible to click on a protein (or select groups of proteins) in the t-SNE map. Note that gene and protein names are displayed when 
                                  your mouse hovers over these points. You can also zoom into these plots and download images by clicking on the -", icon("camera"), "- icon. 
                                  Our database currently covers ~2,500 chicken chromatin proteins, so it is possible that your protein of interest may not have been detected in our experiments.
                                  For more information, please have a look at our instruction video in the Help section."),
                                icon = icon("info-circle"), color = "light-blue"))
                
             ),
        
        
      tabItem(tabName = "T2",
              
              fluidRow(
                column(3, 
                       fluidRow(
                         column(12,
                                box( solidHeader = TRUE, background = "navy", width = NULL,
                                     title = "Set cluster granularity",
                                     radioGroupButtons( inputId = "clustering_height", justified = TRUE,
                                                        # label = "Cluster size",
                                                        choices = c("coarse" = "hier_cl_averag_1",
                                                                    "medium" = "hier_cl_averag_2",
                                                                    "fine" = "hier_cl_averag_3"),
                                                        selected = "hier_cl_averag_1",
                                                        individual = TRUE,
                                                        checkIcon = list(yes = tags$i(class = "fa fa-circle", style = "color: steelblue"), 
                                                                         no = tags$i(class = "fa fa-circle-o", style = "color: steelblue")))
                                     ))),
                       fluidRow(
                         column(12,
                                box( solidHeader = TRUE, background = "navy", width = NULL,
                                     title = "Select clusters to display",
                                     selectizeInput(inputId = "cluster_to_highlight", 
                                                    label = tags$p( "Option A: by ID"),
                                                    choices = unique( DT9[, hier_cl_averag_1 ] ),
                                                    options = list( placeholder = "Select clusters by their ID"),
                                                    multiple = TRUE),
                                     selectizeInput(inputId = "protein_to_highlight_in_cluster", 
                                                    label = "Option B: by protein",
                                                    choices = DT9[ !is.na( hier_cl_averag_1 ), plot_label ], 
                                                    options = list( placeholder = "Search by gene name, protein name or UniProt ID"),
                                                    multiple = TRUE)
                                     )))),
                
                column(width = 3,
                       box( solidHeader = TRUE, background = "navy", width = NULL, align = "center", withSpinner( plotlyOutput("tSNE_plot_cluster")))),
                
                conditionalPanel(condition = "output.show_cluster_plots", 
                                 column( width = 6, 
                                         box( solidHeader = TRUE, background = "navy", width = NULL, align = "center", withSpinner( plotlyOutput("cluster_lineplot")))))
                ),
              
              conditionalPanel(condition = "output.show_cluster_plots",    # Show only when some clusters are selected (this output is calculated on the server side)
                               fluidRow(
                                 box( solidHeader = TRUE,  width = 6, title = tags$b("Gene Ontology enrichment in selected clusters"),
                                      DT::dataTableOutput("GO_CC_table") ),
                                 box( solidHeader = TRUE,  width = 6, title = tags$b("Proteins belonging to selected clusters (click to show proteins in plot)"), 
                                      DT::dataTableOutput("proteins_in_cl_table") ))),
              
              column(12, 
                     infoBox( "Info", width = NULL, value = "How to explore this hierarchical clustering analysis of the proteomics timecourse",
                              subtitle = tags$p(
                              
                              "Firstly, select the `granularity` of the clustering analysis. This setting corresponds to the height at which the clustering dendrogram is cut. 
                              You can choose to display few, but relatively large clusters, many small clusters, or something in between. Secondly, use the dropdown menu options
                              to select which clusters to display. Clusters can be selected in two ways: by their ID or through the proteins that they contain. Once selected, clusters
                              will be shown on the t-SNE map and the corresponding line plots. The GO terms enriched in the selected clusters will be shown in the table at the bottom
                              left, whereas the proteins assigned to each cluster can be found in the bottom right table. Click on those proteins to make them appear in the line plots.
                              You can also zoom into these plots and download images by clicking on the -", icon("camera"), "- icon. For more information, please have a look at our
                              instruction video in the Help section."),
                              
                              icon = icon("info-circle"), color = "light-blue"))),
      
      
      tabItem( tabName = "T3", 
        fluidRow( column(6, align = "center", 
                         tags$h3(tags$b("Download section"))),
                  column(6)),
        
        fluidRow( column(6, align = "center", 
                         tags$h4(tags$i("Click on the buttons below to download the data used to populate this app"))),
                  column(6)),
              
        fluidRow( column(6, align = "center",
                         br(),
                         downloadButton(outputId = "data_download_button_1", label = "Chromatin proteomics data (incl. cluster assignments)", 
                                        style = "color: white; background-color: #2C7DA5; border-color: #1f5f80; font-size: 125%; padding: 12px" )),
                  column(6)),
        
        fluidRow( column(6, align = "center",
                         br(),
                         downloadButton(outputId = "data_download_button_2", label = "Gene Ontology enrichment analysis of clusters",
                                        style = "color: white; background-color: #2C7DA5; border-color: #1f5f80; font-size: 125%; padding: 12px" )),
                  column(6)))
      )))


#### The server side ####

server <- function(input, output, session) {
  
  #### Landing page ####
  
  # If user clicks on first action button, change to tab for individual proteins
  observeEvent( input$bttn_switch_X, {
    updateTabItems(session, "tabs", "T1")
  })
  
  # If user clicks on second action button, change to tab for help page
  observeEvent( input$bttn_switch_Y, {
    updateTabItems(session, "tabs", "HowToPage")
  })
  
  
  
  
  
  
  
  # If user clicks on first action button, change to tab for individual proteins
  observeEvent( input$bttn_switch_1, {
    updateTabItems(session, "tabs", "T1")
  })
  
  # If user clicks on second action button, change to tab for cluster analysis
  observeEvent( input$bttn_switch_2, {
    updateTabItems(session, "tabs", "T2")
  })
  
  
  #### Tab 1 ####
  
  # Create an object that will hold the selected proteins as a reactive value
  rv_prots <- reactiveValues()  
  
  # Reactive updating of selected proteins after user presses add_trace action button
  observeEvent( input$add_trace , {                     # When add_trace button is triggered
    rv_prots$highlight_1 <- input$highlight_1           # Transfer the currently selected proteins to the reactive object
    rv_prots$col_highlight_1 <- input$col_highlight_1   # And the same for the selected colour
    })
  
  # Reactive updating of selected proteins after user presses delete_trace action button
  observeEvent( input$delete_trace , {                  # When delete_trace button is triggered
    rv_prots$highlight_1 <- NULL                        # Set the selected proteins to NULL
    updateSelectizeInput( session, inputId = "highlight_1", selected = character(0))  # And also update the selected proteins in the drop down menu
    })
  
  # Select data by clicking on individual points in the t-SNE map
  observe({
    t_SNE_click_data <- event_data("plotly_click", source = "t_SNE_plot")   # If user clicks on protein in t-SNE plot, record which point that is (note: the key variable in the t_SNE plot corresponds to that proteins ID)
    rv_prots$highlight_1 <- t_SNE_click_data$key                            # Update the selected proteins accordingly
    rv_prots$col_highlight_1 <- isolate( input$col_highlight_1)             # Also update the currently selected colour, but *isolate* it because plots should only be re-drawn if user clicks add_trace
    updateSelectizeInput( session, inputId = "highlight_1", selected = t_SNE_click_data$key )  # And also update the selected proteins in the drop down menu
    })
   
  # Select data by selecting multiple points in the t-SNE map
  observe({
    t_SNE_selected_data <- event_data("plotly_selected", source = "t_SNE_plot")   # If user selects proteins in t-SNE plot, record which points these are (note: the key variable in the t_SNE plot corresponds to that proteins ID)
    rv_prots$highlight_1 <- t_SNE_selected_data$key                               # The rest is the same as for individual points
    rv_prots$col_highlight_1 <- isolate( input$col_highlight_1)            
    updateSelectizeInput( session, inputId = "highlight_1", selected = t_SNE_selected_data$key )
  })
  
  # Create the plotly line plot output
  output$lineplot<-renderPlotly({
  
    p <- ggplot( mDT9[ replicate == input$dataset ],
                 aes(x = timepoint, y = value, group = plot_label, text = hovertext, label = Simple_gene_name))+
          geom_line(alpha = 0.05, size = 0.25)+
          { if( input$show_gene_names ) coord_cartesian(xlim = c(1,8)) }+                                                                                # If this selection is enabled, I need to make space to add the gene names
          { if( input$showICP ) geom_line( data = mDT9[ replicate == input$dataset & human_ICP > 0.7 ], size = 0.25, alpha = 0.3, colour = "magenta") }+
          { if( input$showMedian ) geom_line( data = mDT9[ replicate == input$dataset, .(value = median(value, na.rm = TRUE )), by = timepoint], aes( group = 1, label = NULL, text = NULL), size = 0.5, colour = "seagreen3") }+
          { if( input$showHistones ) geom_line( data = mDT9[ replicate == input$dataset & core_histone ], size = 0.25, colour = "deepskyblue2") }+
          { if( !is.null(rv_prots$highlight_1) ) geom_line(data = mDT9[ plot_label %in% rv_prots$highlight_1 & replicate == input$dataset ], colour = rv_prots$col_highlight_1, size = 0.5) }+
          { if( !is.null(rv_prots$highlight_1) & input$show_gene_names) geom_text(data = mDT9[ plot_label %in% rv_prots$highlight_1 & replicate == input$dataset & timepoint == 7 ], colour = rv_prots$col_highlight_1, size = 3, aes(x = 7.1)) }+
          scale_x_continuous( breaks = 1:7, labels = c("G2", "0 min", "5 min", "10 min", "15 min", "20 min", "25 min"))+
          scale_y_continuous( limits = c(-2.75, 4 ), breaks = seq(-6,6,2))+
          ylab("Fold-change vs. G2 [log2]")+
          theme(line = element_line(size = 0.5), panel.border = element_rect(colour = "black", fill = NA),
                panel.grid = element_blank(), axis.text=element_text(size=8), axis.title.x=element_blank(),
                axis.title.y = element_text( size=9, face = "bold"))

    ggplotly(p, tooltip = "text") %>%
      style(textposition = "right") %>%     # This corrects for the fact that ggplotly does not understand hjust in geom_text
      config(displaylogo = FALSE,
             modeBarButtonsToRemove = c("zoomIn2d", "zoomOut2d", "pan2d", "autoScale2d", "hoverCompareCartesian", "toggleSpikelines"),
             toImageButtonOptions = list( width = 1400, height = 900 ))
    
  }) %>% bindCache( input$dataset, input$show_gene_names,  input$showICP, input$showMedian , input$showHistones, rv_prots$highlight_1, rv_prots$col_highlight_1 )
  
  # Create the plotly t-SNE output
  output$tSNE_plot<-renderPlotly({
    
    p_tSNE <- ggplot(data = DT9, 
                     aes(x = tSNE_dim_1, y = tSNE_dim_2, text = hovertext, label = Simple_gene_name, key = plot_label))+
                geom_point(alpha = 0.3, size = 0.8)+
                { if( input$showICP ) geom_point(data = DT9[ human_ICP > 0.7 ], size = 1, alpha = 0.5, colour = "magenta") }+
                { if( input$showHistones ) geom_point(data = DT9[ core_histone == TRUE ], size = 1.5, alpha = 0.7, colour = "deepskyblue2") }+
                { if( !is.null(rv_prots$highlight_1) ) geom_point(data = DT9[ plot_label %in% rv_prots$highlight_1 ], colour = rv_prots$col_highlight_1, size = 1) }+
                { if( !is.null(rv_prots$highlight_1) & input$show_gene_names) geom_text(data = DT9[ plot_label %in% rv_prots$highlight_1 ], colour = rv_prots$col_highlight_1, size = 3, nudge_y = 2) }+
                xlab("tSNE dimension 1")+
                ylab("tSNE dimension 2")+
                theme(line = element_line(size = 0.5), panel.border = element_rect(colour = "black", fill = NA), axis.ticks = element_blank(),
                       panel.grid = element_blank(), axis.text=element_blank(), axis.title = element_text(size=9, face = "bold"))
    
    ggplotly(p_tSNE, source = "t_SNE_plot", tooltip = "text") %>% 
      config(displaylogo = FALSE,
             modeBarButtonsToRemove = c("zoomIn2d", "zoomOut2d", "pan2d", "autoScale2d", "hoverCompareCartesian", "toggleSpikelines"),
             toImageButtonOptions = list( width = 1400, height = 1100 ))  
    
  })
  
  
  #### Tab 2 ####
  
  # I need two reactive values that control this page: the chosen clustering height and the chosen cluster IDs
  rv <- reactiveValues()  # Creates an object that will hold these values

  # Reactive updating of clustering height: if the user changes the clustering height via picker input, update the clustering height accordingly
  observeEvent(input$clustering_height, {
    rv$current_cl_height <- input$clustering_height    # Update the clustering cut-off
    rv$current_clusters <- NULL                        # And re-set the clusters to display
    })
  
  ## Reactive updating of clusters shown ##
  
    # (1/2) If user selects clusters to show via Option A (cluster ID dropdown menu), update the corresponding reactive value accordingly
    observeEvent(input$cluster_to_highlight, { 
      rv$current_clusters <- input$cluster_to_highlight                                                    # Updates the chosen cluster values for plotting
      updateSelectizeInput( session, inputId = "protein_to_highlight_in_cluster",  selected = character(0) )  # And re-set Option B, so that previously selected proteins become deselected
      rv$selected_proteins_dropdown <- NULL                                                                # And clear any protein selected via option B previously
      })
    
    # (2/2) If user selects cluster to show via Option B (cluster ID inferred from protein selection), update the corresponding reactive value accordingly
    observeEvent(input$protein_to_highlight_in_cluster, {
      cl_containing_these_prots <- DT9[ plot_label %in% input$protein_to_highlight_in_cluster, unique( get(input$clustering_height) ) ]
      rv$current_clusters <- cl_containing_these_prots                                           # Update the reactive value that will affect the plot
      updateSelectizeInput( session, inputId = "cluster_to_highlight", selected = character(0)  )   # At the same time, we need to re-set the menu for Option A
      rv$selected_proteins_dropdown <- DT9[ plot_label %in% input$protein_to_highlight_in_cluster, Majority_protein_IDs ]   # Also create a new reactive variable holding these proteins, so they can be added automatically to the line plot
      })

  # Clearing of plot when everything is de-selected from the dropdown menu options. Note that de-selection in the dropdown menu changes the input values of
  # the dropdown menus to NULL. This doesn't trigger an observeEvent, therefore I force the clearance with the following function.
  observe({
    if( is.null(input$cluster_to_highlight) & 
        is.null(input$protein_to_highlight_in_cluster) 
        ) {
      rv$current_clusters <- NULL      # Remove all selected clusters if inputs are NULL and signal didn't come through level buttons
      }
    })
    
    
  ## Dropdown menu behaviour ## 
    
    # If users changes the clustering cut-off we need to update the choices of clusters and proteins that are
    # available for selection from the dropdown menus. Note that this doesn't actually change the value / selection of these inputs, so it doesn't trigger a reactive event
    observeEvent(rv$current_cl_height, {
      if(            rv$current_cl_height == "hier_cl_averag_1" ) {      # If user selects the first cut-off
        choices <- DT9[, sort(unique(hier_cl_averag_1)) ]                # These will be the choices they should have to selected clusters by ID
        choices_prots <- DT9[ !is.na( hier_cl_averag_1 ), plot_label ]   # These should be the choices they should have to select clusters by protein
        
      } else if(     rv$current_cl_height == "hier_cl_averag_2") {
        choices <- DT9[, sort(unique(hier_cl_averag_2)) ]
        choices_prots <- DT9[ !is.na( hier_cl_averag_2 ), plot_label ]
        
      } else {
        choices <- DT9[, sort(unique(hier_cl_averag_3)) ]
        choices_prots <- DT9[ !is.na( hier_cl_averag_3 ), plot_label ]
      }
      
      updateSelectizeInput( session, inputId = "cluster_to_highlight", choices = choices )    # If cluster height changes, update dropdown menu choices for cluster IDs
      updateSelectizeInput( session, inputId = "protein_to_highlight_in_cluster",             # If cluster height changes, update dropdown menu choices for selection via protein IDs
                         choices = choices_prots )
    })
    
    
  ## Create the plots ##
    
    # Show the selected clusters in a line plot
    output$cluster_lineplot<-renderPlotly({
      
    p2 <- ggplot(mDT9[ replicate %in% c("repA", "repB") ],
                 aes(x = timepoint, y = value, group = Majority_protein_IDs, colour = as.factor(get( rv$current_cl_height )),
                     text = hovertext, label = plot_label))+
          facet_wrap(facets = "replicate", labeller = as_labeller( c("repA" = "replicate A", "repB" = "replicate B") ))+
          scale_x_continuous( breaks = 1:7, labels = c("G2", "0 min", "5 min", "10 min", "15 min", "20 min", "25 min"))+
          scale_y_continuous( limits = c(-2.75, 4 ), breaks = seq(-6,6,2))+
          { if( !is.null(rv$current_clusters)) geom_line( data = mDT9[ replicate %in% c("repA", "repB") ][ get( rv$current_cl_height ) %in% rv$current_clusters ], size = 0.25, alpha = 0.5)}+                     # If clusters are selected, plot them
          { if( !is.null(rv$selected_proteins_dropdown)) geom_line( data = mDT9[ replicate %in% c("repA", "repB")][ Majority_protein_IDs %in% rv$selected_proteins_dropdown ], size = 0.5, colour = "navy")}+     # Highlights the protein selected through the dropdown menu
          { if( !is.null(input$proteins_in_cl_table_rows_selected)) geom_line( data = mDT9[ replicate %in% c("repA", "repB")][ Majority_protein_IDs %in% selected_proteins() ], size = 0.35, colour = "black")}+   # If the user has clicked on a protein in the table, add that protein as a line ( the selected_protein() is a reactive expression containing the respective protein IDs)
          scale_colour_discrete()+
          ylab("Fold-change vs. G2 [log2]")+
          theme(line = element_line(size = 0.5), panel.border = element_rect(colour = "black", fill = NA), title = element_text(size = 9, face = "bold"), 
                panel.grid = element_blank(), axis.text=element_text(size=8),  axis.title.x=element_blank(), axis.title.y=element_text(size=9, face = "bold"), strip.background = element_blank(),  
                legend.position = "right", legend.title = element_blank())
    
    ggplotly(p2, tooltip = "hovertext") %>%
      config(displaylogo = FALSE,
             modeBarButtonsToRemove = c("zoomIn2d", "zoomOut2d", "pan2d", "autoScale2d", "hoverCompareCartesian", "toggleSpikelines"),
             toImageButtonOptions = list( width = 1400, height = 900 )) %>%
      layout( legend = list( bgcolor = "white", bordercolor = "white", borderwidth = 1, y = 0.5, yanchor = "middle"))
    
    })
  
    
    # Show the t-SNE plot
    
    # Output the base plotly plot without any highlights
    output$tSNE_plot_cluster <- renderPlotly({
      
      p_tSNE_cluster <- ggplot(data = DT9, aes(x = tSNE_dim_1, y = tSNE_dim_2, text = hovertext, label = plot_label))+
        geom_point(alpha = 0.3, size = 0.8)+
        { if( !is.null(rv$current_clusters)) geom_point( data = DT9[ get( rv$current_cl_height ) %in% rv$current_clusters ], 
                                                         aes(colour = as.factor(get( rv$current_cl_height ))), size = 1, alpha = 0.6)}+                      # If clusters are selected, plot them
        xlab("tSNE dimension 1")+
        ylab("tSNE dimension 2")+
        xlim(-55, 49)+
        ylim(-53, 81)+
        theme(line = element_line(size = 0.5), panel.border = element_rect(colour = "black", fill = NA), axis.ticks = element_blank(),
              panel.grid = element_blank(), axis.text=element_blank(), axis.title = element_text(size=9, face = "bold"), legend.position = "none")
      
      ggplotly(p_tSNE_cluster, tooltip = "hovertext") %>%
        config(displaylogo = FALSE,
               modeBarButtonsToRemove = c("zoomIn2d", "zoomOut2d", "pan2d", "autoScale2d", "hoverCompareCartesian", "toggleSpikelines"),
               toImageButtonOptions = list( width = 1400, height = 1100 ))  
      
    })
    
  
  ## Create the output tables ##
    
  # Show a table with the enriched GO CC terms
  output$GO_CC_table <- DT::renderDataTable(
    DT::datatable( cl_GO[ cut_height == rv$current_cl_height & `Cluster ID` %in% rv$current_clusters, 
                          .(`Cluster ID`, Aspect, `GO ID`, `GO Name`, `p-value<br>-log10`) ],
                   options = list(pageLength = 10, columnDefs = list(list(className = 'dt-center', targets = "_all")), order = list(list(4, 'desc'))),
                   escape = FALSE, rownames = FALSE, selection = "single") %>%
    formatStyle("p-value<br>-log10", background = styleColorBar(cl_GO$`p-value<br>-log10`, 'lightblue')))

  # Show a table with the protein IDs in each cluster
  output$proteins_in_cl_table <- DT::renderDataTable(
        DT::datatable( DT9[ get(rv$current_cl_height) %in% rv$current_clusters,
                            .(`Protein ID` = Main_Uniprot_ID,
                              `Protein name` = Simple_protein_name, 
                              `Gene name` = Simple_gene_name,
                              `Cluster ID` = get(rv$current_cl_height)) ],
                   options = list(pageLength = 10, columnDefs = list(list(className = 'dt-center', targets = "_all", color = "black"))),
                   rownames = FALSE ))
  
  # Create a reactive object that contains the IDs of the proteins selected by the user through clicking in the proteins_in_cl_table output table
  # The object input$proteins_in_cl_table_rows_selected is an internally created vector with the row indices of that table
  selected_proteins <- reactive({
    DT9[ get(rv$current_cl_height) %in% rv$current_clusters ][ input$proteins_in_cl_table_rows_selected , Majority_protein_IDs ]
  })
  
  ## Various other stuff ##
  
    # Create a variable that is TRUE when some clusters are selected, so that the respective plots can be shown conditionally
    output$show_cluster_plots <- reactive({ !is.null( rv$current_clusters ) })
    outputOptions(output, "show_cluster_plots", suspendWhenHidden = FALSE)  # This needs to be set so that the browser knows about what's happening
  
  #### Download tab ####
    
  # Create the first download 
  output$data_download_button_1 <- downloadHandler(
    filename = function() { paste("mitoChEP-data-", Sys.Date(), ".csv", sep="") },
    content = function(file) { fwrite( DT9_download_version, file) }
  )
    
  # Create the second download 
  output$data_download_button_2 <- downloadHandler(
    filename = function() { paste("mitoChEP-cluster-GO-data-", Sys.Date(), ".csv", sep="") },
      content = function(file) { fwrite( cl_GO, file) }
    )
    
}

shinyApp(ui, server)