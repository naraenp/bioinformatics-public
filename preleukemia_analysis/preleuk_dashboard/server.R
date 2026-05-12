#R Server

#load libraries
library(shiny)
library(ggplot2)
library(dplyr)

#load shiny data
umap_data <- read.csv("data/umap_data_for_shiny.csv", row.names = 1)

#get the unique cell types
cell_type_choices <- c("All", unique(umap_data$celltype))

#server function
function(input, output, session) {
  
  #input selection
  updateSelectInput(session, "celltype_select", choices = cell_type_choices, selected = "All")
  
  #reactive for plot data
  umap_filtered_data <- reactive({
    if ("All" %in% input$celltype_select) {
      umap_data %>% mutate(highlight = "selected")
    } else {
      umap_data %>% mutate(highlight = ifelse(celltype %in% input$celltype_select, "selected", "unselected"))
    }
  })
  
  #plot renderer
  output$umapPlot <- renderPlot({
    ggplot(umap_filtered_data(), aes(x = UMAP_1, y = UMAP_2, color = celltype, alpha = highlight)) +
      geom_point(size = 1.5) +
      scale_alpha_manual(values = c("selected" = 1.0, "unselected" = 0.1), guide = "none") +
      labs(title = "Interactive UMAP of Cell Clusters", color = "Cell Type") +
      theme_minimal() +
      theme(legend.position = "bottom", text = element_text(size = 14))
  })
  
  #renders the image
  output$analysisPlot <- renderImage({
    
    #makes sure the plots are pointing to the right places
    plot_path <- switch(input$plot_select,
                        "UMAP Integrated" = "plots/UMAP_Integrated.png",
                        "Cell Abundance" = "plots/cellabundance_plot.png",
                        "UMAP Cell Type" = "plots/umapumap_celltype.png",
                        "Macrostates" = "plots/macrostates.png",
                        "Fate Probabilities" = "plots/fate_probabilities.png",
                        "Metabolic Pathway Activity" = "plots/invlogp_metabolicpathways.png",
                        "Pseudotime Gene Dynamics" = "plots/umapdpt_pseudotime.png",
                        "PLPS Survival Plot" = "plots/PLPS_survival_plot.png",
                        "Stem11 Survival Plot" = "plots/Stem11_survival_plot.png"
    )
    
    #file name list
    list(src = plot_path,
         contentType = 'image/png',
         width = 600,
         height = 500,
         alt = "selected plot")
  }, deleteFile = FALSE)
}