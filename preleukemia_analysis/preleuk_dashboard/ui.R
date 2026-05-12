#R UI

#load shiny
library(shiny)

#type of interface
fluidPage(
  #app title
  titlePanel("Pre-Leukemic Stem Cell Analysis Dashboard"),
  
  #main layout with different sections
  navbarPage("Project Sections",
             
            #first tap with umap
             tabPanel("Interactive UMAP",
                      sidebarLayout(
                        #sidebar with controls
                        sidebarPanel(
                          h4("UMAP Controls"),
                          #dropdown
                          selectInput("celltype_select", 
                                      "Highlight Cell Type:", 
                                      choices = NULL,
                                      selected = "All",
                                      multiple = TRUE)
                        ),
                        #showing umap
                        mainPanel(
                          h4("UMAP of Hematopoietic Lineages"),
                          plotOutput("umapPlot", height = "600px")
                        )
                      )
             ),
             
              #second tab with all the other plots
             tabPanel("Analysis Plots",
                      sidebarLayout(
                        #sidebar with controls
                        sidebarPanel(
                          h4("Plot Selection"),
                          #drop down
                          selectInput("plot_select", 
                                      "Choose a Plot to Display:",
                                      choices = c(
                                                  "UMAP Integrated",
                                                  "Cell Abundance",
                                                  "UMAP Cell Type",
                                                  "Macrostates",
                                                  "Fate Probabilities",
                                                  "Metabolic Pathway Activity",
                                                  "Pseudotime Gene Dynamics",
                                                  "PLPS Survival Plot",
                                                  "Stem11 Survival Plot"
                                                  ))
                        ),
                        #main panel
                        mainPanel(
                          h4("Saved Analysis Figures"),
                          #shows image
                          imageOutput("analysisPlot")
                        )
                      )
             )
  )
)