#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
# Adapted from:
# Grolemund, G. (2015). Learn Shiny - Video Tutorials. URL:https://shiny.rstudio.com/tutorial/

# Define UI for application that draws a histogram
ui <- fluidPage(

  # Page title
  titlePanel("rnaDualSeq: differentially expressed genes within host-pathogen infection studies"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(

      tags$p("This is a simple Shiny App that is part of the rndDualSeq in R."),
      # br() element to introduce extra vertical spacing ----
      br(),

      tags$b("Description: The purpose of rndDualSeq is to facilitate package
          can the read counts from genes of either host or pathogen and reports
          the differentially expressed genes visually via volcano plots.
          It takes .rda files and perform analysis."),

      # br() element to introduce extra vertical spacing ----
      br(),
      br(),

      # input
      tags$p("Instructions: Below input read counts of host or pathogen with
            the phenotype file.
            Then press 'Run'. Navigate the right to explore the results."),

      # br() element to introduce extra vertical spacing ----
      br(),

      # inputs
      #shinyalert::useShinyalert(force = TRUE),  # Set up shinyalert

      uiOutput("tab1"),
      actionButton(inputId = "data1",
                   label = "Example Dataset Details"),
      fileInput(inputId = "file1",
                label = "Select a host dRNA-seq tab-delimited input file with read counts to visualize. File should be in .rda format",
                accept = c(".rda")),
      fileInput(inputId = "file2",
                label = "Select a pathogen dRNA-seq tab-delimited input file with read counts to visualize. File should be in .rda format",
                accept = c(".rda")),
      fileInput(inputId = "file3",
                label = "Select a phenotype file that contains names of the columns headers corresponding
                          to read counts file, and groups of time/condition. File should be in .rda format",
                accept = c(".rda")),
      selectInput("visMethod", "Choose a host or pathogen:",
                  choices = c("Host", "Pathogen")),

      # br() element to introduce extra vertical spacing ----
      br(),

      # actionButton
      actionButton(inputId = "button1",
                   label = "RUN")
    ),

    # Main panel for displaying outputs ----
    mainPanel(
      # Output: Tabset w/ plot, summary, and table ----
      plotOutput('volc')

    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  dataList <- eventReactive(eventExpr = input$button1, {
    print("test test")
    output$volc <- renderPlot({

      inFile1 <- input$file1
      inFile2 <- input$file2
      inFile3 <- input$file3
      n <- input$n

      if (is.null(inFile1) || is.null(inFile2) || is.na(n)) {
        return(NULL)
      }
      output$before <-
      host_counts <- get(load(inFile1$datapath))
      pathogen_counts <- get(load(inFile2$datapath))
      phenotype <- get(load(inFile3$datapath))
      norm <- rnaDualSeq::norm_TMM(host_counts, phenotype)
      de <- rnaDualSeq::identifyDE(norm, phenotype)
      return(rnaDualSeq::volcanoPlot("Host", de))

    })


  })
  # output$volc <- renderPlot({
  #
  #   inFile1 <- input$file1
  #   inFile2 <- input$file2
  #   inFile3 <- input$file3
  #   n <- input$n
  #
  #   if (is.null(inFile1) || is.null(inFile2) || is.na(n)) {
  #     return(NULL)
  #   }
  #
  #   host_counts <- get(load(inFile1$datapath))
  #   pathogen_counts <- get(load(inFile2$datapath))
  #   phenotype <- get(load(inFile3$datapath))
  #   return (rnaDualSeq::norm_TMM(host_counts, phenotype))
  # })
}

# Run the application
shinyApp(ui = ui, server = server)
