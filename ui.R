



# library(shinyWidgets)
library(shiny)
# library(thematic)
library(shinythemes)
library(plotly)

# functions for UI --------------------------------------------------------








# UI ----------------------------------------------------------------------

# navbarPage
# tabBox

ui <- fluidPage(
  
  shinythemes::themeSelector(),
  
  titlePanel('SLE Bulk RNA-Seq data analysis and visualization',
             windowTitle = 'SLE expression'
             ),
  
  sidebarLayout(
    sidebarPanel(
      textInput(inputId = 'geneSymbol', 
                label = 'Input Gene Symbol',
                value = 'IRF4'
                ),
      width = 3
    ),
    
    mainPanel(h2('DEGs across 33 GSE datasets'),
              h4('note that: not all 33 datasets exists!'),
              # tabsetPanel(
              #   tabPanel("Plot", plotOutput("plot1")),
              #   tabPanel("Table", tableOutput("table1"))
              # )
              )
  ),
  
  fluidRow(
    column(8,
           # plotOutput("plot1"),
           plotlyOutput('plot1_plotly'),
           offset = 3)
  ),
  
  hr(),
  
  fluidRow(
    column(3,
           h5('Select different GSE'),
           selectInput(inputId='GSEnum', label = 'GSE_num',
                       choices = chioce1
                       )
           ),
    column(8,
           tabsetPanel(
             tabPanel('Plot', plotOutput('plot2')),
             tabPanel('Summary', textOutput('text1'))
           )
           # plotOutput('plot2')
           )
  ),
  
  hr(),
  
  fluidRow(
    column(10,
           h4('Expression by Tissue / Cell Type'),
           plotOutput('tissue_plot'), offset = 2,
           plotOutput('cell_plot')
    )
  )
  
)












