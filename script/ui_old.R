

library(shiny)

# UI section
shinyUI(fluidPage(
  titlePanel("RNA-Seq data collection and visualization"),
  #using fluidRow to mage the page
  fluidRow(
    column(2, wellPanel(
      textInput("text",label = h4("Gene Name"), value = "MYC"),
      helpText("Please enter the gene name, eg. MYC"),
      submitButton("Submit",icon("Submit"))
    )
    ),
    column(10,
           h3("Bulk RNA-Seq"),
           h4("Forest Plot"),
           plotOutput("plot27"))
  ), 
  fluidRow(
    column(4,plotOutput("plot31"),offset = 2)
  ),
  fluidRow(
    column(5,h4("Expression by Tissue / Cell Type"),offset = 2)
  ),
  fluidRow(
    column(9,plotOutput("plot32"),offset = 2)
  ),
  fluidRow(
    column(9,plotOutput("plot33"),offset = 2)
  ),
  fluidRow(
    column(5,h4("Cell Types"),offset = 2)
  ),
  fluidRow(
    column(4,plotOutput("plot1"),offset = 2),
    column(6,plotOutput("plot28"))
  ),#get the side bar & first 2 plots in the first row
  fluidRow(
    column(4,plotOutput("plot2"),offset = 2),
    column(6,plotOutput("plot29"))
  ), # get the 3rd  & 4th plots in the second row
  fluidRow(
    column(4,plotOutput("plot10"),offset = 2),
    column(6,plotOutput("plot30"))
  ),
  fluidRow(
    column(5,plotOutput("plot7"),offset = 2),
    column(5,plotOutput("plot8"))
  ),
  fluidRow(
    column(5,plotOutput("plot6"),offset = 2),
    column(5,plotOutput("plot4"))
  ),
  fluidRow(
    column(5,plotOutput("plot9"),offset = 2),
    column(5,plotOutput("plot3"))
  ),
  fluidRow(
    column(5,plotOutput("plot5"),offset = 2)
  ),
  fluidRow(
    column(5,h4("Blood"),offset = 2)
  ),
  fluidRow(
    column(5,plotOutput("plot16"),offset = 2),
    column(5,plotOutput("plot13"))
  ),
  fluidRow(
    column(5,plotOutput("plot15"),offset = 2),
    column(5,plotOutput("plot20"))
  ),
  fluidRow(
    column(5,plotOutput("plot17"),offset = 2),
    column(5,plotOutput("plot21"))
  ),
  fluidRow(
    column(5,h3("scRNA-Seq"),offset = 2)
  ),
  fluidRow(
    column(5,h4("Combined data from GSE162577, GSE96587_batch1 & batch 2 (removed stimulated cells from GSM2560249)"),offset = 2)
  ),
  fluidRow(
    column(10,plotOutput("plot35"),offset = 2),
  ),
  fluidRow(
    column(7,plotOutput("plot34"),offset = 2),
  ),
  # fluidRow(
  #   column(5,h4("GSE135779_child"),offset = 2)
  # ),
  # fluidRow(
  #   column(10,plotOutput("plot37"),offset = 2),
  # ),
  # fluidRow(
  #   column(7,plotOutput("plot36"),offset = 2),
  # )
))

