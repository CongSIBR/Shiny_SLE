# library packages --------------------------------------------------------

library(shiny)
library(bslib)

library(plotly)



# functions  --------------------------------------------------------------



# UI ----------------------------------------------------------------------


ui <- navbarPage("Bulk RNASeq Expression",
  # theme = shinythemes::themeSelector(),
  # theme = shinythemes::shinytheme("united"),
  theme = bslib::bs_theme(bootswatch = 'united',
                          # 'navbar_bg' = "#25443B"
                          ) |>
    bslib::bs_add_rules(
      rules = ".navbar.navbar-default {
                        background-color: #320638 !important;
                    }"
    ),

  windowTitle = "GEO Shiny",
  
  tabPanel(
    "SLE",

    # copy from ui_demo.R
    sidebarLayout(
      sidebarPanel(
        textInput(
          inputId = "geneSymbol",
          label = "Input Gene Symbol: ",
          value = "IRF4"
        ),
        width = 3
      ),
      
      mainPanel(
        h2("DEGs across 33 GSE datasets"),
        h4("note that: not all gene exists in all 33 datasets!"),
        # tabsetPanel(
        #   tabPanel("Plot", plotOutput("plot1")),
        #   tabPanel("Table", tableOutput("table1"))
        # )
      )
    ),
    
    fluidRow(
      column(8,
        # plotOutput("plot1"),
        # plotlyOutput("plot1_plotly"),
        tabsetPanel(tabPanel('Plot', plotlyOutput("plot1_plotly"))),
        offset = 3
      )
    ),
    
    hr(),
    
    fluidRow(column(8,
      h2("Single GSE Expression"),
      offset = 3
    )),
    
    tags$br(),
    
    fluidRow(
      column(
        3,
        h5("Select different GSE"),
        selectInput(
          inputId = "GSEnum", label = "GSE_num",
          choices = chioce1
        )
      ),
      
      column(
        8,
        tabsetPanel(
          tabPanel("Plot", plotOutput("plot2", 
                                      # width = '60%'
                                      width = 800,
                                      height = 600
                                      )),
          tabPanel("Summary", uiOutput("text1"))
        )
        # plotOutput('plot2')
      )
    ),
    
    hr(),
    
    fluidRow(
      column(3, textOutput('Here')), 
      column(7,
        h4("Expression in Different Tissues"),
        plotOutput("tissue_plot"),
        # offset = 2,
        # plotOutput("cell_plot")
      )
    ),
    
    fluidRow(
      column(3, textOutput('Here. ')),
      column(7,
             h4("Expression in Different Cell Type"),
             # plotOutput("tissue_plot"),
             # offset = 2,
             plotOutput("cell_plot")
      )
    )
    
    
  ),


  # UI for IBD

  tabPanel(
    "IBD",
    sidebarLayout(
      sidebarPanel(
        textInput(
          inputId = "gs_IBD",
          label = "Input Gene Symbol",
          value = "IRF4"
        ),
        selectInput(
          inputId = "GSE_IBD",
          label = "IBD GSE number",
          choices =
            c(
              "GSE59071", "GSE23597", "GSE73661",
              "GSE87466"
            )
        ),
        width = 3
      ),
      mainPanel(h3('IBD gene Expression'),
                tabsetPanel(
        tabPanel("Plot", plotOutput("IBD_gene_expr")),
        tabPanel("PCA")
      ))
    )
  ),
  navbarMenu(
    "More",
    tabPanel("Others")
  )
)




