# library packages --------------------------------------------------------

library(shiny)
library(bslib)

library(plotly)

# 不错的框架
# library(golem)
# library(rhino)
# require(periscope)


# functions  --------------------------------------------------------------



# UI ----------------------------------------------------------------------


ui <- navbarPage(
  "SIBR Expression Explorer",
  # theme = shinythemes::themeSelector(),
  # theme = shinythemes::shinytheme("united"),
  theme = bslib::bs_theme(
    bootswatch = "united",
    # 'navbar_bg' = "#25443B"
  ) |>
    bslib::bs_add_rules(
      rules = ".navbar.navbar-default {
                        background-color: #320638 !important;
                    }"
    ),
  
  windowTitle = "SIBR Shiny",
  
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
        actionButton(
          inputId = "click1",
          label = "GO!"
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
        tabsetPanel(tabPanel("Plot", plotlyOutput("plot1_plotly"))),
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
      column(3, textOutput("Here")),
      column(
        7,
        h4("Expression in Different Tissues"),
        plotOutput("tissue_plot"),
        # offset = 2,
        # plotOutput("cell_plot")
      )
    ),
    fluidRow(
      column(3, textOutput("Here. ")),
      column(
        7,
        h4("Expression in Different Cell Type"),
        # plotOutput("tissue_plot"),
        # offset = 2,
        plotOutput("cell_plot")
      )
    )
  ),


  # UI for single cell RNA expression

  tabPanel(
    title = 'scRNA',
    
    sidebarLayout(
      sidebarPanel(
        textInput(
          inputId = "GeneName",
          label = "Input Gene Symbol: ",
          value = "CAMK4"
        ),
        
        actionButton(
          inputId = "click2",
          label = "Input Gene!"
        ),
        
        tags$br(),
        
        selectInput(
          inputId = "singleGSE",
          label = "Choose GSE:",
          choices = chioce_IBD
        ),
        width = 3
      ),
      
      mainPanel(
        h3("Gene Expression in IBD studies"),
        # tabsetPanel(
        #   tabPanel("Plot", plotOutput("plot1")),
        #   tabPanel("Table", tableOutput("table1"))
        # )
        column(
          8,
          tabsetPanel(
            tabPanel("Plot", plotlyOutput("scPlot1",
                                        # width = '60%'
                                        width = 800,
                                        height = 600
            )),
            tabPanel("CellCount", DT::dataTableOutput("scTable1")),
            tabPanel("StudyMeta", reactable::reactableOutput("scST1"))
          )
        )
      )
    ),
    
    hr(),
    
    fluidRow(column(8,
                    h3("Gene Expression in SLE studies"),
                    offset = 3
    )),
    
    tags$br(),
    
    fluidRow(
      column(
        3,
        h5("Select different GSE"),
        selectInput(
          inputId = "singleGSE2", 
          label = "Choose GSE:",
          choices = chioce_SLE
        ),
        
        actionButton(
          inputId = "click3",
          label = "Show SLE!"
        ),
      ),
      column(
        8,
        tabsetPanel(
          tabPanel("Plot", plotlyOutput("scPlot2",
                                      # width = '60%'
                                      width = 800,
                                      height = 600
          )),
          tabPanel("CellCount", DT::dataTableOutput("scTable2")),
          tabPanel("StudyMeta", reactable::reactableOutput("scST2"))
        )
      )
    ),
    
    hr(),
    
    fluidRow(column(8,
                    h3("Gene Expression in COPD studies"),
                    offset = 3
    )),
    
    tags$br(),
    
    fluidRow(
      column(
        3,
        h5("Select different GSE"),
        selectInput(
          inputId = "singleGSE3", 
          label = "Choose GSE:",
          choices = chioce_COPD
        ),
        
        actionButton(
          inputId = "click4",
          label = "Show COPD!"
        ),
      ),
      column(
        8,
        tabsetPanel(
          tabPanel("Plot", plotlyOutput("scPlot3",
                                      # width = '60%'
                                      width = 800,
                                      height = 600
          )),
          tabPanel("CellCount", DT::dataTableOutput("scTable3")),
          tabPanel("StudyMeta", reactable::reactableOutput("scST3"))
        )
      )
    ),
    
    hr(),
    
    fluidRow(column(8,
                    h3("Gene Expression in Psoriasis studies"),
                    offset = 3
    )),
    
    tags$br(),
    
    fluidRow(
      column(
        3,
        h5("Select different GSE"),
        selectInput(
          inputId = "singleGSE4", 
          label = "Choose GSE:",
          choices = chioce_Psoriasis
        ),
        
        actionButton(
          inputId = "click5",
          label = "Show Psoriasis!"
        ),
      ),
      column(
        8,
        tabsetPanel(
          tabPanel("Plot", plotlyOutput("scPlot4",
                                      # width = '60%'
                                      width = 800,
                                      height = 600
          )),
          tabPanel("CellCount", DT::dataTableOutput("scTable4")),
          tabPanel("StudyMeta", reactable::reactableOutput("scST4"))
        )
      )
    ),
    
    
    hr(),
    
    fluidRow(column(8,
                    h3("Gene Expression in Other studies"),
                    offset = 3
    )),
    
    tags$br(),
    
    fluidRow(
      column(
        3,
        h5("Select different GSE"),
        selectInput(
          inputId = "singleGSE5", 
          label = "Choose GSE:",
          choices = chioce_Others
        ),
        
        actionButton(
          inputId = "click6",
          label = "Show Others!"
        ),
      ),
      column(
        8,
        tabsetPanel(
          tabPanel("Plot", plotlyOutput("scPlot5",
                                      # width = '60%'
                                      width = 800,
                                      height = 600
          )),
          tabPanel("CellCount", DT::dataTableOutput("scTable5")),
          tabPanel("StudyMeta", reactable::reactableOutput("scST5"))
        )
      )
    ),
    
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
      mainPanel(
        h3("IBD gene Expression"),
        tabsetPanel(
          tabPanel("Plot", plotOutput("IBD_gene_expr")),
          tabPanel("PCA")
        )
      )
    )
  ),


  # UI for HPA
  tabPanel(
    "HPA",
    sidebarLayout(
      sidebarPanel(
        textInput(
          inputId = "hpa_trans_gene",
          label = "Input Gene Symbol",
          value = "IRF4"
        ),
        downloadButton(
          outputId = "download1",
          label = "Download",
          icon = shiny::icon("download")
        )
      ),
      mainPanel(
        h3("HPA immunecells isoform Expression"),
        tabsetPanel(
          tabPanel("BoxPot", plotOutput("hpa_isoform_immuncells_p1")),
          tabPanel("Heatmap", plotOutput("hpa_isoform_immuncells_p2")),
          tabPanel(
            "BoxPot_Fliter",
            plotOutput("hpa_isoform_immuncells_p3")
          )
        )
      )
    )
  ),
  navbarMenu(
    "More",
    tabPanel("Others")
  )
)


ui <- tagList(
  shiny.info::version(position = "bottom right"),
  shiny.info::powered_by("DDS-SIBR", link = "https://www.sanofi.com/"),
  shiny.info::busy(),
  ui
)


