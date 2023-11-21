




# library packages --------------------------------------------------------


library(tidyverse)

library(shiny)
library(plotly)


# load pandaomics data ----------------------------------------------------

# 最为重要的地方在于数据格式的整理和统一
# 各个数据集间命名格式的统一是极为重要的



# load DESeq2 data --------------------------------------------------------

# should place in global.R
combined_DESeq2_res <- readRDS('./data/combined_DESeq2_res.rds')

# 从以前文件拷过来的
# load('./data/pandaomics_meta_expr_data.Rdata')



# some function -----------------------------------------------------------

# functions for global use



# server ------------------------------------------------------------------

server <- function(input, output, session){
  
  # var used in server
  # dataInput is a function
  dataInput <- reactive({
    gene <- input$geneSymbol
    geneSym <- toupper(gene)
    geneSym
  })
  
  
  # below for render* script
  
  output$plot1 <- renderPlot({
    gene <- input$geneSymbol
    geneSym <- toupper(gene)
    
    df <- combined_DESeq2_res %>% 
      dplyr::filter(Symbol == geneSym) %>% 
      mutate(lower = log2FoldChange - 1.96 * lfcSE,
             up = log2FoldChange + 1.96 * lfcSE
             ) %>% 
      arrange(log2FoldChange)
    
    ggplot(df, aes(x=log2FoldChange, 
                   y = fct_reorder(GSE_name, padj)
                   )) +
      geom_point(aes(color = padj, size = padj)) +
      theme_bw() +
      guides(size= 'none') +
      scale_color_gradientn(na.value = 'gray50',
                            limits = c(0, 0.05),
                            # colours = paletteer::paletteer_c()
                            colours = RColorBrewer::brewer.pal(
                              n = 11, name = 'RdBu'
                            ),
                            breaks=c(0,0.025,0.05)
                            ) +
      scale_size(trans = 'reverse') +
      theme(axis.text = element_text(face = 'bold',size=10))
    
  })
  
  # plotly plot for plot1
  
  output$plot1_plotly <- renderPlotly({
    gene <- input$geneSymbol
    geneSym <- toupper(gene)
    
    df <- combined_DESeq2_res %>% 
      dplyr::filter(Symbol == geneSym) %>% 
      mutate(lower = log2FoldChange - 1.96 * lfcSE,
             up = log2FoldChange + 1.96 * lfcSE
      )
    
    p1 <- ggplot(df, aes(x=log2FoldChange, 
                         y = fct_reorder(GSE_name, padj)
    )) +
      geom_point(aes(color = padj, size = padj)) +
      theme_bw() +
      guides(size= 'none') +
      scale_color_gradientn(na.value = 'gray50',
                            limits = c(0, 0.05),
                            # colours = paletteer::paletteer_c()
                            colours = RColorBrewer::brewer.pal(
                              n = 11, name = 'RdBu'
                            ),
                            breaks=c(0,0.025,0.05)
      ) +
      scale_size(trans = 'reverse') +
      theme(axis.text = element_text(face = 'bold',size=8))
    
    plotly::ggplotly(p1)
  })
  
  
  # plot2
  output$plot2 <- renderPlot({
    gse <- input$GSEnum
    gene <- input$geneSymbol
    geneSym <- toupper(gene)
    
    df_expr <- read_csv(glue::glue('./data/{gse}_expression.csv'))
    df_meta <- read_csv(glue::glue('./data/{gse}_metadata.csv')) 
    
    if('Cell Type' %in% names(df_meta)){
      df_meta <- df_meta %>% 
        rename('Cell_Type'='Cell Type')
      
      df_p <- df_expr %>% filter(gene == geneSym) %>% 
        column_to_rownames('gene') %>% 
        t() %>% as.data.frame() %>% 
        rownames_to_column() %>% as_tibble() %>% 
        left_join(df_meta, by = c('rowname'='name'))
      
      ggplot(df_p, aes_string(x='Cell_Type', 
                              y = geneSym,
                              fill = 'Disease'
      )) +
        geom_boxplot() +
        ggtitle(label = gse) +
        theme_bw() +
        theme(legend.position = "bottom")
      
    } else if('Source' %in% names(df_meta)){
      
      df_p <- df_expr %>% filter(gene == geneSym) %>% 
        column_to_rownames('gene') %>% 
        t() %>% as.data.frame() %>% 
        rownames_to_column() %>% as_tibble() %>% 
        left_join(df_meta, by = c('rowname'='name'))
      
      ggplot(df_p, aes_string(x='Source', 
                              y = geneSym,
                              fill = 'Disease'
      )) +
        geom_boxplot() +
        ggtitle(label = gse) +
        theme_bw() +
        theme(legend.position = "bottom")
    }
    
  }
  )
  
  
  # table1
  
  output$table1 <- renderTable(
    gse_annotation
  )
  
  
  # plot for Tissue expression from
  output$tissue_plot <- renderPlot({
    gene <- input$geneSymbol
    geneSym <- toupper(gene)
    
    # 有的基因有两个，这是要注意的
    tissue_rna <- read_csv(
      './data/rna_tissue_exp_matrix.csv'
    ) %>% 
      filter(`Gene.name` == geneSym) %>% 
      column_to_rownames('Gene.name') %>% 
      t() %>% 
      as.data.frame() %>% 
      rownames_to_column() %>% 
      as_tibble()
      
    # annotation for tissue
    bulk_pair <- read_csv("./data/bulk_pair.csv") %>% 
      pivot_longer(cols = everything())
    
    d_p <- tissue_rna %>% left_join(bulk_pair,
                                    by = c('rowname'='name')
                                    ) %>% 
      arrange(desc(.data[[geneSym]]))
    
    ggplot(d_p, aes(x = fct_reorder(rowname, .data[[geneSym]], 
                                    .desc = TRUE), 
                    y = .data[[geneSym]], 
                    fill = value)) +
      geom_bar(stat = "identity", show.legend = TRUE) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(x = "Human tissues", 
           y = glue::glue("{geneSym} Expression(nTPM)"))
    
  })
  
  
  # plot for Cell expression 
  output$cell_plot <- renderPlot({
    gene <- input$geneSymbol
    geneSym <- toupper(gene)
    
    cell_rna <- read_csv(
      './data/rna_sc_exp_matrix.csv'
    )
    # cell_pair <- read_csv(
    #   './data/sc_pair.csv'
    # )
    
    d_p <- cell_rna %>% 
      filter(`Gene.name` == geneSym) %>% 
      pivot_longer(cols = -`Gene.name`) %>% 
      arrange(desc(value))
    
    ggplot(d_p, aes(x = fct_reorder(name, value, 
                                    .desc = TRUE), 
                    y = value, 
                    fill = value)) +
      geom_bar(stat = "identity", show.legend = TRUE) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(x = "Human Cell ", 
           y = glue::glue("{geneSym} Expression(nTPM)"))
    
  })
  
  
  # render summary text for GSE
  output$text1 <- renderText({
    gse <- input$GSEnum
    
    summary_bulk <- readxl::read_excel(
      './data/Summary_for_bulk_RNASeq.xlsx'
    ) %>% 
      dplyr::rename('GSE'='...1')
    
    summary_bulk %>% filter(GSE == gse) %>% 
      pull(Summary)
    
  })
  
}














