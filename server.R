# library packages --------------------------------------------------------


library(tidyverse)

library(shiny)
library(plotly)


# load pandaomics data ----------------------------------------------------

# 最为重要的地方在于数据格式的整理和统一
# 各个数据集间命名格式的统一是极为重要的



# load DESeq2 data --------------------------------------------------------

# should be placed in global.R
combined_DESeq2_res <- readRDS("./data/combined_DESeq2_res.rds")

# 从以前文件拷过来的
# load('./data/pandaomics_meta_expr_data.Rdata')



# some function -----------------------------------------------------------

# functions for global use

thematic::thematic_shiny(font = "auto")


# server ------------------------------------------------------------------

server <- function(input, output, session) {
  # bslib::bs_themer()

  # var used in server
  # dataInput <- reactive({
  #   gene <- input$geneSymbol
  #   geneSym <- toupper(gene)
  #   geneSym
  # })


  # below are render* script

  output$plot1 <- renderPlot(
    {
      gene <- input$geneSymbol
      geneSym <- toupper(gene)

      df <- combined_DESeq2_res %>%
        dplyr::filter(Symbol == {{ geneSym }}) %>%
        mutate(
          lower = log2FoldChange - 1.96 * lfcSE,
          up = log2FoldChange + 1.96 * lfcSE
        ) %>%
        arrange(log2FoldChange)

      ggplot(df, aes(
        x = log2FoldChange,
        y = fct_reorder(GSE_name, padj, .na_rm = FALSE)
      )) +
        geom_point(aes(color = padj, size = padj)) +
        theme_bw() +
        guides(size = "none") +
        scale_color_gradientn(
          na.value = "gray50",
          limits = c(0, 0.05),
          # colours = paletteer::paletteer_c()
          colours = RColorBrewer::brewer.pal(
            n = 11, name = "RdBu"
          ),
          breaks = c(0, 0.025, 0.05)
        ) +
        scale_size(trans = "reverse") +
        theme(axis.text = element_text(face = "bold", size = 10)) +
        ylab("GSE Number")
    },
    width = "auto",
    height = "auto",
    res = 72
  )

  # plotly plot for plot1

  output$plot1_plotly <- renderPlotly({
    gene <- input$geneSymbol
    geneSym <- toupper(gene)

    df <- combined_DESeq2_res %>%
      dplyr::filter(Symbol == geneSym) %>%
      mutate(
        lower = log2FoldChange - 1.96 * lfcSE,
        up = log2FoldChange + 1.96 * lfcSE
      ) %>%
      arrange(log2FoldChange)

    p1 <- ggplot(df, aes(
      x = log2FoldChange,
      y = fct_reorder(GSE_name, padj, .na_rm = FALSE)
    )) +
      geom_point(aes(color = padj, size = padj)) +
      theme_bw() +
      guides(size = "none") +
      scale_color_gradientn(
        na.value = "gray50",
        limits = c(0, 0.05),
        # colours = paletteer::paletteer_c()
        colours = RColorBrewer::brewer.pal(
          n = 11, name = "RdBu"
        ),
        breaks = c(0, 0.025, 0.05)
      ) +
      scale_size(trans = "reverse") +
      theme(axis.text = element_text(face = "bold", size = 8)) +
      ylab("GSE Number")

    plotly::ggplotly(p1,
                     tooltip = c('x', 'color')
                     )
  })


  # plot2
  output$plot2 <- renderPlot({
    gse <- input$GSEnum
    gene <- input$geneSymbol
    geneSym <- toupper(gene)

    df_expr <- read_csv(glue::glue("./data/{gse}_expression.csv"),
      show_col_types = FALSE
    )
    df_meta <- read_csv(glue::glue("./data/{gse}_metadata.csv"),
      show_col_types = FALSE
    )

    if ("Cell Type" %in% names(df_meta)) {
      df_meta <- df_meta %>%
        rename("Cell_Type" = "Cell Type")

      df_p <- df_expr %>%
        filter(gene == geneSym) %>%
        column_to_rownames("gene") %>%
        t() %>%
        as.data.frame() %>%
        rownames_to_column() %>%
        as_tibble() %>%
        left_join(df_meta, by = c("rowname" = "name"))

      ggplot(df_p, aes_string(
        x = "Cell_Type",
        y = geneSym,
        fill = "Disease"
      )) +
        geom_boxplot() +
        ggtitle(label = gse) +
        theme_bw() +
        ggprism::theme_prism(base_size = 10) +
        paletteer::scale_color_paletteer_d(
          "ggprism::candy_bright"
        ) +
        theme(
          legend.position = "bottom",
          axis.text.x = element_text(angle = 45)
        )
    } else if ("Source" %in% names(df_meta)) {
      df_p <- df_expr %>%
        filter(gene == geneSym) %>%
        column_to_rownames("gene") %>%
        t() %>%
        as.data.frame() %>%
        rownames_to_column() %>%
        as_tibble() %>%
        left_join(df_meta, by = c("rowname" = "name"))

      ggplot(df_p, aes_string(
        x = "Source",
        y = geneSym,
        fill = "Disease"
      )) +
        geom_boxplot() +
        ggtitle(label = gse) +
        theme_bw() +
        theme(legend.position = "bottom") +
        ggprism::theme_prism(base_size = 20) +
        paletteer::scale_color_paletteer_d(
          "ggprism::candy_bright"
        ) +
        theme(legend.position = "bottom")
    }
  })


  # table1 for GSE annotation

  output$table1 <- renderTable(
    gse_annotation
  )


  # plot for Tissue expression from
  output$tissue_plot <- renderPlot({
    gene <- input$geneSymbol
    geneSym <- toupper(gene)

    # 有的基因有两个，这是要注意的
    tissue_rna <- read_csv(
      "./data/rna_tissue_exp_matrix.csv"
    ) %>%
      filter(`Gene.name` == geneSym) %>%
      column_to_rownames("Gene.name") %>%
      t() %>%
      as.data.frame() %>%
      rownames_to_column() %>%
      as_tibble()

    # annotation for tissue
    bulk_pair <- read_csv("./data/bulk_pair.csv") %>%
      pivot_longer(cols = everything())

    d_p <- tissue_rna %>%
      left_join(bulk_pair,
        by = c("rowname" = "name")
      ) %>%
      arrange(desc(.data[[geneSym]])) %>% 
      mutate(rowname = fct_reorder(rowname, .data[[geneSym]],
                                   .desc = TRUE
                                   )
             )

    ggplot(d_p, aes(
      # x = fct_reorder(rowname, .data[[geneSym]],
      #   .desc = TRUE
      # ),
      x = rowname,
      y = .data[[geneSym]],
      fill = value
    )) +
      geom_bar(stat = "identity", show.legend = TRUE) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(
        x = "Human tissues",
        y = glue::glue("{geneSym} Expression(nTPM)")
      )
  })


  # plot for Cell expression
  output$cell_plot <- renderPlot({
    gene <- input$geneSymbol
    geneSym <- toupper(gene)

    cell_rna <- read_csv(
      "./data/rna_sc_exp_matrix.csv"
    )
    # cell_pair <- read_csv(
    #   './data/sc_pair.csv'
    # )

    d_p <- cell_rna %>%
      filter(`Gene.name` == geneSym) %>%
      pivot_longer(cols = -`Gene.name`) %>%
      arrange(desc(value))

    ggplot(d_p, aes(
      x = fct_reorder(name, value,
        .desc = TRUE
      ),
      y = value,
      fill = value
    )) +
      geom_bar(stat = "identity", show.legend = TRUE) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(
        x = "Human Cell ",
        y = glue::glue("{geneSym} Expression(nTPM)")
      )
  })


  # render summary text for GSE
  output$text1 <- renderUI({
    gse <- input$GSEnum

    summary_bulk <- readxl::read_excel(
      "./data/Summary_for_bulk_RNASeq.xlsx"
    )

    summary_text <- summary_bulk %>%
      filter(GSE == gse) %>%
      pull(Summary)

    hlink <- summary_bulk %>%
      filter(GSE == gse) %>%
      pull(PandaOmics)

    # paste0(summary_text, '\n', '\n',hlink)

    HTML(glue::glue(
      '<p> {summary_text} <br> <br> <a href="{hlink}">PandaOmics</a></p>'
    ))
  })



  ## IBD-----------------------------

  data_expr <- reactive({
    gse_num <- input$GSE_IBD

    expSet <- read_csv(
      glue::glue("./data/IBD_Data/{gse_num}_Pandaomics_expression.csv"),
      show_col_types = FALSE
    )

    expSet
  })


  data_meta <- reactive({
    gse_num <- input$GSE_IBD

    metadata <- read_csv(
      glue::glue("./data/IBD_Data/{gse_num}_metadata_pandaomics.csv"),
      show_col_types = FALSE
    )

    metadata
  })

  output$IBD_gene_expr <- renderPlot({
    geneSym <- input$gs_IBD |> toupper()
    gse_num <- input$GSE_IBD

    expSet <- data_expr()

    metadata <- data_meta()

    expSet <- expSet %>% column_to_rownames("gene")
    p <- identical(metadata$name, colnames(expSet))
    if (!p) expSet <- expSet[, match(metadata$name, colnames(expSet))]
    expSet <- log2(expSet + 1)

    top_gene_df <- expSet %>%
      as.data.frame() %>%
      dplyr::filter(rownames(.) == geneSym) %>%
      t() %>%
      data.frame() %>%
      tibble::rownames_to_column("name") %>%
      dplyr::inner_join(dplyr::select(
        metadata,
        name,
        group
      ), by = join_by(name))


    ggplot(
      top_gene_df,
      aes_string(
        x = "group", y = geneSym,
        color = "group"
      )
    ) +
      geom_jitter(width = 0.2, height = 0) +
      geom_boxplot() +
      # theme_classic() +
      ggtitle(label = gse_num) +
      ggprism::theme_prism(base_size = 20) +
      paletteer::scale_color_paletteer_d(
        "ggprism::candy_bright"
      )
  })
}






