# library packages --------------------------------------------------------

suppressPackageStartupMessages(library(tidyverse))

library(shiny)
library(plotly)
# library(reactablefmtr)


# load pandaomics data ----------------------------------------------------

# 最为重要的地方在于数据格式的整理和统一
# 各个数据集间命名格式的统一是极为重要的



# load DESeq2 data --------------------------------------------------------

# should be placed in global.R
# combined_DESeq2_res <- readRDS("./data/combined_DESeq2_res.rds")

# 从以前文件拷过来的
# load('./data/pandaomics_meta_expr_data.Rdata')



# some function -----------------------------------------------------------

# functions for global use

# thematic::thematic_shiny(font = "auto")


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

  gene_event <- eventReactive(input$click1, {
    input$geneSymbol
  })

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
    gene <- gene_event()
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
      ylab("GSE Number") +
      ggtitle(glue::glue("{geneSym} Diff Expression"))

    plotly::ggplotly(p1,
      tooltip = c("x", "color")
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

      ggplot(df_p, aes(
        x = Cell_Type,
        y = .data[[geneSym]],
        fill = Disease
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

      ggplot(df_p, aes(
        x = Source,
        y = .data[[geneSym]],
        fill = Disease
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
      ))

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
      aes(
        x = group, 
        y = .data[[geneSym]],
        color = group
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


  # HPA ---------------------------------------------

  hpa_isoform_immuncells_M <- reactive({
    geneSym <- input$hpa_trans_gene |> toupper()

    hpa_transcript_rna_immunecells %>%
      filter(SYMBOL == geneSym) %>%
      dplyr::select(enstid, starts_with("TPM")) %>%
      column_to_rownames("enstid") -> M

    if (!is.data.frame(M)) {
      validate(paste0("'", geneSym, "' is not a data frame"))
    }

    M
  })

  output$hpa_isoform_immuncells_p1 <- renderPlot({
    geneSym <- input$hpa_trans_gene |> toupper()

    hpa_transcript_rna_immunecells %>%
      filter(SYMBOL == geneSym) %>%
      dplyr::select(enstid, starts_with("TPM")) %>%
      column_to_rownames("enstid") -> M

    M %>%
      rownames_to_column() %>%
      pivot_longer(cols = -rowname) %>%
      distinct() %>%
      separate_wider_delim(
        cols = name, delim = ".",
        names = c("foo", "celltype", "sample"),
        cols_remove = FALSE
      ) %>%
      dplyr::mutate(logvalue = log2(value + 1)) -> df_p


    ggplot(df_p, aes(
      x = celltype, y = logvalue,
      fill = rowname,
      color = rowname
    )) +
      geom_boxplot() +
      ggprism::theme_prism() +
      ggprism::scale_colour_prism(palette = "colors") +
      ggprism::scale_fill_prism(palette = "colors") +
      theme(axis.text.x = element_text(angle = 90)) +
      ylab("logTPM") -> p1

    p1
  })


  output$hpa_isoform_immuncells_p2 <- renderPlot({
    geneSym <- input$hpa_trans_gene |> toupper()

    colors_f <- circlize::colorRamp2(c(-3, 0, 3),
      c("navy", "white", "firebrick3"),
      transparency = 0
    )

    ComplexHeatmap::Heatmap(log2(t(hpa_isoform_immuncells_M()) + 1),
      name = glue::glue("{geneSym}"),
      na_col = "gray",
      border_gp = grid::gpar(col = "black"),
      show_row_dend = T,
      show_column_dend = F,
      cluster_columns = F,
      cluster_rows = T,
      clustering_distance_columns = "euclidean",
      clustering_distance_rows = "euclidean",
      clustering_method_rows = "complete",
      clustering_method_columns = "complete",
      # col = colors_f,
      column_names_rot = 90,
      show_column_names = TRUE,
      show_row_names = TRUE,
      width = unit(5, "cm"),
      row_names_gp = grid::gpar(fontsize = 4),
      # left_annotation = ha_row,
      show_heatmap_legend = TRUE,
      # cell_fun = function(j, i, x, y, width, height, fill) {
      #   grid.text(round(M[i, j], 3), x, y, gp = gpar(fontsize = 4))
      # }
    )
  })


  output$hpa_isoform_immuncells_p3 <- renderPlot({
    M2 <- log2(hpa_isoform_immuncells_M() + 1)

    # M2 <- M2[rowSums(M2) > 1, ]
    keep <- rowSums(M2 >= 1) >= 3
    M2 <- M2[keep, ]

    M2 %>%
      rownames_to_column() %>%
      pivot_longer(cols = -rowname) %>%
      distinct() %>%
      separate_wider_delim(
        cols = name, delim = ".",
        names = c("foo", "celltype", "sample"),
        cols_remove = FALSE
      ) -> df_p2

    ggplot(df_p2, aes(
      x = celltype, y = value,
      fill = rowname,
      color = rowname
    )) +
      geom_boxplot() +
      ggprism::theme_prism() +
      ggprism::scale_colour_prism(palette = "colors") +
      ggprism::scale_fill_prism(palette = "colors") +
      theme(axis.text.x = element_text(angle = 90)) +
      ylab("logTPM") -> p2

    p2
  })


  output$download1 <- downloadHandler(
    filename = function() {
      paste0(input$hpa_trans_gene, ".csv")
    },
    content = function(file) {
      write.csv(hpa_isoform_immuncells_M(), file)
    }
  )


  # scRNA ---------------------------------------------

  # Sys.setenv("VROOM_CONNECTION_SIZE"=131072*20)

  make_gse_plot <- function(scgene, scGSE) {
    p_data <- vroom::vroom(glue::glue("./data/talk2data/{scGSE}_sc_meanCI.csv"),
      delim = ",",
      show_col_types = FALSE
    )
    condition1 <- names(p_data)[[2]]
    celltype <- names(p_data)[[1]]
    genemean <- paste0(scgene, "_mean")
    genelow <- paste0(scgene, "_low")
    geneup <- paste0(scgene, "_up")

    p_data <- p_data %>%
      dplyr::select(1:2, starts_with({{ scgene }})) %>%
      dplyr::rename(
        "Condition" = all_of(condition1),
        "CellType" = all_of(celltype)
      )

    scpp <- ggplot(p_data, aes(
      x = CellType,
      y = .data[[genemean]],
      fill = Condition
    )) +
      geom_col(position = "dodge") +
      geom_errorbar(
        aes(
          x = CellType,
          ymin = .data[[genelow]],
          ymax = .data[[geneup]]
        ),
        width = 0.4, colour = "orange",
        alpha = 0.9, linewidth = 1,
        position = position_dodge(.9)
      ) +
      paletteer::scale_fill_paletteer_d("ggsci::category20_d3") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90))

    # scpp

    plotly::ggplotly(scpp)
  }


  # IBD
  scgene_event <- eventReactive(input$click2, {
    input$GeneName
  })

  output$scPlot1 <- renderPlotly({
    scgene <- scgene_event() |> toupper()
    scGSE <- input$singleGSE

    make_gse_plot(
      scgene = scgene,
      scGSE = scGSE
    )
  })


  # SLE
  scgene_event3 <- eventReactive(input$click3, {
    input$singleGSE2
  })

  output$scPlot2 <- renderPlotly({
    scgene <- scgene_event() |> toupper()
    scGSE <- scgene_event3()

    make_gse_plot(
      scgene = scgene,
      scGSE = scGSE
    )
  })


  # COPD
  scgene_event4 <- eventReactive(input$click4, {
    input$singleGSE3
  })

  output$scPlot3 <- renderPlotly({
    scgene <- scgene_event() |> toupper()
    scGSE <- scgene_event4()

    make_gse_plot(
      scgene = scgene,
      scGSE = scGSE
    )
  })


  # Psoriasis
  scgene_event5 <- eventReactive(input$click5, {
    input$singleGSE4
  })

  output$scPlot4 <- renderPlotly({
    scgene <- scgene_event() |> toupper()
    scGSE <- scgene_event5()

    make_gse_plot(
      scgene = scgene,
      scGSE = scGSE
    )
  })


  # Others
  scgene_event6 <- eventReactive(input$click6, {
    input$singleGSE5
  })

  output$scPlot5 <- renderPlotly({
    scgene <- scgene_event() |> toupper()
    scGSE <- scgene_event6()

    make_gse_plot(
      scgene = scgene,
      scGSE = scGSE
    )
  })


  read_meta_table <- function(scGSE) {
    dd <- read_csv(
      glue::glue("./data/talk2data/table/{scGSE}_table.csv"),
      show_col_types = FALSE
    )

    dd
  }


  output$scTable1 <- DT::renderDataTable({
    scGSE <- input$singleGSE

    dd <- read_meta_table(scGSE = scGSE)

    # dd %>% DT::datatable()

    dd
  })


  output$scTable2 <- DT::renderDataTable({
    scGSE <- input$singleGSE2

    dd <- read_meta_table(scGSE = scGSE)

    # dd %>% DT::datatable()

    dd
  })

  output$scTable3 <- DT::renderDataTable({
    scGSE <- input$singleGSE3

    dd <- read_meta_table(scGSE = scGSE)

    # dd %>% DT::datatable()

    dd
  })

  output$scTable4 <- DT::renderDataTable({
    scGSE <- input$singleGSE4

    dd <- read_meta_table(scGSE = scGSE)

    # dd %>% DT::datatable()

    dd
  })


  output$scTable5 <- DT::renderDataTable({
    scGSE <- input$singleGSE5

    dd <- read_meta_table(scGSE = scGSE)

    # dd %>% DT::datatable()

    dd
  })


  read_summary_table <- function(gse) {
    df <- st %>%
      filter(GSE_Num == {{ gse }}) %>%
      mutate(across(everything(), as.character)) %>%
      pivot_longer(
        cols = everything(),
        values_to = "Value"
      )

    reactable::reactable(df,
      highlight = TRUE,
      searchable = TRUE,
      wrap = TRUE,
      columns = list(
        Value = reactable::colDef(
          cell = function(value) {
            if (value == df$Value[4]) {
              div(
                style = list(
                  "max-height" = "80px",
                  "overflow" = "hidden"
                ),
                # value,
                a(
                  href = "#",
                  onclick = "this.parentNode.style.maxHeight = null; this.style.display = 'none'",
                  "Show more"
                ),
                value
              )
            } else {
              value
            }
          }, minWidth = 200
        ),
        name = reactable::colDef(minWidth = 100)
      )
    )
  }


  output$scST1 <- reactable::renderReactable({
    scGSE <- input$singleGSE

    dd <- read_summary_table(gse = scGSE)

    # dd %>% DT::datatable()

    dd
  })

  output$scST2 <- reactable::renderReactable({
    scGSE <- input$singleGSE2

    dd <- read_summary_table(gse = scGSE)

    # dd %>% DT::datatable()

    dd
  })

  output$scST3 <- reactable::renderReactable({
    scGSE <- input$singleGSE3

    dd <- read_summary_table(gse = scGSE)

    # dd %>% DT::datatable()

    dd
  })

  output$scST4 <- reactable::renderReactable({
    scGSE <- input$singleGSE4

    dd <- read_summary_table(gse = scGSE)

    # dd %>% DT::datatable()

    dd
  })

  output$scST5 <- reactable::renderReactable({
    scGSE <- input$singleGSE5

    dd <- read_summary_table(gse = scGSE)

    # dd %>% DT::datatable()

    dd
  })

  ## mail-----------------------------

  observeEvent(input$send, {
    mail_res <- sendmailR::sendmail(
      from = "",
      to = as.character(input$user_mail),
      subject = "Gene Expression",
      msg = mime_part("get it"),
      engineopts = list(username = "发送邮箱", password = "授权码"),
      control = list(smtpServer = "smtp.163.com:465", verbose = TRUE)
    )

    if (!is.null(mail_res)) {
      output$mailout <- renderText("Mail Sended!")
    } else {
      output$mailout <- renderText("Mail Send fail!")
    }
  })
}
