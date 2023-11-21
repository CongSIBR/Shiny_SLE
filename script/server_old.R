## load the required packages
library(shiny)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(Seurat)
library(patchwork)



# Data from PandaOmics ----------------------------------------------------


# setwd("/Users/yutong/Documents/Intern/Visualization/")
## Read in data
# meta1 <- read.csv("data/GSE131525_metadata.csv")
# expr1 <- read.csv("data/GSE131525_expression.csv")
# meta2 <- read.csv("data/GSE149050_metadata.csv")
# expr2 <- read.csv("data/GSE149050_expression.csv")
# meta3 <- read.csv("data/GSE117836_metadata.csv")
# expr3 <- read.csv("data/GSE117836_expression.csv")
# meta4 <- read.csv("data/GSE139358_metadata.csv")
# expr4 <- read.csv("data/GSE139358_expression.csv")
# meta5 <- read.csv("data/GSE164457_metadata.csv")
# expr5 <- read.csv("data/GSE164457_expression.csv")
# meta6 <- read.csv("data/GSE92387_metadata.csv")
# expr6 <- read.csv("data/GSE92387_expression.csv")
# meta7 <- read.csv("data/GSE110999_metadata.csv")
# expr7 <- read.csv("data/GSE110999_expression.csv")
# meta8 <- read.csv("data/GSE163121_metadata.csv")
# expr8 <- read.csv("data/GSE163121_expression.csv")
# meta9 <- read.csv("data/GSE136731_metadata.csv")
# expr9 <- read.csv("data/GSE136731_expression.csv")
# meta10 <- read.csv("data/GSE118254_metadata.csv")
# expr10 <- read.csv("data/GSE118254_expression.csv")
# meta13 <- read.csv("data/GSE122459_metadata.csv")
# expr13 <- read.csv("data/GSE122459_expression.csv")
# meta15 <- read.csv("data/GSE112087_metadata.csv")
# expr15 <- read.csv("data/GSE112087_expression.csv")
# meta16 <- read.csv("data/GSE162828_metadata.csv")
# expr16 <- read.csv("data/GSE162828_expression.csv")
# meta17 <- read.csv("data/GSE110685_metadata.csv")
# expr17 <- read.csv("data/GSE110685_expression.csv")
# meta20 <- read.csv("data/GSE80183_metadata.csv")
# expr20 <- read.csv("data/GSE80183_expression.csv")
# meta21 <- read.csv("data/GSE72509_metadata.csv")
# expr21 <- read.csv("data/GSE72509_expression.csv")

# save(list = ls(), './data/pandaomics_meta_expr_data.Rdata')

load('./data/pandaomics_meta_expr_data.Rdata')


# DESeq2 ------------------------------------------------------------------

# Read in DESeq2 output information
# DE1 <- read.delim("data/DESeq2/DESeq2_DEGs_lfcShrink_GSE131525_All.txt") |>
#   dplyr::select(Symbol, baseMean, log2FoldChange, lfcSE, pvalue, padj)
# DE1_1 <- read.delim("data/DESeq2/DESeq2_DEGs_lfcShrink_GSE131525_b.txt") |>
#   dplyr::select(Symbol, baseMean, log2FoldChange, lfcSE, pvalue, padj)
# DE1_2 <- read.delim("data/DESeq2/DESeq2_DEGs_lfcShrink_GSE131525_cd4.txt") |>
#   dplyr::select(Symbol, baseMean, log2FoldChange, lfcSE, pvalue, padj)
# DE1_3 <- read.delim("data/DESeq2/DESeq2_DEGs_lfcShrink_GSE131525_cd8.txt") |>
#   dplyr::select(Symbol, baseMean, log2FoldChange, lfcSE, pvalue, padj)
# DE1_4 <- read.delim("data/DESeq2/DESeq2_DEGs_lfcShrink_GSE131525_monocytes.txt") |>
#   dplyr::select(Symbol, baseMean, log2FoldChange, lfcSE, pvalue, padj)
# 
# DE2 <- read.delim("data/DESeq2/DESeq2_DEGs_lfcShrink_GSE72509.txt") |>
#   dplyr::select(Symbol, baseMean, log2FoldChange, lfcSE, pvalue, padj)
# DE3 <- read.delim("data/DESeq2/DESeq2_DEGs_lfcShrink_GSE92387.txt") |>
#   dplyr::select(Symbol, baseMean, log2FoldChange, lfcSE, pvalue, padj)
# DE4 <- read.delim("data/DESeq2/DESeq2_DEGs_lfcShrink_GSE110999.txt") |>
#   dplyr::select(Symbol, baseMean, log2FoldChange, lfcSE, pvalue, padj)
# DE5 <- read.delim("data/DESeq2/DESeq2_DEGs_lfcShrink_GSE112087.txt") |>
#   dplyr::select(Symbol, baseMean, log2FoldChange, lfcSE, pvalue, padj)
# 
# DE6 <- read.delim("data/DESeq2/DESeq2_DEGs_lfcShrink_GSE118254.txt") |>
#   dplyr::select(Symbol, baseMean, log2FoldChange, lfcSE, pvalue, padj)
# DE6_1 <- read.delim("data/DESeq2/DESeq2_DEGs_lfcShrink_GSE118254_activated_naive_b_cells_(an).txt") |>
#   dplyr::select(Symbol, baseMean, log2FoldChange, lfcSE, pvalue, padj)
# DE6_2 <- read.delim("data/DESeq2/DESeq2_DEGs_lfcShrink_GSE118254_antigen_secreting_cell_(asc).txt") |>
#   dplyr::select(Symbol, baseMean, log2FoldChange, lfcSE, pvalue, padj)
# DE6_3 <- read.delim("data/DESeq2/DESeq2_DEGs_lfcShrink_GSE118254_double_negative_b_cells_(dn2).txt") |>
#   dplyr::select(Symbol, baseMean, log2FoldChange, lfcSE, pvalue, padj)
# DE6_4 <- read.delim("data/DESeq2/DESeq2_DEGs_lfcShrink_GSE118254_resting_naive_b_cells_(rn).txt") |>
#   dplyr::select(Symbol, baseMean, log2FoldChange, lfcSE, pvalue, padj)
# DE6_5 <- read.delim("data/DESeq2/DESeq2_DEGs_lfcShrink_GSE118254_switched_memory_b_cells_(sm).txt") |>
#   dplyr::select(Symbol, baseMean, log2FoldChange, lfcSE, pvalue, padj)
# DE6_6 <- read.delim("data/DESeq2/DESeq2_DEGs_lfcShrink_GSE118254_transitional_3_b_cells_(t3).txt") |>
#   dplyr::select(Symbol, baseMean, log2FoldChange, lfcSE, pvalue, padj)
# 
# DE7 <- read.delim("data/DESeq2/DESeq2_DEGs_lfcShrink_GSE136731.txt") |>
#   dplyr::select(Symbol, baseMean, log2FoldChange, lfcSE, pvalue, padj)
# DE8 <- read.delim("data/DESeq2/DESeq2_DEGs_lfcShrink_GSE139358.txt") |>
#   dplyr::select(Symbol, baseMean, log2FoldChange, lfcSE, pvalue, padj)
# DE9 <- read.delim("data/DESeq2/DESeq2_DEGs_lfcShrink_GSE149050_All.txt") |>
#   dplyr::select(Symbol, baseMean, log2FoldChange, lfcSE, pvalue, padj)
# DE9_1 <- read.delim("data/DESeq2/DESeq2_DEGs_lfcShrink_GSE149050_b cells.txt") |>
#   dplyr::select(Symbol, baseMean, log2FoldChange, lfcSE, pvalue, padj)
# DE9_2 <- read.delim("data/DESeq2/DESeq2_DEGs_lfcShrink_GSE149050_cdc.txt") |>
#   dplyr::select(Symbol, baseMean, log2FoldChange, lfcSE, pvalue, padj)
# DE9_3 <- read.delim("data/DESeq2/DESeq2_DEGs_lfcShrink_GSE149050_cmo.txt") |>
#   dplyr::select(Symbol, baseMean, log2FoldChange, lfcSE, pvalue, padj)
# DE9_4 <- read.delim("data/DESeq2/DESeq2_DEGs_lfcShrink_GSE149050_pdc.txt") |>
#   dplyr::select(Symbol, baseMean, log2FoldChange, lfcSE, pvalue, padj)
# DE9_5 <- read.delim("data/DESeq2/DESeq2_DEGs_lfcShrink_GSE149050_pmn.txt") |>
#   dplyr::select(Symbol, baseMean, log2FoldChange, lfcSE, pvalue, padj)
# DE9_6 <- read.delim("data/DESeq2/DESeq2_DEGs_lfcShrink_GSE149050_t cells.txt") |>
#   dplyr::select(Symbol, baseMean, log2FoldChange, lfcSE, pvalue, padj)
# DE10 <- read.delim("data/DESeq2/DESeq2_DEGs_lfcShrink_GSE162828.txt") |>
#   dplyr::select(Symbol, baseMean, log2FoldChange, lfcSE, pvalue, padj)
# DE11 <- read.delim("data/DESeq2/DESeq2_DEGs_lfcShrink_GSE110685.txt") |>
#   dplyr::select(Symbol, baseMean, log2FoldChange, lfcSE, pvalue, padj)


# save(list = ls(), file = './data/all_DESeq2_results.Rdata')

load('./data/all_DESeq2_results.Rdata')


# HPA and so on -----------------------------------------------------------


# Read in human protein atlas data
bulk <- read.csv("rna_tissue_exp_matrix.csv")
single <- read.csv("rna_sc_exp_matrix.csv")
bulk_pair <- read.csv("bulk_pair.csv")
sc_pair <- read.csv("sc_pair.csv")


# single cell -------------------------------------------------------------


# Read in single cell data
sc.combined <- readRDS("sc_combined.rds")
sc_exp <- as.matrix(sc.combined[["RNA"]]@counts)
sc_meta <- sc.combined@meta.data
# sc1<-readRDS("GSE135779_child.rds")



# Server ------------------------------------------------------------------


## Server Section
shinyServer(function(input, output, session) {
  # Plot for GSE131525
  output$plot1 <- renderPlot({
    gene <- toupper(input$text) # Get the gene name interested
    Gen1 <- data.frame(t(expr1[expr1$gene == gene, -1]))
    colnames(Gen1) <- c("Expression")
    mode(Gen1$Expression) <- "numeric"
    Gen1$name <- rownames(Gen1) # get the expression information and manage the expression data
    Genetable1 <- merge(Gen1, meta1, by = "name") # merge the expression & meta data
    p1 <- ggplot(Genetable1, aes(x = Cell.Type, y = Expression, fill = Disease)) +
      geom_boxplot() +
      labs(y = "Expression", title = "GSE131525") +
      theme(legend.position = "bottom")
    p1 # give the plot with ggplot2
  })
  # Plot for GSE149050
  output$plot2 <- renderPlot({
    gene <- toupper(input$text)
    Gen2 <- data.frame(t(expr2[expr2$gene == gene, -1]))
    colnames(Gen2) <- c("Expression")
    mode(Gen2$Expression) <- "numeric"
    Gen2$name <- rownames(Gen2)
    Genetable2 <- merge(Gen2, meta2, by = "name")
    p2 <- ggplot(Genetable2, aes(x = Cell.Type, y = Expression, fill = Disease.Subtype)) +
      geom_boxplot() +
      labs(y = "Expression", title = "GSE149050") +
      theme(legend.position = "bottom")
    p2
  })
  # Plot for GSE117836
  output$plot3 <- renderPlot({
    gene <- toupper(input$text)
    Gen3 <- data.frame(t(expr3[expr3$gene == gene, -1]))
    colnames(Gen3) <- c("Expression")
    mode(Gen3$Expression) <- "numeric"
    Gen3$name <- rownames(Gen3)
    Genetable3 <- merge(Gen3, meta3, by = "name")
    p3 <- ggplot(Genetable3, aes(x = Cell.Type, y = Expression, fill = Disease.Subtype)) +
      geom_boxplot() +
      labs(y = "Expression", title = "GSE117836")
    p3
  })
  # Plot for GSE139358
  output$plot4 <- renderPlot({
    gene <- toupper(input$text)
    Gen4 <- data.frame(t(expr4[expr4$gene == gene, -1]))
    colnames(Gen4) <- c("Expression")
    mode(Gen4$Expression) <- "numeric"
    Gen4$name <- rownames(Gen4)
    Genetable4 <- merge(Gen4, meta4, by = "name")
    p4 <- ggplot(Genetable4, aes(x = Cell.Type, y = Expression, fill = Disease)) +
      geom_boxplot() +
      labs(y = "Expression", title = "GSE139358")
    p4
  })
  # Plot for GSE164457
  output$plot5 <- renderPlot({
    gene <- toupper(input$text)
    Gen5 <- data.frame(t(expr5[expr5$gene == gene, -1]))
    colnames(Gen5) <- c("Expression")
    mode(Gen5$Expression) <- "numeric"
    Gen5$name <- rownames(Gen5)
    Genetable5 <- merge(Gen5, meta5, by = "name")
    p5 <- ggplot(Genetable5, aes(x = Cell.Type, y = Expression, fill = Disease)) +
      geom_boxplot() +
      labs(y = "Expression", title = "GSE164457")
    p5
  })
  # Plot for GSE92387
  output$plot6 <- renderPlot({
    gene <- toupper(input$text)
    Gen6 <- data.frame(t(expr6[expr6$gene == gene, -1]))
    colnames(Gen6) <- c("Expression")
    mode(Gen6$Expression) <- "numeric"
    Gen6$name <- rownames(Gen6)
    Genetable6 <- merge(Gen6, meta6, by = "name")
    p6 <- ggplot(Genetable6, aes(x = Cell.Type, y = Expression, fill = Disease)) +
      geom_boxplot() +
      labs(y = "Expression", title = "GSE92387")
    p6
  })
  # Plot for GSE110999
  output$plot7 <- renderPlot({
    gene <- toupper(input$text)
    Gen7 <- data.frame(t(expr7[expr7$gene == gene, -1]))
    colnames(Gen7) <- c("Expression")
    mode(Gen7$Expression) <- "numeric"
    Gen7$name <- rownames(Gen7)
    Genetable7 <- merge(Gen7, meta7, by = "name")
    p7 <- ggplot(Genetable7, aes(x = Cell.Type, y = Expression, fill = Disease)) +
      geom_boxplot() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(y = "Expression", title = "GSE110999")
    p7
  })
  # Plot for GSE163121
  output$plot8 <- renderPlot({
    gene <- toupper(input$text)
    Gen8 <- data.frame(t(expr8[expr8$gene == gene, -1]))
    colnames(Gen8) <- c("Expression")
    mode(Gen8$Expression) <- "numeric"
    Gen8$name <- rownames(Gen8)
    Genetable8 <- merge(Gen8, meta8, by = "name")
    p8 <- ggplot(Genetable8, aes(x = Cell.Type, y = Expression, fill = Disease)) +
      geom_boxplot() +
      labs(y = "Expression", title = "GSE163121")
    p8
  })
  # Plot for GSE136731
  output$plot9 <- renderPlot({
    gene <- toupper(input$text)
    Gen9 <- data.frame(t(expr9[expr9$gene == gene, -1]))
    colnames(Gen9) <- c("Expression")
    mode(Gen9$Expression) <- "numeric"
    Gen9$name <- rownames(Gen9)
    Genetable9 <- merge(Gen9, meta9, by = "name")
    p9 <- ggplot(Genetable9, aes(x = Cell.Type, y = Expression, fill = Disease)) +
      geom_boxplot() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(y = "Expression", title = "GSE136731")
    p9
  })
  # Plot for GSE118254
  output$plot10 <- renderPlot({
    gene <- toupper(input$text)
    Gen10 <- data.frame(t(expr10[expr10$gene == gene, -1]))
    colnames(Gen10) <- c("Expression")
    mode(Gen10$Expression) <- "numeric"
    Gen10$name <- rownames(Gen10)
    Genetable10 <- merge(Gen10, meta10, by = "name")
    p10 <- ggplot(Genetable10, aes(x = Cell.Type, y = Expression, fill = Disease)) +
      geom_boxplot() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom") +
      labs(y = "Expression", title = "GSE118254")
    p10
  })
  # Plot for GSE122459
  output$plot13 <- renderPlot({
    gene <- toupper(input$text)
    Gen13 <- data.frame(t(expr13[expr13$gene == gene, -1]))
    colnames(Gen13) <- c("Expression")
    mode(Gen13$Expression) <- "numeric"
    Gen13$name <- rownames(Gen13)
    Genetable13 <- merge(Gen13, meta13, by = "name")
    p13 <- ggplot(Genetable13, aes(x = Source, y = Expression, fill = Disease)) +
      geom_boxplot() +
      labs(y = "Expression", title = "GSE122459")
    p13
  })
  # Plot for GSE112087
  output$plot15 <- renderPlot({
    gene <- toupper(input$text)
    Gen15 <- data.frame(t(expr15[expr15$gene == gene, -1]))
    colnames(Gen15) <- c("Expression")
    mode(Gen15$Expression) <- "numeric"
    Gen15$name <- rownames(Gen15)
    Genetable15 <- merge(Gen15, meta15, by = "name")
    p15 <- ggplot(Genetable15, aes(x = Source, y = Expression, fill = Disease)) +
      geom_boxplot() +
      labs(y = "Expression", title = "GSE112087")
    p15
  })
  # Plot for GSE162828
  output$plot16 <- renderPlot({
    gene <- toupper(input$text)
    Gen16 <- data.frame(t(expr16[expr16$gene == gene, -1]))
    colnames(Gen16) <- c("Expression")
    mode(Gen16$Expression) <- "numeric"
    Gen16$name <- rownames(Gen16)
    Genetable16 <- merge(Gen16, meta16, by = "name")
    p16 <- ggplot(Genetable16, aes(x = Source, y = Expression, fill = Disease)) +
      geom_boxplot() +
      labs(y = "Expression", title = "GSE162828")
    p16
  })
  # Plot for GSE110685
  output$plot17 <- renderPlot({
    gene <- toupper(input$text)
    Gen17 <- data.frame(t(expr17[expr17$gene == gene, -1]))
    colnames(Gen17) <- c("Expression")
    mode(Gen17$Expression) <- "numeric"
    Gen17$name <- rownames(Gen17)
    Genetable17 <- merge(Gen17, meta17, by = "name")
    p17 <- ggplot(Genetable17, aes(x = Source, y = Expression, fill = Disease)) +
      geom_boxplot() +
      labs(y = "Expression", title = "GSE110685")
    p17
  })
  
  # Plot for GSE80183
  output$plot20 <- renderPlot({
    gene <- toupper(input$text)
    Gen20 <- data.frame(t(expr20[expr20$gene == gene, -1]))
    colnames(Gen20) <- c("Expression")
    mode(Gen20$Expression) <- "numeric"
    Gen20$name <- rownames(Gen20)
    Genetable20 <- merge(Gen20, meta20, by = "name")
    p20 <- ggplot(Genetable20, aes(x = Source, y = Expression, fill = Disease.Subtype)) +
      geom_boxplot() +
      labs(y = "Expression", title = "GSE80183")
    p20
  })
  # Plot for GSE72509
  output$plot21 <- renderPlot({
    gene <- toupper(input$text)
    Gen21 <- data.frame(t(expr21[expr21$gene == gene, -1]))
    colnames(Gen21) <- c("Expression")
    mode(Gen21$Expression) <- "numeric"
    Gen21$name <- rownames(Gen21)
    Genetable21 <- merge(Gen21, meta21, by = "name")
    p21 <- ggplot(Genetable21, aes(x = Source, y = Expression, fill = Disease)) +
      geom_boxplot() +
      labs(y = "Expression", title = "GSE72509")
    p21
  })
  
  # Forest Plot
  output$plot27 <- renderPlot({
    gene <- input$text
    gene <- toupper(gene)
    DE <- rbind(
      DE1[DE1$Symbol == gene, ], DE1_1[DE1_1$Symbol == gene, ], DE1_2[DE1_2$Symbol == gene, ], 
      DE1_3[DE1_3$Symbol == gene, ],
      DE1_4[DE1_4$Symbol == gene, ], DE2[DE2$Symbol == gene, ], DE3[DE3$Symbol == gene, ], DE4[DE4$Symbol == gene, ],
      DE5[DE5$Symbol == gene, ], DE6[DE6$Symbol == gene, ], DE6_1[DE6_1$Symbol == gene, ], DE6_2[DE6_2$Symbol == gene, ],
      DE6_3[DE6_3$Symbol == gene, ], DE6_4[DE6_4$Symbol == gene, ], DE6_5[DE6_5$Symbol == gene, ], DE6_6[DE6_6$Symbol == gene, ],
      DE7[DE7$Symbol == gene, ], DE8[DE8$Symbol == gene, ], DE9[DE9$Symbol == gene, ], DE9_1[DE9_1$Symbol == gene, ],
      DE9_2[DE9_2$Symbol == gene, ], DE9_3[DE9_3$Symbol == gene, ], DE9_4[DE9_4$Symbol == gene, ], DE9_5[DE9_5$Symbol == gene, ],
      DE9_6[DE9_6$Symbol == gene, ], DE10[DE10$Symbol == gene, ], DE11[DE11$Symbol == gene, ]
    )
    rownames(DE) <- c(
      "GSE131525", "GSE131525_1", "GSE131525_2", "GSE131525_3", "GSE131525_4", "GSE72509",
      "GSE92387", "GSE110999", "GSE112087", "GSE118254", "GSE118254_1", "GSE118254_2",
      "GSE118254_3", "GSE118254_4", "GSE118254_5", "GSE118254_6", "GSE136731", "GSE139358",
      "GSE149050", "GSE149050_1", "GSE149050_2", "GSE149050_3", "GSE149050_4", "GSE149050_5",
      "GSE149050_6", "GSE162828", "GSE110685"
    )
    tissue <- c(
      "pbmcs", "pbmcs", "pbmcs", "pbmcs", "pbmcs", "whole blood", "pbmcs", "b cells", "blood",
      "circulating b cells", "circulating b cells", "circulating b cells", "circulating b cells",
      "circulating b cells", "circulating b cells", "circulating b cells",
      "pbmcs", "human whole blood", "pbmcs", "pbmcs", "pbmcs", "pbmcs", "pbmcs", "pbmcs", "pbmcs",
      "pbmcs", "human whole blood"
    )
    cell <- c("——", "b", "cd4", "cd8", "monocytes", "——", "——", "——", "——", "——", 
              "an", "asc", "dn2", "rn", "sm", "t3", "——", "——", "——", "b", "cdc", 
              "cmo", "pdc", "pmn", "t cells", "——", "——")
    df <- data.frame(
      study = rownames(DE), index = 1:27, effect = DE$log2FoldChange,
      lower = DE$log2FoldChange - 1.96 * DE$lfcSE,
      upper = DE$log2FoldChange + 1.96 * DE$lfcSE,
      logP = -log10(DE$pvalue),
      label = c(paste("P value = ", signif(DE$pvalue, 3), " FDR = ", signif(DE$padj, 3))),
      linety = factor(ifelse(DE$pvalue < 0.05, "solid", "longdash"), levels = c("solid", "longdash"))
    )
    
    forest <- ggplot(data = df, aes(y = index, x = effect, xmin = lower, xmax = upper, color = study)) +
      geom_point(aes(size = logP)) +
      geom_errorbarh(height = .1, aes(linetype = linety)) +
      scale_y_continuous(name = "", breaks = 1:nrow(df), labels = df$study) +
      scale_linetype_manual(values = c("solid", "longdash")) +
      geom_vline(aes(xintercept = 0), colour = "grey", linetype = "dashed") +
      theme_bw() +
      labs(x = "Log2 Fold Change", size = "-log10(P)", title = gene) +
      geom_label_repel(aes(label = label, color = study), point.padding = 10) +
      guides(color = "none", linetype = "none") +
      theme(legend.position = c(0.95, 0.11), legend.background = element_rect(fill = rgb(1, 1, 1, alpha = 0.001), color = NA))
    
    t1 <- ggplot(df) +
      geom_text(aes(y = index, x = 1, label = tissue), vjust = 0) +
      ggtitle("Tissue") +
      xlab("  ") +
      theme_classic(base_size = 13) +
      theme(
        axis.line = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(color = "white"),
        plot.title = element_text(size = 12, hjust = 0.5)
      )
    t2 <- ggplot(df) +
      geom_text(aes(y = index, x = 1, label = cell), vjust = 0) +
      ggtitle("Cell Type") +
      xlab("  ") +
      theme_classic(base_size = 13) +
      theme(
        axis.line = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(color = "white"),
        plot.title = element_text(size = 12, hjust = 0.5)
      )
    mlayout <- rbind(
      c(1, 1, 1, 1, 2, 3),
      c(1, 1, 1, 1, 2, 3),
      c(1, 1, 1, 1, 2, 3),
      c(1, 1, 1, 1, 2, 3),
      c(1, 1, 1, 1, 2, 3)
    )
    p <- grid.arrange(forest, t1, t2, layout_matrix = mlayout)
    p # rinder[ plot
  })
  
  output$plot28 <- renderPlot({
    gene <- input$text
    gene <- toupper(gene)
    # GSE131525
    DE <- rbind(
      DE1[DE1$Symbol == gene, ], DE1_1[DE1_1$Symbol == gene, ], DE1_2[DE1_2$Symbol == gene, ], 
      DE1_3[DE1_3$Symbol == gene, ],
      DE1_4[DE1_4$Symbol == gene, ]
    )
    rownames(DE) <- c("GSE131525", "GSE131525_1", "GSE131525_2", "GSE131525_3", "GSE131525_4")
    tissue <- c("pbmcs", "pbmcs", "pbmcs", "pbmcs", "pbmcs")
    cell <- c("——", "b", "cd4", "cd8", "monocytes")
    df <- data.frame(
      study = rownames(DE), index = 1:5, effect = DE$log2FoldChange,
      lower = DE$log2FoldChange - 1.96 * DE$lfcSE,
      upper = DE$log2FoldChange + 1.96 * DE$lfcSE,
      logP = -log10(DE$pvalue),
      label = c(paste("P value = ", signif(DE$pvalue, 3), " FDR = ", signif(DE$padj, 3))),
      linety = factor(ifelse(DE$pvalue < 0.05, "solid", "longdash"), levels = c("solid", "longdash"))
    )
    
    forest <- ggplot(data = df, aes(y = index, x = effect, xmin = lower, xmax = upper, color = study)) +
      geom_point(aes(size = logP)) +
      geom_errorbarh(height = .1, aes(linetype = linety)) +
      scale_y_continuous(name = "", breaks = 1:nrow(df), labels = df$study) +
      scale_linetype_manual(values = c("solid", "longdash")) +
      geom_vline(aes(xintercept = 0), colour = "grey", linetype = "dashed") +
      theme_bw() +
      labs(x = "Log2 Fold Change", size = "-log10(P)", title = gene) +
      geom_label_repel(aes(label = label, color = study), point.padding = 10) +
      guides(color = "none", linetype = "none") +
      theme(legend.position = c(0.95, 0.1), legend.background = element_rect(fill = rgb(1, 1, 1, 
                                                                                        alpha = 0.001), 
                                                                             color = NA))
    
    t1 <- ggplot(df) +
      geom_text(aes(y = index, x = 1, label = tissue), vjust = 0) +
      ggtitle("Tissue") +
      xlab("  ") +
      theme_classic(base_size = 13) +
      theme(
        axis.line = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(color = "white"),
        plot.title = element_text(size = 12, hjust = 0.5)
      )
    t2 <- ggplot(df) +
      geom_text(aes(y = index, x = 1, label = cell), vjust = 0) +
      ggtitle("Cell Type") +
      xlab("  ") +
      theme_classic(base_size = 13) +
      theme(
        axis.line = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(color = "white"),
        plot.title = element_text(size = 12, hjust = 0.5)
      )
    mlayout <- rbind(
      c(1, 1, 1, 1, 2, 3),
      c(1, 1, 1, 1, 2, 3),
      c(1, 1, 1, 1, 2, 3)
    )
    p <- grid.arrange(forest, t1, t2, layout_matrix = mlayout)
    p
  })
  output$plot29 <- renderPlot({
    gene <- input$text
    gene <- toupper(gene)
    # GSE149050
    DE <- rbind(DE9[DE9$Symbol == gene, ], DE9_1[DE9_1$Symbol == gene, ], DE9_2[DE9_2$Symbol == gene, ], 
                DE9_3[DE9_3$Symbol == gene, ], DE9_4[DE9$Symbol == gene, ], DE9_5[DE9_5$Symbol == gene, ], 
                DE9_6[DE9_6$Symbol == gene, ])
    rownames(DE) <- c(
      "GSE149050", "GSE149050_1", "GSE149050_2", "GSE149050_3", "GSE149050_4", "GSE149050_5",
      "GSE149050_6"
    )
    tissue <- c("pbmcs", "pbmcs", "pbmcs", "pbmcs", "pbmcs", "pbmcs", "pbmcs")
    cell <- c("——", "b", "cdc", "cmo", "pdc", "pmn", "t cells")
    df <- data.frame(
      study = rownames(DE), index = 1:7, effect = DE$log2FoldChange,
      lower = DE$log2FoldChange - 1.96 * DE$lfcSE,
      upper = DE$log2FoldChange + 1.96 * DE$lfcSE,
      logP = -log10(DE$pvalue),
      label = c(paste("P value = ", signif(DE$pvalue, 3), " FDR = ", signif(DE$padj, 3))),
      linety = factor(ifelse(DE$pvalue < 0.05, "solid", "longdash"), levels = c("solid", "longdash"))
    )
    
    forest <- ggplot(data = df, aes(y = index, x = effect, xmin = lower, xmax = upper, color = study)) +
      geom_point(aes(size = logP)) +
      geom_errorbarh(height = .1, aes(linetype = linety)) +
      scale_y_continuous(name = "", breaks = 1:nrow(df), labels = df$study) +
      scale_linetype_manual(values = c("solid", "longdash")) +
      geom_vline(aes(xintercept = 0), colour = "grey", linetype = "dashed") +
      theme_bw() +
      labs(x = "Log2 Fold Change", size = "-log10(P)", title = gene) +
      geom_label_repel(aes(label = label, color = study), point.padding = 10) +
      guides(color = "none", linetype = "none") +
      theme(legend.position = c(0.95, 0.1), legend.background = element_rect(fill = rgb(1, 1, 1, alpha = 0.001), 
                                                                             color = NA))
    
    t1 <- ggplot(df) +
      geom_text(aes(y = index, x = 1, label = tissue), vjust = 0) +
      ggtitle("Tissue") +
      xlab("  ") +
      theme_classic(base_size = 13) +
      theme(
        axis.line = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(color = "white"),
        plot.title = element_text(size = 12, hjust = 0.5)
      )
    t2 <- ggplot(df) +
      geom_text(aes(y = index, x = 1, label = cell), vjust = 0) +
      ggtitle("Cell Type") +
      xlab("  ") +
      theme_classic(base_size = 13) +
      theme(
        axis.line = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(color = "white"),
        plot.title = element_text(size = 12, hjust = 0.5)
      )
    mlayout <- rbind(
      c(1, 1, 1, 1, 2, 3),
      c(1, 1, 1, 1, 2, 3),
      c(1, 1, 1, 1, 2, 3)
    )
    p <- grid.arrange(forest, t1, t2, layout_matrix = mlayout)
    p
  })
  output$plot30 <- renderPlot({
    gene <- input$text
    gene <- toupper(gene)
    ## GSE118254
    DE <- rbind(DE6[DE6$Symbol == gene, ], DE6_1[DE6_1$Symbol == gene, ], DE6_2[DE6_2$Symbol == gene, ], 
                DE6_3[DE6_3$Symbol == gene, ], DE6_4[DE6_4$Symbol == gene, ], DE6_5[DE6_5$Symbol == gene, ], 
                DE6_6[DE6_6$Symbol == gene, ])
    rownames(DE) <- c("GSE118254", "GSE118254_1", "GSE118254_2", "GSE118254_3", "GSE118254_4", "GSE118254_5", 
                      "GSE118254_6")
    tissue <- c(
      "circulating b cells", "circulating b cells", "circulating b cells", "circulating b cells",
      "circulating b cells", "circulating b cells", "circulating b cells"
    )
    cell <- c("——", "an", "asc", "dn2", "rn", "sm", "t3")
    df <- data.frame(
      study = rownames(DE), index = 1:7, effect = DE$log2FoldChange,
      lower = DE$log2FoldChange - 1.96 * DE$lfcSE,
      upper = DE$log2FoldChange + 1.96 * DE$lfcSE,
      logP = -log10(DE$pvalue),
      label = c(paste("P value = ", signif(DE$pvalue, 3), " FDR = ", signif(DE$padj, 3))),
      linety = factor(ifelse(DE$pvalue < 0.05, "solid", "longdash"), levels = c("solid", "longdash"))
    )
    
    forest <- ggplot(data = df, aes(y = index, x = effect, xmin = lower, xmax = upper, color = study)) +
      geom_point(aes(size = logP)) +
      geom_errorbarh(height = .1, aes(linetype = linety)) +
      scale_y_continuous(name = "", breaks = 1:nrow(df), labels = df$study) +
      scale_linetype_manual(values = c("solid", "longdash")) +
      geom_vline(aes(xintercept = 0), colour = "grey", linetype = "dashed") +
      theme_bw() +
      labs(x = "Log2 Fold Change", size = "-log10(P)", title = gene) +
      geom_label_repel(aes(label = label, color = study), point.padding = 10) +
      guides(color = "none", linetype = "none") +
      theme(legend.position = c(0.95, 0.1), legend.background = element_rect(fill = rgb(1, 1, 1, alpha = 0.001), 
                                                                             color = NA))
    
    t1 <- ggplot(df) +
      geom_text(aes(y = index, x = 1, label = tissue), vjust = 0) +
      ggtitle("Tissue") +
      xlab("  ") +
      theme_classic(base_size = 13) +
      theme(
        axis.line = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(color = "white"),
        plot.title = element_text(size = 12, hjust = 0.5)
      )
    t2 <- ggplot(df) +
      geom_text(aes(y = index, x = 1, label = cell), vjust = 0) +
      ggtitle("Cell Type") +
      xlab("  ") +
      theme_classic(base_size = 13) +
      theme(
        axis.line = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(color = "white"),
        plot.title = element_text(size = 12, hjust = 0.5)
      )
    mlayout <- rbind(
      c(1, 1, 1, 1, 2, 3),
      c(1, 1, 1, 1, 2, 3),
      c(1, 1, 1, 1, 2, 3)
    )
    p <- grid.arrange(forest, t1, t2, layout_matrix = mlayout)
    p
  })
  output$plot31 <- renderPlot({
    gene <- input$text
    gene <- toupper(gene)
    # b cells
    DE <- rbind(DE1_1[DE1_1$Symbol == gene, ], DE6[DE6$Symbol == gene, ], DE9_1[DE9_1$Symbol == gene, ], 
                DE4[DE4$Symbol == gene, ], DE3[DE3$Symbol == gene, ])
    rownames(DE) <- c("GSE131525", "GSE118254", "GSE149050", "GSE110999", "GSE92387")
    
    df <- data.frame(
      study = rownames(DE), index = 1:5, effect = DE$log2FoldChange,
      lower = DE$log2FoldChange - 1.96 * DE$lfcSE,
      upper = DE$log2FoldChange + 1.96 * DE$lfcSE,
      logP = -log10(DE$pvalue),
      label = c(paste("P value = ", signif(DE$pvalue, 3), " FDR = ", signif(DE$padj, 3))),
      linety = factor(ifelse(DE$pvalue < 0.05, "solid", "longdash"), levels = c("solid", "longdash"))
    )
    
    forest <- ggplot(data = df, aes(y = index, x = effect, xmin = lower, xmax = upper, color = study)) +
      geom_point(aes(size = logP)) +
      geom_errorbarh(height = .1, aes(linetype = linety)) +
      scale_y_continuous(name = "", breaks = 1:nrow(df), labels = df$study) +
      scale_linetype_manual(values = c("solid", "longdash")) +
      geom_vline(aes(xintercept = 0), colour = "grey", linetype = "dashed") +
      theme_bw() +
      labs(x = "Log2 Fold Change", size = "-log10(P)", title = paste(gene, "from b cells")) +
      geom_label_repel(aes(label = label, color = study), point.padding = 10) +
      guides(color = "none", linetype = "none") +
      theme(legend.position = c(0.95, 0.1), legend.background = element_rect(fill = rgb(1, 1, 1, alpha = 0.001), 
                                                                             color = NA))
    forest
  })
  
  # Expression by Tissue
  output$plot32 <- renderPlot({
    gene <- input$text # read in the interested gene
    gene <- toupper(gene) # make the gene in upper case by default
    bulk_gene <- bulk[bulk$Gene.name == gene, ]
    bar_tissue <- as.data.frame(cbind(t(bulk_gene[, -1]), t(bulk_pair)))
    colnames(bar_tissue) <- c("Expression", "Organ")
    bar_tissue$Expression <- as.numeric(bar_tissue$Expression)
    bar_tissue$Tissue <- rownames(bar_tissue)
    Tissue_plt <- ggplot(bar_tissue, aes(x = reorder(Tissue, -Expression), y = Expression, fill = Organ)) +
      geom_bar(stat = "identity") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(x = "Human tissues", y = "Expression(nTPM)")
    Tissue_plt
  })
  
  # Expression by Cell
  output$plot33 <- renderPlot({
    gene <- input$text
    gene <- toupper(gene)
    single_gene <- single[single$Gene.name == gene, ]
    bar_cell <- as.data.frame(cbind(t(single_gene[, -1]), t(sc_pair)))
    colnames(bar_cell) <- c("Expression", "Tissue")
    bar_cell$Expression <- as.numeric(bar_cell$Expression)
    bar_cell$Cell <- rownames(bar_cell)
    Cell_plt <- ggplot(bar_cell, aes(x = reorder(Cell, -Expression), y = Expression, fill = Tissue)) +
      geom_bar(stat = "identity") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(x = "Human cell types", y = "Expression(nTPM)")
    Cell_plt
  })
  # violin plot for the combined data
  output$plot34 <- renderPlot({
    gene <- input$text
    gene <- toupper(gene)
    sc_gene <- sc_exp[rownames(sc_exp) == gene, ]
    sc_gene <- as.data.frame(sc_gene)
    sc_meta <- sc_meta[, c("orig.ident", "Author.s.cell.type", "Condition")]
    table2 <- merge(sc_gene, sc_meta, by = "row.names", all = T)
    colnames(table2) <- c("Row.names", "Expression", "Sample", "Cell.Type", "Condition")
    p <- ggplot(table2, aes(x = Cell.Type, y = Expression, fill = Condition)) +
      geom_violin() +
      labs(y = "Expression", x = "Cell Type") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    p
  })
  # umap for the combined data
  output$plot35 <- renderPlot({
    p1 <- DimPlot(sc.combined, reduction = "umap", group.by = "Condition") + labs(title = "Disease Status")
    p2 <- DimPlot(sc.combined, reduction = "umap", group.by = "Author.s.cell.type") + labs(title = "Cell Types")
    p1 + p2
  })
  # #violin plot for GSE135779_child
  # output$plot36 <- renderPlot({
  #   gene=input$text
  #   gene<-toupper(gene)
  #   sc1_exp<-as.matrix(sc1[["RNA"]]@counts)
  #   sc1_meta<-sc1@meta.data
  #   sc1_gene<-sc1_exp[rownames(sc1_exp)==gene,]
  #   sc1_gene<-as.data.frame(sc1_gene)
  #   sc1_meta<-sc1_meta[,c("orig.ident","Author.s.cell.type","Condition")]
  #   table1<-merge(sc1_gene,sc1_meta, by = "row.names", all = T)
  #   colnames(table1)<-c("Row.names","Expression","Sample","Cell.Type","Condition")
  #   p<-ggplot(table1, aes(x=Cell.Type,y=Expression,fill=Condition))+
  #     geom_violin()+
  #     labs(y="Expression", x= "Cell Type")+
  #     theme(axis.text.x = element_text(angle = 45, hjust = 1))
  #   p
  # })
  # #umap for GSE135779_child
  # output$plot37 <- renderPlot({
  #   p1 <- DimPlot(sc1, reduction = "umap", group.by = "Condition",raster=FALSE)
  #   p2 <- DimPlot(sc1, reduction = "umap", group.by = "Author.s.cell.type",raster=FALSE)
  #   p1+p2
  # })
})
