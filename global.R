


VERSION <- "0.0.2"

# read file global --------------------------------------------------------



library(stringr)

combined_DESeq2_res <- readRDS('./data/combined_DESeq2_res.rds')


hpa_transcript_rna_immunecells <- readr::read_csv(
  './data/HPA/hpa_transcript_rna_immunecells.csv',
  show_col_types = FALSE
)

# 从以前文件拷过来的
# load('./data/pandaomics_meta_expr_data.Rdata')



chioce1 <- unique(combined_DESeq2_res$GSE_name)
chioce1 <- chioce1[str_length(chioce1) < 10]



gse_annotation <-
  tibble::tibble(
    GSE = c(
      "GSE131525", "GSE131525_b", "GSE131525_cd4", "GSE131525_cd8", "GSE131525_monocytes", "GSE72509",
      "GSE92387", "GSE110999", "GSE112087", "GSE118254", 
      "GSE118254_activated_naive_b_cells_(an)", "GSE118254_antigen_secreting_cell_(asc)",
      "GSE118254_double_negative_b_cells_(dn2)", "GSE118254_resting_naive_b_cells_(rn)", 
      "GSE118254_switched_memory_b_cells_(sm)", "GSE118254_transitional_3_b_cells_(t3)", "GSE136731", "GSE139358",
      "GSE149050", "GSE149050_b_cells", "GSE149050_cdc", "GSE149050_cmo", 
      "GSE149050_pdc", "GSE149050_5",
      "GSE149050_t_cells", "GSE162828", "GSE110685"
    ),
    tissue = c(
      "pbmcs", "pbmcs", "pbmcs", "pbmcs", "pbmcs", "whole blood", "pbmcs", "b cells", "blood",
      "circulating b cells", "circulating b cells", "circulating b cells", "circulating b cells",
      "circulating b cells", "circulating b cells", "circulating b cells",
      "pbmcs", "human whole blood", "pbmcs", "pbmcs", "pbmcs", "pbmcs", "pbmcs", "pbmcs", "pbmcs",
      "pbmcs", "human whole blood"
    ),
    cell = c(
      "——", "b", "cd4", "cd8", "monocytes", "——", "——", "——", "——", "——",
      "an", "asc", "dn2", "rn", "sm", "t3", "——", "——", "——", "b", "cdc",
      "cmo", "pdc", "pmn", "t cells", "——", "——"
    )
  )



# scRNA GSE number

chioce_IBD <- c('PMID34497389', 'GSE164985', 'GSE202052', 'GSE125527')
chioce_SLE <- c('GSE174188', 'GSE163121')
chioce_Psoriasis <- c('PMID35958578', 'GSE194315')
chioce_COPD <- c('GSE168191', 'GSE136831', 'GSE173896')
chioce_Others <- c('PMID35383111', 'GSE128033', 'GSE195452', 'GSE200815')


st <- readxl::read_excel("./data/talk2data/scRNA_metadata.xlsx"
                         )


# function for sever ------------------------------------------------------





# functions for UI --------------------------------------------------------







