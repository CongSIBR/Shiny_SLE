





library(tidyverse)


# tidy DESeq2 data --------------------------------------------------------



all_deseq2 <- list.files(path = './data/DESeq2',
           pattern = 'DESeq2_DEGs_lfcShrink_.*.txt',
           full.names = TRUE
           )


# str_match('./data/DESeq2/DESeq2_DEGs_lfcShrink_GSE110685.txt', 
#            pattern = 'lfcShrink_\\s*(.*?)\\s*.txt')


make_file <- function(f){
  
  f_name <- str_match(f, 
                      pattern = 'lfcShrink_\\s*(.*?)\\s*.txt')[,2]
  
  print(glue::glue('Current file is {f_name}'))
  
  res <- read_tsv(f, show_col_types = FALSE) %>% 
    dplyr::select(Symbol, baseMean, log2FoldChange, lfcSE, pvalue, padj) %>% 
    dplyr::mutate(GSE_name = f_name)
  
  problems(res)
  
  res
}


combined_data <- map(all_deseq2, make_file) %>% bind_rows()


combined_DESeq2_res <- combined_data


# saveRDS(combined_DESeq2_res, file = './data/combined_DESeq2_res.rds')


readRDS('./data/combined_DESeq2_res.rds')



# tidy PandaOmics data ----------------------------------------------------

all_pandaomics_meta <- list.files(path = './data/',
                         pattern = 'GSE.*?_metadata.csv',
                         full.names = TRUE
                         )


all_pandaomics_expr <- list.files(path = './data/',
                                  pattern = 'GSE.*?_expression.csv',
                                  full.names = TRUE
)




# additional meta data ----------------------------------------------------


gse_meta <- tibble(
  GSE = c(
    "GSE131525", "GSE131525_1", "GSE131525_2", "GSE131525_3", "GSE131525_4", "GSE72509",
    "GSE92387", "GSE110999", "GSE112087", "GSE118254", "GSE118254_1", "GSE118254_2",
    "GSE118254_3", "GSE118254_4", "GSE118254_5", "GSE118254_6", "GSE136731", "GSE139358",
    "GSE149050", "GSE149050_1", "GSE149050_2", "GSE149050_3", "GSE149050_4", "GSE149050_5",
    "GSE149050_6", "GSE162828", "GSE110685"
  ),
  tissue = c(
    "pbmcs", "pbmcs", "pbmcs", "pbmcs", "pbmcs", "whole blood", "pbmcs", "b cells", "blood",
    "circulating b cells", "circulating b cells", "circulating b cells", "circulating b cells",
    "circulating b cells", "circulating b cells", "circulating b cells",
    "pbmcs", "human whole blood", "pbmcs", "pbmcs", "pbmcs", "pbmcs", "pbmcs", "pbmcs", "pbmcs",
    "pbmcs", "human whole blood"
  ),
  cell = c("——", "b", "cd4", "cd8", "monocytes", "——", "——", "——", "——", "——", 
            "an", "asc", "dn2", "rn", "sm", "t3", "——", "——", "——", "b", "cdc", 
            "cmo", "pdc", "pmn", "t cells", "——", "——")
)









