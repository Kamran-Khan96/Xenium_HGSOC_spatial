library(tidyverse)
library(tidybulk)
library(sva)

patient_list <- read_csv("tma5_6_metadata.csv") %>%
  pull(sample) %>% unique()

log <- read_csv("tma5_tumour_decoupler_coarse_mean_bulk.csv") %>%
  dplyr::filter(cell_type == "Tumour") %>%
  dplyr::select(sample, symbol, count) %>%
  distinct() %>% mutate(TMA = "TMA5") %>%
  bind_rows(read_csv("tma6_tumour_decoupler_coarse_mean_bulk.csv") %>%
              dplyr::filter(cell_type == "Tumour") %>%
              dplyr::select(sample, symbol, count) %>%
              distinct() %>% mutate(TMA = "TMA6")) %>%
  dplyr::filter(sample %in% patient_list) %>%
  mutate(log_count = log2(count + 1))

wide_data <- log %>%
  select(symbol, sample, log_count) %>%
  pivot_wider(names_from = sample, values_from = log_count, values_fill = 0) %>%
  column_to_rownames(var = "symbol")

expr_matrix <- as.matrix(wide_data)

batch_info <- log %>%
  select(sample, TMA) %>%
  distinct() %>%
  # Match precisely to the columns of the wide matrix
  slice(match(colnames(expr_matrix), sample)) 

batch_vector <- batch_info$TMA

combat_adjusted_matrix <- ComBat(
  dat = expr_matrix,
  batch = batch_vector,
  par.prior = FALSE
)

metadata_lookup <- log %>% 
  select(sample, symbol, TMA) %>% 
  distinct()

merged_tma_adjusted <- as.data.frame(combat_adjusted_matrix) %>%
  rownames_to_column(var = "symbol") %>%
  pivot_longer(cols = -symbol, names_to = "sample", values_to = "adjusted_log_count") %>%
  left_join(metadata_lookup, by = c("sample", "symbol")) %>%
  group_by(symbol) %>%
  mutate( 
    adjusted_log_count = as.numeric(scale(adjusted_log_count)) 
  )
