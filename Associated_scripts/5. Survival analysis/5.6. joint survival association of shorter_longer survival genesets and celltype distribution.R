library(tidyverse)
library(tidybulk)
library(patchwork)
library(ggpubr)
library(survminer)
library(survival)
library(foreach)
library(ggsci)
library(stringr)
library(scales)
library(lubridate)
library(sva)
library(forcats)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ggnewscale)


my_theme <- theme_bw() +
  theme(text=element_text(size=12),
        legend.position = "bottom",
        axis.line = element_line(size=0.1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0), size = 12),
        axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0), size = 12),
        panel.spacing.x=unit(1, "lines"),
        axis.text.x = element_text(size=12, 
                                   #hjust=1, angle=45
        ),
        axis.text.y = element_text(size=12),
        strip.text.x = element_text(),
        strip.text.y = element_text(),
  )



#------------------------------------------------------------
# Spearman's correlation between survival associated genes and cell-types
#------------------------------------------------------------


frac_ct_rads <- read_csv("tma5_sqpy_coarse_tumour_touch_count_frac_diff_radii_fine_labels.csv")  %>%
  dplyr::filter(tumour_type == "Coarse") %>%
  dplyr::select(sample = PatientID, ct = "neighbour_type", count = fraction,  radius_um) %>%
  complete(sample, ct, fill = list(count = 0)) %>% 
  bind_rows(read_csv("tma6_sqpy_coarse_tumour_touch_count_frac_diff_radii_fine_labels.csv")  %>%
              dplyr::filter(tumour_type == "Coarse") %>%
              dplyr::select(sample = PatientID, ct = "neighbour_type", count = fraction, radius_um) %>%
              complete(sample, ct, fill = list(count = 0))) %>%
  dplyr::filter(!radius_um %in% c(3, 5, 75)) %>%
  filter(sample %in% patient_list)


gene_ct_corr_func <- function(rad){
  
  merged_hr_gene_ct <- merged_tma_adjusted %>%
    dplyr::filter(symbol %in% hr_genes) %>%
    dplyr::select(-TMA) %>%
    pivot_wider(names_from = symbol, values_from = adjusted_log_count) %>%
    inner_join(frac_ct_rads %>%
                 filter(radius_um == rad) %>%
                 filter(ct %in% (all_results %>%
                                   filter(radius == rad & P_tma5_6 < 0.05) %>%
                                   pull(cell_type) %>%
                                   unique()
                 )
                 ) %>%
                 dplyr::select(sample, ct, count) %>%
                 pivot_wider(names_from = ct, values_from = count)) %>%
    dplyr::select(-sample)
  
  
  
  actual_genes <- intersect(colnames(merged_hr_gene_ct), hr_genes)
  actual_cts   <- setdiff(colnames(merged_hr_gene_ct), actual_genes)
  
  cor_matrix <- cor(merged_hr_gene_ct, method = "spearman")
  res1       <- corrplot::cor.mtest(merged_hr_gene_ct, conf.level = 0.95, method = "spearman")
  p_matrix   <- res1$p
  rownames(p_matrix) <- rownames(cor_matrix)
  colnames(p_matrix) <- colnames(cor_matrix)
  
  r_sub <- cor_matrix[actual_genes, actual_cts, drop = FALSE]
  p_sub <- p_matrix[actual_genes, actual_cts, drop = FALSE]
  
  tidy_cors_data <- as.data.frame(as.table(r_sub)) %>%
    dplyr::rename(gene = Var1, cell_type = Var2, r = Freq) %>%
    left_join(
      as.data.frame(as.table(p_sub)) %>%
        dplyr::rename(gene = Var1, cell_type = Var2, p = Freq),
      by = c("gene", "cell_type")
    ) %>%
    mutate(
      gene = as.character(gene),
      cell_type = as.character(cell_type)
    ) %>%
    mutate(ct_radius = rad)
}


target_radii <- c(10, 25, 50, 100, 250, 500)

master_ct_gene_corr_df <- target_radii %>%
  purrr::map_df(~ gene_ct_corr_func(rad = .x)) %>%
  dplyr::arrange(p) %>%
  mutate(adj_spearman_p = p.adjust(p, method = "BH"))



#----------------------------------
# plotting correlation heatmap


custom_ct_order <- c(
  "B_Plasma cells IGHG+",
  "T cells CXCR4+",
  "T cells IL32+",
  "T cells CXCL13+",
  "T cells CCL5+",
  "Monocyte_Macrophages CXCL9+",
  "T cells Treg",
  "B_Plasma cells IGHM+",
  "Mesothelial cells",
  "Fibroblast - EGR1+",
  "Fibroblast - IGFBP5+"
)

p <- master_ct_gene_corr_df %>%
  filter(gene %in% en_hr_genes & ct_radius == 10) %>%
  mutate(
    is_sig = ifelse(p < 0.05, TRUE, FALSE),
    label_text = ifelse(is_sig, round(r, digits = 2), ""),
    plot_r = ifelse(is_sig, r, NA) 
  ) %>% 
  as_tibble() %>%
  mutate(gene_hr_type = ifelse(gene %in% good_hr_genes, "Longer-survival", "Shorter-survival")) %>%
  mutate(cell_type = factor(cell_type, levels = custom_ct_order)) %>%
  ggplot(aes(y = gene, x = cell_type, fill = plot_r)) + 
  geom_tile(color = "grey92", linewidth = 0.2) + 
  geom_text(aes(label = label_text), color = "black", size = 3) +
  scale_fill_gradient2(
    low = "#007C7A",
    mid = "white",       
    high = "#E66101",  
    midpoint = 0,        
    name = "Corr",
    na.value = "white"   
  ) +
  facet_grid(gene_hr_type ~ ., scales = "free", space = "free") + 
  my_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank(),
    legend.position = "right",
    plot.margin = margin(10, 10, 10, 10)
  )



#-------------------------------------------
# joint Cox regression between survival associated genes and cell-types (Figure 4B)
#-------------------------------------------

gene_df <- merged_tma_adjusted %>%
  dplyr::filter(symbol %in% hr_genes) %>%
  dplyr::select(-TMA) %>%
  tidyr::pivot_wider(names_from = symbol, values_from = adjusted_log_count)


gene_ct_pairwise_multivariate <- function(rad) {
  
  pairs_df <- master_ct_gene_corr_df %>%
    dplyr::filter(ct_radius == rad) %>%
    dplyr::select(gene, cell_type)
  
  if(nrow(pairs_df) == 0) return(NULL)
  
  
  ct_df <- frac_ct_rads %>%
    dplyr::filter(radius_um == rad & ct %in% unique(pairs_df$cell_type)) %>%
    dplyr::select(sample, ct, count) %>%
    tidyr::pivot_wider(names_from = ct, values_from = count, values_fill = 0)
  
  
  df_model <- ova_tma_clin %>%
    dplyr::filter(sample %in% patient_list) %>%
    inner_join(gene_df, by = "sample") %>%
    inner_join(ct_df, by = "sample")
  
  
  results <- pairs_df %>%
    split(seq(nrow(.))) %>%
    map_df(function(pair) {
      
      t_gene <- pair$gene
      t_ct   <- pair$cell_type
      
      
      if(!(t_gene %in% colnames(df_model)) || !(t_ct %in% colnames(df_model))) return(NULL)
      
      
      df_subset <- df_model %>%
        mutate(
          scaled_gene = as.numeric(scale(.data[[t_gene]])),
          scaled_ct   = as.numeric(scale(.data[[t_ct]]))
        )
      
      
      fit <- tryCatch({
        coxph(Surv(total_living_days, vital_status) ~ scaled_gene + scaled_ct, data = df_subset)
      }, error = function(e) return(NULL))
      
      if (is.null(fit)) return(NULL)
      
      
      tidy_fit <- broom::tidy(fit, exponentiate = TRUE)
      
      gene_res <- tidy_fit %>% dplyr::filter(term == "scaled_gene")
      ct_res   <- tidy_fit %>% dplyr::filter(term == "scaled_ct")
      
      data.frame(
        radius_um = rad,
        gene      = t_gene,
        cell_type = t_ct,
        HR_gene   = round(gene_res$estimate, 3),
        P_gene    = gene_res$p.value,
        HR_ct     = round(ct_res$estimate, 3),
        P_ct      = ct_res$p.value
      )
    })
  
  return(results)
}

target_radii <- c(10, 25, 50, 100, 250, 500)

multivariate_pairwise_results <- target_radii %>%
  purrr::map_df(~ gene_ct_pairwise_multivariate(rad = .x)) %>%
  dplyr::arrange(radius_um, P_gene, P_ct)



#--------------------------------------------------------------------------
## plotting heatmap

good_hr_genes <- c("TNFRSF11B", "LTF", "PTN", "CD44", "CHI3L1")

#==============================================================================

plot_df1 <- multivariate_pairwise_results %>%
  dplyr::filter(radius_um == 10 & gene %in% en_hr_genes) %>%
  mutate(
    sig_stars = case_when(
      P_gene < 0.001 ~ "***", P_gene < 0.01  ~ "**", P_gene < 0.05  ~ "*", TRUE ~ ""
    ),
    plot_HR = ifelse(P_gene < 0.05, HR_gene, NA),
    gene_hr_type = ifelse(gene %in% good_hr_genes, "Longer-survival", "Shorter-survival"),
    cell_type = factor(cell_type, levels = custom_ct_order) 
  )

p1 <- ggplot(plot_df1, aes(x = cell_type, y = gene, fill = plot_HR)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sig_stars), color = "black", size = 5, vjust = 0.75) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 1, 
    name = "Gene HR\n(Adjusted)", na.value = "white"
  ) +
  facet_grid(gene_hr_type ~ ., scales = "free_y", space = "free_y") +
  my_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(
    x = "Adjusted for (Covariate Cell Type)",
    y = "Evaluated Gene"
  )

# ==============================================================================

plot_df2 <- multivariate_pairwise_results %>%
  dplyr::filter(radius_um == 10 & gene %in% en_hr_genes) %>%
  mutate(
    sig_stars = case_when(
      P_ct < 0.001 ~ "***", P_ct < 0.01  ~ "**", P_ct < 0.05  ~ "*", TRUE ~ ""
    ),
    plot_HR = ifelse(P_ct < 0.05, HR_ct, NA),
    gene_hr_type = ifelse(gene %in% good_hr_genes, "Longer-survival", "Shorter-survival"),
    cell_type = factor(cell_type, levels = rev(custom_ct_order))
  )

p2 <- ggplot(plot_df2, aes(x = gene, y = cell_type, fill = plot_HR)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sig_stars), color = "black", size = 5, vjust = 0.75) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 1, 
    name = "Cell Type HR\n(Adjusted)", na.value = "white"
  ) +
  facet_grid(. ~ gene_hr_type, scales = "free_x", space = "free_x") +
  my_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(
    x = "Adjusted for (Covariate Gene)",
    y = "Evaluated Cell Type"
  )














#--------------------------------------------
# Combination KM curves on Fig4C-E
#--------------------------------------------

# Helper function to extract, dichotomize, fit Cox models, annotate all 6 comparisons on plot, and build KM curves

run_spatial_combo_analysis <- function(gene_name, cell_type_name, plot_legend_title) {
  
  gene_df <- bulk_exp_data %>%
    filter(symbol == gene_name) %>%
    { if ("expression" %in% colnames(.)) dplyr::rename(., expr_val = expression)
      else if ("log2_count" %in% colnames(.)) dplyr::rename(., expr_val = log2_count)
      else dplyr::rename(., expr_val = 3) } %>% 
    mutate(gene_status = ifelse(expr_val > median(expr_val, na.rm = TRUE), "High", "Low")) %>%
    select(sample, gene_status)
  
  spatial_df <- frac_ct_rads %>%
    filter(radius_um == "10", ct == cell_type_name) %>%
    select(sample, ct, count) %>%
    group_by(ct) %>%
    mutate(ct_status = ifelse(count > median(count, na.rm = TRUE), "High", "Low")) %>%
    ungroup() %>%
    select(sample, ct_status)
  
  analysis_df <- ova_tma_clin %>%
    filter(sample %in% patient_list) %>%
    left_join(gene_df, by = "sample") %>%
    left_join(spatial_df, by = "sample") %>%
    filter(!is.na(gene_status), !is.na(ct_status)) %>%
    mutate(fraction = factor(
      paste(gene_status, ct_status, sep = "/"),
      levels = c("Low/Low", "Low/High", "High/Low", "High/High")
    ))
  
 
  cox_A <- summary(coxph(Surv(total_living_days, vital_status) ~ fraction, data = analysis_df))
  
  df_refLH <- analysis_df %>% 
    mutate(frac_LH = factor(fraction, levels = c("Low/High", "High/Low", "High/High", "Low/Low")))
  cox_B <- summary(coxph(Surv(total_living_days, vital_status) ~ frac_LH, data = df_refLH))
  
  df_refHL <- analysis_df %>% 
    mutate(frac_HL = factor(fraction, levels = c("High/Low", "High/High", "Low/Low", "Low/High")))
  cox_C <- summary(coxph(Surv(total_living_days, vital_status) ~ frac_HL, data = df_refHL))
  
  get_stats <- function(cox_obj, param_name) {
    hr <- round(cox_obj$conf.int[param_name, "exp(coef)"], 2)
    p_val <- format_raw_p(cox_obj$coefficients[param_name, "Pr(>|z|)"])
    list(hr = hr, p = p_val)
  }
  
  c1 <- get_stats(cox_A, "fractionLow/High")   # L/H vs L/L
  c2 <- get_stats(cox_A, "fractionHigh/Low")   # H/L vs L/L
  c3 <- get_stats(cox_A, "fractionHigh/High")  # H/H vs L/L
  c4 <- get_stats(cox_B, "frac_LHHigh/Low")    # H/L vs L/H
  c5 <- get_stats(cox_B, "frac_LHHigh/High")   # H/H vs L/H
  c6 <- get_stats(cox_C, "frac_HLHigh/High")   # H/H vs H/L
  
  label_block_6 <- paste0(
    "HR (L/H vs L/L): ", c1$hr, " (p", c1$p, ")\n",
    "HR (H/L vs L/L): ", c2$hr, " (p", c2$p, ")\n",
    "HR (H/H vs L/L): ", c3$hr, " (p", c3$p, ")\n",
    "HR (H/L vs L/H): ", c4$hr, " (p", c4$p, ")\n",
    "HR (H/H vs L/H): ", c5$hr, " (p", c5$p, ")\n",
    "HR (H/H vs H/L): ", c6$hr, " (p", c6$p, ")"
  )
  
  fit <- survfit(Surv(total_living_days, vital_status) ~ fraction, data = analysis_df)
  names(fit$strata) <- gsub("fraction=", "", names(fit$strata))
  fit$strata <- fit$strata[intersect(c("Low/Low", "Low/High", "High/Low", "High/High"), names(fit$strata))]
  
  surv_p <- ggsurvplot(
    fit, data = analysis_df, risk.table = TRUE, pval = TRUE,
    legend.title = plot_legend_title, xlab = "Days since diagnosis",
    legend.labs = c("Low/Low", "Low/High", "High/Low", "High/High"),
    xlim = c(0, 3650), break.time.by = 730,
    palette = "jco", ggtheme = my_theme, legend = "bottom"
  )
  
  surv_p$plot <- surv_p$plot + 
    annotate("text", x = Inf, y = Inf, label = label_block_6, 
             hjust = 1.05, vjust = 1.1, size = 2.3, color = "black")
  
  panel <- surv_p$plot / surv_p$table + plot_layout(heights = c(3, 0.7))
  return(panel)
}

# ==============================================================================

panel1 <- run_spatial_combo_analysis(
  gene_name = "CHI3L1", 
  cell_type_name = "T cells CCL5+", 
  plot_legend_title = "CHI3L1 / CCL5+ T cells"
)

panel2 <- run_spatial_combo_analysis(
  gene_name = "COL4A5", 
  cell_type_name = "Fibroblast - IGFBP5+", 
  plot_legend_title = "COL4A5 / IGFBP5+ Fib"
)

panel3 <- run_spatial_combo_analysis(
  gene_name = "CHI3L1", 
  cell_type_name = "Fibroblast - EGR1+", 
  plot_legend_title = "CHI3L1 / EGR1+ Fib"
)

# ==============================================================================

spatial_composite_3panel <- (panel1 | panel2 | panel3)

ggsave(
  filename = "spatial_cell_combinations_radius10_km_3panel_composite.pdf",
  plot = spatial_composite_3panel,
  width = 12,
  height = 5,
  dpi = 1200
)


