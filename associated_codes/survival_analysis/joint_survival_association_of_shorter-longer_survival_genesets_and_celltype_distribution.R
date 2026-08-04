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
library(broom)



#-----------------------------------------------------------------
# load the master table with all celltype distribution abundances

frac_ct_rads <- read_csv("tma5_sqpy_coarse_tumour_touch_count_frac_diff_radii_fine_labels.csv")  %>%
  dplyr::filter(tumour_type == "Coarse") %>%
  dplyr::select(sample = PatientID, ct = "neighbour_type", count = fraction,  radius_um) %>%
  complete(sample, ct, fill = list(count = 0)) %>% 
  bind_rows(read_csv("tma6_sqpy_coarse_tumour_touch_count_frac_diff_radii_fine_labels.csv")  %>%
              dplyr::filter(tumour_type == "Coarse") %>%
              dplyr::select(sample = PatientID, ct = "neighbour_type", count = fraction, radius_um) %>%
              complete(sample, ct, fill = list(count = 0)))




#------------------------------------------------------------------------------------------------------------
# Combination KM curves for shorter-/longer-survival gene sets with celltype abundances around 10 um of tumour cells 

processed_km_df <- read_csv("en_1se_min_genes_ssGSEA_scores_shorter_survival_longer_survival_all_modules.csv") %>%
  dplyr::filter(module %in% c("shorter_survival_EN_1se", "longer_survival_EN_1se")) %>%
  pivot_wider(names_from = module, values_from = module_score) %>%
  inner_join(
    frac_ct_rads %>%
      filter(radius_um == 10) %>%
      filter(ct %in% (all_results %>%
                        filter(radius == 10 & P_tma5_6 < 0.05) %>%
                        pull(cell_type) %>%
                        unique()
      )
      ) %>%
      dplyr::select(sample, ct, count) %>%
      pivot_wider(names_from = ct, values_from = count),
    by = "sample"
  ) %>%
  mutate(
    Fib_EGR1_Status = ifelse(`Fibroblast - EGR1+` > median(`Fibroblast - EGR1+`, na.rm = TRUE), "High", "Low"),
    T_CCL5_Status   = ifelse(`T cells CCL5+` > median(`T cells CCL5+`, na.rm = TRUE), "High", "Low"),
    shorter_survival_Status    = ifelse(shorter_survival_EN_1se > median(shorter_survival_EN_1se, na.rm = TRUE), "High", "Low"),
    longer_survival_Status  = ifelse(longer_survival_EN_1se > median(longer_survival_EN_1se, na.rm = TRUE), "High", "Low")
  ) %>%
  mutate(
    shorter_survival_vs_Fib      = paste0(shorter_survival_Status, "/", Fib_EGR1_Status),
    longer_survival_vs_T   = paste0(longer_survival_Status, "/", T_CCL5_Status),
    longer_survival_vs_Fib = paste0(longer_survival_Status, "/", Fib_EGR1_Status),
    shorter_survival_vs_T        = paste0(shorter_survival_Status, "/", T_CCL5_Status)
  ) %>%
  select(sample, shorter_survival_vs_Fib, longer_survival_vs_T, longer_survival_vs_Fib, shorter_survival_vs_T)

x <- processed_km_df %>%
  pivot_longer(-sample, names_to = "combs", values_to = "fraction") %>%
  inner_join(ova_tma_clin, by = "sample") %>%
  as.data.frame()

x$fraction <- factor(x$fraction, levels = c("Low/Low", "Low/High", "High/Low", "High/High"))


combination_names = c("shorter_survival_vs_Fib", "longer_survival_vs_T", "longer_survival_vs_Fib", "shorter_survival_vs_T")

plot_list <- map(combination_names, function(comb_name) {
  
  sub_df <- x %>% filter(combs == comb_name)
  
  fit_cox <- coxph(Surv(total_living_days, vital_status) ~ fraction, data = sub_df)
  sum_cox <- summary(fit_cox)
  coefs   <- sum_cox$coefficients
  cis     <- sum_cox$conf.int
  
  hr_lh   <- round(cis["fractionLow/High", "exp(coef)"], 2)
  p_lh    <- format_raw_p(coefs["fractionLow/High", "Pr(>|z|)"])
  
  hr_hl   <- round(cis["fractionHigh/Low", "exp(coef)"], 2)
  p_hl    <- format_raw_p(coefs["fractionHigh/Low", "Pr(>|z|)"])
  
  hr_hh   <- round(cis["fractionHigh/High", "exp(coef)"], 2)
  p_hh    <- format_raw_p(coefs["fractionHigh/High", "Pr(>|z|)"])
  
  label_block <- paste0(
    "HR (Low/High vs. Low/Low): ", hr_lh, " (p", p_lh, ")\n",
    "HR (High/Low vs. Low/Low): ", hr_hl, " (p", p_hl, ")\n",
    "HR (High/High vs. Low/Low): ", hr_hh, " (p", p_hh, ")"
  )
  
  km_fit <- survfit(Surv(total_living_days, vital_status) ~ fraction, data = sub_df)
  
  p_km <- ggsurvplot(
    km_fit,
    data = sub_df,
    pval = TRUE,
    legend.labs = c("Low/Low", "Low/High", "High/Low", "High/High"),
    xlab = "Days since diagnosis",
    title = gsub("_", " ", comb_name), 
    risk.table = TRUE,
    palette = "jco",
    ggtheme = my_theme,
    xlim = c(0, 3650),
    break.time.by = 730,
    legend = "none" 
  )
  
  p_km$plot <- p_km$plot + 
    annotate(
      "text", 
      x = Inf, y = Inf, 
      label = label_block, 
      hjust = 1.05, vjust = 1.15, 
      size = 3.0, 
      color = "black"
    )
  
  return(p_km)
})



g1 <- plot_list[[1]]$plot / plot_list[[1]]$table + plot_layout(heights = c(3, 0.75))
g2 <- plot_list[[2]]$plot / plot_list[[2]]$table + plot_layout(heights = c(3, 0.75))
g3 <- plot_list[[3]]$plot / plot_list[[3]]$table + plot_layout(heights = c(3, 0.75))
g4 <- plot_list[[4]]$plot / plot_list[[4]]$table + plot_layout(heights = c(3, 0.75))

final_row_grid <- (g4 | g1) / (g2 | g3) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

#-----------------------------------------------------------------------------
# HR and P-values of the other combinations

radii_to_test <- c(10, 25, 50, 100, 250, 500)

hr_summary_df <- map_dfr(radii_to_test, function(rad) {
  
  sig_cts <- all_results %>%
    filter(radius == rad & P_tma5_6 < 0.05) %>%
    pull(cell_type) %>%
    unique()
  
  if(length(sig_cts) == 0) return(NULL)
  
  radius_counts <- frac_ct_rads %>%
    filter(radius_um == rad & ct %in% sig_cts) %>%
    dplyr::select(sample, ct, count) %>%
    pivot_wider(names_from = ct, values_from = count)
  
  merged_hr_gene_ct <- read_csv("en_1se_min_genes_ssGSEA_scores_shorter_survival_longer_survival_all_modules.csv") %>%
    dplyr::filter(module %in% c("shorter_survival_EN_1se", "longer_survival_EN_1se")) %>%
    pivot_wider(names_from = module, values_from = module_score) %>%
    inner_join(radius_counts, by = "sample") %>%
    inner_join(ova_tma_clin, by = "sample")
  
  map_dfr(sig_cts, function(current_ct) {
    
    med_ct      <- median(merged_hr_gene_ct[[current_ct]], na.rm = TRUE)
    med_shorter_survival   <- median(merged_hr_gene_ct$shorter_survival_EN_1se, na.rm = TRUE)
    med_longer_survival <- median(merged_hr_gene_ct$longer_survival_EN_1se, na.rm = TRUE)
    
    sub_df <- merged_hr_gene_ct %>%
      mutate(
        CT_Status      = ifelse(.data[[current_ct]] > med_ct, "High", "Low"),
        shorter_survival_Status    = ifelse(shorter_survival_EN_1se > med_shorter_survival, "High", "Low"),
        longer_survival_Status  = ifelse(longer_survival_EN_1se > med_longer_survival, "High", "Low")
      ) %>%
      mutate(
        shorter_survival_vs_CT     = factor(paste0(shorter_survival_Status, "/", CT_Status), 
                                 levels = c("Low/Low", "Low/High", "High/Low", "High/High")),
        longer_survival_vs_CT   = factor(paste0(longer_survival_Status, "/", CT_Status), 
                                 levels = c("Low/Low", "Low/High", "High/Low", "High/High"))
      )
    
    fit_shorter_survival <- coxph(Surv(total_living_days, vital_status) ~ shorter_survival_vs_CT, data = sub_df)
    stats_shorter_survival <- broom::tidy(fit_shorter_survival, exponentiate = TRUE, conf.int = TRUE) %>%
      mutate(Comparison_Type = "shorter_survival_vs_CellType")
    
    fit_longer_survival <- coxph(Surv(total_living_days, vital_status) ~ longer_survival_vs_CT, data = sub_df)
    stats_longer_survival <- broom::tidy(fit_longer_survival, exponentiate = TRUE, conf.int = TRUE) %>%
      mutate(Comparison_Type = "longer_survival_vs_CellType")
    
    bind_rows(stats_shorter_survival, stats_longer_survival) %>%
      mutate(
        Radius = rad,
        Cell_Type = current_ct
      ) %>%
      dplyr::select(
        Radius, 
        Cell_Type, 
        Comparison_Type, 
        Group_Comparison = term, 
        Hazard_Ratio = estimate, 
        CI_Lower = conf.low, 
        CI_Upper = conf.high, 
        p_value = p.value
      ) %>%
      mutate(Group_Comparison = gsub("shorter_survival_vs_CT|longer_survival_vs_CT", "", Group_Comparison))
  })
})


polished_hr_table <- hr_summary_df %>%
  mutate(
    Module = ifelse(Comparison_Type == "shorter_survival_vs_CellType", "shorter_survival", "longer_survival"),
    `Formula Structure` = "Module/CellType" 
  ) %>%
  
  pivot_wider(
    id_cols = c(Radius, Cell_Type, Module, `Formula Structure`),
    names_from = Group_Comparison,
    values_from = c(Hazard_Ratio, p_value)
  ) %>%
  
  dplyr::select(
    Radius,
    `Cell Type` = Cell_Type,
    Module,
    `Formula Structure`,
    
    `HR (Low/High)`   = `Hazard_Ratio_Low/High`,
    `P (Low/High)`    = `p_value_Low/High`,
    
    `HR (High/Low)`   = `Hazard_Ratio_High/Low`,
    `P (High/Low)`    = `p_value_High/Low`,
    
    `HR (High/High)`  = `Hazard_Ratio_High/High`,
    `P (High/High)`   = `p_value_High/High`
  ) %>%
  
  mutate(
    across(starts_with("HR"), ~ round(.x, digits = 2))
  ) %>%
  arrange(Radius, `Cell Type`, desc(Module))

