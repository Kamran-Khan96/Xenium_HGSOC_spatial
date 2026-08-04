library(tidyverse)
library(survminer)
library(survival)
library(GSVA)
library(patchwork)

## the df generated from the Elastic-net prioritisation will be used here

en_filt <- df_all %>%
  unite(Model_Lambda, c(Model, Lambda), sep = "-") %>%
  mutate(label_text = ifelse(abs(Coefficient) > 0, round(Coefficient, digits = 4), "")) %>%
  group_by(Gene) %>%
  filter(any(abs(Coefficient) > 0)) %>%
  ungroup()

#----------------------------------------------------------------
## Construct the shorter-survival and longer-survival gene sets

gene_sets <- list(
  
  shorter_survival_EN_1se = en_filt %>% filter(Coefficient > 0 & Model_Lambda == "Elastic Net-1se") %>% pull(Gene) %>% unique(),
  longer_survival_EN_1se = en_filt %>% filter(Coefficient < 0 & Model_Lambda == "Elastic Net-1se") %>% pull(Gene) %>% unique()
  
  )


#--------------------------------------
# prepare the adjust the count data

coarse_tum_wide <- read_csv("all_tumour_merged_tma5_6_combat_adjusted_log2_counts.csv") %>%
  dplyr::select(sample, symbol, count = adjusted_log_count) %>%
  pivot_wider(names_from = sample, values_from = count) %>%
  drop_na(symbol) %>% 
  column_to_rownames("symbol")

coarse_tum_matrix <- as.matrix(coarse_tum_wide)


#----------------------------------
# Run ssGSEA

params <- ssgseaParam(
  exprData = coarse_tum_matrix,
  geneSets = gene_sets,
  normalize = FALSE
)

ssgsea_scores <- gsva(params)

ssgsea_scores_df_en <- as.data.frame(ssgsea_scores) %>%
  rownames_to_column(var = "module") %>%
  as_tibble() %>%
  pivot_longer(-module, names_to = "sample", values_to = "module_score")


#---------------------------------------------------
# KM curves

processed_base <- ssgsea_scores_df_en %>%
  dplyr::filter(module %in% c("shorter_survival_EN_1se", "longer_survival_EN_1se")) %>%
  pivot_wider(names_from = module, values_from = module_score) %>%
  inner_join(ova_tma_clin, by = "sample") %>%
  mutate(
    shorter_survival_Status = ifelse(shorter_survival_EN_1se > median(shorter_survival_EN_1se, na.rm = TRUE), "High", "Low"),
    longer_survival_Status = ifelse(longer_survival_EN_1se > median(longer_survival_EN_1se, na.rm = TRUE), "High", "Low")
  ) %>%
  mutate(
    Combined_Modules = paste0(shorter_survival_Status, "/", longer_survival_Status)
  )

processed_base$shorter_survival_Status   <- factor(processed_base$shorter_survival_Status, levels = c("Low", "High"))
processed_base$longer_survival_Status <- factor(processed_base$longer_survival_Status, levels = c("Low", "High"))
processed_base$Combined_Modules <- factor(
  processed_base$Combined_Modules, 
  levels = c("Low/Low", "Low/High", "High/Low", "High/High")
)

format_raw_p <- function(p) {
  if (is.na(p)) return("=NA")
  if (p < 0.001) {
    return("<0.001")
  } else {
    return(paste0("=", round(p, 3)))
  }
}

fit_shorter_survival <- survfit(Surv(total_living_days, vital_status) ~ shorter_survival_Status, data = processed_base)
cox_shorter_survival <- summary(coxph(Surv(total_living_days, vital_status) ~ shorter_survival_Status, data = processed_base))
hr_r      <- round(cox_shorter_survival$conf.int["shorter_survival_StatusHigh", "exp(coef)"], 2)
p_r       <- format_raw_p(cox_shorter_survival$coefficients["shorter_survival_StatusHigh", "Pr(>|z|)"])
label_shorter_survival <- paste0("HR (High vs Low): ", hr_r, " (p", p_r, ")")

fit_longer_survival <- survfit(Surv(total_living_days, vital_status) ~ longer_survival_Status, data = processed_base)
cox_longer_survival <- summary(coxph(Surv(total_living_days, vital_status) ~ longer_survival_Status, data = processed_base))
hr_p        <- round(cox_longer_survival$conf.int["longer_survival_StatusHigh", "exp(coef)"], 2)
p_p         <- format_raw_p(cox_longer_survival$coefficients["longer_survival_StatusHigh", "Pr(>|z|)"])
label_longer_survival <- paste0("HR (High vs Low): ", hr_p, " (p", p_p, ")")

fit_combined <- survfit(Surv(total_living_days, vital_status) ~ Combined_Modules, data = processed_base)
cox_comb     <- summary(coxph(Surv(total_living_days, vital_status) ~ Combined_Modules, data = processed_base))

hr_lh <- round(cox_comb$conf.int["Combined_ModulesLow/High", "exp(coef)"], 2)
p_lh  <- format_raw_p(cox_comb$coefficients["Combined_ModulesLow/High", "Pr(>|z|)"])

hr_hl <- round(cox_comb$conf.int["Combined_ModulesHigh/Low", "exp(coef)"], 2)
p_hl  <- format_raw_p(cox_comb$coefficients["Combined_ModulesHigh/Low", "Pr(>|z|)"])

hr_hh <- round(cox_comb$conf.int["Combined_ModulesHigh/High", "exp(coef)"], 2)
p_hh  <- format_raw_p(cox_comb$coefficients["Combined_ModulesHigh/High", "Pr(>|z|)"])

processed_base_alt <- processed_base %>%
  mutate(Combined_Modules_LH_Base = factor(Combined_Modules, levels = c("Low/High", "High/Low", "Low/Low", "High/High")))
cox_comb_B   <- summary(coxph(Surv(total_living_days, vital_status) ~ Combined_Modules_LH_Base, data = processed_base_alt))

hr_hl_vs_lh <- round(cox_comb_B$conf.int["Combined_Modules_LH_BaseHigh/Low", "exp(coef)"], 2)
p_hl_vs_lh  <- format_raw_p(cox_comb_B$coefficients["Combined_Modules_LH_BaseHigh/Low", "Pr(>|z|)"])


label_combined <- paste0(
  "HR (Low/High vs Low/Low): ", hr_lh, " (p", p_lh, ")\n",
  "HR (High/Low vs Low/Low): ", hr_hl, " (p", p_hl, ")\n",
  "HR (High/Low vs Low/High): ", hr_hl_vs_lh, " (p", p_hl_vs_lh, ")"
)

names(fit_shorter_survival$strata)    <- gsub("shorter_survival_Status=", "", names(fit_shorter_survival$strata))
names(fit_longer_survival$strata)  <- gsub("longer_survival_Status=", "", names(fit_longer_survival$strata))
names(fit_combined$strata) <- gsub("Combined_Modules=", "", names(fit_combined$strata))
fit_combined$strata        <- fit_combined$strata[c("Low/Low", "Low/High", "High/Low", "High/High")]


p_shorter_survival <- ggsurvplot(
  fit_shorter_survival, data = processed_base, risk.table = TRUE, pval = TRUE, palette = "aaas",
  legend.labs = c("Low", "High"), xlab = "Days since diagnosis",
  title = "shorter_survival Module (Independent)", ggtheme = my_theme,
  xlim = c(0, 3650), break.time.by = 730, legend = "bottom"
)
p_shorter_survival$plot <- p_shorter_survival$plot + 
  annotate("text", x = Inf, y = Inf, label = label_shorter_survival, hjust = 1.05, vjust = 1.15, size = 3.0, color = "black")

p_longer_survival <- ggsurvplot(
  fit_longer_survival, data = processed_base, risk.table = TRUE, pval = TRUE, palette = "aaas",
  legend.labs = c("Low", "High"), xlab = "Days since diagnosis",
  title = "longer_survivalive Module (Independent)", ggtheme = my_theme,
  xlim = c(0, 3650), break.time.by = 730, legend = "bottom"
)
p_longer_survival$plot <- p_longer_survival$plot + 
  annotate("text", x = Inf, y = Inf, label = label_longer_survival, hjust = 1.05, vjust = 1.15, size = 3.0, color = "black")

p_combined <- ggsurvplot(
  fit_combined, data = processed_base, risk.table = TRUE, pval = TRUE, palette = "jco",
  legend.labs = c("Low/Low", "Low/High", "High/Low", "High/High"),
  xlab = "Days since diagnosis", title = "Combined Module Interaction", ggtheme = my_theme,
  xlim = c(0, 3650), break.time.by = 730, legend = "bottom"
)
p_combined$plot <- p_combined$plot + 
  annotate("text", x = Inf, y = Inf, label = label_combined, hjust = 1.05, vjust = 1.15, size = 3.0, color = "black")

panel_shorter_survival    <- p_shorter_survival$plot / p_shorter_survival$table + plot_layout(heights = c(3, 0.8))
panel_longer_survival  <- p_longer_survival$plot / p_longer_survival$table + plot_layout(heights = c(3, 0.8))
panel_combined <- p_combined$plot / p_combined$table + plot_layout(heights = c(3, 0.8))

final_horizontal_row <- (
  panel_shorter_survival | 
    panel_longer_survival | 
    panel_combined
)
