library(dplyr)
library(tidyr)
library(survminer)
library(survival)
library(patchwork)

## coefficients after the elastic-net prioritisation from the 61 survival-associated genes

coef_df <- df_all %>%
  select(symbol = Gene, Coefficient, Selected, Model, Alpha, Lambda) %>%
  filter(Lambda == "1se")

# -------------------------------------------------------------------------
# calculation of the unified scoring


adj_tum_bulk <- merged_tma_adjusted %>%
  dplyr::filter(symbol %in% c("CTPS2", "IGF1R", "MECOM", "MSH2", "MSH6", "NOTCH3", "OAS1", "RNF43", "TFB1M", "ADA", "DLK1", "RBPJ", "SOGA1", 
                              "ADIRF", "APBA2", "CACNA1C", "CCND1", "CD44", "CD9", "CD99", "CDH3", "CXCL12", "ENG", "F2RL1", "FAM120B", "FGFR3", 
                              "ITGB5", "L1CAM", "MLH1", "MTSS1", "PDGFRB", "PPARG", "RORA", "TNFRSF12A", "VCAN", "KLRG1", "MGAM", "TMEM45A", "PLIN1", 
                              "TFF3", "CHI3L1", "COL4A5", "COMP", "EPOR", "ISG20", "PROX1", "RASGRP4", "TNFRSF19", "ADIPOQ", "CD207", "FGF18", "FOXP3", 
                              "ICAM2", "LAG3", "LTF", "PTN", "TCF15", "TERT", "TNFRSF11B", "VWF", "MPO"))

patient_unified_score_df <- adj_tum_bulk %>%
  filter(symbol %in% coef_df$symbol) %>%
  
  group_by(symbol) %>%
  mutate(z_scaled_counts = as.numeric(scale(adjusted_log_count))) %>%
  ungroup() %>%
  
  left_join(coef_df, by = "symbol") %>%
  
  mutate(weighted_score = z_scaled_counts * Coefficient) %>%

  group_by(sample) %>%
  summarise(
    unified_score_Score = sum(weighted_score, na.rm = TRUE),
    .groups = "drop"
  )

# -------------------------------------------------------------------------
# KM plot with the unified scores

km_data <- patient_unified_score_df %>%
  mutate(
    unified_score_Group = ifelse(unified_score_Score > median(unified_score_Score, na.rm = TRUE), "High unified_score", "Low unified_score"),
    unified_score_Group = factor(unified_score_Group, levels = c("Low unified_score", "High unified_score"))
  ) %>%
  inner_join(ova_tma_clin, by = "sample")

# Fit KM and Cox models
km_fit   <- survfit(Surv(total_living_days, vital_status) ~ unified_score_Group, data = km_data)
cox_unified_score <- summary(coxph(Surv(total_living_days, vital_status) ~ unified_score_Group, data = km_data))


hr_unified_score    <- round(cox_unified_score$conf.int["unified_score_GroupHigh unified_score", "exp(coef)"], 2)
lower_unified_score <- round(cox_unified_score$conf.int["unified_score_GroupHigh unified_score", "lower .95"], 2)
upper_unified_score <- round(cox_unified_score$conf.int["unified_score_GroupHigh unified_score", "upper .95"], 2)
p_unified_score     <- format_raw_p(cox_unified_score$coefficients["unified_score_GroupHigh unified_score", "Pr(>|z|)"])


label_unified_unified_score <- paste0("HR: ", hr_unified_score, " (95% CI: ", lower_unified_score, "-", upper_unified_score, ", p ", p_unified_score, ")")

names(km_fit$strata) <- gsub("unified_score_Group=", "", names(km_fit$strata))

# -------------------------------------------------------------------------
# PLotting

p_unified <- ggsurvplot(
  km_fit, 
  data = km_data, 
  unified_score.table = TRUE, 
  pval = TRUE, 
  palette = "aaas",                 
  legend.labs = c("Low unified_score", "High unified_score"), 
  xlab = "Days since diagnosis",
  ggtheme = my_theme,
  xlim = c(0, 3650), 
  break.time.by = 730, 
  legend = "bottom"
)

p_unified$plot <- p_unified$plot + 
  annotate(
    "text", x = Inf, y = Inf, 
    label = label_unified_unified_score, 
    hjust = 1.05, vjust = 1.5, 
    size = 3, color = "black"
  )

panel_unified <- p_unified$plot / p_unified$table + plot_layout(heights = c(3, 0.5))
panel_unified



#============================
# timeROC assessment with the unified score


library(timeROC)

# -------------------------------------------------------------------------
# Compute Time-Dependent ROC for 1, 2, and 5 Years

roc_data <- km_data %>%
  filter(!is.na(total_living_days), !is.na(vital_status), !is.na(unified_score_Score))

# Target time horizons in days
target_times <- c(365, 730, 1825)

time_roc_fit <- timeROC(
  T = roc_data$total_living_days,
  delta = roc_data$vital_status,
  marker = roc_data$unified_score_Score,
  cause = 1,
  times = target_times, 
  iid = TRUE            
)

ci_matrix <- confint(time_roc_fit)$CI_AUC

get_auc_label <- function(index, label_text) {
  auc  <- round(time_roc_fit$AUC[index], 2)
  low  <- round(ci_matrix[index, 1]/100, 2)
  high <- round(ci_matrix[index, 2]/100, 2)
  return(paste0(label_text, " (AUC = ", auc, " [95% CI: ", low, "-", high, "])"))
}

label_1yr  <- get_auc_label(1, "1 Year")
label_2yr  <- get_auc_label(2, "2 Years")
label_5yr  <- get_auc_label(3, "5 Years")

# -------------------------------------------------------------------------
# plotting

df_roc <- rbind(
  data.frame(False_Positive = time_roc_fit$FP[, 1], True_Positive = time_roc_fit$TP[, 1], Timepoint = label_1yr),
  data.frame(False_Positive = time_roc_fit$FP[, 2], True_Positive = time_roc_fit$TP[, 2], Timepoint = label_2yr),
  data.frame(False_Positive = time_roc_fit$FP[, 3], True_Positive = time_roc_fit$TP[, 3], Timepoint = label_5yr)
)

df_roc$Timepoint <- factor(df_roc$Timepoint, levels = c(label_1yr, label_2yr, label_5yr))

p_roc_soothing <- ggplot(df_roc, aes(x = False_Positive, y = True_Positive, color = Timepoint)) +
  geom_line(size = 1.1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey70") +
  
  scale_color_npg() + 
  
  labs(
    x = "1 - Specificity (False Positive Rate)",
    y = "Sensitivity (True Positive Rate)",
  ) +
  scale_x_continuous(limits = c(0, 1), expand = c(0.01, 0.01)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0.01, 0.01)) +
  my_theme +
  theme(
    legend.position = c(0.62, 0.22),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.text = element_text(size = 8.5)
  )

final_layout <- panel_unified | p_roc_soothing
