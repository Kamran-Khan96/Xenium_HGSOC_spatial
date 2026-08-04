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


hr_function_radii <- function(rad){
  
  # load each tma specific tumour-neighbouring spatial cell-type distribution
  
  tma5 <- read_csv("sqpy_niche/tma5_sqpy_coarse_tumour_touch_count_frac_diff_radii_fine_labels.csv")  %>%
    dplyr::filter(tumour_type == "Coarse", radius_um == rad) %>%
    select(sample = PatientID, ct = "neighbour_type", count = fraction) %>%
    complete(sample, ct, fill = list(count = 0)) %>%
    dplyr::filter(sample %in% patient_list) %>%
    mutate(TMA = "TMA5")

  tma6 <- read_csv("tma6_sqpy_coarse_tumour_touch_count_frac_diff_radii_fine_labels.csv")  %>%
    dplyr::filter(tumour_type == "Coarse", radius_um == rad) %>%
    select(sample = PatientID, ct = "neighbour_type", count = fraction) %>%
    complete(sample, ct, fill = list(count = 0)) %>%
    dplyr::filter(sample %in% patient_list) %>%
    mutate(TMA = "TMA6")

  # merge
  tma5_6 <- tma5 %>% bind_rows(tma6) %>%
    dplyr::filter(sample %in% patient_list)
  
  cts <- tma5_6 %>% select(ct) %>% pull(ct) %>% unique()
  
  
  # perform univariate Cox regression using quartile split
  
  hr_tma5_6 <- map_dfr(cts, function(g){
    
    x <- ova_tma_clin %>%
      dplyr::filter(sample %in% patient_list) %>%
      left_join(
        tma5_6 %>%
          filter(ct == g) %>%
          select(sample, count, ct),
        by = "sample"
      ) %>%
      mutate(count  = replace_na(count, 0),
             ct = replace_na(ct, g)) %>%
      mutate(tile = ntile(count, 4)) %>%
      dplyr::filter(tile %in% c(1,4)) %>%
      mutate(fraction = ifelse(tile == 4, "High", "Low")) %>%
      as.data.frame()
    
    x$fraction <- factor(x$fraction, levels = c("Low", "High"), ordered = TRUE)
    
    if (nlevels(x$fraction) > 1) {
      summary_hr <- coxph(Surv(total_living_days, vital_status) ~ fraction , data = x) %>%
        summary()
      
      data.frame(cell_type = g,
                 HR_tma5_6     = summary_hr$coefficients[1,2],
                 P_tma5_6      = summary_hr$coefficients[1,5]
      )
    } else {
      NULL
    }
    
  }) %>% arrange(P_tma5_6)
  
  write_csv(hr_tma5_6, paste0("HR_top25_v_bot25_tma5+6_sqpy_tumour+prolif_neigh_frac_", rad,"um.csv"))
  
}

radii = c(10, 25, 50, 100, 250, 500)

lapply(radii, hr_function_radii)




#----------------------------------------------------------
# plot


folder_path <- ""
base_pattern <- "HR_top25_v_bot25_tma5\\+6_sqpy_tumour\\+prolif_neigh_frac_"

files <- list.files(folder_path, pattern = paste0(base_pattern, ".*\\.csv$"), full.names = TRUE)

all_results <- map_df(files, function(f) {
  radius_val <- str_extract(basename(f), "(?<=_)\\d+(?=um)")
  
  read_csv(f) %>%
    mutate(radius = as.numeric(radius_val))
}) %>% 
  dplyr::select(cell_type, HR_tma5_6, P_tma5_6, radius)

sig_cell_types <- all_results %>%
  filter(P_tma5_6 < 0.05) %>%
  pull(cell_type) %>%
  unique()

plot_data <- all_results %>%
  filter(cell_type %in% sig_cell_types) %>%
  mutate(sig_label = case_when(
    P_tma5_6 < 0.001 ~ "***",
    P_tma5_6 < 0.01  ~ "**",
    P_tma5_6 < 0.05  ~ "*",
    TRUE          ~ ""
  ))


lineage_order <- c(
  "B_Plasma cells IGHM+", "B_Plasma cells IGHG+",
  "T cell effector", "T cells Treg", "T cells CXCL13+", "T cells CCL5+", "T cells CXCR4+", "T cells IL32+",
  "Dendritic cells", "Monocyte_Macrophages CXCL9+",
  "Mesothelial cells",
  "Fibroblast - EGR1+", "Fibroblast - IGFBP5+"
)

plot_data_final <- plot_data %>%
  mutate(cell_type = factor(cell_type, levels = lineage_order)) %>%
  filter(!is.na(cell_type))

p4 <- ggplot(plot_data_final, aes(x = factor(radius), y = HR_tma5_6, group = cell_type)) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "grey50", size = 0.6) +
  geom_line(color = "grey30", alpha = 0.5, size = 0.7,
  ) +
  geom_point(aes(fill = HR_tma5_6 > 1), shape = 21, size = 3, color = "black") +
  geom_text(aes(label = sig_label), vjust = -0.8, size = 4.5, fontface = "bold") +
  facet_wrap(~cell_type, nrow = 3, ncol = 5, scales = "free_y") +
  scale_fill_manual(values = c("TRUE" = "#B2182B", "FALSE" = "#2166AC"), guide = "none") + 
  my_theme + 
  labs(
    x = expression(Radius~(mu*m)), 
    y = "Hazard Ratio",
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold.italic", size = 9),
    panel.spacing.x = unit(1, "lines") 
  )




#-----------------------------------------------------------
# Combination KM curves in Figure 3D-E

#----------------------
# Panel 3D: T cells CXCL13+ and Monocyte_Macrophages CXCL9+

format_raw_p <- function(p) {
  if (is.na(p)) return("=NA")
  if (p < 0.001) {
    return("<0.001")
  } else {
    return(paste0("=", round(p, 3)))
  }
}

x1 <- ova_tma_clin %>%
  dplyr::filter(sample %in% patient_list) %>%
  left_join(
    frac_ct_rads %>%
      dplyr::filter(radius_um == "10") %>%
      dplyr::filter(ct %in% c("T cells CXCL13+", "Monocyte_Macrophages CXCL9+")) %>%
      dplyr::select(sample, ct, count) %>%
      group_by(ct) %>%
      mutate(fraction = ifelse(count > median(count), "High", "Low")) %>%
      ungroup() %>%
      select(-count) %>%
      tidyr::pivot_wider(names_from = ct, values_from = fraction) %>%
      mutate(fraction = interaction(.[["T cells CXCL13+"]], .[["Monocyte_Macrophages CXCL9+"]], sep = "/")),
    by = "sample"
  )

x1$fraction <- factor(x1$fraction, levels = c("Low/Low", "Low/High", "High/Low", "High/High"))

cox_1 <- summary(coxph(Surv(total_living_days, vital_status) ~ fraction, data = x1))
hr_lh_1 <- round(cox_1$conf.int["fractionLow/High", "exp(coef)"], 2)
p_lh_1  <- format_raw_p(cox_1$coefficients["fractionLow/High", "Pr(>|z|)"])
hr_hl_1 <- round(cox_1$conf.int["fractionHigh/Low", "exp(coef)"], 2)
p_hl_1  <- format_raw_p(cox_1$coefficients["fractionHigh/Low", "Pr(>|z|)"])
hr_hh_1 <- round(cox_1$conf.int["fractionHigh/High", "exp(coef)"], 2)
p_hh_1  <- format_raw_p(cox_1$coefficients["fractionHigh/High", "Pr(>|z|)"])

label_block_1 <- paste0(
  "HR (Low/High vs. Low/Low): ", hr_lh_1, " (p", p_lh_1, ")\n",
  "HR (High/Low vs. Low/Low): ", hr_hl_1, " (p", p_hl_1, ")\n",
  "HR (High/High vs. Low/Low): ", hr_hh_1, " (p", p_hh_1, ")"
)

fit1 <- survfit(Surv(total_living_days, vital_status) ~ fraction, data = x1)
names(fit1$strata) <- gsub("fraction=", "", names(fit1$strata))
fit1$strata        <- fit1$strata[c("Low/Low", "Low/High", "High/Low", "High/High")]

surv1 <- ggsurvplot(
  fit1, data = x1, risk.table = TRUE, pval = TRUE,
  legend.title = "CXCL13 / CXCL9", xlab = "Days since diagnosis",
  xlim = c(0, 3650), break.time.by = 730,
  palette = "jco", ggtheme = my_theme, legend = "bottom"
)
surv1$plot <- surv1$plot + 
  annotate("text", x = Inf, y = Inf, label = label_block_1, hjust = 1.05, vjust = 1.15, size = 3.2, color = "black")


# ==============================================================================
# Panel 3E: Il32+ T cells & IGFBP5+ Fibroblasts

x2 <- ova_tma_clin %>%
  dplyr::filter(sample %in% patient_list) %>%
  left_join(
    frac_ct_rads %>%
      dplyr::filter(radius_um == "10") %>%
      dplyr::filter(ct %in% c("T cells IL32+", "Fibroblast - IGFBP5+")) %>%
      dplyr::select(sample, ct, count) %>%
      group_by(ct) %>%
      mutate(fraction = ifelse(count > median(count), "High", "Low")) %>%
      ungroup() %>%
      select(-count) %>%
      tidyr::pivot_wider(names_from = ct, values_from = fraction) %>%
      mutate(fraction = interaction(.[["T cells IL32+"]], .[["Fibroblast - IGFBP5+"]], sep = "/")),
    by = "sample"
  )


x2$fraction <- factor(x2$fraction, levels = c("Low/Low", "Low/High", "High/Low", "High/High"))
cox_2_A     <- summary(coxph(Surv(total_living_days, vital_status) ~ fraction, data = x2))

hr_lh_2 <- round(cox_2_A$conf.int["fractionLow/High", "exp(coef)"], 2)
p_lh_2  <- format_raw_p(cox_2_A$coefficients["fractionLow/High", "Pr(>|z|)"])

hr_hl_2 <- round(cox_2_A$conf.int["fractionHigh/Low", "exp(coef)"], 2)
p_hl_2  <- format_raw_p(cox_2_A$coefficients["fractionHigh/Low", "Pr(>|z|)"])


x2_alt   <- x2 %>% mutate(fraction_alt = factor(fraction, levels = c("Low/High", "High/Low", "Low/Low", "High/High")))
cox_2_B  <- summary(coxph(Surv(total_living_days, vital_status) ~ fraction_alt, data = x2_alt))

hr_hl_vs_lh_2 <- round(cox_2_B$conf.int["fraction_altHigh/Low", "exp(coef)"], 2)
p_hl_vs_lh_2  <- format_raw_p(cox_2_B$coefficients["fraction_altHigh/Low", "Pr(>|z|)"])


label_block_2 <- paste0(
  "HR (Low/High vs. Low/Low): ", hr_lh_2, " (p", p_lh_2, ")\n",
  "HR (High/Low vs. Low/Low): ", hr_hl_2, " (p", p_hl_2, ")\n",
  "HR (High/Low vs. Low/High): ", hr_hl_vs_lh_2, " (p", p_hl_vs_lh_2, ")"
)

fit2 <- survfit(Surv(total_living_days, vital_status) ~ fraction, data = x2)
names(fit2$strata) <- gsub("fraction=", "", names(fit2$strata))
fit2$strata        <- fit2$strata[c("Low/Low", "Low/High", "High/Low", "High/High")]

surv2 <- ggsurvplot(
  fit2, data = x2, risk.table = TRUE, pval = TRUE,
  legend.title = "IL32 / IGFBP5", xlab = "Days since diagnosis",
  legend.labs = c("Low/Low", "Low/High", "High/Low", "High/High"),
  xlim = c(0, 3650), break.time.by = 730, 
  palette = "jco", ggtheme = my_theme, legend = "bottom"
)
surv2$plot <- surv2$plot + 
  annotate("text", x = Inf, y = Inf, label = label_block_2, hjust = 1.05, vjust = 1.15, size = 3.2, color = "black")


# ==============================================================================
# generate final plot


panel1 <- surv1$plot / surv1$table + plot_layout(heights = c(3, 0.7))
panel2 <- surv2$plot / surv2$table + plot_layout(heights = c(3, 0.7))

spatial_composite <- (panel1 | panel2)



