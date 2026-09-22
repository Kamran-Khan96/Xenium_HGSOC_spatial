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
library(forestmodel)

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
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        strip.text.x = element_text(),
        strip.text.y = element_text(),
  )

#----------------------
# Clinical data
#----------------------

patient_list <- read_csv("tma5_6_mean_bulk_adj_metadata.csv") %>%
  pull(sample) %>% unique()

ova_tma_clin <- read_csv("OV_TMA_clin_KK_16-02-2026.csv") %>%
  select(sample = OVDB_ID, Pathdate, "dod_Aug 2020") %>% distinct() %>% 
  mutate(
    Pathdate = dmy(Pathdate),
    dod = dmy(`dod_Aug 2020`),
    
    # Vital status: 1 = death, 0 = censored
    vital_status = if_else(is.na(dod), 0L, 1L),
    
    # Censoring date for those alive
    censor_date = dmy("01/08/2020"),
    
    # Survival time in days
    total_living_days = case_when(
      vital_status == 1 ~ as.numeric(dod - Pathdate),
      vital_status == 0 ~ as.numeric(censor_date - Pathdate)
    )
  ) %>%
  # filter out samples with Metastasis/post-treatment
  dplyr::filter(!sample %in%
                  c("OV88",
                    "OV182",
                    "OV213")) %>%
  dplyr::filter(sample %in% patient_list)


#----------------------------------------
# load the adjusted bulk tumour counts
#-----------------------------------------

# merged_tma_adjusted - the dataframe after Combat adjustment


genes <- merged_tma_adjusted %>% select(symbol) %>% pull(symbol) %>% unique()


hr_tma5_6 <- map_dfr(genes, function(g){
  
  x <- ova_tma_clin %>%
    filter(sample %in% patient_list) %>%
    left_join(
      merged_tma_adjusted  %>%
        filter(symbol == g) %>%
        select(sample, count = adjusted_log_count, symbol, TMA),
      by = "sample"
    ) %>%
    mutate(symbol = replace_na(symbol, g)) %>%
    mutate(count = as.numeric(scale(count))) %>%
    rename(!!g := count) %>% 
    as.data.frame() 
  
  
  summary_hr <- coxph(Surv(total_living_days, vital_status) ~ count, data = x) %>%
    summary()
  
  data.frame(
    symbol = g,
    HR     = summary_hr$coefficients[1,2],
    P     = summary_hr$coefficients[1,5]
  )
  
}) %>% 
  arrange(P)


hr_genes <- hr_tma5_6 %>%
  dplyr::filter(P < 0.05 & !symbol %in% c("CFD", "CAPS")) %>%
  pull(symbol)

#---------------------------------------------
# Forest plot
#---------------------------------------------

screened_models_list <- map(genes, function(g) {
  x <- ova_tma_clin %>%
    filter(sample %in% patient_list) %>%
    left_join(
      merged_tma_adjusted %>%
        filter(symbol == g) %>%
        select(sample, count = adjusted_log_count, symbol, TMA),
      by = "sample"
    ) %>%
    mutate(symbol = replace_na(symbol, g)) %>%
    mutate(count = as.numeric(scale(count))) %>% 
    rename(!!g := count) %>% 
    as.data.frame()
  
  f <- as.formula(paste0("Surv(total_living_days, vital_status) ~ `", g, "`"))
  
  fit <- coxph(f, data = x)
  summary_hr <- summary(fit)
  
  list(
    gene = g,
    model = fit,
    p_val = summary_hr$coefficients[1, 5]
  )
})


#-----------------------------------------
# generation of forest plot for the genes with significant HR


selected_models <- map(hr_genes, function(g) {
  match <- keep(screened_models_list, ~ .x$gene == g)[[1]]
  return(match$model)
})

names(selected_models) <- hr_genes
selected_models_sorted <- selected_models[order(map_dbl(selected_models, ~ exp(coef(.x)[1])), decreasing = TRUE)]

custom_panels <- list(
  list(heading = "Gene", display = ~variable, fontface = "italic", width = 0.15),
  list(heading = "", item = "forest", line_x = 1, width = 0.40), 
  list(heading = "HR (95% CI)", display = ~sprintf("%.2f (%.2f - %.2f)", trans(estimate), trans(conf.low), trans(conf.high)), width = 0.30),
  list(heading = "P", display = ~format.pval(p.value, digits = 3, eps = 0.001), width = 0.15)
)

forest_plot <- forest_model(
  model_list = selected_models_sorted,
  merge_models = TRUE,
  panels = custom_panels,
  exponentiate = TRUE, 
  limits = c(log(0.5), log(1.8)), 
  breaks = log(c(0.5, 1.0, 1.5, 1.8)),
  theme = theme_forest() + 
    theme(
      panel.grid.major.x = element_line(color = "#EEEEEE", linetype = "dotted"),
      axis.text.y = element_blank(),
      plot.margin = margin(t = 30, r = 20, b = 30, l = 20)
    ),
  
  format_options = forest_model_format_options(
    shape = 16,
    color = "#2C3E50", 
    banded = TRUE,
    text_size = 3.5        
  )
)

forest_plot <- forest_plot + 
  geom_vline(xintercept = 0, linetype = "solid", color = "#95A5A6", linewidth = 0.5)



