library(tidyverse)
library(glmnet)
library(survival)

patient_list <- read_csv("tma5_6_mean_bulk_adj_metadata.csv") %>%
  pull(sample) %>% unique()

ova_tma_clin <- read_csv("OV_TMA_clin_KK_16-02-2026.csv") %>%
  select(sample = OVDB_ID, Pathdate, "dod_Aug 2020") %>% distinct() %>% 
  mutate(
    # Convert dates safely
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


#---------------------------------------
# load and keep counts for only the survival-associated genes

adj_tum_bulk <- merged_tma_adjusted %>%
  dplyr::filter(symbol %in% c("CTPS2", "IGF1R", "MECOM", "MSH2", "MSH6", "NOTCH3", "OAS1", "RNF43", "TFB1M", "ADA", "DLK1", "RBPJ", "SOGA1", 
                              "ADIRF", "APBA2", "CACNA1C", "CCND1", "CD44", "CD9", "CD99", "CDH3", "CXCL12", "ENG", "F2RL1", "FAM120B", "FGFR3", 
                              "ITGB5", "L1CAM", "MLH1", "MTSS1", "PDGFRB", "PPARG", "RORA", "TNFRSF12A", "VCAN", "KLRG1", "MGAM", "TMEM45A", "PLIN1", 
                              "TFF3", "CHI3L1", "COL4A5", "COMP", "EPOR", "ISG20", "PROX1", "RASGRP4", "TNFRSF19", "ADIPOQ", "CD207", "FGF18", "FOXP3", 
                              "ICAM2", "LAG3", "LTF", "PTN", "TCF15", "TERT", "TNFRSF11B", "VWF", "MPO")) ## using the 61 HR associated genes



#---------------------------------
# prepare the data
X_combat <- adj_tum_bulk %>% 
  dplyr::select(-TMA) %>%
  pivot_wider(names_from = sample, values_from = adjusted_log_count, values_fn = mean) %>% 
  column_to_rownames(var = "symbol") %>% 
  as.matrix()


clin <- ova_tma_clin %>%
  mutate(
    status = vital_status,
    time = total_living_days
  )

common_samples <- intersect(colnames(X_combat), clin$sample)

X <- X_combat[, common_samples]
clin2 <- clin %>% filter(sample %in% common_samples)


clin2 <- clin2[match(common_samples, clin2$sample), ]
clin2 <- clin2 %>%
  mutate(time = ifelse(total_living_days <= 0, 0.00001, total_living_days))

y <- Surv(time = clin2$time, event = clin2$status)


x <- t(X) 

x_scaled <- scale(x)


#--------------------------------------------
# apply elastic-net

fit_en <- cv.glmnet(
  x_scaled, y,
  family = "cox",
  alpha = 0.5
)


#--------------------------------------------
# extract the coefficients

extract_coefs <- function(fit, alpha_value, lambda_type) {
  
  cm <- as.matrix(coef(fit, s = lambda_type))
  
  df <- cm %>%
    as.data.frame() %>%
    rownames_to_column("Gene")
  
  if (ncol(df) == 1) {
    names(df)[2] <- "Coefficient"
  } else {
    names(df)[2] <- "Coefficient"
  }
  
  df %>%
    mutate(
      Alpha = alpha_value,
      Lambda = lambda_type,
      Selected = Coefficient != 0
    )
}



df_all <- bind_rows(
  extract_coefs(fit_en, alpha_value = 0.5, lambda_type = "lambda.min"),
  extract_coefs(fit_en, alpha_value = 0.5, lambda_type = "lambda.1se")
)


df_all <- df_all %>%
  mutate(
    Model = case_when(
      Alpha == 0.5 ~ "Elastic Net"
    ),
    Lambda = ifelse(Lambda == "lambda.min", "min", "1se")
  ) %>%
  select(Gene, Coefficient, Selected, Model, Alpha, Lambda)

