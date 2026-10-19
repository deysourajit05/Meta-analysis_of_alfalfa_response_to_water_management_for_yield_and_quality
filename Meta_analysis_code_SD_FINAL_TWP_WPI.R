#Code used to clear the console
rm(list = ls())
cat("\014")
if (!is.null(dev.list())) dev.off()


library(mice) #For imputation
library(revtools)
library(dplyr)
library(metafor)
library(purrr)
library(ggplot2)
library(tidyr) # For pivot_longer
library(patchwork) # For combining plots
library(forcats)   # For ordering factors in plots
library(cowplot)
library(meta)
library(maps)
library(geodata)
library(terra)
library(sf)
library(ggspatial)
library(ggpmisc)

#This dataset is prepared using Zotero by just removing the duplicated files as two datasources was used
bib <- revtools::read_bibliography("Meta_analysis_cleaned_dataset.ris")

# Screening the abstracts based on abstracts and title using Revtools
revtools::screen_abstracts(bib)

# Loading raw data
meta <- read.csv("Metaanalysis_dataset_SD_Test.csv", header=T)
meta
# 
# #Run multiple imputation
# ###############################################################

# the function for imputation
impute_sd <- function(data, mean_var, sd_var, new_name, m = 1000, seed = 123) {
  
  df_tmp <- data %>%
    dplyr::select(all_of(c(mean_var, sd_var)))
  
  imp <- mice(
    df_tmp,
    m = m,
    method = c("", "pmm"),  # do not impute mean, impute SD only
    seed = seed,
    printFlag = FALSE
  )
  
  sd_imp <- mice::complete(imp, action = 1)[[sd_var]]
  
  data %>%
    mutate(!!new_name := sd_imp)
}

#Define what needs to be imputed
impute_specs <- tibble::tibble(
  mean_var = c(
    "Yield_kg_trt",
    "Yield_kg_ctrl",
    "CP_trt",
    "CP_ctrl",
    "WPI_kg.m3_trt",
    "TWP_kg.m3_trt",
    "WPI_kg.m3_ctrl",
    "TWP_kg.m3_ctrl"),
  
  sd_var = c(
    "Yield_kg_SD_trt",
    "Yield_kg_SD_ctrl",
    "CP_SD_trt",
    "CP_SD_ctrl",
    "WPI_kg.m3_SD_trt",
    "TWP_kg.m3_SD_trt",
    "WPI_kg.m3_SD_ctrl",
    "TWP_kg.m3_SD_ctrl"),
  
  new_name = c(
    "Yield_SD_trt_imp",
    "Yield_SD_ctrl_imp",
    "CP_SD_trt_imp",
    "CP_SD_ctrl_imp",
    "WPI_kg.m3_SD_trt_imp",
    "TWP_kg.m3_SD_trt_imp",
    "WPI_kg.m3_SD_ctrl_imp",
    "TWP_kg.m3_SD_ctrl_imp"
  )
)

#Looping this up
meta_imp <- meta

for (i in seq_len(nrow(impute_specs))) {
  meta_imp <- impute_sd(
    data    = meta_imp,
    mean_var = impute_specs$mean_var[i],
    sd_var   = impute_specs$sd_var[i],
    new_name = impute_specs$new_name[i],
    m        = 1000,
    seed     = 123
  )
}

#Data downloading with four new imputed columns
# write.csv(meta_imp, "Metanalysis_dataset_imputed.csv", row.names = FALSE)



# chk <- meta_imp %>%
#   transmute(Study_ID, Location_ID,
#             r_mean_trt  = IWUE_kg.ha.mm_trt     / WPI_kg.m3_trt,
#             r_mean_ctrl = IWUE_kg.ha.mm_ctrl    / WPI_kg.m3_ctrl,
#             r_sd_trt    = IWUE_kg.ha.mm_SD_trt  / WPI_kg.m3_SD_trt,
#             r_sd_ctrl   = IWUE_kg.ha.mm_SD_ctrl / WPI_kg.m3_SD_ctrl)
# 
# # Rows where any ratio is not 10 (allowing rounding)
# chk %>% filter(if_any(starts_with("r_"), ~ !is.na(.x) & abs(.x - 10) > 0.05))
# 
# # Rows that have a value in one version but not the other
# meta_imp %>% filter(is.na(IWUE_kg.ha.mm_trt) != is.na(WPI_kg.m3_trt) |
#                       is.na(IWUE_kg.ha.mm_SD_trt) != is.na(WPI_kg.m3_SD_trt))
# 
# # WPI vs IWUE
# chk %>% filter(if_any(starts_with("r_"), ~ !is.na(.x) & abs(.x - 10) > 0.05))
# 
# # TWP vs TWP
# meta_imp %>%
#   transmute(Study_ID, Location_ID,
#             r_mean_trt  = TWP_kg.ha.mm_trt     / TWP_kg.m3_trt,
#             r_mean_ctrl = TWP_kg.ha.mm_ctrl    / TWP_kg.m3_ctrl,
#             r_sd_trt    = TWP_kg.ha.mm_SD_trt  / TWP_kg.m3_SD_trt,
#             r_sd_ctrl   = TWP_kg.ha.mm_SD_ctrl / TWP_kg.m3_SD_ctrl) %>%
#   filter(if_any(starts_with("r_"), ~ !is.na(.x) & abs(.x - 10) > 0.05))
# 
# meta_imp %>%
#   mutate(row_id = row_number(),
#          ratio  = TWP_kg.ha.mm_trt / TWP_kg.m3_trt) %>%
#   filter(Study_ID == 14, abs(ratio - 10) > 0.05) %>%
#   select(row_id, Study_ID, Location_ID, Year, Trt_type,
#          TWP_kg.ha.mm_trt, TWP_kg.m3_trt, ratio)

#################################################################################

meta_es <- meta_imp %>%
  dplyr::mutate(
    #Log-response and variance for Yield
    Yield_LnR = log(Yield_kg_trt / Yield_kg_ctrl),
    
    Yield_Variance = (Yield_SD_ctrl_imp^2 / (Rep_ctrl * Yield_kg_ctrl^2)) +
      (Yield_SD_trt_imp^2  / (Rep_trt  * Yield_kg_trt^2)),
    
    #Log-response and variance for CP
    CP_LnR = log(CP_trt / CP_ctrl),
    
    CP_Variance = (CP_SD_ctrl_imp^2 / (Rep_ctrl * CP_ctrl^2)) +
      (CP_SD_trt_imp^2  / (Rep_trt  * CP_trt^2)),
    
    #Log-response and variance for WPI
    WPI_LnR = log(WPI_kg.m3_trt / WPI_kg.m3_ctrl),
    
    WPI_Variance = (WPI_kg.m3_SD_ctrl_imp^2 / (Rep_ctrl * WPI_kg.m3_ctrl^2)) +
      (WPI_kg.m3_SD_trt_imp^2  / (Rep_trt  * WPI_kg.m3_trt^2)),
    
    #Log-response and variance for TWP
    TWP_LnR = log(TWP_kg.m3_trt / TWP_kg.m3_ctrl),
    
    TWP_Variance = (TWP_kg.m3_SD_ctrl_imp^2 / (Rep_ctrl * TWP_kg.m3_ctrl^2)) +
      (TWP_kg.m3_SD_trt_imp^2  / (Rep_trt  * TWP_kg.m3_trt^2)),
    
    #weights
    Yield_weight = 1/Yield_Variance,
    CP_weight = 1/CP_Variance,
    WPI_weight = 1/WPI_Variance,
    TWP_weight = 1/TWP_Variance
  )

#####################################################################################

#Making the categories for Moderators with continuous values
# For Degree of Deficit irrigation
meta_es <- meta_es %>%
  mutate(
    Degree_of_irri_deficit_grp = case_when(
      Degree_of_irri_deficit < 0.35 ~ "Severe (<0.35)",
      Degree_of_irri_deficit >= 0.35 & Degree_of_irri_deficit <= 0.70 ~ "Moderate (0.35-0.70)",
      Degree_of_irri_deficit > 0.70 ~ "Slight (>0.70)",
      TRUE ~ NA_character_
    )
  )

# For Climate regions
meta_es <- meta_es %>%
  mutate(
    MAP_mm_grp = case_when(
      MAP_mm < 250 ~ "<250 mm",
      MAP_mm >= 250 & MAP_mm <= 500 ~ "250–500 mm",
      MAP_mm > 500 ~ ">500 mm",
      TRUE ~ NA_character_
    )
  )

meta_es$MAP_mm_grp <- factor(
  meta_es$MAP_mm_grp,
  levels = c("<250 mm", "250–500 mm", ">500 mm")
)

# For mean annual temperature
meta_es <- meta_es %>%
  mutate(
    MAT_degreeC_grp = case_when(
      MAT_degreeC <= 10 ~ "≤10 °C",
      MAT_degreeC > 10 ~ ">10 °C",
      TRUE ~ NA_character_
    )
  )

meta_es$MAT_degreeC_grp <- factor(
  meta_es$MAT_degreeC_grp,
  levels = c("≤10 °C", ">10 °C")
)


# For Soil Bulk Density
meta_es <- meta_es %>%
  mutate(
    Soil_BD_grp = case_when(
      Soil_BD_g.cm.3 < 1.4 ~ "<1.4 g cm^-3",
      Soil_BD_g.cm.3 >= 1.4 & Soil_BD_g.cm.3 <= 1.5 ~ "1.4–1.5 g cm^-3",
      Soil_BD_g.cm.3 > 1.5 ~ ">1.5 g cm^-3",
      TRUE ~ NA_character_
    )
  )

meta_es$Soil_BD_grp <- factor(
  meta_es$Soil_BD_grp,
  levels = c("<1.4 g cm^-3", "1.4–1.5 g cm^-3", ">1.5 g cm^-3")
)

# For Soil pH
meta_es <- meta_es %>%
  mutate(
    Soil_pH_grp = case_when(
      Soil_pH <= 7 ~ "≤7",
      Soil_pH > 7 ~ ">7",
      TRUE ~ NA_character_
    )
  )

meta_es$Soil_pH_grp <- factor(
  meta_es$Soil_pH_grp,
  levels = c("≤7", ">7")
)


# For Soil Organic Matter
meta_es <- meta_es %>%
  mutate(
    Soil_OM_grp = case_when(
      Soil_OM_. <= 1.0 ~ "≤1 %",
      Soil_OM_. > 1.0 ~ ">1 %",
      TRUE ~ NA_character_
    )
  )

meta_es$Soil_OM_grp <- factor(
  meta_es$Soil_OM_grp,
  levels = c("≤1 %", ">1 %")
)

# For Winterhardiness
meta_es <- meta_es %>%
  mutate(
    Winterhardiness_grp = case_when(
      Winterhardiness %in% c("4a", "4b") ~ "4a–4b",
      Winterhardiness %in% c("5a", "5b", "6a", "6b", "7a", "7b") ~ "5a–7b",
      Winterhardiness %in% c("8a", "8b", "9a", "9b") ~ "8a–9b",
      TRUE ~ NA_character_
    )
  )

meta_es$Winterhardiness_grp <- factor(
  meta_es$Winterhardiness_grp,
  levels = c("4a–4b", "5a–7b", "8a–9b")
)

# For lat and long
meta_es <- meta_es %>%
  mutate(
    Latitude_grp = case_when(
      Latitude < 35 ~ "Southern region",
      Latitude >= 35 & Latitude <= 45 ~ "Mid region",
      Latitude > 45 ~ "Northern region",
      TRUE ~ NA_character_
    ),
    # Create Longitude groups
    Longitude_grp = case_when(
      Longitude > -100 ~ "Humid East",
      Longitude < -100 ~ "Arid West",
      TRUE ~ NA_character_
    )
  )

# For No. of cuttings
meta_es <- meta_es %>%
  mutate(
    Cut_number_grp = case_when(
      No._of_cuttings <= 4 ~ "≤4",
      No._of_cuttings >= 5 & No._of_cuttings <= 6 ~ "5-6",
      No._of_cuttings >= 7 ~ "7",
      TRUE ~ NA_character_
    )
  )

meta_es$Cut_number_grp <- factor(
  meta_es$Cut_number_grp,
  levels = c("≤4", "5-6", "7")
)


##################################################################################################

#Now this is to get the Qm and P value directly for yield, CP, WPI, and TWP all at a time for all the moderators in a loop

outcomes <- c("Yield", "CP", "WPI", "TWP")
moderators <- c("Irrigation_method", "Irri_amt", "Degree_of_irri_deficit_grp", 
                "Irrigation_scheduling", "MAP_mm_grp", "MAT_degreeC_grp", 
                "Texture", "Soil_BD_grp", "Soil_pH_grp", "Soil_OM_grp", 
                "Winterhardiness_grp", "Latitude_grp", "Longitude_grp", "Cut_number_grp")

model_list <- list()
qm_summary <- data.frame()

for (out in outcomes) {
  for (mod in moderators) {
    
    yi_col <- paste0(out, "_LnR")
    V_col  <- paste0(out, "_Variance")
    W_col  <- paste0(out, "_weight")
    
    # Run Model
    res <- try(rma.mv(
      yi = meta_es[[yi_col]],
      V  = meta_es[[V_col]],
      W  = meta_es[[W_col]],
      mods = as.formula(paste("~", mod)),
      random = ~ 1 | Study_ID / Location_ID,
      data = meta_es
    ), silent = TRUE)
    
    if (!inherits(res, "try-error")) {
      model_list[[out]][[mod]] <- res
      
      # Extract QM, df, and p-value
      temp_res <- data.frame(
        Outcome = out,
        Moderator = mod,
        QM_stat = round(res$QM, 2),
        df = res$p - 1, # Number of coefficients minus intercept
        p_val = round(res$QMp, 3)
      )
      qm_summary <- rbind(qm_summary, temp_res)
    }
  }
}

# View results
print(qm_summary)


##################################################################################


meta_es$Study_ID <- as.factor(meta_es$Study_ID)
meta_es$Location_ID <- as.factor(meta_es$Location_ID)
meta_es$Irrigation_method <- as.factor(meta_es$Irrigation_method)
meta_es$Texture <- as.factor(meta_es$Texture)
meta_es$No._of_cuttings <- as.factor(meta_es$No._of_cuttings)


#####################################################################################
model_nested <- rma.mv(
  yi = Yield_LnR,
  V = Yield_Variance,
  random = ~ 1 | Study_ID / Location_ID, # Nesting: Location within Study
  data = meta_es
)

summary(model_nested)

#########################################################################################
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(metafor)
library(cowplot)

# -------------------------------------------------------------------------
# 1. Helper Function to Calculate Total I^2 for rma.mv models
# -------------------------------------------------------------------------
# Standard rma.mv doesn't give a single I^2. This calculates the total 
# proportion of variance attributed to heterogeneity (Study + Location).
calc_i2 <- function(model) {
  n <- model$k
  p <- model$p
  # Typical sampling variance (W)
  wi <- 1/diag(model$V)
  s2w <- (n - p) / (sum(wi) - sum(wi^2)/sum(wi))
  # Sum of all random effect variance components
  total_sigma2 <- sum(model$sigma2)
  i2 <- (total_sigma2 / (total_sigma2 + s2w)) * 100
  return(round(i2, 1))
}

# -------------------------------------------------------------------------
# 2. Main Function for Modeled Frequency Panels
# -------------------------------------------------------------------------
generate_modeled_freq_panel <- function(df, outcome_var, x_label, model_obj) {
  # Filter data for specific outcome
  df_sub <- df %>% filter(Outcome_Var == outcome_var & !is.na(Effect_Size_LnR))
  
  # --- Extract Stats from the Model ---
  model_estimate <- as.numeric(coef(model_obj)) # Modeled Mean (LnR)
  model_p        <- model_obj$pval              # Modeled P-value
  n_obs          <- model_obj$k                 # Sample Size
  i2_val         <- calc_i2(model_obj)          # Heterogeneity (I^2)
  
  # Back-transform Estimate to %
  xc_val <- (exp(model_estimate) - 1) * 100
  
  # Histogram math for curve scaling
  bin_width <- (max(df_sub$Effect_Size_LnR, na.rm=T) - min(df_sub$Effect_Size_LnR, na.rm=T)) / 20
  raw_sd    <- sd(df_sub$Effect_Size_LnR, na.rm=T)
  
  ggplot(df_sub, aes(x = Effect_Size_LnR)) + 
    geom_histogram(bins = 20, color = "white", fill = "gray60", alpha = 0.7) +
    
    # Normal curve centered at the MODEL ESTIMATE
    stat_function(fun = function(x) dnorm(x, mean = model_estimate, sd = raw_sd) * n_obs * bin_width, 
                  color = "black", size = 1) +
    
    geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 0.7) +
    
    # --- ANNOTATIONS ---
    # 1. xc line (Model Mean)
    annotate("text", x = -Inf, y = Inf, label = paste0("x[c] == ", round(xc_val, 2), "*'%'"),
             hjust = -0.2, vjust = 1.5, size = 3.5, parse = TRUE) +
    
    # # 2. R2 line 
    # annotate("text", x = -Inf, y = Inf, label = "R^2 == 0.93",
    #          hjust = -0.25, vjust = 3.2, size = 3.5, parse = TRUE) +
    
    # 3. I2 line (New Heterogeneity metric)
    annotate("text", x = -Inf, y = Inf, label = paste0("I^2 == ", i2_val, "*'%'"),
             hjust = -0.25, vjust = 2.8, size = 3.5, parse = TRUE) +
    
    # 4. P-value (Model P-value)
    annotate("text", x = -Inf, y = Inf, 
             label = ifelse(model_p < 0.001, "P < 0.001", paste0("P == ", round(model_p, 3))),
             hjust = -0.2, vjust = 5.7, size = 3.5, parse = TRUE) +
    
    # 5. Sample size (n)
    annotate("text", x = -Inf, y = Inf, label = paste0("n == ", n_obs),
             hjust = -0.3, vjust = 7.3, size = 3.5, parse = TRUE) +
    
    theme_bw() +
    labs(y = NULL, x = x_label) + # Y label removed here; added globally later
    theme(
      panel.grid = element_blank(),
      axis.title.x = element_text(size = 9, face = "bold"),
      axis.text = element_text(size = 8)
    )
}

# -------------------------------------------------------------------------
# 3. Fit Models and Generate Plot
# -------------------------------------------------------------------------

# First, ensure your data is pivoted
plot_data_freq <- meta_es %>%
  pivot_longer(
    cols = c(Yield_LnR, CP_LnR, WPI_LnR, TWP_LnR),
    names_to = "Outcome_Var", 
    values_to = "Effect_Size_LnR"
  )

# Fit your four multilevel models (Adjust variable names as per your dataset)
m1 <- rma.mv(Yield_LnR, Yield_Variance, random = ~ 1 | Study_ID/Location_ID, data = meta_es)
m2 <- rma.mv(CP_LnR,    CP_Variance,    random = ~ 1 | Study_ID/Location_ID, data = meta_es)
m3 <- rma.mv(WPI_LnR,  WPI_Variance,  random = ~ 1 | Study_ID/Location_ID, data = meta_es)
m4 <- rma.mv(TWP_LnR,   TWP_Variance,   random = ~ 1 | Study_ID/Location_ID, data = meta_es)

# Generate individual panels
pf1 <- generate_modeled_freq_panel(plot_data_freq, "Yield_LnR", "Effect size of forage yield (LRR)", m1)
pf2 <- generate_modeled_freq_panel(plot_data_freq, "CP_LnR",    "Effect size of crude protein (LRR)", m2)
pf3 <- generate_modeled_freq_panel(plot_data_freq, "WPI_LnR", "Effect size of irrigation water productivity (LRR)", m3)
pf4 <- generate_modeled_freq_panel(plot_data_freq, "TWP_LnR", "Effect size of total water productivity (LRR)", m4)

# Combine using patchwork
final_grid <- (pf1 | pf2) / (pf3 | pf4) + 
  plot_annotation(tag_levels = 'a', tag_suffix = ')') & 
  theme(plot.tag = element_text(size = 11, face = "bold"))

# Final Touch: Add the shared Y-axis label using cowplot
final_plot <- ggdraw(final_grid) +
  draw_label("Frequency", x = 0.008, y = 0.5, angle = 90, vjust = 0.5, size = 11, fontface = "bold")

print(final_plot)

ggsave("Frequency_plot_Grid.png", final_plot, width = 8, height = 8, dpi = 300)




# List of your models
# 1. Re-run the base nested models to get the SEs
res_yield <- rma.mv(yi = Yield_LnR, V = Yield_Variance, 
                    random = ~ 1 | Study_ID/Location_ID, data = meta_es)

res_cp    <- rma.mv(yi = CP_LnR, V = CP_Variance, 
                    random = ~ 1 | Study_ID/Location_ID, data = meta_es)

res_wpi  <- rma.mv(yi = WPI_LnR, V = WPI_Variance, 
                   random = ~ 1 | Study_ID/Location_ID, data = meta_es)

res_TWP   <- rma.mv(yi = TWP_LnR, V = TWP_Variance, 
                    random = ~ 1 | Study_ID/Location_ID, data = meta_es)

# 2. Extract and convert to Percentage Scale
get_stats <- function(model) {
  est <- as.numeric(model$beta)
  se  <- as.numeric(model$se)
  
  # xc is the percentage change
  xc <- (exp(est) - 1) * 100
  # SE for percentage change (Delta Method approximation)
  xc_se <- exp(est) * se * 100
  
  return(round(xc_se, 2))
}

# 3. Print the SEs for your "+-" labels
cat("Yield SE (%): ", get_stats(res_yield), "\n")
cat("CP SE (%):    ", get_stats(res_cp), "\n")
cat("WPI SE (%):  ", get_stats(res_wpi), "\n")
cat("TWP SE (%):   ", get_stats(res_TWP), "\n")

########################################### Q-Q plot ##########################################################################

# 1. Individual Panel Function (No changes needed)
generate_qq_panel <- function(df, outcome_var, specific_x_label) {
  df_sub <- df %>% filter(Outcome_Var == outcome_var & !is.na(Effect_Size_LnR))
  
  ggplot(df_sub, aes(sample = Effect_Size_LnR)) +
    stat_qq(size = 2, alpha = 0.6) + 
    stat_qq_line(color = "red", linetype = "dashed", linewidth = 0.8) +
    # Y is NULL because we add it with cowplot later
    labs(x = specific_x_label, y = NULL) + 
    theme_bw() +
    theme(
      panel.grid = element_blank(), 
      axis.title.x = element_text(size = 9, face = "bold"),
      axis.text = element_text(size = 8, color = "black"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
    )
}

# 2. Create the panels
pq1 <- generate_qq_panel(plot_data_freq, "Yield_LnR", "Normal quantiles for yield (LRR)")
pq2 <- generate_qq_panel(plot_data_freq, "CP_LnR",    "Normal quantiles for crude protein (LRR)")
pq3 <- generate_qq_panel(plot_data_freq, "WPI_LnR", expression(bold("Normal quantiles for ") * bold(WP)[bold(I)] * bold(" (LRR)")))
pq4 <- generate_qq_panel(plot_data_freq, "TWP_LnR", "Normal quantiles for TWP (LRR)")

# 3. Create the 2x2 Grid using Patchwork
# This handles the a, b, c, d tags perfectly
qq_grid <- (pq1 | pq2) / (pq3 | pq4) + 
  plot_annotation(
    tag_levels = 'a', 
    tag_suffix = ')'
  ) & 
  theme(
    plot.tag.position = c(0.03, 1.02),
    plot.tag = element_text(size = 11, face = "bold")
  )

# 4. Add the Shared Y-Axis Label using Cowplot
# This puts the label on the far left without "stealing" a letter tag
final_qq_plot <- cowplot::ggdraw(qq_grid) +
  cowplot::draw_label("Observed quantiles (LRR)", 
                      x = 0.01, y = 0.5, angle = 90, 
                      vjust = 0.5, size = 11, fontface = "bold")

# 5. Print the final result
print(final_qq_plot)

ggsave("QQ_plot_Grid.png", final_qq_plot, width = 8, height = 8, dpi = 300)

#################################### Robust Variance Estimation (RVE)##########################################################################
library(metafor)
library(clubSandwich)

# Define outcomes and their corresponding variance columns
outcomes_vars <- list(
  Yield = c("Yield_LnR", "Yield_Variance"),
  CP    = c("CP_LnR",    "CP_Variance"),
  WPI  = c("WPI_LnR",  "WPI_Variance"),
  TWP   = c("TWP_LnR",   "TWP_Variance")
)

# List to store results
robust_results <- list()

for (name in names(outcomes_vars)) {
  yi_col <- outcomes_vars[[name]][1]
  vi_col <- outcomes_vars[[name]][2]
  
  # 1. Fit the Nested Model (Location within Study)
  # Using your structure: random = ~ 1 | Study_ID / Location_ID
  res <- rma.mv(
    yi = meta_es[[yi_col]],
    V  = meta_es[[vi_col]],
    random = ~ 1 | Study_ID / Location_ID,
    data = meta_es
  )
  
  # 2. Apply Robust Variance Estimation (Clustered by Study_ID)
  # clubSandwich = TRUE is the gold standard for your n=193 dataset
  res_robust <- robust(res, cluster = meta_es$Study_ID, clubSandwich = TRUE)
  
  # Store result
  robust_results[[name]] <- res_robust
}

# Example: View Robust results for Yield
print(robust_results[["Yield"]])

summary_table <- data.frame()

for (name in names(robust_results)) {
  m <- robust_results[[name]]
  
  summary_table <- rbind(summary_table, data.frame(
    Outcome = name,
    Estimate = round(m$beta, 3),
    CI_Lower = round(m$ci.lb, 3),
    CI_Upper = round(m$ci.ub, 3),
    P_Value_Robust = round(m$pval, 4)
  ))
}

print(summary_table)

####################################### REML V/S ROBUST ##############################################################################

# 1. Define the variables and their corresponding variance columns
outcomes <- list(
  Yield = c("Yield_LnR", "Yield_Variance"),
  CP    = c("CP_LnR",    "CP_Variance"),
  WPI  = c("WPI_LnR",  "WPI_Variance"),
  TWP   = c("TWP_LnR",   "TWP_Variance")
)

# 2. Loop through each outcome to compare REML vs. Robust (RVE) results
for (name in names(outcomes)) {
  
  # Fit the nested REML model (Study_ID / Location_ID)
  m_reml <- rma.mv(yi = meta_es[[outcomes[[name]][1]]], 
                   V  = meta_es[[outcomes[[name]][2]]], 
                   random = ~ 1 | Study_ID/Location_ID, 
                   data = meta_es)
  
  # Apply Cluster-Robust Variance Estimation (RVE) with clubSandwich
  m_robust <- robust(m_reml, cluster = meta_es$Study_ID, clubSandwich = TRUE)
  
  # Extract standard errors and calculate CI widths
  reml_se <- m_reml$se
  rob_se  <- m_robust$se
  reml_w  <- 2 * 1.96 * reml_se
  rob_w   <- m_robust$ci.ub - m_robust$ci.lb
  
  # Print the comparison results
  cat(paste0("--- ", name, " ---\n"))
  cat("REML Width:  ", round(reml_w, 3), "\n")
  cat("Robust Width:", round(rob_w, 3), "\n")
  cat("SE Inflation:", round(((rob_se - reml_se) / reml_se) * 100, 2), "%\n\n")
}

################################################################################################
# 
# library(ggplot2)
# 
# # 1. Prepare the data from your RVE results
# # (Back-transforming LnR to Percent Change)
# forest_data <- data.frame(
#   Outcome = c("Yield", "Crude Protein (CP)", "WPI", "TWP"),
#   Estimate_Pct = c(-29.0, 0.7, 25.2, 8.3),
#   Lower_CI = c(-38.8, -1.3, 3.6, -12.8),
#   Upper_CI = c(-17.6, 2.6, 51.4, 34.6),
#   P_Value = c(0.0002, 0.485, 0.024, 0.443)
# )
# 
# # Add a significance label for the plot
# forest_data$Significance <- ifelse(forest_data$P_Value < 0.05, "Significant", "Non-Significant")
# 
# # 2. Create the Forest Plot
# ggplot(forest_data, aes(x = Outcome, y = Estimate_Pct)) +
#   # Add a dashed line at zero (No Effect)
#   geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
#   # Add the error bars (Confidence Intervals)
#   geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI, color = Significance), 
#                 width = 0.2, linewidth = 1.2) +
#   # Add the point estimate (The Mean)
#   geom_point(aes(fill = Significance), shape = 21, size = 5, color = "black") +
#   # Flip coordinates to have outcomes on the left (Traditional Forest Plot style)
#   coord_flip() +
#   # Labels and Theme
#   labs(
#     x = NULL, 
#     y = "Percentage Change (%)",
#     subtitle = "Effect size of treatment on forage outcomes"
#   ) +
#   scale_color_manual(values = c("Significant" = "red3", "Non-Significant" = "black")) +
#   scale_fill_manual(values = c("Significant" = "red", "Non-Significant" = "white")) +
#   theme_bw() +
#   theme(
#     panel.grid.major.y = element_blank(),
#     axis.title.x = element_text(face = "bold", size = 11),
#     axis.text = element_text(size = 10, color = "black", face = "bold"),
#     legend.position = "none" # Hide legend to keep it clean
#   ) +
#   # Optional: Add text labels showing the actual % on the plot
#   geom_text(aes(label = paste0(round(Estimate_Pct, 1), "%")), 
#             vjust = -1.5, size = 3.5, fontface = "bold")
# 

######################################### PUBLICATION BIASS#######################################
# Modified Egger's Test for Nested Data
# 1. Define outcomes and their variance columns
outcomes_list <- list(
  Yield = c("Yield_LnR", "Yield_Variance"),
  CP    = c("CP_LnR",    "CP_Variance"),
  WPI  = c("WPI_LnR",  "WPI_Variance"),
  TWP   = c("TWP_LnR",   "TWP_Variance")
)

# 2. Run the Modified Egger's Test for each
bias_results <- data.frame()

for (name in names(outcomes_list)) {
  yi_col <- outcomes_list[[name]][1]
  vi_col <- outcomes_list[[name]][2]
  
  # Fit the nested model with SE as a moderator (sqrt of variance)
  m_egger <- rma.mv(
    yi = meta_es[[yi_col]], 
    V  = meta_es[[vi_col]], 
    mods = ~ sqrt(meta_es[[vi_col]]), # The SE moderator
    random = ~ 1 | Study_ID/Location_ID, 
    data = meta_es
  )
  
  # Apply Robust Variance Estimation (Clustered by Study_ID)
  m_robust <- robust(m_egger, cluster = meta_es$Study_ID, clubSandwich = TRUE)
  
  # Extract the p-value for the moderator (Publication Bias Test)
  # In metafor mods models, the moderator is usually the 2nd coefficient
  p_val_bias <- m_robust$pval[2] 
  
  bias_results <- rbind(bias_results, data.frame(
    Outcome = name,
    Egger_P_Value = round(p_val_bias, 4),
    Interpretation = ifelse(p_val_bias > 0.05, "No Evidence of Bias", "Potential Bias Found")
  ))
}

print(bias_results)

library(grid)
# 1. Prepare your data for the plot (Use your model results here)
# We calculate the 'funnel' bounds (95% CI) manually for ggplot
generate_funnel_data <- function(model_obj, outcome_label, p_val_bias) {
  df <- data.frame(
    yi = model_obj$yi,
    se = sqrt(model_obj$vi),
    Outcome = outcome_label
  )
  
  # Calculate the funnel cone
  max_se <- max(df$se)
  se_seq <- seq(0, max_se, length.out = 100)
  ll95 <- model_obj$beta - 1.96 * se_seq
  ul95 <- model_obj$beta + 1.96 * se_seq
  
  cone <- data.frame(se = se_seq, lower = as.numeric(ll95), upper = as.numeric(ul95))
  
  list(data = df, cone = cone, p = p_val_bias, n = model_obj$k, estimate = as.numeric(model_obj$beta))
}

# 2. Function to create each individual panel

plot_single_funnel <- function(obj_list, x_lab, tag) {
  ggplot(obj_list$data, aes(x = yi, y = se)) +
    # The Funnel Cone
    geom_line(data = obj_list$cone, aes(x = lower, y = se), linetype = "dashed", color = "gray40") +
    geom_line(data = obj_list$cone, aes(x = upper, y = se), linetype = "dashed", color = "gray40") +
    # Individual Studies
    geom_point(alpha = 0.6, shape = 21, fill = "gray70", size = 1.5) +
    # Vertical line at mean effect
    geom_vline(xintercept = obj_list$estimate, color = "red", linetype = "solid") +
    
    # --- FIXED ANNOTATION BLOCK ---
    # We use x = -Inf and y = -Inf because scale_y_reverse flips the coordinate system
    annotate("text", x = -Inf, y = -Inf, 
             label = paste0("P = ", format(obj_list$p, nsmall = 3), "\nn = ", obj_list$n),
             hjust = -0.2, vjust = 3.2, size = 3.5, fontface = "bold", color = "black") +
    
    scale_y_reverse() + # Reversed Y-axis (Standard Error 0 at top)
    labs(x = x_lab, y = NULL, tag = tag) +
    theme_bw() +
    theme(
      plot.tag = element_text(face = "bold", size = 12),
      plot.tag.position = c(0.005, 1.02), # Moves a, b, c, d to the very top left
      panel.grid = element_blank(),
      axis.title.x = element_text(size = 9, face = "bold"),
      axis.text = element_text(size = 8)
    )
}

# Update p_val_bias with the exact values from your bias_results table
d1 <- generate_funnel_data(m1, "Yield", 0.107)
d2 <- generate_funnel_data(m2, "CP", 0.320)
d3 <- generate_funnel_data(m3, "WPI", 0.472)
d4 <- generate_funnel_data(m4, "TWP", 0.058)

# 2. Re-generate the 4 panels
p1 <- plot_single_funnel(d1, "Effect size of forage yield (LRR)", "a)")
p2 <- plot_single_funnel(d2, "Effect size of crude protein (LRR)", "b)")
p3 <- plot_single_funnel(d3, expression(bold("Effect size of ") * bold(WP)[bold(I)] * bold(" (LRR)")), "c)")
p4 <- plot_single_funnel(d4, "Effect size of TWP (LRR)", "d)")

# 3. Re-combine with the shared Y-axis label workaround
main_grid <- (p1 | p2) / (p3 | p4)
y_label <- grid::textGrob("Standard Error", rot = 90, gp = gpar(fontsize = 11, fontface = "bold"))

final_funnel_plot <- wrap_elements(panel = y_label) + main_grid + 
  plot_layout(widths = c(1, 40)) 

print(final_funnel_plot)

ggsave("Funnel_grid.png", final_funnel_plot, width = 8, height = 8, dpi = 300)

############################ Trim and fill ##############################################

# =====================================================================
# TRIM-AND-FILL (study-level aggregated effect sizes)
# =====================================================================
library(metafor)
library(dplyr)

outcomes_list <- list(
  Yield = c("Yield_LnR", "Yield_Variance"),
  CP    = c("CP_LnR",    "CP_Variance"),
  WPI   = c("WPI_LnR",   "WPI_Variance"),
  TWP   = c("TWP_LnR",   "TWP_Variance")
)

run_trimfill <- function(rho) {
  results <- data.frame()
  models  <- list()
  
  for (name in names(outcomes_list)) {
    yi_col <- outcomes_list[[name]][1]
    vi_col <- outcomes_list[[name]][2]
    
    dat <- meta_es %>%
      filter(!is.na(.data[[yi_col]]), !is.na(.data[[vi_col]])) %>%
      transmute(Study_ID, yi = .data[[yi_col]], vi = .data[[vi_col]])
    
    dat <- escalc(measure = "GEN", yi = yi, vi = vi, data = dat)
    
    # One aggregated effect size per study (assumed within-study correlation = rho)
    agg <- aggregate(dat, cluster = dat$Study_ID, struct = "CS", rho = rho)
    
    res <- rma(yi, vi, data = agg, method = "REML")
    tf  <- trimfill(res)
    models[[name]] <- tf
    
    results <- rbind(results, data.frame(
      Outcome      = name,
      rho          = rho,
      k_studies    = res$k,
      Filled_k0    = tf$k0,
      Side_filled  = tf$side,
      Pct_original = round((exp(as.numeric(coef(res))) - 1) * 100, 2),
      P_original   = round(res$pval, 4),
      Pct_adjusted = round((exp(as.numeric(coef(tf)))  - 1) * 100, 2),
      CI_adj_lower = round((exp(tf$ci.lb) - 1) * 100, 2),
      CI_adj_upper = round((exp(tf$ci.ub) - 1) * 100, 2),
      P_adjusted   = round(tf$pval, 4)
    ))
  }
  list(results = results, models = models)
}

# --- Main analysis (rho = 0.5) ---
tf_main <- run_trimfill(rho = 0.5)
print(tf_main$results)

# --- Sensitivity to the assumed correlation ---
tf_sens <- bind_rows(lapply(c(0.3, 0.5, 0.8), function(r) run_trimfill(r)$results))
print(tf_sens %>% select(Outcome, rho, Filled_k0, Pct_original, Pct_adjusted, P_adjusted))

# --- Funnel plots with filled studies (open circles = imputed) ---
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
titles <- c(Yield = "a) Yield", CP = "b) Crude protein", WPI = "c) WPI", TWP = "d) TWP")
for (name in names(tf_main$models)) {
  funnel(tf_main$models[[name]], main = titles[name], xlab = "LRR (study-level)")
}
par(mfrow = c(1, 1))

# --- Save results ---
write.csv(tf_main$results, "TrimFill_results.csv", row.names = FALSE)
write.csv(tf_sens,         "TrimFill_rho_sensitivity.csv", row.names = FALSE)

library(ggplot2)
library(patchwork)
library(grid)

# --- 1. Prepare data from a trimfill object ---
make_tf_data <- function(tf) {
  filled <- if (is.null(tf$fill)) rep(FALSE, length(tf$yi)) else tf$fill
  df <- data.frame(yi = as.numeric(tf$yi), se = sqrt(tf$vi), Imputed = filled)
  
  est    <- as.numeric(coef(tf))
  se_seq <- seq(0, max(df$se), length.out = 100)
  cone   <- data.frame(se = se_seq, lower = est - 1.96 * se_seq, upper = est + 1.96 * se_seq)
  
  list(data = df, cone = cone, estimate = est)
}

# --- 2. One panel ---
plot_tf_funnel <- function(tf, res_row, x_lab, tag) {
  d <- make_tf_data(tf)
  
  lab <- paste0("Imputed studies = ", res_row$Filled_k0,
                "\nk = ", res_row$k_studies, " studies")
  if (res_row$Filled_k0 > 0) {
    lab <- paste0(lab, "\nAdjusted: ", sprintf("%+.1f", res_row$Pct_adjusted), "%")
  }
  
  ggplot(d$data, aes(x = yi, y = se)) +
    geom_line(data = d$cone, aes(x = lower, y = se), linetype = "dashed", color = "gray40") +
    geom_line(data = d$cone, aes(x = upper, y = se), linetype = "dashed", color = "gray40") +
    geom_vline(xintercept = d$estimate, color = "red", linewidth = 0.7) +
    geom_point(data = subset(d$data, !Imputed), shape = 21, fill = "gray70", size = 2) +
    geom_point(data = subset(d$data,  Imputed), shape = 21, fill = "white", color = "black", size = 2) +
    annotate("text", x = -Inf, y = -Inf, label = lab,
             hjust = -0.1, vjust = 1.4, size = 3.2, fontface = "bold", lineheight = 0.95) +
    scale_y_reverse() +
    labs(x = x_lab, y = NULL, tag = tag) +
    theme_bw() +
    theme(plot.tag = element_text(face = "bold", size = 12),
          plot.tag.position = c(0.005, 1.02),
          panel.grid = element_blank(),
          axis.title.x = element_text(size = 9, face = "bold"),
          axis.text = element_text(size = 8))
}

# --- 3. Panels ---
r <- tf_main$results
t1 <- plot_tf_funnel(tf_main$models$Yield, r[r$Outcome == "Yield", ],
                     "Effect size of forage yield (LRR)", "a)")
t2 <- plot_tf_funnel(tf_main$models$CP,    r[r$Outcome == "CP", ],
                     "Effect size of crude protein (LRR)", "b)")
t3 <- plot_tf_funnel(tf_main$models$WPI,   r[r$Outcome == "WPI", ],
                     expression(bold("Effect size of ") * bold(WP)[bold(I)] * bold(" (LRR)")), "c)")
t4 <- plot_tf_funnel(tf_main$models$TWP,   r[r$Outcome == "TWP", ],
                     "Effect size of TWP (LRR)", "d)")

# --- 4. Combine with shared y-axis label ---
tf_grid  <- (t1 | t2) / (t3 | t4)
y_label  <- textGrob("Standard Error", rot = 90, gp = gpar(fontsize = 11, fontface = "bold"))
final_tf_plot <- wrap_elements(panel = y_label) + tf_grid + plot_layout(widths = c(1, 40))

print(final_tf_plot)

# --- 5. Save ---
ggsave("Trim_and_fill_grid.png", final_tf_plot, width = 8, height = 8, dpi = 300)

##########################################################################################

# --- SECTION A: Run All Meta-Analysis Models ---

outcomes <- c("Yield", "CP", "WPI", "TWP")
moderators <- c("Irrigation_method", "Irri_amt", "Degree_of_irri_deficit_grp", 
                "Irrigation_scheduling", "MAP_mm_grp", "MAT_degreeC_grp", 
                "Texture", "Soil_BD_grp", "Soil_pH_grp", "Soil_OM_grp", 
                "Winterhardiness_grp", "Latitude_grp", "Longitude_grp")

model_list <- list()
for (out in outcomes) {
  model_list[[out]] <- list()
  for (mod in moderators) {
    yi_col <- paste0(out, "_LnR")
    V_col  <- paste0(out, "_Variance")
    W_col  <- paste0(out, "_weight")
    
    # Construct the formula *without* an intercept to get means for all levels
    form_no_intercept <- as.formula(paste("~ 0 +", mod))
    
    res <- try(rma.mv(
      yi = meta_es[[yi_col]],
      V  = meta_es[[V_col]],
      W  = meta_es[[W_col]],
      mods = form_no_intercept, 
      random = ~ 1 | Study_ID / Location_ID,
      data = meta_es
    ), silent = TRUE)
    
    if (!inherits(res, "try-error")) {
      model_list[[out]][[mod]] <- res
    }
  }
}


# --- SECTION B: Extract Data from Models for Plotting ---

# Function to safely extract subgroup means and CIs (for models WITHOUT intercept)
extract_subgroup_data_with_n <- function(model, moderator_name, outcome_name) {
  if (inherits(model, "try-error") || is.null(model)) return(NULL)
  
  coef_names <- names(coef(model))
  subgroups <- gsub(moderator_name, "", coef_names, fixed = TRUE) %>% trimws()
  
  # Calculate N (observations) for each subgroup from the raw data
  # We filter the raw data for each level to count rows
  counts <- map_dbl(subgroups, function(lvl) {
    sum(meta_es[[moderator_name]] == lvl, na.rm = TRUE)
  })
  
  # Predict means and CIs
  new_data <- diag(length(coef_names))
  colnames(new_data) <- coef_names
  preds <- predict(model, newmods = new_data)
  
  # --- CONVERT TO PERCENTAGE ---
  # Formula: % Change = (exp(LnR) - 1) * 100
  est_percent <- (exp(as.numeric(preds$pred)) - 1) * 100
  lo_percent  <- (exp(as.numeric(preds$ci.lb)) - 1) * 100
  hi_percent  <- (exp(as.numeric(preds$ci.ub)) - 1) * 100
  
  data.frame(
    Outcome = outcome_name,
    Moderator = moderator_name,
    Subgroup = subgroups,
    Estimate = est_percent,
    CI_Lower = lo_percent,
    CI_Upper = hi_percent,
    N_Obs = counts,
    Label = paste0(subgroups, " (n=", counts, ")") 
  )
}

all_subgroup_data <- map_dfr(names(model_list), function(out_name) {
  map_dfr(names(model_list[[out_name]]), function(mod_name) {
    extract_subgroup_data_with_n(model_list[[out_name]][[mod_name]], mod_name, out_name)
  })
})

############################### FOR ONLY IRRIGATION RELATED PLOT ##########################################
# Filter for the 4 specific moderators to block together
target_mods <- c("Irrigation_method", "Irri_amt", "Degree_of_irri_deficit_grp", "Irrigation_scheduling")
plot_data <- all_subgroup_data %>% 
  filter(Moderator %in% target_mods)

# Adjust these based on your specific results to ensure all error bars are visible
yield_lims <- c(-80, 80)
cp_lims    <- c(-15, 15)
wpi_lims  <- c(-150, 150)
TWP_lims   <- c(-150, 150)


# Left side = Current name in data, Right side = New name for plot
mod_names <- c(
  "Irrigation_method"          = "Irrigation method",
  "Irri_amt"                   = "Irrigation amount",
  "Degree_of_irri_deficit_grp" = "Degree of deficit",
  "Irrigation_scheduling"      = "Irrigation scheduling"
)


generate_stacked_plot <- function(df, outcome_name, x_label, x_lims) {
  df_sub <- df %>% filter(Outcome == outcome_name)
  
  ggplot(df_sub, aes(y = Label, x = Estimate)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "blue") +
    geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0.2, color = "black") +
    geom_point(shape = 15, size = 3) +
    # Added labeller here to rename the moderators
    facet_grid(Moderator ~ ., 
               scales = "free_y", 
               space = "free_y", 
               switch = "y",
               labeller = as_labeller(mod_names)) + 
    theme_bw() +
    labs(x = x_label, y = "") +
    xlim(x_lims) +
    theme(
      strip.text.y.left = element_text(angle = 0, face = "bold", size = 9),
      strip.background = element_rect(fill = "gray95"),
      strip.placement = "outside",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 8)
    )
}


# Define the desired order for all subgroup levels (Bottom to Top)
desired_order <- c(
  # Degree_of_irri_deficit_grp
  "Slight (>0.70)", "Moderate (0.35-0.70)", "Severe (<0.35)", 
  
  # MAT
  ">10 °C", "≤10 °C", 
  
  # Irri_amt
  "Over-irrigation", "Limited irrigation", "Zero-irrigation", 
  
  # Irrigation_method
  "Surface irrigation", "Sprinkler irrigation", "Subsurface drip irrigation", 
  
  # Irrigation_scheduling
  "Other", "Supplemental irrigation based", "Soil moisture based", "Evapotranspiration based"
)


# 2. Apply the order to the Subgroup and sync the Label
plot_data <- plot_data %>%
  # Force Subgroup to follow the manual list
  mutate(Subgroup = fct_relevel(as.character(Subgroup), desired_order)) %>%
  # Sort the dataframe by the new factor order
  arrange(Subgroup) %>%
  # Ensure the Label (which includes the 'n=') follows the Subgroup order
  mutate(Label = fct_inorder(as.character(Label)))


# Using the same generate_stacked_plot function as before
p1 <- generate_stacked_plot(plot_data, "Yield", "Effect on forage yield (%)", yield_lims)
p2 <- generate_stacked_plot(plot_data, "CP",    "Effect on crude protein (%)",    cp_lims)
p3 <- generate_stacked_plot(plot_data, "WPI",  "Effect on irrigation water productivity (%)",  wpi_lims)
p4 <- generate_stacked_plot(plot_data, "TWP",   "Effect on total water productivity (%)",   TWP_lims)

# This makes the labels from panels (a) and (c) apply to the whole row
p2_clean <- p2 + theme(
  strip.background = element_blank(),
  strip.text.y.left = element_blank()
)

p4_clean <- p4 + theme(
  strip.background = element_blank(),
  strip.text.y.left = element_blank()
)

# Panel (a) and (b) in row 1, Panel (c) and (d) in row 2
irrigation_grid <- (p1 | p2_clean) / (p3 | p4_clean) + 
  plot_annotation(
    tag_levels = 'a', 
    tag_suffix = ')'
  )
print(irrigation_grid)
ggsave("Irrigation_Grid.png", irrigation_grid, width = 12, height = 8, dpi = 300)



################ FOR ONLY CLIMATE< SOIL AND WINTERHARDINESS RELATED PLOT #################################

# --- 1. Filter specifically for the Environmental/Soil moderators ---
# We use a unique name 'plot_data_weather' to avoid overwriting issues
target_mods <- c("MAP_mm_grp", "MAT_degreeC_grp", "Texture", "Soil_BD_grp",
                 "Soil_pH_grp", "Soil_OM_grp", 
                 "Winterhardiness_grp", "Latitude_grp", "Longitude_grp")
plot_data <- all_subgroup_data %>% 
  filter(Moderator %in% target_mods)

# Define the limits for each outcome (adjust as needed to fix "Removed 1 row" warnings)
yield_lims <- c(-100, 50)
cp_lims    <- c(-20, 20)
wpi_lims  <- c(-100, 250)
TWP_lims   <- c(-100, 250)

# Create a "translation" vector for your moderator names
mod_names <- c(
  "MAP_mm_grp"          = "Mean annual\nprecipitation",
  "MAT_degreeC_grp"     = "Mean annual\ntemperature",
  "Texture"             = "Soil texture",
  "Soil_BD_grp"         = "Soil bulk density",
  "Soil_pH_grp"         = "Soil pH",
  "Soil_OM_grp"         = "Soil\norganic matter",
  "Winterhardiness_grp" = "Winterhardiness",
  "Latitude_grp" = "Latitude",
  "Longitude_grp" = "Longitude"
)

# 2. Updated function with the labeller
generate_stacked_plot <- function(df, outcome_name, x_label, x_lims) {
  df_sub <- df %>% filter(Outcome == outcome_name)
  
  ggplot(df_sub, aes(y = Label, x = Estimate)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "blue") +
    geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0.2, color = "black") +
    geom_point(shape = 15, size = 3) +
    # Added labeller here to rename the moderators
    facet_grid(Moderator ~ ., 
               scales = "free_y", 
               space = "free_y", 
               switch = "y",
               labeller = as_labeller(mod_names)) + 
    theme_bw() +
    labs(x = x_label, y = "") +
    xlim(x_lims) +
    theme(
      strip.text.y.left = element_text(angle = 0, face = "bold", size = 9),
      strip.background = element_rect(fill = "gray95"),
      strip.placement = "outside",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 8)
    )
}

# Define the desired order for all subgroup levels (Bottom to Top)
desired_order <- c(
  #Latitude
  "Southern region", "Mid region", "Northern region",
  
  #Longitude
  "Arid west", "Humid east",
  
  # MAP
  ">500 mm", "250–500 mm",  "<250 mm", 
  
  # MAT
  ">10 °C", "≤10 °C", 
  
  # Soil Bulk Density
  ">1.5 g cm^-3", "1.4–1.5 g cm^-3",  "<1.4 g cm^-3", 
  
  # Soil Organic Matter
  ">1 %", "≤1 %", 
  
  # Soil pH 
  ">7", "≤7", 
  
  # Soil Texture
  "Coarse", "Medium", "Fine", 
  
  # Winterhardiness
  "8a–9b","5a–7b",  "4a–4b"
) 



meta_es <- meta_es %>%
  mutate(
    Latitude_grp = case_when(
      Latitude < 35 ~ "Southern region",
      Latitude >= 35 & Latitude <= 45 ~ "Mid region",
      Latitude > 45 ~ "Northern region",
      TRUE ~ NA_character_
    ),
    # Create Longitude groups
    Longitude_grp = case_when(
      Longitude > -100 ~ "Humid east",
      Longitude < -100 ~ "Arid west",
      TRUE ~ NA_character_
    )
  )

# 2. Apply the order to the Subgroup and sync the Label
plot_data <- plot_data %>%
  # Force Subgroup to follow the manual list
  mutate(Subgroup = fct_relevel(as.character(Subgroup), desired_order)) %>%
  # Sort the dataframe by the new factor order
  arrange(Subgroup) %>%
  # Ensure the Label (which includes the 'n=') follows the Subgroup order
  mutate(Label = fct_inorder(as.character(Label)))



# Using the same generate_stacked_plot function as before
p1 <- generate_stacked_plot(plot_data, "Yield", "Effect on forage yield (%)", yield_lims)
p2 <- generate_stacked_plot(plot_data, "CP",    "Effect on crude protein (%)",    cp_lims)
p3 <- generate_stacked_plot(plot_data, "WPI",  "Effect on irrigation water productivity (%)",  wpi_lims)
p4 <- generate_stacked_plot(plot_data, "TWP",   "Effect on total water productivity (%)",   TWP_lims)

# --- 2. Remove moderator labels from the right-hand panels (b and d) ---
# This makes the labels from panels (a) and (c) apply to the whole row
p2_clean <- p2 + theme(
  strip.background = element_blank(),
  strip.text.y.left = element_blank()
)

p4_clean <- p4 + theme(
  strip.background = element_blank(),
  strip.text.y.left = element_blank()
)

# Panel (a) and (b) in row 1, Panel (c) and (d) in row 2
weather_grid <- (p1 | p2_clean) / (p3 | p4_clean) + 
  plot_annotation(
    tag_levels = 'a', 
    tag_suffix = ')'
  )
print(weather_grid)

ggsave("Weather_Soil_Grid.png", weather_grid, width = 12, height = 10, dpi = 300)

###################################################################################################################

# --- Function to extract SE in percentage scale ---
extract_se_percentage <- function(model_list, outcomes, target_mods) {
  se_summary <- data.frame()
  
  for (out in outcomes) {
    for (mod in target_mods) {
      m <- model_list[[out]][[mod]]
      
      if (!is.null(m) && !inherits(m, "try-error")) {
        # 1. Get coefficients and log-scale SE
        coef_names <- names(coef(m))
        subgroups  <- gsub(mod, "", coef_names, fixed = TRUE)
        log_se     <- m$se  # Standard errors from the model
        log_beta   <- as.numeric(m$beta)
        
        # 2. Convert to Percentage Scale SE (Delta Method)
        # Formula: SE_percent = exp(beta) * SE_log * 100
        se_percent <- exp(log_beta) * log_se * 100
        
        temp <- data.frame(
          Outcome = out,
          Moderator = mod,
          Subgroup = subgroups,
          Estimate_Pct = round((exp(log_beta) - 1) * 100, 2),
          SE_PlusMinus = round(se_percent, 2)
        )
        se_summary <- rbind(se_summary, temp)
      }
    }
  }
  return(se_summary)
}

# --- Execute and Print ---
# target_mods <- c("Irrigation_method", "Irri_amt", "Degree_of_irri_deficit_grp", "Irrigation_scheduling")
target_mods <- c("MAP_mm_grp", "MAT_degreeC_grp", "Texture", "Soil_BD_grp",
                 "Soil_pH_grp", "Soil_OM_grp",
                 "Winterhardiness_grp", "Latitude_grp", "Longitude_grp")
se_results <- extract_se_percentage(model_list, outcomes, target_mods)

# This prints a table you can use to fill in your "±" values
print(se_results)


####################### Building Regression models #################################################################

# --- 1. Filter for Irrigated Data (No NAs, No Zeroes) ---
meta_es_irrigated <- meta_es %>% 
  filter(!is.na(Irri_amt_mm_trt) & Irri_amt_mm_trt > 0)

responses <- c("Yield_LnR", "CP_LnR")
moderators <- c("Irrigation_method", "Irri_amt", "Degree_of_irri_deficit", 
                "Irrigation_scheduling", "MAP_mm", "MAT_degreeC", 
                "Texture", "Soil_BD_g.cm.3", "Soil_pH", "Soil_OM_.", 
                "Winterhardiness", "Latitude", "Longitude")

individual_results <- list()

# --- 2. Run Stepwise Models ---
for (resp in responses) {
  temp_data <- meta_es_irrigated[, c(resp, moderators)]
  data_clean <- na.omit(temp_data)
  
  # Check for moderators with more than one level
  valid_mods <- moderators[sapply(data_clean[moderators], function(x) length(unique(x)) > 1)]
  
  if (length(valid_mods) > 0) {
    f_full <- as.formula(paste(resp, "~", paste(valid_mods, collapse = " + ")))
    f_null <- as.formula(paste(resp, "~ 1"))
    
    m_full <- lm(f_full, data = data_clean)
    m_null <- lm(f_null, data = data_clean)
    
    individual_results[[paste0(resp, "_back")]] <- step(m_full, direction = "backward", trace = 0)
    individual_results[[paste0(resp, "_fwd")]]  <- step(m_null, 
                                                        scope = list(lower = f_null, upper = f_full), 
                                                        direction = "forward", trace = 0)
    message("Successfully modeled irrigated: ", resp)
  }
}

# --- 3. Extract Summaries into Table ---
ind_summary_list <- list() # Initialize the missing list

for (model_name in names(individual_results)) {
  current_model <- individual_results[[model_name]]
  selected_vars <- labels(terms(current_model))
  if(length(selected_vars) == 0) selected_vars <- "None (Intercept only)"
  
  ind_summary_list[[model_name]] <- data.frame(
    Model_Name = model_name,
    Moderators_Selected = paste(selected_vars, collapse = ", "),
    Adj_R2 = round(summary(current_model)$adj.r.squared, 3),
    AIC = round(AIC(current_model), 2),
    BIC = round(BIC(current_model), 2),
    R2 = round(summary(current_model)$r.squared, 3),
    n = length(residuals(current_model))
  )
}

# --- 4. Combine and Display ---
summary_table_final <- do.call(rbind, ind_summary_list)
rownames(summary_table_final) <- NULL
print(summary_table_final)

# Extract the coefficients from your forward model
yield_coef <- coef(individual_results$Yield_LnR_fwd)
CP_coef <- coef(individual_results$CP_LnR_back)
print(yield_coef)
print(CP_coef)

# For the individual coefficients (p-values for each category)
summary(individual_results$Yield_LnR_fwd)
summary(individual_results$CP_LnR_back)

# For the overall model significance (The F-statistic p-value)
summary(individual_results$Yield_LnR_fwd)$fstatistic

###########################################################################################################################


# # --- STEP 1: Reshaping with Unique Row IDs for Total Water Analysis ---
# plot_data_total <- meta_es %>%
#   mutate(Row_ID = row_number()) %>% 
#   select(Row_ID, Study_ID, MAP_mm_grp, 
#          Irri_amt_mm_trt, Irri_amt_mm_Ctrl, 
#          Yield_kg_trt, Yield_kg_ctrl, 
#          CP_trt, CP_ctrl, 
#          WPI_kg.m3_trt, WPI_kg.m3_ctrl) %>%
#   
#   pivot_longer(cols = c(Irri_amt_mm_trt, Irri_amt_mm_Ctrl),
#                names_to = "Type_X", values_to = "Total_Water") %>%
#   pivot_longer(cols = c(Yield_kg_trt, Yield_kg_ctrl),
#                names_to = "Type_Y1", values_to = "Yield_Value") %>%
#   pivot_longer(cols = c(CP_trt, CP_ctrl),
#                names_to = "Type_Y2", values_to = "CP_Value") %>%
#   pivot_longer(cols = c(WPI_kg.m3_trt, WPI_kg.m3_ctrl),
#                names_to = "Type_Y3", values_to = "WPI_Value") %>%
#   
#   filter(
#     gsub(".*_", "", tolower(Type_X)) == gsub(".*_", "", tolower(Type_Y1)),
#     gsub(".*_", "", tolower(Type_X)) == gsub(".*_", "", tolower(Type_Y2)),
#     gsub(".*_", "", tolower(Type_X)) == gsub(".*_", "", tolower(Type_Y3))
#   ) %>%
#   select(-Type_X, -Type_Y1, -Type_Y2, -Type_Y3)
# 
# # --- STEP 2: Pre-compute N per facet (numeric, not label) ---
# n_counts_total <- plot_data_total %>%
#   group_by(MAP_mm_grp) %>%
#   summarise(n_val = sum(!is.na(Total_Water)), .groups = "drop")
# 
# plot_data_total <- left_join(plot_data_total, n_counts_total, by = "MAP_mm_grp")
# 
# # --- STEP 3: Fixed Function for Box-Based Alignment ---
# generate_total_water_reg_plot <- function(df, y_var, y_label, y_lims,
#                                           stats_y_pos = "top", # Use "top" or "bottom"
#                                           show_x_title = FALSE) {
#   
#   # Determine if we are bottom-aligning for the N label
#   is_bottom <- stats_y_pos == "bottom"
#   
#   ggplot(df, aes(x = Total_Water, y = .data[[y_var]])) +
#     geom_point(alpha = 0.4, size = 1.5, shape = 21, fill = "gray80") +
#     geom_smooth(method = "lm", formula = y ~ x + I(x^2), color = "black", size = 0.8) +
#     
#     # --- 1. Regression Equation (Anchored to BOX, not points) ---
#     stat_poly_eq(
#       formula = y ~ x + I(x^2),
#       aes(label = paste(after_stat(eq.label), after_stat(rr.label), 
#                         after_stat(p.value.label), sep = "*', '*")),
#       parse = TRUE,
#       label.x = "right",    # Keywords force it to the physical edge
#       label.y = stats_y_pos, # Use "bottom" for CP
#       size = 3,
#       hjust = 1.2          # Nudges it slightly left so it doesn't hit the border
#     ) +
#     
#     # --- 2. N Label (Stays in corner) ---
#     geom_text(
#       data = distinct(df, MAP_mm_grp, n_val),
#       aes(label = paste0("(n = ", n_val, ")")),
#       x = Inf, 
#       y = if(is_bottom) -Inf else Inf,
#       hjust = 1.1, 
#       vjust = if(is_bottom) -3.5 else 2.3, 
#       size = 3,
#       inherit.aes = FALSE
#     ) +
#     
#     coord_cartesian(xlim = c(0, 2000), ylim = y_lims, clip = "off") +
#     facet_wrap(~ MAP_mm_grp, scales = "free_y") +
#     theme_bw() +
#     labs(y = NULL, x = if (show_x_title) "Irrigation water (mm)" else NULL) +
#     theme(
#       panel.grid.major = element_blank(), # Removes major grid lines
#       panel.grid.minor = element_blank(), # Removes minor grid lines
#       strip.background = element_rect(fill = "gray95"),
#       strip.text = element_text(face = "bold")
#     )
# }
# 
# 
# # --- STEP 4: Generate plots ---
# 
# # 1. Yield – top (Using "top" keyword)
# p1_yield <- generate_total_water_reg_plot(
#   plot_data_total,
#   "Yield_Value",
#   "Yield (kg ha⁻¹)",
#   c(0, 30000),
#   stats_y_pos = "top"  
# )
# 
# # 2. CP – middle (Using "bottom" keyword)
# p2_cp <- generate_total_water_reg_plot(
#   plot_data_total,
#   "CP_Value",
#   "Crude Protein (%)",
#   c(0, 30),
#   stats_y_pos = "bottom" 
# ) +
#   theme(strip.text = element_blank(),
#         strip.background = element_blank())
# 
# # 3. WPI – bottom (Using "top" keyword)
# p3_WPI <- generate_total_water_reg_plot(
#   plot_data_total,
#   "WPI_Value",
#   "WPI (kg ha⁻¹ mm⁻¹)",
#   c(0, 180),
#   stats_y_pos = "top", 
#   show_x_title = TRUE
# ) +
#   theme(strip.text = element_blank(),
#         strip.background = element_blank())
# 
# 
# # --- STEP 5: Combine ---
# final_grid <- (p1_yield / p2_cp / p3_WPI) +
#   plot_annotation(tag_levels = "a", tag_suffix = ")") &
#   theme(
#     plot.tag = element_text(face = "bold"),
#     axis.title.y = element_text(face = "bold", size = 10),
#     # Bold the common X-axis label here
#     axis.title.x = element_text(face = "bold", size = 11), 
#     plot.margin = margin(t = 5, r = 5, b = 5, l = 25) # Increased left margin for safety
#   )
# 
# # Explicitly re-apply labels to ensure they stick
# final_grid[[1]] <- final_grid[[1]] + 
#   ylab("Yield (kg ha⁻¹)") + 
#   theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
# 
# final_grid[[2]] <- final_grid[[2]] + 
#   ylab("Crude Protein (%)") + 
#   theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
# 
# final_grid[[3]] <- final_grid[[3]] + 
#   ylab("Irrigation water productivity (kg m\u207b\u00b3)")
# 
# print(final_grid)


# =====================================================================
# FIGURE 6 – Mixed-effects quadratic regressions (accounts for study clustering)
# =====================================================================
library(lme4)
library(performance)

# --- STEP 1: Long data – one row per unique treatment or control observation ---
trt_obs <- meta_es %>%
  transmute(Study_ID, Location_ID, Year, MAP_mm_grp,
            Total_Water = Irri_amt_mm_trt,
            Yield_Value = Yield_kg_trt,
            CP_Value    = CP_trt,
            WPI_Value   = WPI_kg.m3_trt)

ctrl_obs <- meta_es %>%
  transmute(Study_ID, Location_ID, Year, MAP_mm_grp,
            Total_Water = Irri_amt_mm_Ctrl,
            Yield_Value = Yield_kg_ctrl,
            CP_Value    = CP_ctrl,
            WPI_Value   = WPI_kg.m3_ctrl)

plot_data_total <- bind_rows(trt_obs, ctrl_obs) %>%
  distinct() %>%                                   # shared controls counted only once
  filter(!is.na(Total_Water)) %>%
  mutate(Study_ID    = factor(Study_ID),
         Location_ID = factor(Location_ID))

map_levels <- levels(factor(meta_es$MAP_mm_grp))   # panel order

# --- STEP 2: Fit one mixed quadratic model (one outcome, one MAP group) ---
get_r2 <- function(r, nm) {
  v <- tryCatch(as.numeric(r[[nm]]), error = function(e) NA)
  if (length(v) == 0) NA else v
}

fmt_num <- function(v) {
  if (v == 0 || abs(v) >= 0.001) return(format(signif(v, 3), trim = TRUE, scientific = FALSE))
  e   <- floor(log10(abs(v)))
  m   <- signif(v / 10^e, 3)
  sup <- chartr("-0123456789",
                "\u207b\u2070\u00b9\u00b2\u00b3\u2074\u2075\u2076\u2077\u2078\u2079",
                as.character(e))
  paste0(m, " \u00d7 10", sup)
}

eq_label <- function(cf) {
  sgn <- function(v) if (v < 0) " \u2212 " else " + "
  paste0("y = ", fmt_num(cf[1]),
         sgn(cf[2]), fmt_num(abs(cf[2])), " x",
         sgn(cf[3]), fmt_num(abs(cf[3])), " x\u00b2")
}

fit_mixed_quad <- function(d, y_var) {
  d <- d %>%
    filter(!is.na(.data[[y_var]])) %>%
    mutate(y = .data[[y_var]], x100 = Total_Water / 100)   # x per 100 mm for stable fitting
  
  out <- list(n = nrow(d), k = n_distinct(d$Study_ID),
              loc = n_distinct(d$Location_ID), model = NULL)
  if (out$k < 3) return(out)                                # too few studies for a curve
  
  m_full <- lmer(y ~ x100 + I(x100^2) + (1 | Study_ID/Location_ID), data = d, REML = FALSE)
  m_null <- lmer(y ~ 1                + (1 | Study_ID/Location_ID), data = d, REML = FALSE)
  
  if (isSingular(m_full)) {                       # Location variance ~ 0: drop that level
    m_full <- lmer(y ~ x100 + I(x100^2) + (1 | Study_ID), data = d, REML = FALSE)
    m_null <- lmer(y ~ 1                + (1 | Study_ID), data = d, REML = FALSE)
  }
  
  out$p   <- anova(m_null, m_full)$`Pr(>Chisq)`[2]          # likelihood-ratio test
  r2      <- tryCatch(r2_nakagawa(m_full), error = function(e) NULL)
  out$r2m <- get_r2(r2, "R2_marginal")
  out$r2c <- get_r2(r2, "R2_conditional")
  
  b <- fixef(m_full)
  out$coef_mm <- c(b[1], b[2] / 100, b[3] / 100^2)          # coefficients back in mm units
  
  grid <- data.frame(x100 = seq(min(d$x100), max(d$x100), length.out = 100))
  X    <- model.matrix(~ x100 + I(x100^2), grid)
  se   <- sqrt(diag(X %*% as.matrix(vcov(m_full)) %*% t(X)))
  grid$fit <- as.vector(X %*% b)
  grid$lwr <- grid$fit - 1.96 * se
  grid$lwr <- pmax(grid$lwr, 0)
  grid$upr <- grid$fit + 1.96 * se
  grid$Total_Water <- grid$x100 * 100
  
  out$pred  <- grid
  out$model <- m_full
  out
}

# --- STEP 3: Plot one outcome across the MAP groups ---
generate_total_water_reg_plot <- function(df, y_var, y_label, y_lims,
                                          stats_y_pos = "top",
                                          show_x_title = FALSE) {
  
  fits <- setNames(lapply(map_levels, function(g)
    fit_mixed_quad(filter(df, MAP_mm_grp == g), y_var)), map_levels)
  
  pred_df <- bind_rows(lapply(map_levels, function(g)
    if (!is.null(fits[[g]]$model)) mutate(fits[[g]]$pred, MAP_mm_grp = g)))
  
  lab_df <- bind_rows(lapply(map_levels, function(g) {
    f  <- fits[[g]]
    st <- ifelse(f$k   == 1, "study",    "studies")
    lc <- ifelse(f$loc == 1, "location", "locations")
    lab <- if (is.null(f$model)) {
      sprintf("n = %d (%d %s, %d %s)\nfew studies for a fitted curve",
              f$n, f$k, st, f$loc, lc)
    } else {
      sprintf("%s\nR\u00b2m = %s, R\u00b2c = %s, %s\nn = %d (%d %s, %d %s)",
              eq_label(f$coef_mm),
              ifelse(is.na(f$r2m), "\u2013", sprintf("%.2f", f$r2m)),
              ifelse(is.na(f$r2c), "\u2013", sprintf("%.2f", f$r2c)),
              ifelse(f$p < 0.001, "P < 0.001", sprintf("P = %.3f", f$p)),
              f$n, f$k, st, f$loc, lc)
    }
    data.frame(MAP_mm_grp = g, label = lab)
  }))
  
  pts <- df %>% filter(!is.na(.data[[y_var]]))
  pts$MAP_mm_grp     <- factor(pts$MAP_mm_grp,     levels = map_levels)
  pred_df$MAP_mm_grp <- factor(pred_df$MAP_mm_grp, levels = map_levels)
  lab_df$MAP_mm_grp  <- factor(lab_df$MAP_mm_grp,  levels = map_levels)
  
  is_top <- stats_y_pos == "top"
  
  p <- ggplot() +
    geom_point(data = pts, aes(x = Total_Water, y = .data[[y_var]]),
               alpha = 0.4, size = 1.5, shape = 21, fill = "gray80") +
    geom_ribbon(data = pred_df, aes(x = Total_Water, ymin = lwr, ymax = upr),
                fill = "gray70", alpha = 0.5) +
    geom_line(data = pred_df, aes(x = Total_Water, y = fit),
              color = "black", linewidth = 0.8) +
    geom_text(data = lab_df, aes(label = label),
              x = 1000, y = if (is_top) Inf else -Inf,
              hjust = 0.5, vjust = if (is_top) 1.2 else -0.3,
              size = 2.6, lineheight = 0.95, inherit.aes = FALSE) +
    coord_cartesian(xlim = c(0, 2000), ylim = y_lims) +
    facet_wrap(~ MAP_mm_grp) +
    theme_bw() +
    labs(y = y_label, x = if (show_x_title) "Irrigation water (mm)" else NULL) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "gray95"),
          strip.text   = element_text(face = "bold"),
          axis.title.y = element_text(face = "bold", size = 10),
          axis.title.x = element_text(face = "bold", size = 11))
  
  list(plot = p, fits = fits)
}

# --- STEP 4: Generate the three rows ---
res_yield <- generate_total_water_reg_plot(
  plot_data_total, "Yield_Value",
  expression(bold("Yield (kg ha"^-1*")")),
  c(0, 35000), "top")

res_cp <- generate_total_water_reg_plot(
  plot_data_total, "CP_Value",
  "Crude protein (%)",
  c(0, 30), "bottom")

res_wpi <- generate_total_water_reg_plot(
  plot_data_total, "WPI_Value",
  expression(bold("Irrigation water productivity (kg m"^-3*")")),
  c(0, 10), "top", show_x_title = TRUE)

p1_yield <- res_yield$plot +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
p2_cp    <- res_cp$plot +
  theme(strip.text = element_blank(), strip.background = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank())
p3_WPI   <- res_wpi$plot +
  theme(strip.text = element_blank(), strip.background = element_blank())

sapply(list(Yield = res_yield, CP = res_cp, WPI = res_wpi), function(r)
  sapply(r$fits, function(f) if (is.null(f$model)) NA else
    ifelse(length(VarCorr(f$model)) == 1, "Study only", "Study/Location")))


# --- STEP 5: Combine ---
final_grid <- (p1_yield / p2_cp / p3_WPI) +
  plot_annotation(tag_levels = "a", tag_suffix = ")") &
  theme(plot.tag    = element_text(face = "bold"),
        plot.margin = margin(t = 5, r = 5, b = 5, l = 25))

print(final_grid)

# --- STEP 6: Results table (for Supplementary / reviewer reply) ---
fig6_table <- bind_rows(lapply(list(Yield = res_yield, CP = res_cp, WPI = res_wpi), function(r)
  bind_rows(lapply(names(r$fits), function(g) {
    f <- r$fits[[g]]
    data.frame(MAP_group = g, n = f$n, studies = f$k, locations = f$loc,
               b0 = if (is.null(f$model)) NA else signif(f$coef_mm[1], 4),
               b1 = if (is.null(f$model)) NA else signif(f$coef_mm[2], 4),
               b2 = if (is.null(f$model)) NA else signif(f$coef_mm[3], 4),
               R2_marginal    = if (is.null(f$model)) NA else round(f$r2m, 2),
               R2_conditional = if (is.null(f$model)) NA else round(f$r2c, 2),
               P_LRT          = if (is.null(f$model)) NA else signif(f$p, 3))
  }))), .id = "Outcome")
print(fig6_table)

ggsave("Regression_WPI.png", final_grid, width = 15, height = 10, units = "in", dpi = 300)

write.csv(fig6_table, "Fig6_mixed_model_results.csv", row.names = FALSE)



# =====================================================================
# FIGURE 6 – Mixed-effects quadratic regressions (accounts for study clustering)
# Main figure: <250 and 250–500 mm | Supplementary figure: >500 mm
# =====================================================================
library(lme4)
library(performance)
library(dplyr)
library(ggplot2)
library(patchwork)

# --- STEP 1: Long data – one row per unique treatment or control observation ---
trt_obs <- meta_es %>%
  transmute(Study_ID, Location_ID, Year, MAP_mm_grp,
            Total_Water = Irri_amt_mm_trt,
            Yield_Value = Yield_kg_trt,
            CP_Value    = CP_trt,
            WPI_Value   = WPI_kg.m3_trt)

ctrl_obs <- meta_es %>%
  transmute(Study_ID, Location_ID, Year, MAP_mm_grp,
            Total_Water = Irri_amt_mm_Ctrl,
            Yield_Value = Yield_kg_ctrl,
            CP_Value    = CP_ctrl,
            WPI_Value   = WPI_kg.m3_ctrl)

plot_data_total <- bind_rows(trt_obs, ctrl_obs) %>%
  distinct() %>%                                   # shared controls counted only once
  filter(!is.na(Total_Water)) %>%
  mutate(Study_ID    = factor(Study_ID),
         Location_ID = factor(Location_ID))

map_levels  <- levels(factor(meta_es$MAP_mm_grp))   # all groups, in order
main_groups <- map_levels[1:2]                      # <250 mm, 250–500 mm
supp_groups <- map_levels[3]                        # >500 mm

# --- STEP 2: Fit one mixed quadratic model (one outcome, one MAP group) ---
get_r2 <- function(r, nm) {
  v <- tryCatch(as.numeric(r[[nm]]), error = function(e) NA)
  if (length(v) == 0) NA else v
}
fmt_r2 <- function(v) {
  if (is.na(v)) return("\u2013")
  if (v < 0.005) return("< 0.01")
  sprintf("%.2f", v)
}

fmt_num <- function(v) {
  if (v == 0 || abs(v) >= 0.001) return(format(signif(v, 3), trim = TRUE, scientific = FALSE))
  e   <- floor(log10(abs(v)))
  m   <- signif(v / 10^e, 3)
  sup <- chartr("-0123456789",
                "\u207b\u2070\u00b9\u00b2\u00b3\u2074\u2075\u2076\u2077\u2078\u2079",
                as.character(e))
  paste0(m, " \u00d7 10", sup)
}

eq_label <- function(cf) {
  sgn <- function(v) if (v < 0) " \u2212 " else " + "
  paste0("y = ", fmt_num(cf[1]),
         sgn(cf[2]), fmt_num(abs(cf[2])), " x",
         sgn(cf[3]), fmt_num(abs(cf[3])), " x\u00b2")
}

fit_mixed_quad <- function(d, y_var) {
  d <- d %>%
    filter(!is.na(.data[[y_var]])) %>%
    mutate(y = .data[[y_var]], x100 = Total_Water / 100)   # x per 100 mm for stable fitting
  
  out <- list(n = nrow(d), k = n_distinct(d$Study_ID),
              loc = n_distinct(d$Location_ID), model = NULL)
  if (out$k < 3) return(out)                                # too few studies for a curve
  
  m_full <- lmer(y ~ x100 + I(x100^2) + (1 | Study_ID/Location_ID), data = d, REML = FALSE)
  m_null <- lmer(y ~ 1                + (1 | Study_ID/Location_ID), data = d, REML = FALSE)
  
  if (isSingular(m_full)) {                       # Location variance ~ 0: drop that level
    m_full <- lmer(y ~ x100 + I(x100^2) + (1 | Study_ID), data = d, REML = FALSE)
    m_null <- lmer(y ~ 1                + (1 | Study_ID), data = d, REML = FALSE)
  }
  
  out$p   <- anova(m_null, m_full)$`Pr(>Chisq)`[2]          # likelihood-ratio test
  r2      <- tryCatch(r2_nakagawa(m_full), error = function(e) NULL)
  out$r2m <- get_r2(r2, "R2_marginal")
  out$r2c <- get_r2(r2, "R2_conditional")
  
  b <- fixef(m_full)
  out$coef_mm <- c(b[1], b[2] / 100, b[3] / 100^2)          # coefficients back in mm units
  
  grid <- data.frame(x100 = seq(min(d$x100), max(d$x100), length.out = 100))
  X    <- model.matrix(~ x100 + I(x100^2), grid)
  se   <- sqrt(diag(X %*% as.matrix(vcov(m_full)) %*% t(X)))
  grid$fit <- as.vector(X %*% b)
  grid$lwr <- pmax(grid$fit - 1.96 * se, 0)
  grid$upr <- grid$fit + 1.96 * se
  grid$Total_Water <- grid$x100 * 100
  
  out$pred  <- grid
  out$model <- m_full
  out
}

# --- STEP 3: Plot one outcome for the chosen MAP groups ---
generate_total_water_reg_plot <- function(df, y_var, y_label, y_lims,
                                          groups,                    # which MAP groups to draw
                                          x_lims = c(0, 2000),
                                          stats_y_pos = "top",
                                          show_x_title = FALSE) {
  
  df <- df %>% filter(MAP_mm_grp %in% groups)
  
  fits <- setNames(lapply(groups, function(g)
    fit_mixed_quad(filter(df, MAP_mm_grp == g), y_var)), groups)
  
  pred_df <- bind_rows(lapply(groups, function(g)
    if (!is.null(fits[[g]]$model)) mutate(fits[[g]]$pred, MAP_mm_grp = g)))
  
  lab_df <- bind_rows(lapply(groups, function(g) {
    f  <- fits[[g]]
    st <- ifelse(f$k   == 1, "study",    "studies")
    lc <- ifelse(f$loc == 1, "location", "locations")
    lab <- if (is.null(f$model)) {
      sprintf("n = %d (%d %s, %d %s)", f$n, f$k, st, f$loc, lc)
    } else {
      sprintf("%s\nR\u00b2m = %s, R\u00b2c = %s, %s\nn = %d (%d %s, %d %s)",
              eq_label(f$coef_mm),
              fmt_r2(f$r2m),
              fmt_r2(f$r2c),
              ifelse(f$p < 0.001, "P < 0.001", sprintf("P = %.3f", f$p)),
              f$n, f$k, st, f$loc, lc)
    }
    data.frame(MAP_mm_grp = g, label = lab)
  }))
  
  pts <- df %>% filter(!is.na(.data[[y_var]]))
  pts$MAP_mm_grp    <- factor(pts$MAP_mm_grp,    levels = groups)
  lab_df$MAP_mm_grp <- factor(lab_df$MAP_mm_grp, levels = groups)
  
  is_top <- stats_y_pos == "top"
  
  p <- ggplot() +
    geom_point(data = pts, aes(x = Total_Water, y = .data[[y_var]]),
               alpha = 0.4, size = 1.5, shape = 21, fill = "gray80")
  
  if (nrow(pred_df) > 0) {                         # only draw curves where a model was fitted
    pred_df$MAP_mm_grp <- factor(pred_df$MAP_mm_grp, levels = groups)
    p <- p +
      geom_ribbon(data = pred_df, aes(x = Total_Water, ymin = lwr, ymax = upr),
                  fill = "gray70", alpha = 0.5) +
      geom_line(data = pred_df, aes(x = Total_Water, y = fit),
                color = "black", linewidth = 0.8)
  }
  
  p <- p +
    geom_text(data = lab_df, aes(label = label),
              x = mean(x_lims), y = if (is_top) Inf else -Inf,
              hjust = 0.5, vjust = if (is_top) 1.2 else -0.3,
              size = 2.6, lineheight = 0.95, inherit.aes = FALSE) +
    coord_cartesian(xlim = x_lims, ylim = y_lims) +
    facet_wrap(~ MAP_mm_grp) +
    theme_bw() +
    labs(y = y_label, x = if (show_x_title) "Irrigation water (mm)" else NULL) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "gray95"),
          strip.text   = element_text(face = "bold"),
          axis.title.y = element_text(face = "bold", size = 10),
          axis.title.x = element_text(face = "bold", size = 11))
  
  list(plot = p, fits = fits)
}

# --- STEP 4: Build a 3-row figure (yield, CP, WPI) for a set of groups ---
make_fig6 <- function(groups, x_lims = c(0, 2000)) {
  res_yield <- generate_total_water_reg_plot(
    plot_data_total, "Yield_Value", expression(bold("Yield (kg ha"^-1*")")),
    c(0, 35000), groups, x_lims, "top")
  res_cp <- generate_total_water_reg_plot(
    plot_data_total, "CP_Value", "Crude protein (%)",
    c(0, 30), groups, x_lims, "bottom")
  res_wpi <- generate_total_water_reg_plot(
    plot_data_total, "WPI_Value", expression(bold("Irrigation water productivity (kg m"^-3*")")),
    c(0, 10), groups, x_lims, "top", show_x_title = TRUE)
  
  p1 <- res_yield$plot + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  p2 <- res_cp$plot + theme(strip.text = element_blank(), strip.background = element_blank(),
                            axis.text.x = element_blank(), axis.ticks.x = element_blank())
  p3 <- res_wpi$plot + theme(strip.text = element_blank(), strip.background = element_blank())
  
  grid <- (p1 / p2 / p3) +
    plot_annotation(tag_levels = "a", tag_suffix = ")") &
    theme(plot.tag    = element_text(face = "bold"),
          plot.margin = margin(t = 5, r = 5, b = 5, l = 25))
  
  list(grid = grid, res = list(Yield = res_yield, CP = res_cp, WPI = res_wpi))
}

# --- STEP 5a: Main figure (<250 and 250–500 mm) ---
fig6_main <- make_fig6(main_groups, x_lims = c(0, 2000))
print(fig6_main$grid)
ggsave("Fig6_main.png", fig6_main$grid, width = 9, height = 9, units = "in", dpi = 300)

# --- STEP 5b: Supplementary figure (>500 mm) ---
fig6_supp <- make_fig6(supp_groups, x_lims = c(0, 400))      # data span ~0–300 mm
print(fig6_supp$grid)
ggsave("Fig6_supplementary.png", fig6_supp$grid, width = 5, height = 9, units = "in", dpi = 300)

# Which random structure was used in each fitted model
sapply(fig6_main$res, function(r)
  sapply(r$fits, function(f) if (is.null(f$model)) NA else
    ifelse(length(VarCorr(f$model)) == 1, "Study only", "Study/Location")))

# --- STEP 6: Results table (main-figure groups; supplementary has no fitted models) ---
fig6_table <- bind_rows(lapply(c(fig6_main$res), function(r)
  bind_rows(lapply(names(r$fits), function(g) {
    f <- r$fits[[g]]
    data.frame(MAP_group = g, n = f$n, studies = f$k, locations = f$loc,
               b0 = if (is.null(f$model)) NA else signif(f$coef_mm[1], 4),
               b1 = if (is.null(f$model)) NA else signif(f$coef_mm[2], 4),
               b2 = if (is.null(f$model)) NA else signif(f$coef_mm[3], 4),
               R2_marginal    = if (is.null(f$model)) NA else round(f$r2m, 2),
               R2_conditional = if (is.null(f$model)) NA else round(f$r2c, 2),
               P_LRT          = if (is.null(f$model)) NA else signif(f$p, 3))
  }))), .id = "Outcome")
print(fig6_table)

##################################### GEtting thresholds ######################################################

library(metafor); library(dplyr)

meta_es <- meta_es %>% mutate(rel_irr = Irri_amt_mm_trt / Irri_amt_mm_Ctrl)
zones    <- levels(factor(meta_es$MAP_mm_grp))[1:2]
deficits <- c(0.10, 0.20, 0.30, 0.40, 0.50)      # water deficit relative to own full irrigation

pred_quad <- function(zone) {
  d <- meta_es %>%
    filter(MAP_mm_grp == zone, !is.na(Yield_LnR), !is.na(Yield_Variance),
           is.finite(rel_irr), rel_irr <= 1) %>%
    mutate(dr = rel_irr - 1)
  # Quadratic through relative yield = 1 at full irrigation (no intercept)
  m <- rma.mv(Yield_LnR, Yield_Variance, mods = ~ dr + I(dr^2) - 1,
              random = ~ 1 | Study_ID/Location_ID, data = d)
  nd <- cbind(-deficits, deficits^2)
  p  <- predict(m, newmods = nd)
  data.frame(Zone = zone,
             water_deficit   = paste0(deficits * 100, "%"),
             yield_retained  = round(exp(p$pred)  * 100, 1),
             CI_low          = round(exp(p$ci.lb) * 100, 1),
             CI_up           = round(exp(p$ci.ub) * 100, 1),
             WPI_change_pct  = round((exp(p$pred) / (1 - deficits) - 1) * 100, 1),
             quad_term_P     = signif(m$pval[2], 3),
             k = m$k, studies = n_distinct(d$Study_ID))
}

deficit_table <- bind_rows(lapply(zones, pred_quad))
print(deficit_table)

##############################################################################################

library(pbkrtest)

kr_test <- function(f) {
  if (is.null(f$model)) return(c(F = NA, ddf = NA, P = NA))
  mf <- model.frame(f$model)                       # data the model was fitted on
  m1 <- update(f$model, REML = TRUE, data = mf)
  m0 <- update(m1, . ~ . - x100 - I(x100^2))
  s  <- KRmodcomp(m1, m0)$stats
  c(F = round(s$Fstat, 2), ddf = round(s$ddf, 1), P = signif(s$p.value, 3))
}

kr_table <- bind_rows(
  lapply(list(Yield = res_yield, CP = res_cp, WPI = res_wpi), function(r)
    bind_rows(lapply(names(r$fits), function(g)
      data.frame(MAP_group = g, t(kr_test(r$fits[[g]])))))),
  .id = "Outcome")

print(kr_table)


balance <- thr_c %>%
  mutate(r_star         = 1 / (2 * slope),                                   # optimum fraction of full irrigation
         irr_mm_balance = round(r_star * avg_full_irr_mm),
         rel_yield_pct  = round(exp(slope * (r_star - 1)) * 100, 1),         # yield as % of full-irrigation yield
         WPI_change_pct = round((exp(slope * (r_star - 1)) / r_star - 1) * 100, 1),
         r_star         = round(r_star, 2)) %>%
  select(Zone, slope, r_star, irr_mm_balance, rel_yield_pct, WPI_change_pct)

print(balance)

################################################################################

# Calculate average baseline values for Full Irrigation (Control)
baseline_stats <- meta_es %>%
  distinct(Study_ID, Location_ID, Year, MAP_mm_grp,
           Irri_amt_mm_Ctrl, Yield_kg_ctrl, WPI_kg.m3_ctrl, CP_ctrl) %>%   # each control once
  group_by(MAP_mm_grp) %>%
  summarize(
    avg_yield_ctrl = mean(Yield_kg_ctrl, na.rm = TRUE),
    avg_water_ctrl = mean(Irri_amt_mm_Ctrl, na.rm = TRUE),
    avg_wpi_ctrl   = mean(WPI_kg.m3_ctrl, na.rm = TRUE),
    avg_cp_ctrl    = mean(CP_ctrl, na.rm = TRUE),
    n_studies      = n_distinct(Study_ID),
    n_controls     = n()
  ) %>%
  arrange(MAP_mm_grp)

print(baseline_stats)

#################################################################################

meta_es <- meta_es %>%
  mutate(rel_irr = Irri_amt_mm_trt / Irri_amt_mm_Ctrl)

meta_es %>%
  filter(is.finite(rel_irr), rel_irr <= 1,
         MAP_mm_grp %in% levels(factor(MAP_mm_grp))[1:2]) %>%
  mutate(band = cut(rel_irr, c(-Inf, 0.4, 0.6, 0.85, 1),
                    labels = c("<40%", "40–60%", "60–85%", "85–100%"))) %>%
  count(MAP_mm_grp, band)


# =====================================================================
# YIELD AND WPI AT REDUCED IRRIGATION LEVELS + BALANCE POINT
# (within-study meta-regression, arid and semi-arid zones)
# =====================================================================
library(metafor)
library(dplyr)

# --- 1. Relative irrigation (treatment / full-irrigation control) ---
meta_es <- meta_es %>%
  mutate(rel_irr = Irri_amt_mm_trt / Irri_amt_mm_Ctrl)

zones    <- levels(factor(meta_es$MAP_mm_grp))[1:2]   # <250 mm, 250–500 mm
levels_r <- c(0.75, 0.50)                              # water levels to predict

# --- 2. Average full irrigation per zone (each control counted once) ---
full_irr <- meta_es %>%
  distinct(Study_ID, Location_ID, Year, MAP_mm_grp,
           Irri_amt_mm_Ctrl, Yield_kg_ctrl, WPI_kg.m3_ctrl, CP_ctrl) %>%
  group_by(MAP_mm_grp) %>%
  summarise(full_irr_mm = mean(Irri_amt_mm_Ctrl, na.rm = TRUE), .groups = "drop") %>%
  mutate(Zone = as.character(MAP_mm_grp)) %>%
  select(Zone, full_irr_mm)

# --- 3. Fit model and predict for one zone ---
pred_zone <- function(zone) {
  d <- meta_es %>%
    filter(MAP_mm_grp == zone, !is.na(Yield_LnR), !is.na(Yield_Variance),
           is.finite(rel_irr), rel_irr <= 1)                 # deficit treatments only
  
  # Relative yield constrained to 1 at full irrigation
  m <- rma.mv(Yield_LnR, Yield_Variance, mods = ~ I(rel_irr - 1) - 1,
              random = ~ 1 | Study_ID/Location_ID, data = d)
  b <- as.numeric(coef(m))
  
  # Predictions at fixed water levels
  p <- predict(m, newmods = levels_r - 1)
  pred_tab <- data.frame(
    Zone         = zone,
    level        = paste0(round(levels_r * 100), "% of full"),
    rel_irr      = levels_r,
    yield_pct    = round(exp(p$pred)  * 100, 1),
    yield_CI_low = round(exp(p$ci.lb) * 100, 1),
    yield_CI_up  = round(exp(p$ci.ub) * 100, 1),
    WPI_change   = round((exp(p$pred) / levels_r - 1) * 100, 1))
  
  # Balance point: maximizes relative yield x relative WPI (r* = 1 / 2b)
  r_star <- 1 / (2 * b)
  pb <- predict(m, newmods = r_star - 1)
  bal_tab <- data.frame(
    Zone         = zone,
    level        = "Balance point",
    rel_irr      = round(r_star, 2),
    yield_pct    = round(exp(pb$pred)  * 100, 1),
    yield_CI_low = round(exp(pb$ci.lb) * 100, 1),
    yield_CI_up  = round(exp(pb$ci.ub) * 100, 1),
    WPI_change   = round((exp(pb$pred) / r_star - 1) * 100, 1))
  
  out <- bind_rows(pred_tab, bal_tab)
  out$slope   <- round(b, 3)
  out$slope_P <- signif(m$pval, 3)
  out$k       <- m$k
  out$studies <- n_distinct(d$Study_ID)
  out
}

# --- 4. Run for both zones and add irrigation in mm ---
water_levels <- bind_rows(lapply(zones, pred_zone)) %>%
  left_join(full_irr, by = "Zone") %>%
  mutate(irr_mm      = round(rel_irr * full_irr_mm),
         full_irr_mm = round(full_irr_mm)) %>%
  select(Zone, level, rel_irr, irr_mm, full_irr_mm,
         yield_pct, yield_CI_low, yield_CI_up, WPI_change,
         slope, slope_P, k, studies)

print(water_levels)

# --- 5. Save for the manuscript / supplement ---
write.csv(water_levels, "Irrigation_levels_yield_WPI.csv", row.names = FALSE)





###################################################################################

class(meta_es$Degree_of_irri_deficit)          # should be "numeric"
summary(meta_es$Degree_of_irri_deficit)        # values like 0–1 (or 0–100 if in %)
table(meta_es$Degree_of_irri_deficit > 1, useNA = "ifany")

# Compare with the calculated ratio
meta_es %>%
  select(Study_ID, Degree_of_irri_deficit, rel_irr) %>%
  filter(abs(Degree_of_irri_deficit - rel_irr) > 0.1) %>%
  print(n = 30)


meta_es %>%
  filter(is.na(Degree_of_irri_deficit)) %>%
  summarise(n_NA          = n(),
            rainfed_trt   = sum(Irri_amt_mm_trt == 0, na.rm = TRUE),
            missing_irr   = sum(is.na(Irri_amt_mm_trt) | is.na(Irri_amt_mm_Ctrl)),
            other         = n() - rainfed_trt - missing_irr)


ctrl <- list(optimizer = "optim", optmethod = "Nelder-Mead")

# 1. Relative irrigation (includes rainfed treatments as 0)
meta_es <- meta_es %>%
  mutate(rel_irr = Irri_amt_mm_trt / Irri_amt_mm_Ctrl)

# 2. Sensitivity: all rows vs deficit only (excluding over-irrigation, rel_irr > 1)
outcomes <- list(Yield = c("Yield_LnR", "Yield_Variance"),
                 CP    = c("CP_LnR",    "CP_Variance"),
                 WPI   = c("WPI_LnR",   "WPI_Variance"),
                 TWP   = c("TWP_LnR",   "TWP_Variance"))

sens <- bind_rows(lapply(names(outcomes), function(nm) {
  yi <- outcomes[[nm]][1]; vi <- outcomes[[nm]][2]
  fit <- function(d) rma.mv(d[[yi]], d[[vi]], random = ~ 1 | Study_ID/Location_ID,
                            data = d, control = ctrl)
  d_all <- meta_es %>% filter(!is.na(.data[[yi]]), !is.na(.data[[vi]]))
  d_def <- d_all %>% filter(is.na(rel_irr) | rel_irr <= 1)
  m_all <- fit(d_all); m_def <- fit(d_def)
  data.frame(Outcome = nm,
             k_all = m_all$k, pct_all = round((exp(coef(m_all)) - 1) * 100, 1), P_all = signif(m_all$pval, 3),
             k_def = m_def$k, pct_def = round((exp(coef(m_def)) - 1) * 100, 1), P_def = signif(m_def$pval, 3))
}))
print(sens)

# 3. Threshold: % of full irrigation giving 90% of full-irrigation yield (deficit rows only)
zones <- levels(factor(meta_es$MAP_mm_grp))[1:2]

thr_zone_c <- function(zone, target = 0.90) {
  d <- meta_es %>%
    filter(MAP_mm_grp == zone, !is.na(Yield_LnR), !is.na(Yield_Variance),
           is.finite(rel_irr), rel_irr <= 1)
  m <- rma.mv(Yield_LnR, Yield_Variance, mods = ~ I(rel_irr - 1) - 1,
              random = ~ 1 | Study_ID/Location_ID, data = d, control = ctrl)
  b <- coef(m)
  data.frame(Zone = zone, n = m$k, studies = n_distinct(d$Study_ID),
             slope = round(b, 3), slope_P = signif(m$pval, 3),
             min_rel = round(min(d$rel_irr), 2), max_rel = round(max(d$rel_irr), 2),
             rel_irr_for_90pct_yield = round(1 + log(target) / b, 2))
}

full_irr <- meta_es %>%
  distinct(Study_ID, Location_ID, Year, MAP_mm_grp, Irri_amt_mm_Ctrl) %>%
  group_by(MAP_mm_grp) %>%
  summarise(avg_full_irr_mm = mean(Irri_amt_mm_Ctrl, na.rm = TRUE), .groups = "drop") %>%
  mutate(Zone = as.character(MAP_mm_grp)) %>% select(Zone, avg_full_irr_mm)

thr_c <- bind_rows(lapply(zones, thr_zone_c)) %>%
  left_join(full_irr, by = "Zone") %>%
  mutate(irr_mm_for_90pct_yield = round(rel_irr_for_90pct_yield * avg_full_irr_mm))
print(thr_c)
# =====================================================================
# Within-study irrigation thresholds (yield) + CP check
# =====================================================================
library(metafor)
library(dplyr)

ctrl <- list(optimizer = "optim", optmethod = "Nelder-Mead")

# Relative irrigation: treatment irrigation as a fraction of its own control (full) irrigation
meta_es <- meta_es %>%
  mutate(rel_irr = Irri_amt_mm_trt / Irri_amt_mm_Ctrl)

zones <- levels(factor(meta_es$MAP_mm_grp))[1:2]   # arid and semi-arid only

# --- 1. Yield threshold: fraction of full irrigation giving 90% of full-irrigation yield ---
thr_zone <- function(zone, target = 0.90) {
  d <- meta_es %>%
    filter(MAP_mm_grp == zone, !is.na(Yield_LnR), !is.na(Yield_Variance),
           is.finite(rel_irr))
  m <- rma.mv(Yield_LnR, Yield_Variance, mods = ~ rel_irr,
              random = ~ 1 | Study_ID/Location_ID, data = d, control = ctrl)
  b <- coef(m)
  data.frame(Zone = zone,
             n = m$k,
             studies = n_distinct(d$Study_ID),
             intercept = round(b[1], 3),
             slope = round(b[2], 3),
             slope_P = signif(m$pval[2], 3),
             rel_irr_for_90pct_yield = round((log(target) - b[1]) / b[2], 2))
}

thr_within <- bind_rows(lapply(zones, thr_zone))

# Convert the fraction into mm using average full irrigation of each zone
full_irr <- meta_es %>%
  distinct(Study_ID, Location_ID, Year, MAP_mm_grp, Irri_amt_mm_Ctrl) %>%
  group_by(MAP_mm_grp) %>%
  summarise(avg_full_irr_mm = mean(Irri_amt_mm_Ctrl, na.rm = TRUE), .groups = "drop") %>%
  mutate(Zone = as.character(MAP_mm_grp)) %>%
  select(Zone, avg_full_irr_mm)

thr_within <- thr_within %>%
  left_join(full_irr, by = "Zone") %>%
  mutate(irr_mm_for_90pct_yield = round(rel_irr_for_90pct_yield * avg_full_irr_mm))

print(thr_within)

# --- 2. CP check: does CP change with the degree of irrigation deficit? ---
cp_zone <- function(zone) {
  d <- meta_es %>%
    filter(MAP_mm_grp == zone, !is.na(CP_LnR), !is.na(CP_Variance), is.finite(rel_irr))
  m <- rma.mv(CP_LnR, CP_Variance, mods = ~ rel_irr,
              random = ~ 1 | Study_ID/Location_ID, data = d, control = ctrl)
  data.frame(Zone = zone, n = m$k, studies = n_distinct(d$Study_ID),
             slope = round(coef(m)[2], 4), slope_P = signif(m$pval[2], 3))
}

cp_within <- bind_rows(lapply(zones, cp_zone))
print(cp_within)

# Quick look at the range of relative irrigation in each zone
meta_es %>%
  filter(is.finite(rel_irr), MAP_mm_grp %in% zones) %>%
  group_by(MAP_mm_grp) %>%
  summarise(min_rel = round(min(rel_irr), 2), max_rel = round(max(rel_irr), 2), n = n())


# Which studies have treatment irrigation > control irrigation?
meta_es %>% filter(rel_irr > 1) %>% count(MAP_mm_grp, Study_ID)

meta_es %>%
  filter(rel_irr > 1) %>%
  select(Study_ID, Location_ID, Year, Trt_type, Control_type,
         Irri_amt_mm_trt, Irri_amt_mm_Ctrl, Yield_kg_trt, Yield_kg_ctrl, Yield_LnR) %>%
  print(n = 40)

##########################################################################################################################

# # --- Sensitivity analysis Setup ---
# outcomes <- c("Yield", "CP", "WPI", "TWP") 
# model_list_overall <- list()
# for (out in outcomes) {
#   model_list_overall[[out]] <- rma.mv(
#     yi = meta_es[[paste0(out, "_LnR")]],
#     V = meta_es[[paste0(out, "_Variance")]],
#     random = ~ 1 | Study_ID / Location_ID,
#     data = meta_es
#   )
# }
# 
# 
# # Study-Level Jackknife Function (MODIFIED for Percentage Change)
# run_jackknife_study_percent <- function(model, outcome_name, data_full) {
#   unique_study_ids <- unique(model$mf.r[]) %>% unlist() %>% as.character()
#   num_studies <- length(unique_study_ids)
#   
#   # Initialize variables in LnR format first
#   estimates_lnr <- rep(NA, num_studies)
#   ci_lower_lnr  <- rep(NA, num_studies)
#   ci_upper_lnr  <- rep(NA, num_studies)
#   
#   for (i in 1:num_studies) {
#     current_id_to_remove <- unique_study_ids[i]
#     
#     subset_data <- data_full[data_full$Study_ID != current_id_to_remove, ]
#     
#     # Run the model from scratch on the new subset_data 
#     res.loo <- rma.mv(
#       yi = subset_data[[paste0(outcome_name, "_LnR")]],
#       V = subset_data[[paste0(outcome_name, "_Variance")]],
#       random = ~ 1 | Location_ID, 
#       data = subset_data
#     )
#     
#     estimates_lnr[i] <- coef(res.loo)
#     ci_lower_lnr[i]  <- res.loo$ci.lb
#     ci_upper_lnr[i]  <- res.loo$ci.ub
#   }
#   
#   # --- Conversion to Percentage Change ---
#   # Apply the formula: (exp(LnR) - 1) * 100
#   estimates_percent <- (exp(estimates_lnr) - 1) * 100
#   ci_lower_percent  <- (exp(ci_lower_lnr) - 1) * 100
#   ci_upper_percent  <- (exp(ci_upper_lnr) - 1) * 100
#   
#   data.frame(
#     Outcome = outcome_name,
#     Study_Removed = unique_study_ids,
#     Estimate_Adj  = estimates_percent, # Renaming back to Estimate_Adj for plot function compatibility
#     CI_Lower_Adj  = ci_lower_percent,
#     CI_Upper_Adj  = ci_upper_percent
#   )
# }
# 
# 
# # Run the study-level jackknife for all specified outcomes and combine (using the new function)
# all_jackknife_results_study <- map_dfr(outcomes, ~run_jackknife_study_percent(model_list_overall[[.x]], .x, meta_es))
# 
# 
# # Jackknife Plotting Function
# plot_jackknife_results_grid <- function(df, outcome_name, x_label) {
#   df_sub <- df %>% filter(Outcome == outcome_name)
#   
#   # Calculate original overall mean and 95% CIs (from your existing list)
#   original_mean_lnr <- coef(model_list_overall[[outcome_name]]) 
#   original_mean <- (exp(original_mean_lnr) - 1) * 100
#   original_ci_lb <- (exp(model_list_overall[[outcome_name]]$ci.lb) - 1) * 100
#   original_ci_ub <- (exp(model_list_overall[[outcome_name]]$ci.ub) - 1) * 100
#   
#   ggplot(df_sub, aes(y = reorder(Study_Removed, Estimate_Adj), x = Estimate_Adj)) +
#     geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.5) + 
#     # Overall mean reference lines
#     geom_vline(xintercept = original_mean, linetype = "solid", color = "red", alpha = 0.6) + 
#     geom_vline(xintercept = c(original_ci_lb, original_ci_ub), linetype = "dashed", color = "red", alpha = 0.5) +
#     # Jackknife points
#     geom_errorbarh(aes(xmin = CI_Lower_Adj, xmax = CI_Upper_Adj), height = 0, color = "gray30") + 
#     geom_point(size = 1.2) +
#     theme_bw() +
#     labs(y = "", x = x_label) + 
#     xlim(-40, 60) + 
#     theme(
#       panel.grid.major = element_blank(), # Removes major grid lines
#       panel.grid.minor = element_blank(), # Removes minor grid lines
#       axis.text.y = element_text(size = 6),
#       panel.grid.major.y = element_blank(),
#       axis.title.x = element_text(size = 9, face = "bold")
#     )
# }
# 
# 
# # Generate individual panels
# pj1 <- plot_jackknife_results_grid(all_jackknife_results_study, "Yield", "Percentage change of yield")
# pj2 <- plot_jackknife_results_grid(all_jackknife_results_study, "CP",    "Percentage change of crude protein")
# pj3 <- plot_jackknife_results_grid(all_jackknife_results_study, "WPI",  "Percentage change of irrigation water productivity")
# pj4 <- plot_jackknife_results_grid(all_jackknife_results_study, "TWP",   "Percentage change of total water productivity")
# 
# # --- 2. Create the 2x2 Layout ---
# # We use the original pj1, pj2, pj3, pj4 (none are "cleaned")
# jackknife_grid <- (pj1 | pj2) / (pj3 | pj4) + 
#   plot_annotation(
#     tag_levels = 'a',
#     tag_suffix = ')'
#   ) +
#   # widths = c(1, 1) ensures equal spacing for the Y-axis labels in both columns
#   plot_layout(guides = 'collect', widths = c(1, 1)) &
#   theme(
#     plot.tag.position = c(0.03, 1.02),
#     plot.tag = element_text(size = 11, face = "bold"),
#     # Ensure Y axis text is visible for all
#     axis.text.y = element_text(size = 6) 
#   )
# 
# # --- 3. Add common Y-axis label using cowplot ---
# final_sensitivity_plot <- cowplot::ggdraw(jackknife_grid) +
#   cowplot::draw_label("Study ID removed", x = 0.005, y = 0.5, angle = 90, 
#                       vjust = 0.5, size = 11, fontface = "bold")
# 
# print(final_sensitivity_plot)


# =====================================================================
# SENSITIVITY ANALYSIS – Leave-one-study-out (jackknife)
# =====================================================================
library(metafor)
library(dplyr)
library(purrr)
library(ggplot2)
library(patchwork)
library(cowplot)

ctrl <- list(optimizer = "optim", optmethod = "Nelder-Mead")

# --- 1. Overall models (same structure as main analysis) ---
outcomes <- c("Yield", "CP", "WPI", "TWP")
model_list_overall <- list()
for (out in outcomes) {
  model_list_overall[[out]] <- rma.mv(
    yi      = meta_es[[paste0(out, "_LnR")]],
    V       = meta_es[[paste0(out, "_Variance")]],
    random  = ~ 1 | Study_ID / Location_ID,
    data    = meta_es,
    control = ctrl
  )
}

# --- 2. Study-level jackknife (leave one STUDY out at a time) ---
run_jackknife_study_percent <- function(outcome_name, data_full) {
  yi_col <- paste0(outcome_name, "_LnR")
  vi_col <- paste0(outcome_name, "_Variance")
  
  # Only the studies that contribute to this outcome
  used <- !is.na(data_full[[yi_col]]) & !is.na(data_full[[vi_col]])
  unique_study_ids <- sort(unique(as.character(data_full$Study_ID[used])))
  
  bind_rows(lapply(unique_study_ids, function(id) {
    subset_data <- data_full[as.character(data_full$Study_ID) != id, ]
    res.loo <- rma.mv(
      yi      = subset_data[[yi_col]],
      V       = subset_data[[vi_col]],
      random  = ~ 1 | Study_ID / Location_ID,   # same as main model
      data    = subset_data,
      control = ctrl
    )
    data.frame(Outcome       = outcome_name,
               Study_Removed = id,
               Estimate_Adj  = (exp(as.numeric(coef(res.loo))) - 1) * 100,
               CI_Lower_Adj  = (exp(res.loo$ci.lb) - 1) * 100,
               CI_Upper_Adj  = (exp(res.loo$ci.ub) - 1) * 100)
  }))
}

all_jackknife_results_study <- map_dfr(outcomes, ~ run_jackknife_study_percent(.x, meta_es))

# Check: number of studies per outcome (should be 16 or fewer)
print(all_jackknife_results_study %>% count(Outcome))

# --- 3. Plotting function ---
plot_jackknife_results_grid <- function(df, outcome_name, x_label) {
  df_sub <- df %>% filter(Outcome == outcome_name)
  
  m <- model_list_overall[[outcome_name]]
  original_mean  <- (exp(as.numeric(coef(m))) - 1) * 100
  original_ci_lb <- (exp(m$ci.lb) - 1) * 100
  original_ci_ub <- (exp(m$ci.ub) - 1) * 100
  
  ggplot(df_sub, aes(y = reorder(Study_Removed, Estimate_Adj), x = Estimate_Adj)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
    geom_vline(xintercept = original_mean, linetype = "solid", color = "red", alpha = 0.6) +
    geom_vline(xintercept = c(original_ci_lb, original_ci_ub),
               linetype = "dashed", color = "red", alpha = 0.5) +
    geom_errorbarh(aes(xmin = CI_Lower_Adj, xmax = CI_Upper_Adj), height = 0, color = "gray30") +
    geom_point(size = 1.2) +
    coord_cartesian(xlim = c(-40, 60)) +          # zooms without deleting error bars
    theme_bw() +
    labs(y = "", x = x_label) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.y  = element_text(size = 6),
          axis.title.x = element_text(size = 9, face = "bold"))
}

# --- 4. Panels ---
pj1 <- plot_jackknife_results_grid(all_jackknife_results_study, "Yield", "Percentage change of yield")
pj2 <- plot_jackknife_results_grid(all_jackknife_results_study, "CP",    "Percentage change of crude protein")
pj3 <- plot_jackknife_results_grid(all_jackknife_results_study, "WPI",   "Percentage change of irrigation water productivity")
pj4 <- plot_jackknife_results_grid(all_jackknife_results_study, "TWP",   "Percentage change of total water productivity")

# --- 5. 2x2 layout ---
jackknife_grid <- (pj1 | pj2) / (pj3 | pj4) +
  plot_annotation(tag_levels = "a", tag_suffix = ")") +
  plot_layout(widths = c(1, 1)) &
  theme(plot.tag.position = c(0.03, 1.02),
        plot.tag    = element_text(size = 11, face = "bold"),
        axis.text.y = element_text(size = 6))

# --- 6. Shared y-axis label ---
final_sensitivity_plot <- ggdraw(jackknife_grid) +
  draw_label("Study ID removed", x = 0.005, y = 0.5, angle = 90,
             vjust = 0.5, size = 11, fontface = "bold")

print(final_sensitivity_plot)

ggsave("Sensitivity_Grid.png", final_sensitivity_plot, width = 12, height = 8, dpi = 300)

all_jackknife_results_study %>%
  group_by(Outcome) %>%
  summarise(n_studies = n(),
            min_est = round(min(Estimate_Adj), 1),
            max_est = round(max(Estimate_Adj), 1),
            all_CI_exclude_0 = all(CI_Lower_Adj > 0 | CI_Upper_Adj < 0))

############################################################################################

# Checking that if the imputed data modified the results using Leave one out analysis

meta_modi <- meta %>%
  dplyr::mutate(
    #Log-response and variance for Yield
    Yield_LnR = log(Yield_kg_trt / Yield_kg_ctrl),
    
    Yield_Variance = (Yield_kg_SD_ctrl^2 / (Rep_ctrl * Yield_kg_ctrl^2)) +
      (Yield_kg_SD_trt^2  / (Rep_trt  * Yield_kg_trt^2)),
    
    #Log-response and variance for CP
    CP_LnR = log(CP_trt / CP_ctrl),
    
    CP_Variance = (CP_SD_ctrl^2 / (Rep_ctrl * CP_ctrl^2)) +
      (CP_SD_trt^2  / (Rep_trt  * CP_trt^2)),
    
    #Log-response and variance for WPI
    WPI_LnR = log(WPI_kg.m3_trt / WPI_kg.m3_ctrl),
    
    WPI_Variance = (WPI_kg.m3_SD_ctrl^2 / (Rep_ctrl * WPI_kg.m3_ctrl^2)) +
      (WPI_kg.m3_SD_trt^2  / (Rep_trt  * WPI_kg.m3_trt^2)),
    
    #Log-response and variance for TWP
    TWP_LnR = log(TWP_kg.m3_trt / TWP_kg.m3_ctrl),
    
    TWP_Variance = (TWP_kg.m3_SD_ctrl^2 / (Rep_ctrl * TWP_kg.m3_ctrl^2)) +
      (TWP_kg.m3_SD_trt^2  / (Rep_trt  * TWP_kg.m3_trt^2)),
    
    #weights
    Yield_weight = 1/Yield_Variance,
    CP_weight = 1/CP_Variance,
    WPI_weight = 1/WPI_Variance,
    TWP_weight = 1/TWP_Variance
  )

# --- Sensitivity Analysis: Imputed vs. Observed SDs ---
comparison_list <- list()
outcomes <- c("Yield", "CP", "WPI", "TWP")

for (out in outcomes) {
  # 1. Main Model (Imputed + Original - from meta_es)
  # Filter to ensure rows have data for this specific outcome
  data_full <- meta_es[!is.na(meta_es[[paste0(out, "_LnR")]]), ]
  
  mod_full <- rma.mv(
    yi = data_full[[paste0(out, "_LnR")]],
    V = data_full[[paste0(out, "_Variance")]],
    random = ~ 1 | Study_ID / Location_ID,
    data = data_full
  )
  
  # 2. Observed Only Model (Original SDs only - from meta)
  # Filter to ensure rows have data for this specific outcome in 'meta'
  data_obs <- meta_modi[!is.na(meta_modi[[paste0(out, "_LnR")]]), ]
  
  mod_obs <- rma.mv(
    yi = data_obs[[paste0(out, "_LnR")]],
    V = data_obs[[paste0(out, "_Variance")]],
    random = ~ 1 | Study_ID / Location_ID,
    data = data_obs
  )
  
  # Helper to format as Percentage Change: (exp(x)-1)*100
  fmt_res <- function(m) {
    est <- (exp(coef(m)) - 1) * 100
    lo  <- (exp(m$ci.lb) - 1) * 100
    hi  <- (exp(m$ci.ub) - 1) * 100
    paste0(round(est, 2), "% (", round(lo, 2), ", ", round(hi, 2), ")")
  }
  
  # Store results in a row
  comparison_list[[out]] <- data.frame(
    Outcome = out,
    Full_k = mod_full$k,
    Full_Estimate_Percent = fmt_res(mod_full),
    Observed_k = mod_obs$k,
    Observed_Estimate_Percent = fmt_res(mod_obs),
    stringsAsFactors = FALSE
  )
}

# Final Table for Supplementary
sensitivity_table <- do.call(rbind, comparison_list)
print(sensitivity_table)
##########################################################################################

# Checking weights of observed and imputed data added in supplementary table2

library(dplyr)
library(ggplot2)
library(tidyr)

outcomes <- c("Yield", "CP", "WPI", "TWP")
all_weight_data <- list()
summary_list <- list()

for (out in outcomes) {
  # 1. Identify the SD column in original data
  actual_col <- colnames(meta_modi)[grep(paste0(out, ".*SD"), colnames(meta_modi), ignore.case = TRUE)][1]
  
  if (!is.na(actual_col)) {
    # 2. Flag imputed rows and filter for this specific outcome
    # We assume meta_es and meta_modi are aligned by row
    temp_es <- meta_es
    temp_es$is_imputed <- is.na(meta_modi[[actual_col]])
    
    # Filter for rows that actually have an Effect Size for this outcome
    dat_sub <- temp_es[!is.na(temp_es[[paste0(out, "_LnR")]]), ]
    
    # 3. Run the model to get weights
    mod <- rma.mv(
      yi = dat_sub[[paste0(out, "_LnR")]], 
      V = dat_sub[[paste0(out, "_Variance")]],
      random = ~ 1 | Study_ID / Location_ID, 
      data = dat_sub
    )
    
    # Extract weights
    dat_sub$Weight_Pct <- weights(mod)
    dat_sub$Outcome_Name <- out
    
    # Store for plotting
    all_weight_data[[out]] <- dat_sub %>% select(Outcome_Name, is_imputed, Weight_Pct)
    
    # 4. Create Summary for this outcome
    summary_list[[out]] <- dat_sub %>%
      group_by(is_imputed) %>%
      summarise(
        Outcome = out,
        N_Studies = n(),
        Total_Weight = sum(Weight_Pct),
        Avg_Weight = mean(Weight_Pct),
        Max_Weight = max(Weight_Pct)
      ) %>%
      mutate(Source = if_else(is_imputed, "Imputed", "Observed"))
  }
}

# --- Final Table ---
final_summary_table <- bind_rows(summary_list) %>%
  select(Outcome, Source, N_Studies, Total_Weight, Avg_Weight, Max_Weight) %>%
  mutate(across(where(is.numeric), ~ round(., 2)))

print("Final Weight Distribution Summary:")
print(final_summary_table)



#############################################################################################

# --- SECTION 1: Data Preparation and Bias Adjustment ---

# Run the base meta-analysis model using the 'meta' package
m.yield <- metagen(
  TE = Yield_LnR,          # Treatment Effect (Log Response Ratio)
  seTE = Yield_SE,         # Standard Error of the Treatment Effect
  data = meta_es,
  studlab = Study_ID / Location_ID,      # Study labels for identification
  comb.random = TRUE       # Use a random-effects model
)

# Apply the trim-and-fill method to estimate and impute missing studies on the left side
m.yield.tf <- trimfill(m.yield, side = "left")


# --- SECTION 2: Compare Overall Pooled Effects (Sensitivity Analysis) ---
# This is the key comparison for the main sensitivity analysis in your paper

# Extract the overall pooled effect from the original data model
original_pooled <- data.frame(
  Dataset = "Original Data (Before Adjustment)",
  Metric = "Yield_Overall_Effect",
  P_value = summary(m.yield)$random$p,        # P-value for the overall random effects model
  Estimate_LnR = summary(m.yield)$TE.random, # The average effect size (LnR)
  stringsAsFactors = FALSE
)

# Extract the overall pooled effect from the bias-adjusted data model
adjusted_pooled <- data.frame(
  Dataset = "Bias Adjusted (Trim-Fill)",
  Metric = "Yield_Overall_Effect",
  P_value = summary(m.yield.tf)$random$p,
  Estimate_LnR = summary(m.yield.tf)$TE.random,
  stringsAsFactors = FALSE
)

# Combine results into a single, clean table for direct comparison
final_comparison_table <- bind_rows(original_pooled, adjusted_pooled)

# Print the comparison table
print(final_comparison_table)


##################################################### Making the MAP ########################################

# --- 1. Get the US Map Data ---
us_states <- map_data("state")

# --- 2. Prepare Your Data (Piped Approach) ---

# Select unique latitude and longitude pairs from your meta_es data frame
# This handles multiple occurrences of the same location by keeping only one row per coordinate pair
unique_locations <- meta_es %>%
  select(Longitude, Latitude, Reference, Study_ID, Location_ID) %>% # Select necessary columns
  distinct() %>%                         # Keep only unique rows
  # Optional: Filter to ensure points are within typical US bounds if you have global data
  filter(Longitude > -130 & Longitude < -65 & Latitude > 20 & Latitude < 50)


# --- 3. Create the Map Visualization ---

ggplot() +
  # Layer 1: Plot the state boundaries (polygons)
  geom_polygon(
    data = us_states, 
    aes(x = long, y = lat, group = group), 
    fill = "grey90",         # Light grey fill
    color = "black",         # Black state borders
    size = 0.3
  ) +
  # Layer 2: Plot your unique data points
  geom_point(
    data = unique_locations, 
    aes(x = Longitude, y = Latitude), # Color points by Study ID
    size = 6,              
    alpha = 0.5,
    show.legend = FALSE # Hide the legend if you have too many studies
  ) +
  # Use the Albers projection for a nice US view and set limits
  coord_map("albers", lat0 = 39, lat1 = 45) + 
  theme_bw() +      # Use a clean black and white theme
  labs(
    title = "",
    x = "Longitude",
    y = "Latitude"
  )

# Use ggsave to save the map as an image file
# ggsave("US_Locations_Map_Unique_Points.png", width = 10, height = 7, dpi = 300)


##################################### FINAL MAP ######### WorldClim 1970-2000 baseline used here ##############################################

study_points_sf <- st_as_sf(meta_es, 
                            coords = c("Longitude", "Latitude"), 
                            crs = 4326)

# Geometry  
sf::sf_use_s2(FALSE) 
us_mask_sf <- st_as_sf(maps::map("state", plot = FALSE, fill = TRUE)) %>%
  st_make_valid() %>% 
  st_union()
us_vect <- vect(us_mask_sf)

# Process Growing Season Precip
prec_global <- geodata::worldclim_global(var = "prec", res = 10, path = tempdir())
gsp_raster <- sum(prec_global[[4:10]]) # April to October
gsp_masked <- crop(gsp_raster, us_vect) %>% mask(us_vect)

gsp_df <- as.data.frame(gsp_masked, xy = TRUE)
colnames(gsp_df) <- c("long", "lat", "prec_sum")


map <- ggplot() +
  # Background Heatmap
  geom_tile(data = gsp_df, aes(x = long, y = lat, fill = prec_sum)) +
  
  # Custom RGB Scale for Precipitation
  scale_fill_gradientn(
    colors = c("#e41a1c", "#4daf4a", "#377eb8"), 
    name = "Mean annual\nprecipitation (mm)",
    na.value = "transparent"
  ) +
  
  # State Borders
  geom_sf(data = st_as_sf(maps::map("state", plot = FALSE, fill = TRUE)), 
          fill = NA, color = "black", size = 0.2) +
  
  geom_sf(data = study_points_sf, aes(shape = "Study sites"), 
          size = 3, fill = "white", color = "black", alpha = 0.8) +
  
  # Define the manual scale for the points
  scale_shape_manual(name = "", values = c("Study sites" = 21)) +
  
  # Scale bar
  annotation_scale(location = "bl", width_hint = 0.2) +
  
  # North arrow
  annotation_north_arrow(location = "tl", which_north = "true", 
                         pad_x = unit(0.2, "in"), pad_y = unit(0.2, "in"),
                         style = north_arrow_fancy_orienteering) +
  
  # Coordinates and Theme
  coord_sf(crs = 4326, xlim = c(-125, -67), ylim = c(25, 50)) + 
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    legend.box = "vertical",           # Stack legends vertically
    legend.spacing.y = unit(0, "cm"),  # Move them closer together
    axis.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 11),
  ) +
  labs(
    title = "",
    x = "Longitude",
    y = "Latitude"  
  )

map

ggsave("Map.png", map, width = 14, height = 8, dpi = 300)

##############################################################################################
############################# ADDRESSING REVIWER'S COMMENTS #################################
# =====================================================================
# CONTROL DEFINITION CHECK
# Reviewer comment: definition of full irrigation differs among studies
# =====================================================================
library(metafor)
library(dplyr)

# --- 1. Group control-definition methods into four categories ---
group_control <- function(x) {
  case_when(
    x %in% c("ET_based", "ET_based+Irri_time", "Check_book")       ~ "ET-based",
    x %in% c("MAD_based", "Soil_moisture_tension_based")                ~ "Soil-moisture-based",
    x == "supplementary_irri"                                           ~ "Supplemental",
    TRUE                                                                ~ "Other (cutting time, experience)")
}

# Check: which original labels ended up in "Other"? (catches new labels or typos)
meta_es %>%
  filter(!is.na(Irri_trt_type)) %>%
  mutate(ctrl_group = group_control(Irri_trt_type)) %>%
  filter(ctrl_group == "Other (timing, supplemental, checkbook)") %>%
  distinct(Irri_trt_type) %>%
  print()

# --- 2. Moderator test + sensitivity analysis for one outcome ---
control_check <- function(outcome, yi_col, vi_col) {
  d <- meta_es %>%
    filter(!is.na(.data[[yi_col]]), !is.na(.data[[vi_col]]), !is.na(Irri_trt_type)) %>%
    mutate(ctrl_group = group_control(Irri_trt_type))
  
  cat("\n=====", outcome, "=====\n")
  print(d %>% group_by(ctrl_group) %>%
          summarise(n_effects = n(), n_studies = n_distinct(Study_ID), .groups = "drop"))
  
  # Moderator test (model-based and cluster-robust)
  m_grp <- rma.mv(d[[yi_col]], d[[vi_col]],
                  mods   = ~ factor(d$ctrl_group),
                  random = ~ 1 | Study_ID/Location_ID,
                  data   = d)
  m_rob <- robust(m_grp, cluster = d$Study_ID, clubSandwich = TRUE)
  
  # Sensitivity: exclude experience-based controls
  d_ex  <- filter(d, Irri_trt_type != "Experience_based")
  m_all <- rma.mv(d[[yi_col]],    d[[vi_col]],    random = ~ 1 | Study_ID/Location_ID, data = d)
  m_ex  <- rma.mv(d_ex[[yi_col]], d_ex[[vi_col]], random = ~ 1 | Study_ID/Location_ID, data = d_ex)
  
  data.frame(
    Outcome           = outcome,
    QM_df             = m_grp$QMdf[1],
    QM                = round(m_grp$QM, 2),
    QM_P              = round(m_grp$QMp, 3),
    Robust_F          = round(m_rob$QM, 2),
    Robust_P          = round(m_rob$QMp, 3),
    k_all             = m_all$k,
    studies_all       = n_distinct(d$Study_ID),
    pct_all           = round((exp(as.numeric(coef(m_all))) - 1) * 100, 1),
    CI_all            = sprintf("%.1f to %.1f", (exp(m_all$ci.lb) - 1) * 100, (exp(m_all$ci.ub) - 1) * 100),
    k_excl_exp        = m_ex$k,
    studies_excl_exp  = n_distinct(d_ex$Study_ID),
    pct_excl_exp      = round((exp(as.numeric(coef(m_ex))) - 1) * 100, 1),
    CI_excl_exp       = sprintf("%.1f to %.1f", (exp(m_ex$ci.lb) - 1) * 100, (exp(m_ex$ci.ub) - 1) * 100)
  )
}

# --- 3. Run for all outcomes ---
control_results <- bind_rows(
  control_check("Yield", "Yield_LnR", "Yield_Variance"),
  control_check("CP",    "CP_LnR",    "CP_Variance"),
  control_check("WPI",   "WPI_LnR",   "WPI_Variance"),
  control_check("TWP",   "TWP_LnR",   "TWP_Variance")
)

print(control_results)

# --- 4. Save for supplementary material / reviewer reply ---
write.csv(control_results, "Control_definition_check.csv", row.names = FALSE)


#################################################################################################################
library(sf); library(dplyr); library(ggplot2); library(tigris)
options(tigris_use_cache = TRUE)

irr <- read.csv("alfalfa_irrigated_county_2022.csv") %>%
  mutate(GEOID    = sprintf("%02d%03d", State.ANSI, County.ANSI),
         withheld = trimws(Value) == "(D)",
         acres    = suppressWarnings(as.numeric(gsub(",", "", Value))))

counties <- counties(cb = TRUE, resolution = "20m") %>%
  filter(!STATEFP %in% c("02", "15", "60", "66", "69", "72", "78")) %>%
  left_join(irr, by = "GEOID")

irr_map <- ggplot() +
  geom_sf(data = counties, fill = "gray95", color = NA) +                       # none reported
  geom_sf(data = filter(counties, withheld), fill = "gray70", color = NA) +    # withheld (D)
  geom_sf(data = filter(counties, !is.na(acres)), aes(fill = acres), color = NA) +
  geom_sf(data = states, fill = NA, color = "gray30", linewidth = 0.3) +
  geom_sf(data = sites, shape = 21, fill = "white", color = "black", size = 2.5) +
  scale_fill_viridis_c(trans = "log10", labels = scales::comma,
                       name = "Irrigated alfalfa\n(acres harvested)", direction = -1) +
  coord_sf(crs = 5070) +
  theme_bw() +
  theme(panel.grid = element_blank())

print(irr_map)
ggsave("Irrigated_alfalfa_sites_map.png", irr_map, width = 9, height = 6, dpi = 300)

################################################################################################################

# =====================================================================
# HETEROGENEITY: where it is (multilevel I2) and what explains it (pseudo-R2)
# Run AFTER your moderator script, so the *_grp columns exist in meta_es
# =====================================================================
library(metafor)
library(dplyr)

# --- 1. Outcomes ---
outcomes <- list(Yield = c("Yield_LnR", "Yield_Variance"),
                 CP    = c("CP_LnR",    "CP_Variance"),
                 WPI   = c("WPI_LnR",   "WPI_Variance"),
                 TWP   = c("TWP_LnR",   "TWP_Variance"))

# --- 2. Your moderators ---
moderators <- c("Irrigation_method", "Irri_amt", "Degree_of_irri_deficit_grp",
                "Irrigation_scheduling", "MAP_mm_grp", "MAT_degreeC_grp",
                "Texture", "Soil_BD_grp", "Soil_pH_grp", "Soil_OM_grp",
                "Winterhardiness_grp", "Latitude_grp", "Longitude_grp")

missing <- setdiff(moderators, names(meta_es))
if (length(missing) > 0) message("Not found in meta_es (run your moderator script first): ",
                                 paste(missing, collapse = ", "))
moderators <- intersect(moderators, names(meta_es))

# --- 3. Multilevel I2: share of variance at study and location level (Nakagawa & Santos 2012) ---
i2_ml <- function(m) {
  W <- diag(1 / m$vi); X <- model.matrix(m)
  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  typical_v <- (m$k - m$p) / sum(diag(P))
  tot <- sum(m$sigma2) + typical_v
  round(100 * c(I2_total = sum(m$sigma2), I2_study = m$sigma2[1], I2_location = m$sigma2[2]) / tot, 1)
}

i2_table <- bind_rows(lapply(names(outcomes), function(o) {
  yi <- outcomes[[o]][1]; vi <- outcomes[[o]][2]
  d  <- meta_es %>% filter(!is.na(.data[[yi]]), !is.na(.data[[vi]]))
  m  <- rma.mv(d[[yi]], d[[vi]], random = ~ 1 | Study_ID/Location_ID, data = d)
  data.frame(Outcome = o, t(i2_ml(m)),
             sigma2_study = round(m$sigma2[1], 4), sigma2_location = round(m$sigma2[2], 4))
}))
print(i2_table)

# --- 4. Moderators: QM and pseudo-R2 (% of heterogeneity explained) ---
mod_R2 <- function(o, mod) {
  yi <- outcomes[[o]][1]; vi <- outcomes[[o]][2]
  d  <- meta_es %>% filter(!is.na(.data[[yi]]), !is.na(.data[[vi]]), !is.na(.data[[mod]]))
  d$grp <- factor(d[[mod]])
  if (nlevels(d$grp) < 2) return(NULL)
  # null model fitted on the same rows, so R2 compares like with like
  m0 <- rma.mv(d[[yi]], d[[vi]], random = ~ 1 | Study_ID/Location_ID, data = d)
  m1 <- rma.mv(d[[yi]], d[[vi]], mods = ~ grp, random = ~ 1 | Study_ID/Location_ID, data = d)
  r2 <- max(0, (sum(m0$sigma2) - sum(m1$sigma2)) / sum(m0$sigma2))
  data.frame(Outcome = o, Moderator = mod, levels = nlevels(d$grp),
             k = m1$k, studies = n_distinct(d$Study_ID),
             QM = round(m1$QM, 2), QM_P = signif(m1$QMp, 3),
             pseudo_R2_pct = round(100 * r2, 1))
}

r2_table <- bind_rows(lapply(names(outcomes), function(o)
  bind_rows(lapply(moderators, function(mod)
    tryCatch(mod_R2(o, mod), error = function(e) {
      message("Skipped ", o, " / ", mod, ": ", conditionMessage(e)); NULL })))))

r2_table <- r2_table %>% arrange(Outcome, desc(pseudo_R2_pct))   # strongest moderators first
print(r2_table, row.names = FALSE)

# --- 5. Save ---
write.csv(i2_table, "Heterogeneity_I2_levels.csv", row.names = FALSE)
write.csv(r2_table, "Heterogeneity_moderators.csv", row.names = FALSE)

