rm(list = ls())
library(tidyverse)
library(brms)
library("emmeans")

source(here::here("code/functions.R"))

# Read files
adsorption_24hr <- read_csv(here::here("cleaned_data/adsorption_24hr.csv"))
adsorption_8hr <- read_csv(here::here("cleaned_data/adsorption_8hr.csv"))

# Compute summary statistics
mean_adsorptions <- bind_rows(
  mean_sd_fcn(adsorption_8hr) %>%
    mutate(study = 8),
  mean_sd_fcn(adsorption_24hr) %>%
    mutate(study = 24)
) 

# Model differences between adsorbents in 8 hr study
data_8hr <- adsorption_8hr %>%
  filter(!sample %in% c("initial", "control"))

# Fit the model
fit_8hr <- lm(log(gu_conc) ~ factor(sample), data = data_8hr)

# intercept free model
fit_8hr <- brms::brm(log(gu_conc) ~ 0 + sample,
               data = data_8hr,
               prior = set_prior("normal(0,10)", class = "b"),                    
               file = here::here("models/8hr_brm.rds"))

# Compare the difference between distribution of mean for three adsorbents vs nylon
nylon_delta <- as_draws_df(fit_8hr) %>% 
  select(ac = b_sampleacrylic_copolymer,
         cn = b_samplecellulose_nitrate,
         gz = b_samplegauze,
         ny = b_samplenylon) %>%
  mutate(across(everything(), ~exp(.x))) %>%
  group_by(row_number()) %>%
  mutate(avg = sum(ac + cn + gz)/3) %>% 
  ungroup() %>%
  mutate(diff = avg - ny) %>%
  summarize(mean = mean(diff),
            q2.5 =  quantile(diff, probs = 0.025),
            q97.5 = quantile(diff, probs = 0.975)) %>%
  mutate(across(where(is.numeric), ~format(., scientific = TRUE, digits=2)))

# Model the data to find differences between adsorbents in 24 hr study
compare_data <- adsorption_24hr %>%
  filter(!sample %in% c("initial", "control"))

# Fit intercept free model
fit <- brm(log(gu_conc) ~ 0 + sample * day,
           data = compare_data,
           prior = set_prior("normal(0,10)", class = "b"),                    
           file = here::here("models/24hr_brm.rds")
)


# Compute differences by drawing from posterior and subtracting
post_draws <- as_draws_df(fit)

cred_ints <- post_draws %>%
  select(ac = b_sampleacrylic_copolymer,
         cn = b_samplecellulose_nitrate,
         g = b_samplegauze) %>%
  mutate(ac_cn = exp(ac) - exp(cn),
         ac_g = exp(ac) - exp(g),
         cn_g = exp(cn) - exp(g)) %>%
  select(-c(ac, cn, g)) %>%
  pivot_longer(cols = everything(), names_to = "ads", values_to = "diff") %>%
  group_by(ads) %>%
  summarize(mean = mean(diff),
            q2.5 =  quantile(diff, probs = 0.025),
            q97.5 = quantile(diff, probs = 0.975)) %>%
  mutate(across(where(is.numeric), ~format(., scientific = TRUE, digits=2)))


