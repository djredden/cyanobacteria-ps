rm(list = ls())
source(here::here("code/functions.R"))

#### Read in the DWTP intake grab sample data 
#### (NOTE: 16S data from Sep 16 and 29 not provided/run)

path <- here::here("raw_data/dwtp_grab_samples/2021-qpcr-results.xlsx")
sheets <- readxl::excel_sheets(path)

# Read and tidy 16S data
cp_16s_curve <- readxl::read_xlsx(path, sheet = sheets[1],
                                  skip = 2) %>% 
  janitor::clean_names() %>% 
  # Remove the standard that was not used in developing the curve
  filter(is.na(notes_2)) %>% 
  mutate(target = "total_cyano") %>% 
  select(target, sample_id, cq = fam_ct) 

cp_tox_curve <- readxl::read_xlsx(path, sheet = sheets[3],
                                  skip = 2) %>% 
  janitor::clean_names() %>% 
  select(sample_id, contains("_ct")) %>% 
  mutate(sample_id = as.numeric(sample_id)) %>% 
  filter(!is.na(sample_id)) %>% 
  rename(mcye_ndaf = "fam_ct",
         cyra = "cy3_ct",
         sxta = "txr_ct") %>% 
  pivot_longer(-sample_id, names_to = "target",
               values_to = "cq")
 
cp_curve_results <- bind_rows(cp_16s_curve, cp_tox_curve) %>% 
  nest(data = -target) %>% 
  mutate(model = map(data, ~ lm(cq ~ log10(sample_id), data =.))) %>% 
  mutate(estimates = map(model, broom::tidy),
         stats = map(model, broom::glance)) %>%
  unnest(c(estimates, stats),
         names_sep = "_") %>% 
  pivot_wider(
    id_cols = c(target, stats_adj.r.squared),
    names_from = estimates_term,
    values_from = estimates_estimate
  ) %>% 
  rename(
    "intercept" = `(Intercept)`,
    "slope" = `log10(sample_id)`,
    "r_squared" = stats_adj.r.squared
  ) %>% 
  # Compute the efficiency of each qPCR target
  mutate(efficiency = ((10^(-1 / slope)) - 1) * 100,
         r_squared = round(r_squared, 6)) 
 



