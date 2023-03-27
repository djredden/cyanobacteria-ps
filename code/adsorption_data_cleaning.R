rm(list = ls())
library(tidyverse)
library(scales)
source(here::here("code/functions.R"))

file_list <- list.files(here::here("raw_data/adsorption_data"))

# Read std curve data
curve_results <- read_csv(here::here("output/curve-results.csv"))

# Read and clean
adsorption_data_raw <- map(file_list, read_ads_data) %>%
  bind_rows() %>% 
  select(sample, content, fluor, assay = biological_set_name, cq) %>% 
  #creates a column for sample location and material based on sample names
  separate(sample, c("lake","dash", "material", "time"), sep = " ") %>% 
  mutate(target = case_when(
    (fluor == "cal orange 560" & assay == "16s rrna") ~ "iac",
    (fluor == "cal orange 560" & assay == "toxins") ~ "cyra",
    (fluor == "fam" & assay == "16s rrna") ~ "total_cyano",
    (fluor == "fam" & assay == "toxins") ~ "mcye_ndaf",
    (fluor == "cal red 610" & assay == "toxins") ~ "sxta"
  )) %>% 
  filter(!target == "iac") %>% 
  inner_join(curve_results, by = c("target")) %>% 
  # Compute the gene copy concentrations per uL in the well
  mutate(gu_well = 10^((cq - intercept) / slope))

adsorption_8hr_test <- adsorption_data_raw %>% 
  #Removes rows that are not part of the sampling program (ie. not from FL), including CP
  filter(lake %in% c("wl", "gl", "cl", "lp")) %>% 
  mutate(sample = case_when(material == "c"~"control",
                            material == "ac"~"acrylic_copolymer",
                            material == "cn"~"cellulose_nitrate",
                            material == "n"~"nylon",
                            material == "g"~"gauze",
                            TRUE~"initial")) %>% 
  mutate(lake = case_when(lake == "wl"~"william",
                          lake == "gl"~"grand",
                          lake == "lp"~"lumsden",
                          lake == "cl"~"chocolate")) %>% 
  mutate(type = case_when(sample == "control"~"volume",
                          sample == "initial"~"volume",
                          TRUE~"area")) %>% 
  select(-dash, -material, -time) %>% 
  # Add column containing passive sampler area (25 cm diameter filter) and grab sample volume (100 mL)
    mutate(area_vol = if_else(type == "area", 
                            pi*(2.5/2)^2,
                            100)
  ) %>% 
  mutate(gu_conc = gu_well * 500 / area_vol) %>% 
  filter(target == "total_cyano") %>% 
  mutate(lod = 9*500/area_vol) 
                       
write_csv(adsorption_8hr_test, here::here("cleaned_data/adsorption_8hr.csv"))
  
# 24hr adsorption data

ads_24_raw <- read_ads_data("CyanoDTec_16S_FL_MONC_Adsorption_Results.xlsx") %>% 
  select(sample, content, fluor, assay = biological_set_name, cq, day, run, id, dilution) %>% 
  separate(sample, c("replicate","dash", "material", "time"), sep = " ") %>% 
  mutate(sample = case_when(material == "c"~"control",
                            material == "ac"~"acrylic_copolymer",
                            material == "cn"~"cellulose_nitrate",
                            material == "n"~"nylon",
                            material == "g"~"gauze",
                            TRUE~"initial")) %>% 
  mutate(type = case_when(sample == "control"~"volume",
                          sample == "initial"~"volume",
                          TRUE~"area")) %>% 
  #creates a column for sample location and material based on sample names
  separate(sample, c("sample","dash", "material", "time"), sep = " ") %>% 
  mutate(target = case_when(
    (fluor == "cal orange 560" & assay == "16s rrna") ~ "iac",
    (fluor == "cal orange 560" & assay == "toxins") ~ "cyra",
    (fluor == "fam" & assay == "16s rrna") ~ "total_cyano",
    (fluor == "fam" & assay == "toxins") ~ "mcye_ndaf",
    (fluor == "cal red 610" & assay == "toxins") ~ "sxta"
  )) %>% 
  filter(!target == "iac") %>% 
  inner_join(curve_results, by = c("target")) %>% 
  # Compute the gene copy concentrations per uL in the well
  mutate(gu_well = 10^((cq - intercept) / slope)) %>% 
  #multiplying initial GU by dilution factor
  mutate(gu_well = ifelse(dilution == 2, gu_well*2, gu_well)) %>% 
  mutate(area_vol = case_when(
    type == "volume" & run == "d" ~ 8,
    type == "volume" ~10,
    type == "area" ~ pi*(2.5/2)^2
  )) %>% 
  mutate(gu_conc = gu_well * 500 / area_vol) %>% 
  filter(target == "total_cyano") %>% 
  mutate(lod = 9 * 500 / area_vol) %>% 
  mutate(day = ifelse(day == 1, "one", "two")) %>% 
  ungroup() %>% 
  mutate(run = case_when(run == "a"~"run 1a",
                         run == "b"~"run 1b",
                         run == "c"~"run 2a",
                         run == "d"~"run 2b")) %>% 
  select(assay, cq, day, run, id, dilution, sample, target,
         gu_well, area_vol, gu_conc, lod)


# Save the cleaned files

write_csv(adsorption_24hr, here::here("cleaned_data/adsorption_24hr.csv"))
write_csv(adsorption_8hr, here::here("cleaned_data/adsorption_8hr.csv"))






