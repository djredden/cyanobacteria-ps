rm(list = ls())
library(tidyverse)
library(scales)
source(here::here("code/functions.R"))

theme_set(theme_bw() +
            theme(
              axis.text = element_text(size = 14, color = "black"),
              axis.title = element_text(size = 16),
              legend.title = element_text(size = 16),
              legend.text = element_text(size = 14),
              strip.text = element_text(size = 14, face = "bold"),
              strip.background = element_rect(fill = "white"),
              panel.spacing.x = unit(0.5, "lines"),
              panel.spacing.y = unit(1.5, "lines")
            ))

#### Read in the DWTP intake grab sample data 
#### (NOTE: 16S data from Sep 16 and 29 not provided/run)

path <- here::here("raw_data/dwtp_grab_samples/2021-qpcr-results.xlsx")
sheets <- readxl::excel_sheets(path)

# Read biorad curve results
biorad <- read_csv(here::here("output/biorad.csv")) %>% 
  select("sample_id" = "starting_quantity_(sq)",
         target, cq) %>% 
  mutate(instrument = "BioRad") 
  

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


# Labels for plot
lab1 <- c(expression(paste("16S rRNA")),
          expression(paste(italic("mcyE/ndaF"))), 
          expression(paste(italic("cyrA"))),
          expression(paste(italic("sxtA")))
          )




curves <- bind_rows(cp_16s_curve, cp_tox_curve) %>% 
  mutate(instrument = "Cepheid") %>% 
  bind_rows(biorad) %>% 
  mutate(label = case_when(
    target == "total_cyano" ~ "16S~rRNA",
    target == "mcye_ndaf" ~ "italic(`mcyE/ndaF`)",
    target == "cyra" ~ "italic(`cyrA`)",
    target == "sxta" ~ "italic(`sxtA`)",
    TRUE ~ target
  )) %>%
  mutate(label = factor(label,
                        levels =
                          c("16S~rRNA",
                            "italic(`mcyE/ndaF`)",
                            "italic(`cyrA`)",
                            "italic(`sxtA`)"
                          )
  )) %>% 
  ggplot(aes(log(sample_id), cq, shape = label)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE,
              color = ggsci::pal_jco("default")(1)) +
  facet_wrap(vars(instrument), scales = "free") +
  scale_shape_discrete(name = "Gene Target", 
                       labels = lab1) +
  
  theme(legend.text.align = 0) +
  labs(y = "Cq",
       x = "Log Gene Concentration")

ggsave(curves,
       filename = here::here("output/standard_curves.tiff"),
       height = 6, width = 10, dpi = 300
)


# NTC checks

field_ntcs <- read_csv(here::here("cleaned_data/field_data_compiled.csv")) %>% 
  filter(content == "ntc") %>% 
  arrange(target)
 




