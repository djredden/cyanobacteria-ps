rm(list = ls())
library(tidyverse)
library(patchwork)
library(ggpubr)
library(scales)
source(here::here("code/functions.R"))

theme_set(theme_bw(base_size = 24) +
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

cleaned_data <- read_csv(here::here("cleaned_data/plot_data.csv")) %>%
  mutate(label = case_when(
    target == "sxtA" ~ "bolditalic(`sxtA`)",
    target == "mcyE/ndaF" ~ "bolditalic(`mcyE/ndaF`)",
    target == "cyrA" ~ "bolditalic(`cyrA`)",
    TRUE ~ target)) %>% 
  mutate(detect = if_else(conc < lod | is.na(conc), FALSE, TRUE))



cleaned_data %>% 
  filter(target == "mcyE/ndaF") %>% 
  group_by(location, type) %>% 
  summarise(prop_detected = 100 * mean(detect)) %>% 
  ungroup() %>% 
  ggplot(aes(location, prop_detected, fill = type)) +
  geom_bar(stat = "identity", 
           position = "dodge") +
  ggsci::scale_fill_jco() +
  theme_classic() +
  theme(legend.position = "bottom") +
  labs(x = NULL,
       y = expression(paste("Positive ", italic("mcyE"), " Detections (%)")),
       fill = NULL)

ggsave(here::here("output/gr_abstract.jpg"), height = 5, width = 4, dpi = 300)
