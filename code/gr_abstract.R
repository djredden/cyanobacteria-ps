rm(list = ls())
library(tidyverse)
library(patchwork)
library(ggpubr)
library(ggtext)
library(scales)
source(here::here("code/functions.R"))

theme_set(theme_classic(base_size = 24) +
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
  filter(target != "16S rRNA") %>% 
  group_by(location, type, target) %>% 
  summarise(prop_detected = 100 * mean(detect)) %>% 
  mutate(target = paste("*", target, "*", sep = ""),
         target = if_else(target == "*mcyE/ndaF*", "*mcyE/<br>ndaF*", target)) %>% 
  ungroup() %>% 
  ggplot(aes(target, prop_detected, fill = type)) +
  geom_bar(stat = "identity", 
           position = "dodge") +
  ggsci::scale_fill_jco() +
  facet_wrap(vars(location)) +
  theme(legend.position = "bottom",
        axis.text.x = element_markdown()) +
  labs(x = NULL,
       y = expression(paste("Positive Gene Target Detections (%)")),
       fill = NULL)

ggsave(here::here("output/gr_abstract.png"), height = 8, width = 8, dpi = 300)
