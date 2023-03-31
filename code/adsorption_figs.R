rm(list = ls())
library(tidyverse)
library(ggtext)
library(scales)
library(ggtext)

source(here::here("code/functions.R"))

theme_set(theme_bw() +
            theme(
              axis.text = element_text(size = 10, color = "black"),
              axis.title = element_text(size = 12),
              legend.title = element_text(size = 11),
              legend.text = element_text(size = 10),
              strip.text = element_text(size = 12, face = "bold"),
              strip.background = element_rect(fill = "white"),
              panel.spacing.x = unit(0.5, "lines"),
              panel.spacing.y = unit(1, "lines")
            ))

# Read files
adsorption_24hr <- read_csv(here::here("cleaned_data/adsorption_24hr.csv")) %>% 
  mutate(sample = case_when(
    sample == "acrylic_copolymer" ~ "Acrylic Copolymer",
    sample == "cellulose_nitrate" ~ "Cellulose Nitrate",
    sample == "gauze" ~ "Gauze",
    TRUE ~ sample
  ))

adsorption_8hr <- read_csv(here::here("cleaned_data/adsorption_8hr.csv")) %>% 
  mutate(sample = case_when(
    sample == "acrylic_copolymer" ~ "Acrylic Copolymer",
    sample == "cellulose_nitrate" ~ "Cellulose Nitrate",
    sample == "gauze" ~ "Gauze",
    TRUE ~ sample
  ))

# Draw boxplot for the 8 hr adsorption study
SI_8hr_adsorption <- adsorption_8hr %>%
  filter(
    !sample %in% c("initial", "control"),
    gu_conc > lod
  ) %>%
  mutate(sample = str_to_title(sample)) %>%
  ggplot(aes(x = sample, y = gu_conc)) +
  stat_summary(fun.data = bp.vals,
               geom = "boxplot",
               fill = ggsci::pal_jco("default")(1)) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  theme_classic() +
  labs(
    x = NULL,
    y = expression(Gene ~ Concentration ~ (GU ~ cm^{
      -2
    }))
  )

ggsave(here::here("output/SI_8hr_adsorption.jpg"), SI_8hr_adsorption,
       height = 6, width = 8, dpi = 300
)

# Draw boxplot for the 24 hr adsorption study
adsorption_24hr %>%
  filter(!sample %in% c("control", "initial")) %>%
  mutate(DOC = if_else(day == "one", 5.12, 10.18)) %>%
  mutate(DOC = factor(DOC)) %>% 
  group_by(day, DOC, sample) %>%
  ggplot(aes(sample, gu_conc, fill = DOC)) +
  stat_summary(fun.data = bp.vals,
               geom = "boxplot") +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  theme_classic(base_size = 18) +
  theme(
    strip.background = element_blank(),
    strip.text = element_blank()
  ) +
  scale_x_discrete(
    breaks = c("Acrylic Copolymer", "Cellulose Nitrate", "Gauze"),
    labels = c(
      "Acrylic <br> Copolymer", "Cellulose <br> Nitrate",
      "Gauze"
    )
  ) +
  ggsci::scale_fill_jco(
    breaks = c(5.12, 10.18),
    labels = c("Day 1", "Day 2")
  ) +
  facet_wrap(~day) +
  theme(axis.text.x = element_markdown()) +
  labs(
    x = NULL,
    y = expression(Gene ~ Concentration ~ (GU ~ cm^{
      -2
    })),
    fill = "Sampling Event"
  )

ggsave(here::here("output/figure_1.tiff"),
       height = 5, width = 11, dpi = 300
)
