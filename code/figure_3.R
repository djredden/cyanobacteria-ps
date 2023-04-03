rm(list = ls())
library(tidyverse)
library(patchwork)
library(ggpubr)
library(scales)
source(here::here("code/functions.R"))

theme_set(theme_bw() +
            theme(
              axis.text = element_text(size = 14, color = "black"),
              axis.title = element_text(size = 16),
              legend.title = element_text(size = 14),
              legend.text = element_text(size = 14),
              strip.text = element_text(size = 16, face = "bold"),
              strip.background = element_rect(fill = "white"),
              panel.spacing.x = unit(0.5, "lines"),
              panel.spacing.y = unit(1, "lines")
            ))

cleaned_data <- read_csv(here::here("cleaned_data/plot_data.csv")) %>% 
  mutate(label = case_when(
    label == "italic(`sxtA`)" ~ "bolditalic(`sxtA`)",
    label == "italic(`mcyE/ndaF`)" ~ "bolditalic(`mcyE/ndaF`)",
    label == "italic(`cyrA`)" ~ "bolditalic(`cyrA`)",
    TRUE ~ label)) %>% 
  mutate(label = factor(label, levels = 
                          c("bolditalic(`mcyE/ndaF`)",
                            "bolditalic(`cyrA`)",
                            "bolditalic(`sxtA`)")))

# Figure 3 - times series of all toxin data faceted by sample type and target

figure_3 <- cleaned_data %>%
  filter(assay == "toxins") %>%
  # Replace NAs with zeros for plotting
  mutate(conc = if_else(is.na(conc), 0, conc)) %>%
  mutate(type = if_else(
    type == "Grab", 
    "bold(Grab~(GU~mL^{-1}))",
    "bold(Passive~(GU~cm^{-2}))"
  )) %>%
  ggplot(aes(x = date, y = conc, fill = location, shape = frequency)) +
  geom_point(position = "jitter", size = 3) +
  geom_hline(aes(yintercept = lod, linetype = factor(lod))) +
  scale_shape_manual(values = c(24, 21)) +
  ggsci::scale_fill_jco() +
  guides(
    fill = guide_legend(override.aes = list(shape = 22, size = 3)),
    shape = guide_legend(override.aes = list(fill = "black", size = 3))
  ) +
  facet_grid(vars(type), vars(label),
    labeller = label_parsed
  ) +
  ylab(NULL) +
  scale_y_log10(
    limits = c(1, 1e6),
    labels = trans_format("log10", math_format(10^.x)),
    expand = c(0, 0)
  ) +
  labs(
    y = "Gene concentration",
    x = "",
    fill = "Location",
    shape = "Sample Frequency",
    linetype = "LOD"
  )

ggsave(figure_3,
  filename = here::here("output/figure_3.tiff"),
  height = 8, width = 12, dpi = 300
)
