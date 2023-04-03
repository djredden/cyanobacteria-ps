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
    target == "sxtA" ~ "bolditalic(`sxtA`)",
    target == "mcyE/ndaF" ~ "bolditalic(`mcyE/ndaF`)",
    target == "cyrA" ~ "bolditalic(`cyrA`)",
    TRUE ~ target)) %>% 
  mutate(detect = if_else(conc < lod | is.na(conc), "< LOD", "> LOD"),
         detect = factor(detect, levels = c("> LOD", "< LOD"))) %>% 
  mutate(conc = if_else(is.na(conc), 1, conc)) %>% 
  mutate(date = lubridate::yday(date),
         date_cont = as.Date(date, origin = "2021-01-01")) %>%
           mutate(location = factor(location,
                           levels = c("Inlet", "Outlet", "Intake")
  )) %>%
  mutate(type = factor(type)) 

### GAM for mcyE
mcye_data <- cleaned_data %>%
  filter(target == "mcyE/ndaF") %>% 
  group_by(location, type) %>%
  mutate(conc = scale(log(conc))) %>% 
  ungroup() 
  
# Fit the model
gam2 <- mgcv::gam(
  conc ~ s(date, bs = "cc", m = 2) +
    s(date, location, bs = "fs", m=1) +
    s(date, type, bs = "fs", m=1),
  data = mcye_data,
  family = "gaussian",
  knots = list(0, 365),
  method = "REML"
)

saveRDS(gam2, file = here::here("models/gam2"))

# Model checks
gam2_checks <- gratia::appraise(gam2)

ggsave(gam2_checks,
       filename = here::here("output/gam2_diagnostics.tiff"),
       height = 8, width = 12, dpi = 300
)

# Create DF to back transform model predictions
scaling_mcye <- cleaned_data %>%
  filter(target == "mcyE/ndaF") %>% 
  group_by(location, type) %>% 
  summarise(
    mean = mean(log(conc)),
    sd = sd(log(conc)),
    n = n()
  )

# Generate predictions along finer grid than original data
pred_data_mcye <- with(
  mcye_data,
  expand.grid(
    date = seq(from = min(date), to = max(date), by = 1),
    location = c("Intake", "Inlet", "Outlet"),
    type = c("Passive", "Grab")
  )
)

# Predict and back transform to natural scale
gam2_pred <- predict(gam2, se.fit = TRUE, newdata = pred_data_mcye)

gam2_predictions <- pred_data_mcye %>%
  mutate(
    fit = gam2_pred$fit,
    se = gam2_pred$se.fit
  ) %>% 
  left_join(scaling_mcye, by = c("location", "type")) %>%
  mutate(
    ymin = fit - 2 * se,
    ymax = fit + 2 * se
  ) %>% 
  mutate(across(
    c(fit, ymin, ymax),
    ~ exp(.x * sd + mean)
  )) %>% 
  mutate(date_cont = as.Date(date, origin = "2021-01-01")) %>%
  mutate(type = if_else(type == "Grab", "Grab~(GU~mL^{-1})",
                        "Passive~(GU~cm^{-2})"
  ))

# Generate plots
p3 <- gam2_predictions %>%
  filter(location != "Intake") %>% 
  ggplot(aes(date_cont, fit)) +
  facet_grid(vars(type), vars(location)
  ) +
  # Plot predictions
  geom_line() +
  geom_ribbon(
    aes(ymin = ymin, ymax = ymax),
    alpha = 0.25,
    fill = ggsci::pal_jco("default")(1)
  ) +
  # Add original data
  geom_point(
    data = cleaned_data %>%
      filter(location != "Intake",
             target == "mcyE/ndaF") %>%
      mutate(detect = if_else(conc > lod, "> LOD", "< LOD"),
             detect = factor(detect, levels = c("> LOD", "< LOD"))) %>% 
      mutate(type = if_else(type == "Grab", "Grab~(GU~mL^{-1})",
                            "Passive~(GU~cm^{-2})"
      )),
    aes(date_cont, conc, color = detect), size = 3) +
  scale_color_manual(values = c("black", "darkgray")) +
  
  coord_panel_ranges(panel_ranges = list(
    list(y = c(0, 4e4)),
    list(y = c(0, 4e4)),
    list(y = c(0, 1e4)),
    list(y = c(0, 1e4))
  )) +
  scale_y_continuous(
    labels = function(x) format(x, scientific = TRUE, digits = 2),
    breaks = function(y) {
      c(
        0, 0.25 * max(y), 0.5 * max(y),
        0.75 * max(y), max(y)
      )
    }
  ) +
  theme(strip.text.y = element_blank(),
        legend.position = "none") +
  labs(
    x = NULL,
    y = "Gene Concentration",
    color = NULL
  )

p4 <- gam2_predictions %>%
  filter(location == "Intake", type == "Passive~(GU~cm^{-2})") %>% 
  mutate(location = "bold(Intake)",
         type = if_else(type == "Passive~(GU~cm^{-2})",
                        "bold(Passive~(GU~cm^{-2}))",
                        "bold(Grab~(GU~mL^{-1}))")) %>% 
  ggplot(aes(date_cont, fit)) +
  facet_grid(vars(type), vars(location),
             labeller = label_parsed
  ) +
  # Plot predictions
  geom_line() +
  geom_ribbon(
    aes(ymin = ymin, ymax = ymax),
    alpha = 0.25,
    fill = ggsci::pal_jco("default")(1)
  ) +
  # Add original data
  geom_point(
    data = cleaned_data %>%
      filter(location == "Intake", target == "mcyE/ndaF") %>% 
      mutate(location = "bold(Intake)",
             type = if_else(type == "Passive",
                            "bold(Passive~(GU~cm^{-2}))",
                            "bold(Grab~(GU~mL^{-1}))")),
    aes(date_cont, conc, color = detect), size = 3) +
  
  scale_color_manual(values = c("black", "darkgray")) +
  coord_panel_ranges(panel_ranges = list(
    list(y = c(0, 1e3)),
    list(y = c(0, 6e4))
  )) +
  scale_y_continuous(
    labels = function(x) format(x, scientific = TRUE, digits = 2),
    breaks = function(y) {
      c(
        0, 0.25 * max(y), 0.5 * max(y),
        0.75 * max(y), max(y)
      )
    }
  ) + 
  labs(
    x = NULL,
    y = NULL,
    color = NULL
  )

p3 + p4 + plot_layout(widths = c(2, 1))

ggsave(here::here("output/figure_4.tiff"), height = 8, width = 12, dpi = 300)
