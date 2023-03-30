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
              legend.title = element_text(size = 16),
              legend.text = element_text(size = 14),
              strip.text = element_text(size = 14, face = "bold"),
              strip.background = element_rect(fill = "white"),
              panel.spacing.x = unit(0.5, "lines"),
              panel.spacing.y = unit(1.5, "lines")
            ))

cleaned_data <- read_csv(here::here("cleaned_data/plot_data.csv")) %>% 
  mutate(date = lubridate::yday(date)) %>%
  mutate(location = factor(location,
                           levels = c("Inlet", "Outlet", "Intake")
  )) %>%
  mutate(type = factor(type),
         date_cont = as.Date(date, origin = "2021-01-01"))

# GAM fit of 16s data

# Prep data for mgcv
test_data <- cleaned_data %>%
  filter(conc > lod, 
         target == "16S rRNA",
         frequency == "Weekly") %>%
  group_by(location, type) %>%
  mutate(conc = scale(log(conc))) %>%
  ungroup() 
  
# Fit the model
gam1 <- mgcv::gam(
  conc ~ s(date, bs = "cc", m = 2) +
    s(date, location, bs = "fs", m=1) +
    s(date, type, bs = "fs", m=1),
  data = test_data,
  family = "gaussian",
  knots = list(0, 365),
  method = "REML"
)

saveRDS(gam1, file = here::here("models/gam1"))

# Model checks
par(mfrow = c(2, 2))
mgcv::gam.check(gam1)
par(mfrow = c(1, 1))

# Create DF to back transform model predictions
scaling_16s <- cleaned_data %>%
  filter(conc > lod, target == "16S rRNA") %>%
  group_by(location, type) %>%
  summarise(
    mean = mean(log(conc)),
    sd = sd(log(conc)),
    n = n()
  )

# Generate predictions along finer grid than original data
pred_data <- with(
  test_data,
  expand.grid(
    date = seq(from = min(date), to = max(date), by = 1),
    location = c("Intake", "Inlet", "Outlet"),
    type = c("Passive", "Grab")
  )
)

# Predict and back transform to natural scale
gam1_pred <- predict(gam1, se.fit = TRUE, newdata = pred_data)

gam1_predictions <- pred_data %>%
  mutate(
    fit = gam1_pred$fit,
    se = gam1_pred$se.fit
  ) %>%
  left_join(scaling_16s, by = c("location", "type")) %>%
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

p1 <- gam1_predictions %>%
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
             target == "16S rRNA") %>%
      mutate(type = if_else(type == "Grab", "Grab~(GU~mL^{-1})",
                            "Passive~(GU~cm^{-2})")),
    aes(date_cont, conc),
    size = 3
  ) +
  coord_panel_ranges(panel_ranges = list(
    list(y = c(0, 6e5)),
    list(y = c(0, 6e5)),
    list(y = c(0, 4e6)),
    list(y = c(0, 4e6))
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
  theme(strip.text.y = element_blank()) +
  labs(
    x = NULL,
    y = "Gene Concentration"
  )

p2 <- gam1_predictions %>% 
  filter(location == "Intake") %>% 
  mutate(type = if_else(type == "Grab~(GU~mL^{-1})", 
                        "bold(Grab~(GU~mL^{-1}))",
                        "bold(Passive~(GU~cm^{-2}))"),
         location = if_else(location == "Intake", "bold(Intake)", as.character(location))) %>% 
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
      filter(location == "Intake",
             target == "16S rRNA") %>%
      mutate(type = if_else(type == "Grab", "bold(Grab~(GU~mL^{-1}))",
                            "bold(Passive~(GU~cm^{-2}))"),
             location = if_else(location == "Intake", "bold(Intake)", as.character(location))),
    aes(date_cont, conc),
    size = 3
  ) +
  coord_panel_ranges(panel_ranges = list(
    list(y = c(0, 2.5e5)),
    list(y = c(0, 4e6))
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
  theme(strip.text = element_text(face = "bold")) +
  labs(
    x = NULL,
    y = NULL
  )

p1 + p2 + plot_layout(widths = c(2, 1))

ggsave(here::here("output/figure_2.tiff"), height = 8, width = 12, dpi = 300)
