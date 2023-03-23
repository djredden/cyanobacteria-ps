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
              strip.text.x = element_text(size = 14, face = "bold.italic")
              strip.background = element_rect(fill = "white"),
              panel.spacing.x = unit(0.5, "lines"),
              panel.spacing.y = unit(1, "lines"))
))

raw_data <- readxl::read_xlsx(here::here("raw_data/FletchersLakeDetections.xlsx")) %>%
  mutate(LOD = case_when(
    Location == "Intake" & Type == "Grab" ~ 45,
    Location %in% c("Inlet", "Outlet") & Type == "Grab" ~ 180,
    Type == "Passive" ~ 260
  )) %>%
  # Toni removed unpaired passive samples (she didn't take grabs until July 13)
  filter(Date > "2021-07-10")

# Time series of all toxin data

time_series_facets <- raw_data %>%
  filter(Assay == "Toxins") %>%
  mutate(hline = case_when(
    Type == "Grab" & Location == "Intake" ~ 180,
    Type == "Grab" & Location != "Intake" ~ 45,
    Type == "Passive" ~ 260,
    TRUE ~ 0
  )) %>%
  mutate(Type = if_else(Type == "Grab", "Grab~(GU~mL^{-1})",
                        "Passive~(GU~cm^{-2})"
  )) %>%
  ggplot(aes(x = Date, y = Conc, fill = Location, shape = Frequency)) +
  geom_point(position = "jitter", size = 3) +
  geom_hline(aes(yintercept = hline, linetype = factor(hline))) +
  scale_shape_manual(values = c(24, 21)) +
  ggsci::scale_fill_jco() +
  guides(
    fill = guide_legend(override.aes = list(shape = 22, size = 3)),
    shape = guide_legend(override.aes = list(fill = "black", size = 3))
  ) +
  facet_grid(Type ~ factor(Target, levels = c("mcyE/ndaF", "cyrA", "sxtA")),
             labeller = label_parsed
  ) +
  ylab(NULL) +
  scale_y_log10(
    limits = c(1, 1000000),
    labels = trans_format("log10", math_format(10^.x)),
    expand = c(0, 0)
  ) +
  labs(
    y = "Gene Concentration",
    x = "",
    fill = "Location",
    shape = "Sample Frequency",
    linetype = "LOD"
  ) 

ggsave(time_series_facets,
       filename = here::here("output/time_series_toxins.tiff"),
       height = 8, width = 12, dpi = 300
)

# Look at correlations in the grab vs passive time series
raw_data %>%
  filter(Conc > LOD, Assay == "16S rRNA") %>%
  group_by(Type, Location) %>%
  mutate(week = lubridate::week(Date)) %>% # sample dates were often off by a day or two at the WTP
  select(-Date) %>%
  pivot_wider(id_cols = c(week, Location), names_from = Type, values_from = Conc, values_fn = mean) %>%
  ggplot(aes(Grab, Passive)) +
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x) +
  stat_cor(aes(label = paste(..rr.label..)), # adds R^2 value
           r.accuracy = 0.01
  ) +
  facet_wrap(~Location, scales = "free")

# GAM fit of 16s data

# Prep data for mgcv
test_data <- raw_data %>%
  filter(Conc > LOD, Assay == "16S rRNA") %>%
  group_by(Location, Type) %>%
  mutate(Conc = scale(Conc)) %>%
  ungroup() %>%
  mutate(Date = lubridate::yday(Date)) %>%
  mutate(Location = factor(Location,
                           levels = c("Inlet", "Outlet", "Intake")
  )) %>%
  mutate(Type = factor(Type))

# Fit the model
gam1 <- mgcv::gam(
  Conc ~ s(Date, bs = "cc", m = 2) +
    s(Date, by = Location, m=1, bs = "tp") +
    s(Date, by = Type, m=1,  bs = "tp"),
  data = test_data,
  family = "gaussian",
  knots = list(0, 365)
)

# Generate predictions and pair with original data
gam1_pred <- predict(gam1, se.fit = TRUE)
gam1_data <- transform(test_data,
                       fit = gam1_pred$fit,
                       se = gam1_pred$se.fit
)

scaling_16s <- raw_data %>%
  filter(Conc > LOD, Assay == "16S rRNA") %>%
  group_by(Location, Type) %>%
  summarise(
    mean = mean(Conc),
    sd = sd(Conc),
    n = n()
  )

plot1_data <- left_join(gam1_data, scaling_16s, by = c("Location", "Type")) %>%
  mutate(across(c(Conc, fit), ~ .x * sd + mean),
         se = (se * sqrt(n) * sd + mean) / sqrt(n)
  ) %>%
  mutate(mean  = exp(mean),
         se = se) %>%
  mutate(Type = if_else(Type == "Grab", "Grab~(GU~mL^{-1})",
                        "Passive~(GU~cm^{-2})"
  )) %>%
  mutate(date_cont = as.Date(Date, origin = "2021-01-01"))


# GAM plot for 16s
p1 <- ggplot(
  data = plot1_data %>%
    filter(Location != "Intake"),
  aes(x = date_cont, y = Conc)) +
  facet_grid(Type ~ Location,
             labeller = label_parsed,
             scales = "free_y") +
  geom_ribbon(
    aes(
      ymin = fit - 2 * se,
      ymax = fit + 2 * se),
    alpha = 0.25,
    fill = ggsci::pal_jco("default")(1)) +
  geom_line(aes(y = fit)) +
  geom_point(size =3) +
  theme(
    strip.text.y = element_blank()
  ) +
  coord_panel_ranges(panel_ranges = list(
    list(y = c(0, 6e5)),
    list(y = c(0, 6e5)),
    list(y = c(0, 4e6)),
    list(y = c(0, 4e6))
  )) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE, digits=2),
                     breaks = function(y) c(0, 0.25*max(y), 0.5*max(y),
                                            0.75*max(y), max(y))) +
  labs(
    x = NULL,
    y = "Gene Concentration"
  )


p2 <- ggplot(
  data = plot1_data %>%
    filter(Location == "Intake"),
  aes(x = date_cont, y = Conc)
) +
  facet_grid(Type ~ Location,
             labeller = label_parsed,
             scales = "free_y"
  ) +
  geom_ribbon(
    aes(
      ymin = fit - 2 * se,
      ymax = fit + 2 * se
    ),
    alpha = 0.25,
    fill = ggsci::pal_jco("default")(1)
  ) +
  geom_line(aes(y = fit)) +
  geom_point(size=3) +
  coord_panel_ranges(panel_ranges = list(
    list(y = c(0, 6e7)),
    list(y = c(0, 4e6))
  )) +
  
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE, digits=2),
                     breaks = function(y) c(0, 0.25*max(y), 0.5*max(y),
                                            0.75*max(y), max(y))) +
  labs(
    x = NULL,
    y = NULL)

p1 + p2 + plot_layout(widths = c(2, 1))

ggsave(here::here("output/16s_gam.tiff"), height = 8, width = 12, dpi = 300)


### GAM for mcyE
mcye_data <- raw_data %>%
  filter(Target == "mcyE/ndaF") %>%
  mutate(detect = if_else(Conc > LOD, "> LOD", "< LOD")) %>%
  mutate(detect = factor(detect, levels = c("> LOD", "< LOD"))) %>%
  group_by(Location, Type) %>%
  mutate(Conc = scale(log(Conc))) %>%
  ungroup() %>%
  mutate(Date = lubridate::yday(Date)) %>%
  mutate(Location = factor(Location,
                           levels = c("Inlet", "Outlet", "Intake")
  )) %>%
  mutate(Type = factor(Type))


# Fit the model
gam2 <- mgcv::gam(
  Conc ~ s(Date, bs = "cc", m = 2) +
    s(Date, by = Location, bs = "re") +
    s(Date, by = Type, bs = "re"),
  data = mcye_data,
  family = "gaussian",
  knots = list(0, 365)
)

# Model checks
par(mfrow = c(2,2))
mgcv::gam.check(gam2)
par(mfrow = c(1,1))

# Generate predictions and pair with original data
gam2_pred <- predict(gam2, se.fit = TRUE)
gam2_data <- transform(mcye_data,
                       fit = gam2_pred$fit,
                       se = gam2_pred$se.fit
)

scaling_mcye <- raw_data %>%
  filter(Target == "mcyE/ndaF") %>%
  group_by(Location, Type) %>%
  summarise(
    mean = mean(Conc),
    sd = sd(Conc),
    n = n()
  )

plot2_data <- left_join(gam2_data, scaling_mcye, by = c("Location", "Type")) %>%
  mutate(across(c(Conc, fit), ~ .x * sd + mean),
         se = (se * sqrt(n) * sd + mean) / sqrt(n)
  ) %>%
  mutate(Type = if_else(Type == "Grab", "Grab~(GU~mL^{-1})",
                        "Passive~(GU~cm^{-2})"
  )) %>%
  mutate(date_cont = as.Date(Date, origin = "2021-01-01"))


# GAM plot for mcyE/ndaF
p3 <- ggplot(
  data = plot2_data %>%
    filter(Location != "Intake"),
  aes(x = date_cont, y = Conc)) +
  geom_point(aes(color = detect), size=3) +
  scale_color_manual(values = c("black", "darkgray")) +
  facet_grid(Type ~ Location,
             labeller = label_parsed) +
  geom_ribbon(
    aes(
      ymin = fit - 2 * se,
      ymax = fit + 2 * se
    ),
    alpha = 0.25,
    fill = ggsci::pal_jco("default")(1)
  ) +
  geom_line(aes(y = fit)) +
  theme(
    strip.text.y = element_blank(),
    legend.position = "none"
  ) +
  coord_panel_ranges(panel_ranges = list(
    list(y = c(-0.5e4, 5.5e4)),
    list(y = c(-0.5e4, 5.5e4)),
    list(y = c(-0.5e4, 5.5e4)),
    list(y = c(-0.5e4, 5.5e4))
  )) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE, digits=2)) +
  labs(
    x = NULL,
    y = "Gene Concentration",
    color = NULL
  )

p4 <- ggplot(
  data = plot2_data %>%
    filter(Location == "Intake"),
  aes(x = date_cont, y = Conc)
) +
  facet_grid(Type ~ Location,
             labeller = label_parsed
  ) +
  geom_ribbon(
    aes(
      ymin = fit - 2 * se,
      ymax = fit + 2 * se
    ),
    alpha = 0.25,
    fill = ggsci::pal_jco("default")(1)
  ) +
  geom_line(aes(y = fit)) +
  geom_point(aes(color = detect), size = 3) +
  scale_color_manual(values = c("black", "darkgray")) +
  coord_panel_ranges(panel_ranges = list(
    list(y = c(-1.0e4, 1.5e5)),
    list(y = c(-0.5e4, 5.5e4))
  )) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE, digits=2)) +
  labs(
    x = NULL,
    y = NULL,
    color = "Detection"
  )

p3 + p4 + plot_layout(widths = c(2, 1))

ggsave(here::here("output/mcye_gam.tiff"), height = 8, width = 12, dpi = 300)

par(mfrow= c(2,2))
mgcv::gam.check(gam1)
par(mfrow = c(1,1))
