rm(list = ls())
library(tidyverse)
library(scales)

adsorption_8hr <- readxl::read_xlsx(here::here("raw_data/8hrAdsorptionDetects_16S.xlsx")) %>%
  rename_with(
    .cols = everything(),
    ~ tolower(gsub(" ", "_", .x))
  ) %>%
  mutate(
    across(
      where(is.character), ~ tolower(.)
    )
  ) %>%
  mutate(
    gu_conc =
      case_when(
        type == "volume" ~ vol,
        TRUE ~ area
      )
  ) %>%
  select(-vol, -area) %>%
  mutate(
    area_cm =
      case_when(
        sample %in% c(
          "acrylic copolymer", "cellulose nitrate",
          "nylon", "gauze"
        ) ~ 4.9087
      )
  ) %>%
  mutate(
    vol_ml =
      case_when(
        !sample %in% c(
          "acrylic copolymer", "nylon",
          "cellulose nitrate", "gauze"
        ) ~ 20
      )
  ) %>%
  mutate(lod = case_when(
    sample %in% c(
      "acrylic copolymer", "nylon",
      "cellulose nitrate", "gauze"
    ) ~ 9 * 500 / area_cm,
    !sample %in% c(
      "acrylic copolymer", "nylon",
      "cellulose nitrate", "gauze"
    ) ~ 9 * 500 / vol_ml
  ))



adsorption_24hr <- readxl::read_xlsx(here::here("raw_data/24hrAdsorptionDetects_16S.xlsx")) %>%
  rename_with(
    .cols = everything(),
    ~ tolower(gsub(" ", "_", .x))
  ) %>%
  mutate(
    across(
      where(is.character), ~ tolower(.)
    )
  ) %>%
  mutate(
    gu_conc =
      case_when(
        type == "volume" ~ vol,
        TRUE ~ area
      )
  ) %>%
  mutate(
    area_cm =
      case_when(
        sample %in% c("acrylic copolymer", "nylon", "cellulose nitrate", "gauze") ~ 4.9087
      )
  ) %>%
  mutate(
    vol_ml =
      case_when(
        !sample %in% c(
          "acrylic copolymer", "nylon",
          "cellulose nitrate", "gauze"
        ) ~ 20
      )
  ) %>%
  mutate(lod = case_when(
    sample %in% c(
      "acrylic copolymer", "nylon",
      "cellulose nitrate", "gauze"
    ) ~ 9 * 500 / area_cm,
    !sample %in% c(
      "acrylic copolymer", "nylon",
      "cellulose nitrate", "gauze"
    ) ~ 9 * 500 / vol_ml
  ))

# Save the cleaned files

write_csv(adsorption_24hr, here::here("cleaned_data/adsorption_24hr.csv"))
write_csv(adsorption_8hr, here::here("cleaned_data/adsorption_8hr.csv"))






