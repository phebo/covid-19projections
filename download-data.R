# Copyright (C) 2020, Phebo Wibbens, Wesley Koo, and Anita McGahan

#   This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.


# This script downloads data from Oxford, Johns Hopkins, and The Economist
# The updated data are written to the input/ directory

library(tidyverse)
outlExDeath <- 0.2 # Threshold for removing data points with very low recorded total deaths compared to expected deaths

dfOxG <- read_csv("https://github.com/OxCGRT/covid-policy-tracker/blob/master/data/OxCGRT_latest.csv?raw=true",
                  col_types = list(RegionName = col_character(), RegionCode = col_character()))
dfJhCaseUS <- read_csv("https://github.com/CSSEGISandData/COVID-19/blob/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv?raw=true")
dfJhCaseG <- read_csv("https://github.com/CSSEGISandData/COVID-19/blob/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv?raw=true")
dfJhDeathUS <- read_csv("https://github.com/CSSEGISandData/COVID-19/blob/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_US.csv?raw=true")
dfJhDeathG <- read_csv("https://github.com/CSSEGISandData/COVID-19/blob/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv?raw=true")
dfEconEd <- read_csv("https://github.com/TheEconomist/covid-19-excess-deaths-tracker/blob/master/output-data/excess-deaths/all_weekly_excess_deaths.csv?raw=true")

# Reformat Johns-Hopkins data
dfJhUS <- bind_rows(dfJhCaseUS %>% mutate(var = "case"),
                    dfJhDeathUS %>% mutate(var = "death")) %>%
  mutate(geo = paste("US", Province_State , sep = " - ")) %>%
  select(-(UID:Combined_Key))
dfJhG <- bind_rows(dfJhCaseG %>% mutate(var = "case"),
                   dfJhDeathG %>% mutate(var = "death")) %>%
  mutate(geo = `Country/Region`) %>%
  select(-(`Province/State`:Long)) 
dfJh <- bind_rows(dfJhUS, dfJhG) %>%
  group_by(geo, var) %>% summarise_all(sum) %>%
  pivot_longer(-c(geo,var), names_to = "date", values_to = "cum") %>%
  mutate(date = as.Date(date, format = "%m/%d/%y")) %>%
  filter(!is.na(date), !is.na(cum))
write_csv(dfJh, "input/jh-database.csv")

# Reformat Economist mortality data
dfEcon <- dfEconEd %>%
  filter(country == region | country == "United States", region != "United States" ) %>%
  mutate(
    geo = case_when(
      country == "United States" ~ paste("US", region, sep = " - "),
      country == "Britain" ~ "United Kingdom",
      TRUE ~ country),
    date = as.Date(week * 7, origin="2019-12-29"),
    total_deaths = round(total_deaths))
stopifnot(all(dfEcon$year == 2020))
dfPop <- dfEcon %>% filter(year == 2020, week == 1) %>% select(geo, population)
dfEcon <- dfEcon %>%
  select(geo, date, deathTot = total_deaths, deathExp = expected_deaths) %>%
  filter(deathTot > deathExp * (1-outlExDeath))
write_csv(dfEcon, "input/econ-database.csv")
write_csv(dfPop, "input/econ-population.csv")

# Reformat Oxford policy data
dfOx <- bind_rows(
  dfOxG %>% filter(is.na(RegionName)) %>%
    select(geo = CountryName, geoCode = CountryCode, date = Date, matches("^C\\d_"), matches("^H[1-3]_")) %>%
    filter(geo != "United States"),
  dfOxG %>% filter(CountryName == "United States", !is.na(RegionName)) %>%
    select(geo = RegionName, geoCode = RegionCode, date = Date, matches("^C\\d_"), matches("^H[1-3]_")) %>%
    mutate(geo = paste("US", geo, sep = " - "), C5_Flag = as.numeric(C5_Flag)))
dfOx <-
  left_join(
    dfOx %>% select(-ends_with("Flag")) %>%
      pivot_longer(-(geo:date), values_to = "level") %>%
      separate(name, c("polCode", "polName"), sep = "_"),
    dfOx %>% select(geo:date, ends_with("Flag")) %>%
      pivot_longer(-(geo:date), names_to = "polCode", values_to = "flag") %>%
      mutate(polCode = substr(polCode, 1, 2))) %>%
  filter(!is.na(level)) %>%
  arrange(geo, polCode, date) %>% distinct(geo, polCode, level, flag, .keep_all = T)  %>%
  mutate(date = as.Date(as.character(date), format = "%Y%m%d"))
write_csv(dfOx, "input/oxford-policy.csv")
