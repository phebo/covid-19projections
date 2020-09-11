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


#### IMPORTANT NOTES ####

# This script estimates Covid-19 projections using a Bayesian model. It requires:
# 1. Data from Johns Hopkins, Oxford, and The economist in this repository; they can be updated using the "download-data.R" script
# 2. The Bayesian estimation package Stan, see: https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started

# It takes several hours to estimate the model on a regular PC

library(RColorBrewer)
library(stats4)
library(tidyverse)
library(rstan)
source("functions.R")

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
time.now <- format(Sys.time(), format='%y%m%d-%H%M%S')
print(time.now)
suppressWarnings(dir.create(file.path("output")))

vGeo <- c("Canada", "Australia", "South Korea", "Brazil", "China", "India", "Singapore", "New Zealand", "Iceland",
          "Russia", "Indonesia", "Pakistan", "Nigeria", "Bangladesh", "Mexico", "Japan", "Taiwan")
dates = c(as.Date("2020-02-01"), Inf)

#### Read and clean data ####

dfJh <- read_csv("input/jh-database.csv")
dfOx <- read_csv("input/oxford-policy.csv")

vGeo <- sort(vGeo)
vDate <- sort(unique(dfJh$date))
vDate <- vDate[as.numeric(vDate) %% 7 == as.numeric(min(dfEcon$date)) %% 7]
vDate <- vDate[vDate >= dates[1] & vDate <= dates[2]]
stopifnot(as.numeric(vDate - lag(vDate))[-1] == 7)

# Epidemiology data
dfE <- dfJh %>% mutate(geo = case_when(geo == "Korea, South" ~ "South Korea",
                                       geo == "Taiwan*" ~ "Taiwan",
                                       TRUE ~ geo)) %>%
  filter(geo %in% vGeo, date %in% vDate) %>%
  arrange(geo, var, date) %>%
  group_by(geo, var) %>% mutate(new = cum - lag(cum, default = 0)) %>% ungroup() %>%
  select(-cum) %>% pivot_wider(names_from = var, values_from = new)

# Policy data; For policy definitions, see https://github.com/OxCGRT/covid-policy-tracker/blob/master/documentation/codebook.md
dfOxSub <- dfOx %>% filter(geo %in% vGeo)
dfPolSel <- dfOxSub %>% distinct(polCode, level) %>% arrange(polCode, level) %>%
  filter(level > 0) %>% rename(polCodeSel = polCode, levelSel = level)
dfP <-
  bind_rows(
    dfPolSel %>% rowwise() %>% mutate(data = list(dfOxSub %>% filter(polCode == polCodeSel, level >= levelSel)), value = 1),
    dfPolSel %>% rowwise() %>% mutate(data = list(dfOxSub %>% filter(polCode == polCodeSel, level < levelSel)), value = 0)) %>%
  unnest(data) %>% mutate(pol = paste(polCode, polName, levelSel, sep = " - ")) %>% select(geo, date, pol, value) %>%
  distinct()
vPolAll <- sort(unique(dfP$pol))
dfP <- full_join(expand_grid(geo = vGeo, pol = vPolAll, date = vDate), dfP) %>% arrange(geo, pol, date) %>%
  group_by(geo, pol) %>% fill(value) %>% ungroup() %>% filter(date %in% vDate) %>%
  mutate(pol2 = pol) %>% separate(pol2, c("polCode", "polName", "level"), sep = " - ") %>%
  mutate(polS = paste(polCode, level, sep = " - "),
         value = ifelse(is.na(value), 0, value))

fNew <- dfE %>% pivot_longer(c(case, death)) %>% filter(value > 0) %>%
  ggplot(aes(x = date, y = value, color = fct_rev(name))) +
  geom_point(size = 0.4) +
  facet_wrap(~ geo, ncol = 5) + ggtitle("New reported events per week (log scale)") + 
  xlab(element_blank()) + ylab(element_blank()) + theme(legend.title = element_blank()) +
  scale_y_continuous(labels = scales::comma, trans="log10", limits = c(1, 1e6))

fPol2 <- dfP %>% mutate(pol = paste(polCode, polName, sep = " - ")) %>% filter(value == 1) %>% group_by(geo, pol, date) %>%
  summarize(level = max(level)) %>%
  ggplot(aes(x = date, y = fct_rev(geo), fill = level)) + geom_tile(height=0.8) +
  facet_grid(~ pol, labeller = label_wrap_gen(8)) +
  scale_fill_brewer(palette = "YlOrRd", name="Policy\nlevel") + 
  ylab(element_blank()) + xlab(element_blank())

pdf(paste0("output/chart-addl.pdf"), width=12, height=8, onefile=T)
  print(fPol2)
  print(fNew)
dev.off()


