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

library(tidyverse)
library(RColorBrewer)

df <- read_csv("model-out.csv") %>%
  mutate(Var = factor(Var, levels  = c("Death", "Reported case", "Infection")))
vGeo2 <- read_csv("input/geo-sel.csv") %>% pull(Geo)

dfLastDate <- df %>% arrange(Geo, Var, Date) %>%
  filter(!is.na(New), !is.na(NewEst), Var == "Reported case", substr(Geo,1,5) != "China",
         substr(Geo, 1, 2) != "US" | Geo == "US - 12 states total") %>%
  group_by(Geo) %>% summarize(Date = last(Date))
dfCo <- bind_rows(
  dfLastDate %>% inner_join(df %>% filter(Var == "Reported case") %>% select(Geo, Date, Var, Cum)),
  dfLastDate %>% inner_join(df %>% filter(Var == "Infection") %>% select(Geo, Date, Var, Cum = CumEst))
  ) %>%
  mutate(Var = fct_recode(Var, `Reported cases` = "Reported case", `Estimated infections` = "Infection"))
vGeo <- dfCo %>% filter(Var == "Estimated infections") %>% top_n(8, Cum) %>% arrange(-Cum) %>% pull(Geo) 

# Reported cases vs actual infections
dfCo %>% filter(Geo %in% vGeo) %>% ggplot(aes(x = Var, y = Cum, fill = reorder(Geo, Cum))) + geom_col() +
  scale_y_continuous(labels = scales::comma) + 
  scale_fill_manual(values = brewer.pal(8,"Set1")) +
  ylab(element_blank()) + xlab(element_blank()) +
  theme(legend.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("output/fig-est-rep.png", width=4, height=5)
dfCo %>% group_by(Var) %>% summarize(sum(Cum))

# Projected cumulative infections
df %>% filter(Geo %in% vGeo, Var == "Infection", between(Date, as.Date("2020-03-01"), as.Date("2020-06-01") )) %>%
  mutate(Geo = fct_rev(factor(Geo, levels = vGeo))) %>%
  ggplot(aes(x=Date, y=CumEst, fill=Geo)) +
  geom_col(width = 1) +
  scale_y_continuous(labels = scales::comma) + 
  scale_fill_manual(values = brewer.pal(8,"Set1"), guide = guide_legend(ncol=2)) +
  ylab(element_blank()) + xlab(element_blank()) +
  theme(legend.title = element_blank(), legend.position = "bottom") +
  ggtitle("Estimated cumulative infections")
ggsave("output/fig-cum-infections.png", width=3.5, height=4)

df %>% filter(Geo %in% vGeo, Date %in% c(as.Date("2020-04-29"), as.Date("2020-05-31")), Var == "Infection") %>%
  select(Geo, Date, CumEst, CumLow, CumHigh) %>% select(-Geo) %>% group_by(Date) %>% summarise_all(sum)
df %>% filter(Geo %in% vGeo, Date %in% c(as.Date("2020-04-27")), Var == "Reported case") %>%
  select(Geo, Date, Cum) %>% select(-Geo) %>% group_by(Date) %>% summarise_all(sum)



# 'Good' geographies
Date1 <- as.Date("2020-03-01")
Date2 <- as.Date("2020-05-01")
vGeo3 <- c("Austria", "US - Florida", "US - Louisiana", "Switzerland", "US - Michigan", "Spain", "Norway", "Germany")
df %>% filter(Date >= as.Date(Date1), Date < as.Date(Date2), Geo %in% vGeo3) %>%
  ggplot(aes(x=Date, y=NewEst, color=Var)) +
  geom_line(aes(y=NewLow), linetype = 2, size = 0.25) + geom_line(aes(y=NewHigh), linetype = 2, size = 0.25) +
  geom_point(aes(y=New), size=0.25) + geom_line(size = 0.25) + facet_wrap(~Geo, ncol = 2) +
  scale_x_date(breaks = "1 month", labels = function(x) substr(format(x, format="%d %b"),2,10)) +
  scale_y_continuous(labels = scales::comma, trans="log10", limits = c(1,3e5)) +
  scale_color_manual(values = brewer.pal(3,"Set1"), guide = guide_legend(reverse = TRUE)) +
  xlab(element_blank()) + ylab("New events per day (logarithmic scale)") +
  theme(legend.title=element_blank(), legend.position="top", axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("output/fig-good.png", width=3.5, height=8, dpi = 600)


# Georgia & Texas impact of lifting orders
require(rstan)
envBase <- envLift <- new.env()
load("output/fit-model-200430-081533.RData", envir = envBase) # With current restrictions in place
load("output/fit-model-200430-081419.RData", envir = envLift) # With restrictions lifted
dfBaseRaw <- with(envBase,
               bind_rows(
                 expand.grid(iter = 1:nIter, Geo = vGeo, Day2=1:nTPred) %>% as_tibble() %>% mutate(Var = "Infection", Log = as.vector(rstan::extract(fit)$lirPred)),
                 expand.grid(iter = 1:nIter, Var = c("Death", "Case"), Geo = vGeo, Day2=1:nTPred) %>% as_tibble() %>% mutate(Log = as.vector(rstan::extract(fit)$lrrPred))
               ) %>% full_join(dfDates %>% select(Geo, Tmax, End))
               %>% mutate(Day = Tmax + Day2) %>% select(-Day2) %>%
                 filter(Geo %in% c("US - Texas", "US - Georgia")))
dfLiftRaw <- with(envLift,
                  bind_rows(
                    expand.grid(iter = 1:nIter, Geo = vGeo, Day2=1:nTPred) %>% as_tibble() %>% mutate(Var = "Infection", Log = as.vector(rstan::extract(fit)$lirPred)),
                    expand.grid(iter = 1:nIter, Var = c("Death", "Case"), Geo = vGeo, Day2=1:nTPred) %>% as_tibble() %>% mutate(Log = as.vector(rstan::extract(fit)$lrrPred))
                  ) %>% full_join(dfDates %>% select(Geo, Tmax, End))
                  %>% mutate(Day = Tmax + Day2) %>% select(-Day2) %>%
                    filter(Geo %in% c("US - Texas", "US - Georgia")))
dfBaseLift <- bind_rows(
  dfBaseRaw %>% group_by(Geo, Var, iter) %>% summarize(Cum = sum(exp(Log))) %>%
    group_by(Geo, Var) %>% summarize(CumEst = median(Cum), CumLow = quantile(Cum, probs=0.025), CumHigh = quantile(Cum, probs=0.975)) %>%
    ungroup() %>% mutate(Case = "Old"),
  dfLiftRaw %>% group_by(Geo, Var, iter) %>% summarize(Cum = sum(exp(Log))) %>%
    group_by(Geo, Var) %>% summarize(CumEst = median(Cum), CumLow = quantile(Cum, probs=0.025), CumHigh = quantile(Cum, probs=0.975)) %>%
    ungroup() %>% mutate(Case = "New"))
dfBaseLift %>% filter(Var == "Death") %>% mutate(Case = factor(Case, levels = c("Old", "New"), labels = c("Stay home order", "Eased policy"))) %>%
  ggplot(aes(x = Case, y = CumEst, fill = Geo)) +
  geom_col() +
  scale_y_continuous(labels = scales::comma) + ylab(element_blank()) + xlab(element_blank()) +
  scale_fill_manual(values = brewer.pal(3,"Set1"), guide = guide_legend(position = "bottom", title = element_blank(), ncol=1)) +
  ggtitle("Projected number of deaths", subtitle = "Total in May and June") +
  theme(legend.position = "bottom")
ggsave("output/fig-texas-georgia.png", width = 3.5, height=4, dpi = 600)
dfBaseLift %>% filter(Var == "Death")
with(
  dfBaseLift %>% filter(Var == "Death") %>% group_by(Case) %>% summarize(CumEst = sum(CumEst), CumLow = sum(CumLow)),
  c(CumEst[Case == "New"] - CumEst[Case == "Old"], CumLow[Case == "New"] - CumLow[Case == "Old"]))
