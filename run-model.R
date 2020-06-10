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
# 1. Johns Hopkins data from Github, see: https://github.com/CSSEGISandData/COVID-19
# 2. The Bayesian estimation package Stan, see: https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started

# It takes one to several hours to estimate the model projections on a regular PC
# For further information see: https://www.covid-19projections.com/
# The first time the model is run useInit should be set to FALSE. It can be set to TRUE in next runs to speed up the model estimation.


library(tidyverse)
library(rstan)
library(RColorBrewer)
library(xtable)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')  # Recommended on a Windows machine
time.now <- format(Sys.time(), format='%y%m%d-%H%M%S')
print(time.now)
set.seed(654321) # To make Bayesian estimates reproducable

filetype <- ".png" # Filetype for output charts
test <- F # Do a test run (much faster)?
useInit <- F # Use initial value from previous run? Needs to be FALSE first time model is used

# Minimum and maximum lags between infection and resp. reported case and reported death:
if(test) {
  lagDeath <- c(15,25)
  lagCase <- c(8,17)
} else {
  lagDeath <- c(10,30)
  lagCase <- c(5,20)
}
PolSel <- c(1,4,5.1,5.2,7,8,9) # Which policies to take into account in analysis
#PolSel <- c(1,7,8,9) # Which policies to take into account in analysis
nTPred <- 60 # Additional prediction days
mortality <- 0.01; mortSig <- 0.5; # Parameters for log-normal distribution; mortSig = 0.5 means 95% interval ~0.3%-2.6%
lr2 <- log(1e3); phi2 <- 1e-4 # Parameters for outlier negative binomial

dfJh1 <- read_csv("../COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv")
dfJh2 <- read_csv("../COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv")
dfJh3 <- read_csv("../COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_us.csv")
dfJh4 <- read_csv("../COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_us.csv")
dfDates <- read_csv("input/data-start-end.csv")
dfPRaw <- read_csv("input/covid-policy.csv")
dfPNames <- read_csv("input/covid-policy-names.csv") %>% select(-Description)
vGeo2 <- read_csv("input/geo-sel.csv") %>% pull(Geo)
dir.create("output", showWarnings = FALSE)

#### Tidy data ####
dfE1 <- dfJh1 %>% select(-c(Lat,Long)) %>%
  pivot_longer(-c(`Province/State`,`Country/Region`), names_to = "Date", values_to = "Cum") %>%
  mutate(Var = "Death")
dfE2 <- dfJh2 %>% select(-c(Lat,Long)) %>%
  pivot_longer(-c(`Province/State`,`Country/Region`), names_to = "Date", values_to = "Cum") %>%
  mutate(Var = "Case")
dfE3 <- dfJh3 %>% select(-c(UID:Admin2,Country_Region:Population)) %>%
  pivot_longer(-Province_State, names_to = "Date", values_to = "Cum") %>%
  group_by(Province_State, Date) %>% summarize(Cum = sum(Cum)) %>% ungroup() %>%
  mutate(Var = "Death", Date = as.Date(Date, format = "%m/%d/%y"), Geo = paste("US -", Province_State)) %>%
  select(-Province_State)
dfE4 <- dfJh4 %>% select(-c(UID:Admin2,Country_Region:Combined_Key)) %>%
  pivot_longer(-Province_State, names_to = "Date", values_to = "Cum") %>%
  group_by(Province_State, Date) %>% summarize(Cum = sum(Cum)) %>% ungroup() %>%
  mutate(Var = "Case", Date = as.Date(Date, format = "%m/%d/%y"), Geo = paste("US -", Province_State)) %>%
  select(-Province_State)

dfE <- bind_rows(dfE1, dfE2) %>%
  mutate(
    Date = as.Date(Date, format = "%m/%d/%y"),
    Geo = case_when(
      `Country/Region` == "China" & `Province/State` == "Hubei" ~ "China - Hubei",
      `Country/Region` == "China" & `Province/State` == "Henan" ~ "China - Henan",
      `Country/Region` == "China" ~ "China - Other mainland",
      `Country/Region` == "Canada" & `Province/State` == "Ontario" ~ "Canada - Ontario",
      `Country/Region` == "Canada" & `Province/State` == "Quebec" ~ "Canada - Quebec",
      `Country/Region` == "Taiwan*" ~ "Taiwan",
      TRUE ~ `Country/Region`)) %>%
  select(-c(`Province/State`,`Country/Region`)) %>% filter(Geo != "US") %>%
  bind_rows(dfE3) %>% bind_rows(dfE4) %>%
  group_by(Geo, Date, Var) %>% summarize(Cum = sum(Cum)) %>% ungroup() %>%
  arrange(Var, Geo, Date) %>%
  group_by(Geo, Var) %>% mutate(
    New = Cum - lag(Cum, default=0),
    change = ifelse(New == 0 & lead(New, default=0) >= 10, floor((New + lead(New, default=0)) / 2), 0),
    Adj = ifelse(New >= 0, New, 0) + change - lag(change, default = 0)) %>% ungroup()

dfCo <- dfE %>% group_by(Geo, Var) %>% summarize(Total = sum(New)) %>% pivot_wider(names_from = Var, values_from=Total)

dfDates <- dfDates %>% mutate(Start = as.Date(Start, format = "%d/%m/%y"), End = as.Date(End, format = "%d/%m/%y"))
vGeo <- dfDates %>% pull(Geo)
dfESub <- inner_join(
  dfDates,
  dfE
) %>%
  filter(Date >= Start, ifelse(!is.na(End), Date <= End, TRUE) ) %>%
  mutate(
    Day = as.integer(Date - Start, units = "days") + 1 + lagDeath[2],
    Geo = factor(Geo, levels = vGeo),
    Var = factor(Var, levels = c("Death", "Case")))
mRep <- xtabs(Adj ~ Var + Geo + Day, dfESub)
lmaxDeath <- log(max(mRep["Death",,]))
vTmax <- dfESub %>% filter(Var == "Case") %>% group_by(Geo) %>% summarize(Tmax = max(Day)) %>% pull(Tmax)
# Check that date range for cases & deaths are the same:
stopifnot(vTmax == dfESub %>% filter(Var == "Death") %>% group_by(Geo) %>% summarize(Tmax = max(Day)) %>% pull(Tmax))
nT <- max(vTmax)
nGeo <- length(vGeo)

# Make all policies unique, and lockdown policies cumulative
stopifnot(anyDuplicated(dfPRaw %>% select(Geo, Policy)) == 0)
dfP <- dfPRaw %>% select(-c(SourceStart, SourceEnd)) %>%
  mutate(PolStart = as.Date(PolStart, format = "%d-%b-%y"), PolEnd = as.Date(PolEnd, format = "%d-%b-%y"))
dfP <- dfP %>%
  left_join(dfDates) %>% mutate(
    DayStart = pmax(1, as.integer(PolStart - Start, units = "days") + 1 + lagDeath[2]),
    DayEnd = pmax(1, as.integer(PolEnd - Start, units = "days") + 1 + lagDeath[2])) %>%
  select(-c(Start, End)) %>%
  left_join(dfPNames) %>% mutate(PolName = paste(Policy, PolName))
dfPFull <- expand_grid(Geo = vGeo, Date = unique(dfE$Date), Policy = unique(dfP$Policy)) %>%
  left_join(bind_rows(
    dfP %>% select(Geo, Policy, Date = PolStart) %>% mutate(Val=1),
    dfP %>% select(Geo, Policy, Date = PolEnd) %>% mutate(Val=0))) %>%
  arrange(Geo, Policy, Date) %>%
  group_by(Geo, Policy) %>% fill(Val) %>% ungroup() %>%
  mutate(Val = ifelse(is.na(Val), 0, Val))  %>%
  left_join(dfPNames) %>% mutate(PolName = paste(Policy, PolName))
stopifnot(with(dfPFull %>% select(-PolName) %>% pivot_wider(names_from = Policy, values_from = Val),
              all(`1` >= `7`, `7` >= `8`, `8` >= `9`))) # Check that policies 1, 7, 8 and 9 are cumulative
dfTmp <- dfPFull %>% select(-PolName) %>% pivot_wider(names_from = Policy, values_from = Val); with(dfTmp, which(!(`1` >= `7` & `7` >= `8` &`8` >= `9`)))

# Select most relevant policies:
nPol <- length(PolSel)
dfP <- dfP %>%
  filter(Policy %in% PolSel) %>%
  left_join(tibble(Policy = PolSel, PolNew = 1:nPol)) %>%
  select(-Policy) %>% rename(Policy = PolNew) %>%
  mutate(Geo = factor(Geo, levels = vGeo))

mPolStart <- xtabs(DayStart ~ Geo + Policy, dfP)
mPolStart[mPolStart == 0] <- nT + 1
mPolEnd <- xtabs(DayEnd ~ Geo + Policy, dfP %>% mutate(Policy = factor(Policy)))
mPolEnd[mPolEnd == 0 | mPolEnd > nT] <- nT + 1

stopifnot(rownames(mPolStart) == vGeo, colnames(mPolStart) == 1:nPol)
stopifnot(rownames(mPolEnd) == vGeo, colnames(mPolEnd) == 1:nPol)
stopifnot(dimnames(mRep)$Day == (lagDeath[2]+1):nT, dimnames(mRep)$Geo == vGeo)

lData <- list(nGeo = nGeo, nT = nT, nTPred = nTPred, nPol = nPol, nK = 2, lagMin = c(lagDeath[1], lagCase[1]), 
              lagMax = c(lagDeath[2], lagCase[2]), vTmax = vTmax, mRep = mRep, mPolStart = mPolStart, mPolEnd = mPolEnd,
              lmaxDeath = lmaxDeath, mortality = mortality, mortSig = mortSig, lr2 = lr2, phi2 = phi2)

m <- stan_model("model.stan")
if(useInit) load("init.RData") else init <- "random"
if(test){
  fit <- sampling(m, data = lData, chains=2, iter = 200, thin = 1, init = init[c(1,2)], control = list(max_treedepth = 10, adapt_delta = 0.8))
} else {
  fit <- sampling(m, data = lData, chains=4, iter = 500, thin = 1, init = init, control = list(max_treedepth = 12, adapt_delta = 0.9))
}
save(list = ls(), file = paste0("output/fit-model-", time.now, ".RData"))

print(summary(fit, pars=c("g0Mu", "g1Mu", "g2Mu", "dgMu", "g0Sig", "gSig", "phi", "p", "muLag", "sigLag", "gDraw"))$summary)
#print(summary(fit, pars=c("dgMu", "dgSig", "g0Mu", "g0Sig", "gMu", "gSig", "phi", "p", "muLag", "sigLag", "gDraw"))$summary)
#print(xtable(summary(fit, pars=c("dgMu", "phi", "p", "muLag", "sigLag"))$summary[,c('50%','2.5%','97.5%','n_eff','Rhat')],
#             digits=c(0,3,3,3,0,2)), type="html", file=paste0("output/table-", time.now, ".html"))

samples <- extract(fit)
nIter <- length(samples$g0Mu)
if(!test) {
  iters <- sample(1:nIter, 4)
  init <- lapply(iters, function(i)
    lapply(samples[c("lir1", "g0", "g0Mu", "g0Sig", "gMu", "gSig", "dg", "dgMu", "dgSig", "phi", "p", "eventFrac", "muLag", "sigLag", "g0Draw", "dgDraw")],
           function(x) if(length(dim(x)) == 3) x[i,,] else if(length(dim(x)) == 2) x[i,] else x[i]))
  save(init, file = "init.RData")
}

apply(samples$dg[,vGeo == "France",],2,function(x) quantile(x, probs=c(0.025,0.5,0.975)))

dfDates <- dfDates %>% mutate(Tmax = vTmax)
dfOutRaw <- bind_rows(
  expand.grid(iter = 1:nIter, Geo = vGeo, Day=1:nT) %>% as_tibble() %>% mutate(Var = "Infection", Log = as.vector(samples$lir)),
  expand.grid(iter = 1:nIter, Var = c("Death", "Case"), Geo = vGeo, Day=(lagDeath[2]+1):nT) %>% as_tibble() %>% mutate(Log = as.vector(samples$lrr))
) %>% full_join(dfDates %>% select(Geo, Tmax, End)) %>% filter(Day <= Tmax)
dfOutPred <- bind_rows(
  expand.grid(iter = 1:nIter, Geo = vGeo, Day2=1:nTPred) %>% as_tibble() %>% mutate(Var = "Infection", Log = as.vector(samples$lirPred)),
  expand.grid(iter = 1:nIter, Var = c("Death", "Case"), Geo = vGeo, Day2=1:nTPred) %>% as_tibble() %>% mutate(Log = as.vector(samples$lrrPred))
) %>% full_join(dfDates %>% select(Geo, Tmax, End)) %>% mutate(Day = Tmax + Day2) %>% select(-Day2)

dfOut <- bind_rows(dfOutRaw, dfOutPred) %>% mutate(New = exp(Log)) %>% arrange(Geo, Var, Day) %>% 
  group_by(Geo, Var, iter) %>% mutate(Cum = cumsum(New)) %>%
  group_by(Geo, Day, Var) %>%
  summarize(NewEst = median(New), NewLow = quantile(New, probs=0.025), NewHigh = quantile(New, probs=0.975),
            CumEst = median(Cum), CumLow = quantile(Cum, probs=0.025), CumHigh = quantile(Cum, probs=0.975)) %>%
  ungroup() %>%
  inner_join(dfDates) %>%
  mutate(Date = Start + Day - 1 - lagDeath[2])
dfOut2 <- left_join(
  dfOut %>% select(Geo, Date, Var, NewEst:CumHigh),
  dfE %>% select(Geo, Date, Var, New, Cum))
dfOutUs <- dfOut2 %>% filter(substr(Geo,1,2) == "US") %>%
  group_by(Date, Var) %>% summarize_at(vars(NewEst:Cum), sum, na.rm=T) %>%
  mutate(Geo = "US - 12 states total")
dfOut2 <- dfOut2 %>%
  bind_rows(dfOutUs) %>%
  bind_rows(dfE %>% filter(Geo %in% c("Japan", "New Zealand", "Singapore", "Taiwan", "Korea, South", "China - Other mainland"))) %>%
  mutate(Var = factor(Var, levels = c("Death", "Case", "Infection"), labels = c("Death", "Reported case", "Infection")),
         Cum = ifelse(Cum == 0, NA, Cum), New = ifelse(New == 0, NA, New)) %>%
  select(-c(change, Adj))
write_csv(dfOut2, "model-out.csv")

muLag <- apply(samples$muLag, 2, median)
sigLag <- apply(samples$sigLag, 2, median)
dfLag <- bind_rows(
  tibble(Lag = lagCase[1]:lagCase[2], w = dnorm(Lag, muLag[2], sigLag[2]), Weight = w/sum(w), Var = "Case"),
  tibble(Lag = lagDeath[1]:lagDeath[2], w = dnorm(Lag, muLag[1], sigLag[1]), Weight = w/sum(w), Var = "Death")
)
dfTest <- t(apply(samples$eventFrac[,2,], 2, quantile, probs=c(0.025,0.5,0.975))) %>% as_tibble() %>%
  mutate(Geo = vGeo)

dfGeoRaw <- bind_rows(
  expand.grid(iter = 1:nIter, Geo = vGeo) %>% as_tibble() %>% mutate(x = expm1(as.vector(samples$g)), Var = "g"),
  expand.grid(iter = 1:nIter, Geo = vGeo) %>% as_tibble() %>% mutate(x = expm1(as.vector(samples$g0)), Var = "g0"))
dfGeo <- dfGeoRaw %>% group_by(Geo, Var) %>% 
  summarize(Estimate = median(x), Low = quantile(x, probs=0.025), High = quantile(x, probs=0.975)) %>% ungroup() %>%
  left_join(dfP %>% left_join(dfDates %>% select(Geo, DataEnd = End) %>% mutate(DataEnd = ifelse(is.na(DataEnd), max(dfE$Date), DataEnd))) %>%
              filter(PolEnd > DataEnd | is.na(PolEnd)) %>% arrange(Geo, PolStart) %>% group_by(Geo) %>%
              summarize(Policy = last(Policy), PolName = substr(last(PolName),3,30))) %>%
  mutate(PolName = ifelse(Var == "g0", "No restrictions (g0)", PolName), Policy = ifelse(Var == "g0", 0, Policy))
dfGeo2 <- dfOutRaw %>% filter(Var == "Infection", is.na(End), Day == Tmax) %>% select(-Var) %>%
  left_join(dfGeoRaw %>% filter(Var == "g") %>% select(iter, Geo, g=x)) %>%
  mutate(
    DaysTo1 = ifelse(Log > 0, ifelse(g < 0,  -Log / log1p(g), Inf), 0),
    DaysTo100 =  ifelse(Log > log(100), ifelse(g < 0,  -(Log - log(100)) / log1p(g), Inf), 0)) %>%
  select(iter, Geo, DaysTo1, DaysTo100) %>% pivot_longer(c(DaysTo1, DaysTo100), names_to = "Var", values_to = "x") %>%
  group_by(Geo, Var) %>% summarize(Estimate = median(x), Low = quantile(x, probs=0.025), High = quantile(x, probs=0.975)) %>% ungroup() %>%
  arrange(Var, Geo)
  #pivot_wider(id_cols = c(Geo, Var), names_from = Var, values_from = Estimate:High)
write_csv(dfGeo2, paste0("output/table-days-", time.now, ".csv"))

dfGRaw <- expand.grid(iter = 1:nIter, Policy = 1:nPol)  %>% as_tibble() %>% mutate(x = expm1(as.vector(samples$gDraw)))
dfG <- dfGRaw %>% group_by(Policy) %>%
  summarize(Estimate = median(x), Low = quantile(x, probs=0.025), High = quantile(x, probs=0.975)) %>% ungroup() %>%
  left_join(dfP %>% group_by(Policy) %>% summarize(PolName = first(PolName))) %>%
  mutate(PolName = ifelse(is.na(PolName), "0 No restrictions", PolName))

dfG2Raw <- bind_rows(
  expand.grid(iter = 1:nIter, Geo = vGeo, Var = "g0") %>% as_tibble() %>% mutate(x = samples$g0),
  expand.grid(iter = 1:nIter, Geo = vGeo, Var = "g1") %>% as_tibble() %>% mutate(x = samples$g1),
  expand.grid(iter = 1:nIter, Geo = vGeo, Var = paste0("dg", 2:nPol))  %>% as_tibble() %>% mutate(x = expm1(as.vector(samples$dg))))
dfG2Raw <- bind_rows(dfG2Raw,
                     dfG2Raw %>% filter(Var != "g0") %>% group_by(iter, Geo) %>% summarize(x = sum(x)) %>% ungroup() %>%
                       mutate(Var = "g2"))
dfG2 <- dfG2Raw %>% group_by(Geo, Var) %>%
  summarize(Est = median(x), Low = quantile(x, probs=0.025), High = quantile(x, probs=0.975))

save(list = ls(), file = paste0("output/fit-model-", time.now, ".RData"))
save(dfOut2, dfGeo, dfGeo2, vGeo2, file="data-app/fit-model-data.Rdata")

#### Charts ####

Date1 <- as.Date("2020-01-15")
Date2 <- as.Date("2020-07-01")

figCum <- dfOut2 %>%
  filter(Date >= as.Date(Date1), Date < as.Date(Date2), Geo %in% vGeo2) %>% mutate(Geo = factor(Geo, levels = vGeo2)) %>%
  ggplot(aes(x=Date, y=CumEst, color=Var)) +
  geom_line(aes(y=CumLow), linetype = 2, size = 0.25) + geom_line(aes(y=CumHigh), linetype = 2, size = 0.25) +
  geom_point(aes(y=Cum), size=0.25) + geom_line(size = 0.25) + facet_wrap(~Geo, ncol = 5) +
  scale_x_date(breaks = c(as.Date("2020-02-01"), as.Date("2020-04-01"), as.Date("2020-06-01")), labels = function(x) substr(format(x, format="%d %b"),2,10),
               date_minor_breaks = "1 month", limits = c(as.Date("2020-01-15"), Date2)) +
  scale_y_continuous(labels = scales::comma, trans="log10", limits = c(1,1e7)) +
  scale_color_manual(values = brewer.pal(3,"Set1"), guide = guide_legend(reverse = TRUE)) +
  ggtitle("Cumulative events (logarithmic scale)", subtitle = "Dot = Reported data, Line = Model estimate, Dash = 95% interval") +
  xlab(element_blank()) + ylab(element_blank()) + theme(legend.title=element_blank(), legend.position="bottom")
ggsave(paste0("fig-cumulative-", time.now, filetype), plot = figCum, path = "output", width = 8, height = 10)

figNew <- dfOut2 %>% 
  filter(Date >= Date1, Date < Date2, Geo %in% vGeo2) %>% mutate(Geo = factor(Geo, levels = vGeo2)) %>%
  ggplot(aes(x=Date, y=NewEst, color=Var)) +
  geom_line(aes(y=NewLow), linetype = 2, size = 0.25) + geom_line(aes(y=NewHigh), linetype = 2, size = 0.25) +
  geom_point(aes(y=New), size=0.25) + geom_line(size = 0.25) + facet_wrap(~Geo, ncol = 5) +
  scale_x_date(breaks = c(as.Date("2020-02-01"), as.Date("2020-04-01"), as.Date("2020-06-01")), labels = function(x) substr(format(x, format="%d %b"),2,10),
               date_minor_breaks = "1 month", limits = c(as.Date("2020-01-15"), Date2)) +
  scale_y_continuous(labels = scales::comma, trans="log10", limits = c(1,1e6)) +
  scale_color_manual(values = brewer.pal(3,"Set1"), guide = guide_legend(reverse = TRUE)) +
  ggtitle("New events per day (logarithmic scale)", subtitle = "Dot = Reported data, Line = Model estimate, Dash = 95% interval") +
  xlab(element_blank()) + ylab(element_blank()) + theme(legend.title=element_blank(), legend.position="bottom")
ggsave(paste0("fig-new-", time.now, filetype), plot = figNew, path = "output", width = 8, height = 10)

# figG <- dfGeo %>% 
#   ggplot(aes(x = factor(Geo, levels = dfGeo %>% filter(Var == "g") %>% arrange(Estimate) %>% pull(Geo)),
#              y = Estimate, ymin = Low, ymax = High, color = fct_reorder(PolName, Policy))) +
#   geom_hline(yintercept = 0, linetype=2) + 
#   geom_pointrange(size=0.4) + 
#   scale_y_continuous(labels = function(x) scales::percent(x, accuracy=1),
#                      breaks = seq(-0.2, 0.6, 0.2), minor_breaks = seq(-0.2, 0.6, 0.1)) +
#   scale_color_manual(values = c("Grey40","blue","Orange","Red"), name = "Policy in place") +
#   ylab("Growth rate of new infections per day (estimate and 95% interval)") +
#   coord_flip() + xlab(element_blank())
# ggsave(paste0("fig-g-", time.now, filetype), plot = figG, path = "output", width = 6, height = 5)

figG2 <- dfGeo %>% filter(Var != "g0", Geo != "China - Hubei") %>%
  ggplot(aes(x = factor(Geo, levels = dfGeo %>% filter(Var == "g") %>% arrange(Estimate) %>% pull(Geo)),
#             color = fct_reorder(PolName, Policy),
             y = Estimate, ymin = Low, ymax = High)) +
  geom_hline(yintercept = 0, linetype=2) + 
  geom_pointrange(size=0.4) + 
  scale_y_continuous(labels = function(x) scales::percent(x, accuracy=1)) +
#  scale_color_manual(values = c("Blue","Orange","Red"), name = "Policy in place",
#                     guide = guide_legend(title.position = "top", ncol = 2)) +
  ggtitle("Growth rate of new infections per day", subtitle = "(Estimate and 95% interval)") +
  coord_flip() + xlab(element_blank()) + ylab(element_blank()) +
  theme(legend.position="bottom", plot.title = element_text(size = 12))
ggsave(paste0("fig-g2-", time.now, filetype), plot = figG2, path = "output", width = 4.2, height = 6)

figPol <- dfPFull %>% filter(Policy %in% PolSel) %>% ggplot(aes(x=Date, y=Val, fill = PolName)) + geom_col(width = 1) + facet_wrap(~Geo) +
  scale_x_date(date_breaks = "1 month", labels = function(x) format(x, format="%b")) +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  scale_fill_manual(name = element_blank(), values = brewer.pal(9,"Set1"), guide = guide_legend(nrow = 3)) +
  theme(legend.title=element_blank(), legend.position="top") +
  xlab(element_blank())
ggsave(paste0("fig-policies-", time.now, filetype), plot = figPol, path = "output", width = 8, height = 11)

figGPol <- dfG %>% ggplot(aes(x = reorder(substr(PolName, 3,30), -Policy), y = Estimate, ymin = Low, ymax = High)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, color="red", linetype=2) + ylab("Infection growth rate per day (Estimate and 95% interval)") +
  coord_flip() + xlab(element_blank())
ggsave(paste0("fig-g-policy-", time.now, filetype), plot = figGPol, path = "output", width = 7, height = 3)

figTest <- dfTest %>%
  ggplot(aes(x=reorder(Geo,-`50%`), y=`50%`,ymin=`2.5%`,ymax=`97.5%`)) + geom_pointrange() + 
  coord_flip() + xlab(element_blank()) + ylab("Fraction of infections tested (median estimate and 95% interval)")
ggsave(paste0("fig-test-", time.now, filetype), plot = figTest, path = "output", width = 7, height = 5)

figLag <- dfLag %>%
  ggplot(aes(x=Lag, y=Weight)) + geom_col() + facet_wrap(~Var, ncol = 1) +
  ylab("Fraction of cases/deaths reported after time lag of getting infected") +
  xlab("Time lag (Days)")
ggsave(paste0("fig-lag-", time.now, filetype), plot = figLag, path = "output", width = 7, height = 5)

figGDetail <- ggplot(dfG2, aes(x=Geo, y=Est, ymin=Low, ymax=High)) + 
  geom_pointrange() + coord_flip() + xlab(element_blank()) + facet_grid(~Var, scales = "free_x")
ggsave(paste0("fig-g3-", time.now, filetype), plot = figGDetail, path = "output", width = 10, height = 5)
