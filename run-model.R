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

# It takes several hours to estimate the model projections on a regular PC


library(RColorBrewer)
library(stats4)
library(tidyverse)
library(rstan)

#### Set options ####
minPop <- 5e6
minDate <- as.Date("2020-02-01")
downloadData <- F
polG1 <- c("C1 - 1", "C2 - 1", "C3 - 1", "C7 - 1")
polExcl <- c("C4 - 1", "C6 - 3", "C8 - 2", "H2 - 1", "H3 - 1")
holidays <-as.Date(c("2020-06-20", "2020-08-28")) # Assume C1 School closure = 3 during summer holidays
lagCaseMax <- 2 # Lag of case reporting (#weeks)
lagDeathMax <- 4 # Lag of deaths (#weeks)
mortMu <- 0.01; mortSig <- 0.5; # Parameters for log-normal distribution; mortSig = 0.5 means 95% interval of */exp(1.96*0.5)=2.6
pOutl <- 1e-3
idgSig <- 0.02 # s.d. of change in idiosyncratic growth rate (AR(2) process)
chains <- 4; iter <- 1000; thin <- 10

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
time.now <- format(Sys.time(), format='%y%m%d-%H%M%S')
print(time.now)
suppressWarnings(dir.create(file.path("output")))


#### Read data ####

dfJh <- read_csv("input/jh-database.csv")
dfEcon <- read_csv("input/econ-database.csv")
dfPop <- read_csv("input/econ-population.csv")
dfOx <- read_csv("input/oxford-policy.csv")

dfGeoAdd <- read_csv("input/geo-add.csv")

#### Clean data ####
dfPop <- bind_rows(dfPop, dfGeoAdd %>% select(-geoJh))
vGeo <- sort(dfPop %>% filter(population >= minPop, geo != "US - New York City") %>% pull(geo))

vDate <- sort(unique(dfJh$date))
vDate <- vDate[as.numeric(vDate) %% 7 == as.numeric(min(dfEcon$date)) %% 7]
vDate <- vDate[vDate >= minDate]
stopifnot(as.numeric(vDate - lag(vDate))[-1] == 7)

# Epidemiology data
dfJh <- dfJh %>% left_join(dfGeoAdd %>% select(geo = geoJh, geoNew = geo)) %>%
  mutate(geo = ifelse(is.na(geoNew), geo, geoNew)) %>% select(-geoNew)
dfE <- dfJh %>% filter(geo %in% vGeo, date %in% vDate) %>%
  arrange(geo, var, date) %>%
  group_by(geo, var) %>% mutate(new = cum - lag(cum, default = 0)) %>% ungroup() %>%
  select(-cum) %>% pivot_wider(names_from = var, values_from = new) %>%
  left_join(dfEcon) %>%
  mutate_at(vars(case:deathExp), ~ ifelse(is.na(.), -1, .)) # Set missing values to -1
mCase <- xtabs(case ~ geo + date, dfE)
mCase[mCase < 0] <- 0
mDeathRep <- xtabs(death ~ geo + date, dfE)
mDeathRep[mDeathRep < 0] <- 0
mDeathTot <- xtabs(deathTot ~ geo + date, dfE)
mDeathExp <- xtabs(deathExp ~ geo + date, dfE)
stopifnot(rownames(mCase) == vGeo, colnames(mCase) == vDate)

# Policy data; For policy definitions, see https://github.com/OxCGRT/covid-policy-tracker/blob/master/documentation/codebook.md
dfOxSub <- dfOx %>% filter(geo %in% vGeo)
dfPolSel <- dfOxSub %>% distinct(polCode, level) %>% arrange(polCode, level) %>% 
  filter(level > 0) %>% rename(polCodeSel = polCode, levelSel = level)
dfPAll <-
  bind_rows(
    dfPolSel %>% rowwise() %>% mutate(data = list(dfOxSub %>% filter(polCode == polCodeSel, level >= levelSel)), value = 1),
    dfPolSel %>% rowwise() %>% mutate(data = list(dfOxSub %>% filter(polCode == polCodeSel, level < levelSel)), value = 0)) %>%
  unnest(data) %>% mutate(pol = paste(polCode, polName, levelSel, sep = " - ")) %>% select(geo, date, pol, value) %>%
  distinct()
vPolAll <- sort(unique(dfPAll$pol))
dfPAll <- full_join(expand_grid(geo = vGeo, pol = vPolAll, date = vDate), dfPAll) %>% arrange(geo, pol, date) %>%
  group_by(geo, pol) %>% fill(value) %>% ungroup() %>% filter(date %in% vDate) %>%
  mutate(pol2 = pol) %>% separate(pol2, c("polCode", "polName", "level"), sep = " - ") %>%
  mutate(polS = paste(polCode, level, sep = " - "),
         value = ifelse(polCode == "C1" & date >= holidays[1] & date <= holidays[2], 1, value),
         value = ifelse(is.na(value), 0, value))
dfPCor <- expand_grid(pol1 = vPolAll, pol2 = vPolAll) %>%
  mutate(cor = as.vector(cor(as.matrix(dfPAll %>% select(geo, pol, date, value) %>% pivot_wider(names_from = pol) %>% select(-c(geo, date))))))
stopifnot(nrow(dfPAll) == length(vGeo) * length(vDate) * length(vPolAll), !is.na(dfPAll$value))

mPol <- dfPAll %>% filter(!polS %in% c(polExcl, polG1)) %>%
  mutate(date = factor(date, levels = as.character(vDate))) %>%
  xtabs(value ~ geo + date + pol, .)
vPol <- dimnames(mPol)$pol
stopifnot(dimnames(mPol)$date == vDate, dimnames(mPol)$geo == vGeo,
          length(vPolAll) == length(vPol) + length(polExcl) + length(polG1))
mPolChange <- matrix(as.numeric(cbind(rep(TRUE, length(vGeo)), apply(mPol[,-length(vDate),] != mPol[,-1,], c(1,2), any))), nrow=length(vGeo))

mPolG1 <- dfPAll %>% filter(polS %in% polG1) %>% group_by(geo, date) %>% summarize(value = sum(value)) %>%
  mutate(value = ifelse(value > 1, 1, 0)) %>%
  xtabs(value ~ geo + date, .)
stopifnot(rownames(mPolG1) == vGeo, colnames(mPolG1) == vDate)

dfTest <- dfOxSub %>% filter(polCode == "H2") %>%
  full_join(expand_grid(geo = vGeo, date = vDate)) %>% arrange(geo, date) %>%
  group_by(geo) %>% fill(level) %>% ungroup() %>% filter(date %in% vDate)
mTest <- xtabs(level ~ geo + date, dfTest)
mTest[mTest == 0] <- 1
stopifnot(rownames(mTest) == vGeo, colnames(mTest) == vDate)

#### Fit model ####

# Estimate outlier distributions:
nll <- function(lLam, lPhi) -sum(dnbinom(y, mu = exp(lLam), size = exp(lPhi), log=T))

y <- as.vector(mCase)
fit <- mle(nll, start = list(lLam = 0, lPhi = 0))
summary(fit)
outlCase <- exp(coef(fit))

y <- as.vector(mDeathRep)
fit <- mle(nll, start = list(lLam = 0, lPhi = 0))
summary(fit)
outlDeath <- exp(coef(fit))

lData = list(mortMu = mortMu, mortSig = mortSig, lagDeathMax = lagDeathMax, lagCaseMax = lagCaseMax,
             nGeo = length(vGeo), nT = length(vDate), nPol = length(vPol), nTest = max(mTest),
             mPol = mPol, mPolChange = mPolChange, mPolG1 = mPolG1, mTest = mTest,
             mCase = mCase, mDeathRep = mDeathRep, mDeathTot = mDeathTot, mDeathExp = mDeathExp,
             outlCase = outlCase, outlDeath = outlDeath, pOutl = pOutl, idgSig = idgSig)

m <- stan_model("model.stan")
fit <- sampling(m, data = lData, chains = chains, iter = iter, thin = thin, control = list(max_treedepth = 12, adapt_delta = 0.9))
#fit <- sampling(m, data = lData, chains = 2, iter = 300)
save(list = ls(), file = paste0("output/image-", time.now, ".RData"))
print(fit, pars = c("deathAdj", "pLagCase", "pLagDeath", "phiCase", "phiDeathRep","phiDeathTot", "idgLam1", "idgLam2", 
                    "lmortality", "fracCaseMu", "fracCaseSig", "fracDeathMu", "fracDeathSig"))


#### Process model output ####
sim <- rstan::extract(fit)
nIter <- length(sim$phiCase)

dfOutRaw <- bind_rows(
  expand_grid(date = vDate, geo = vGeo, iter = 1:nIter) %>%  mutate(name = "infection", value = exp(as.vector(sim$logy))),
  expand_grid(date = vDate[-(1:lagCaseMax)], geo = vGeo, iter = 1:nIter) %>%  mutate(name = "case", value = exp(as.vector(sim$lCaseEst))),
  expand_grid(date = vDate[-(1:lagDeathMax)], geo = vGeo, iter = 1:nIter) %>%  mutate(name = "death", value = exp(as.vector(sim$lDeathEst))),
  expand_grid(date = vDate[-(1:lagDeathMax)], geo = vGeo, iter = 1:nIter) %>%  mutate(name = "deathTot", value = exp(as.vector(sim$lDeathTotEst))),
  expand_grid(pol = vPol, iter = 1:nIter) %>% mutate(name = "dg", value = as.vector(sim$dg)) %>%
    separate(pol, c("polCode", "polName", "level"), sep = " - ") %>%
    mutate(pol = paste(polCode, polName), level = as.integer(level)) %>% select(-c(polCode, polName)),
  expand_grid(level = 1:max(mTest), geo = vGeo, iter = 1:nIter) %>% mutate(name = "fracCase", value = as.vector(sim$fracCase)),
  expand_grid(level = 1:max(mTest), geo = vGeo, iter = 1:nIter) %>% mutate(name = "fracDeath", value = as.vector(sim$fracDeath)),
  expand_grid(geo = vGeo, iter = 1:nIter) %>% mutate(name = "g0", value = as.vector(sim$g0)),
  expand_grid(geo = vGeo, iter = 1:nIter) %>% mutate(name = "g1", value = as.vector(sim$g1)),
  expand_grid(date = vDate[-1], geo = vGeo, iter = 1:nIter) %>% mutate(name = "idg", value = as.vector(sim$idg)),
  expand_grid(date = vDate[-1], geo = vGeo, iter = 1:nIter) %>% mutate(name = "g", value = as.vector(sim$g)))
dfOutRaw <- dfOutRaw %>% bind_rows(
  dfOutRaw %>% filter(name == "dg") %>% mutate(name = "dgCum") %>%
    group_by(iter, pol) %>% mutate(value = cumsum(value)) %>% ungroup(),
  dfOutRaw %>% filter(name == "g") %>% mutate(name = "expg", value = expm1(value) * 100),
  dfOutRaw %>% filter(name == "idg") %>% left_join(dfOutRaw %>% filter(name == "g1") %>% select(geo, iter, g1 = value)) %>%
    mutate(name = "g1+idg", value = value + g1) %>% select(-g1) %>% filter(date >= as.Date("2020-04-01")))
dfOut <- dfOutRaw %>% group_by(name, geo, date, pol, level) %>%
  summarize(estimate = median(value), low = quantile(value, probs=0.025), high = quantile(value, probs=0.975)) %>% ungroup()

dfOutE <- dfOut %>% filter(name %in% c("infection", "case", "death", "deathTot"), date %in% vDate) %>%
  full_join(dfE %>% select(c(geo:deathTot)) %>% pivot_longer(case:deathTot, values_to = "reported")) %>%
  filter(name != "deathTot" | reported > 0) %>% select(-c(pol, level))
dfOutEPop <- left_join(dfOutE, dfPop) %>% mutate_at(vars(estimate:reported), ~ . / population * 1e4)

print(dfOutRaw %>% filter(name %in% c("g0", "g1", "idg")) %>% group_by(iter, name) %>% summarize(value = mean(value)) %>%
        group_by(name) %>% summarize(estimate = median(value), low = quantile(value, probs=0.025), high = quantile(value, probs=0.975))) %>%
  mutate(daily = expm1(estimate / 7))

#### Make charts ####

fGPol <- dfOut %>% filter(name == "dgCum") %>%
  group_by(pol) %>% mutate(estimate = estimate - c(0, estimate[-length(estimate)])) %>%
  left_join(dfOut %>% filter(name == "dgCum") %>% group_by(pol) %>% summarize(levelMax = max(level))) %>%
  mutate(low = ifelse(level == levelMax, low, NA), high = ifelse(level == levelMax, high, NA)) %>%
  ggplot(aes(x = fct_rev(pol), y = estimate, ymin = low, ymax = high, fill = fct_rev(factor(level)))) + geom_col() + geom_errorbar() + coord_flip()  +
  scale_fill_brewer(palette = "RdBu", name="Policy\nlevel") +
  xlab(element_blank()) + ylab(element_blank()) + theme(axis.text.y = element_text(hjust=0)) +
  labs(title = "Reduction of weekly growth rate",
       subtitle = "Bar = estimate; Line = 95% interval")

fPolSum <- dfPAll %>% group_by(polCode, polName, level, date) %>% summarize(frac = sum(value) / n()) %>%
  mutate(pol = paste(polCode, polName, sep = " - ")) %>%
  group_by(pol, date) %>% mutate(frac = frac - c(frac[-1], 0)) %>%
  ggplot(aes(x = date, y = frac, fill = level)) + geom_area() + facet_wrap(~ pol, ncol = 2)  + xlab(element_blank()) + ylab(element_blank()) +
  scale_fill_brewer(palette = "RdBu", name="Policy\nlevel", direction = -1) +
  scale_y_continuous(labels = function(x) scales::percent(x, accuracy=1)) +
  labs(title = "Percentage of jurisdictions with a given policy level in place",
       subtitle = "Source: Oxford coronavirus government response tracker") +
  theme(strip.text = element_text(size = 7.5))

fDash <- bind_rows(dfOut %>% filter(name == "g", date == max(date, na.rm = T)),
                   dfOutEPop %>% filter(name == "infection", date == max(date, na.rm = T))) %>%
  mutate(name = factor(name, levels = c("infection", "g"),
                       labels = c("New infections per 10,000 people", "Growth rate of new infections"))) %>%
  ggplot(aes(x = fct_rev(geo), y = estimate, ymin = low, ymax = high)) + 
  geom_hline(yintercept = 0, color = "grey60", linetype=2) + geom_pointrange() +
  facet_grid(~ name, scales = "free") + coord_flip()  + xlab(element_blank()) + ylab(element_blank()) +
  labs(title = paste0("Current status (weekly rates as of ", format(max(dfOut$date, na.rm = T), "%d %b"),")"),
       subtitle = "Dot = median estimate; Line = 95% interval") +
  theme(axis.text.y = element_text(hjust=0))

fGBase <- dfOut %>% filter(name == "g1+idg") %>%
  ggplot(aes(x = date, y = estimate, ymin = low, ymax = high)) + geom_ribbon(fill="grey70") + geom_line() +
  facet_wrap( ~ geo, ncol = 5) + xlab(element_blank()) + ylab(element_blank()) +
  labs(title = "Weekly base growth rate: variation in policy effectivenss",
       subtitle = "Effect of just recommended (level 1) policies for schools (C1),\nworkplaces (C2), events (C3), and internal movement (C7)") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), strip.text = element_text(size = 7))

pdf(paste0("output/charts-main-", time.now, ".pdf"), width=5.5, height=7, onefile=T)
  print(fGPol)
  print(fPolSum)
  print(fDash)
  print(fGBase)
dev.off()

fNew <- dfOutE %>%
  ggplot(aes(x = date, color = name)) +
  geom_point(aes(y = reported), size = 0.4) + geom_line(aes(y=estimate)) +
  geom_line(aes(y=low), lty=2) + geom_line(aes(y=high), lty=2) +
  facet_wrap(~ geo) + ggtitle("New events per week (log scale)") + 
  labs(subtitle = "Dot = reported; Line = model estimate; Dashed = 95% interval") +
  xlab(element_blank()) + ylab(element_blank()) + 
  scale_y_continuous(labels = scales::comma, trans="log10", limits = c(10, 1e6))

fNewCase <- dfOutEPop %>% filter(name %in% c("infection", "case")) %>%
  ggplot(aes(x = date, color = name)) +
  geom_point(aes(y = reported), size = 0.4) + geom_line(aes(y=estimate)) +
  geom_line(aes(y=low), lty=2) + geom_line(aes(y=high), lty=2)  +
  facet_wrap(~ geo) + xlab(element_blank()) + ylab(element_blank()) + 
  labs(subtitle = "Dot = reported; Line = model estimate; Dashed = 95% interval") +
  ggtitle("New infections & identified cases per week per 10,000 people")

fNewDeath <- dfOutEPop %>% filter(name %in% c("death", "deathTot")) %>%
  ggplot(aes(x = date, color = name)) +
  geom_point(aes(y = reported), size = 0.4) + geom_line(aes(y=estimate)) +
  geom_line(aes(y=low), lty=2) + geom_line(aes(y=high), lty=2) +
  facet_wrap(~ geo) + xlab(element_blank()) + ylab(element_blank()) + 
  labs(subtitle = "Dot = reported; Line = model estimate; Dashed = 95% interval") +
  ggtitle("New Covid and total deaths per week per 10,000 people")

fPol <- dfPAll %>% separate(pol, c("polCode", "polName", "level"), sep = " - ") %>% mutate(pol = paste(polCode, polName)) %>%
  ggplot(aes(x = date, y = value, fill = pol)) + geom_col(width = 7) + facet_wrap(~ geo) +
  scale_x_date(date_breaks = "1 month", labels = function(x) format(x, format="%b")) +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  scale_fill_manual(name = element_blank(), values = brewer.pal(12,"Set3"), guide = guide_legend(nrow = 3)) +
  theme(legend.title=element_blank(), legend.position="top") +
  xlab(element_blank()) + ggtitle("Policy levels recorded in Oxford database")

dfPCorSum <- dfPCor %>% group_by(pol1) %>% summarise(value = sum(cor)) %>% arrange(-value) %>% mutate(pol = pol1) %>%
  separate(pol, c("polCode", "polName", "level"), sep = " - ") %>% mutate(pol2 = paste(polCode, level, sep = " - "))
fPolCor <- dfPCor %>% separate(pol2, c("polCode", "polName", "level"), sep = " - ") %>%
  mutate(pol2 = paste(polCode, level, sep = " - ")) %>%
  ggplot(aes(x = factor(pol2, levels = dfPCorSum$pol2), y = fct_rev(factor(pol1, levels = dfPCorSum$pol1)), fill = cor)) + geom_raster() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab(element_blank()) + ylab(element_blank()) +
  ggtitle("Correlation coefficient between policies")

fGAll <- dfOut %>% filter(name == "g") %>%
  ggplot(aes(x = date, y = estimate, ymin = low, ymax = high)) + geom_hline(yintercept = 0) +
  geom_ribbon(fill="grey70") + geom_line() + facet_wrap( ~ geo)  + xlab(element_blank()) + ylab(element_blank()) +
  ggtitle("Weekly growth rate estimates")

fFrac <- dfOut %>% filter(substr(name, 1, 4) == "frac") %>% select(-c(date, pol)) %>% 
  left_join(dfTest %>% filter(date == max(date)) %>% select(geo, level) %>% mutate(current = T)) %>%
  mutate(name = paste(name, level), current = ifelse(is.na(current), F, current)) %>%
  ggplot(aes(x = fct_rev(geo), y = estimate, ymin = low, ymax = high, color = fct_rev(factor(current)))) + geom_pointrange(size = 0.3) +
  facet_grid(~name) + coord_flip()  + xlab(element_blank()) + ylab(element_blank()) +
  scale_color_brewer(type = "qual", name = "Current testing\npolicy") +
  ggtitle("Fraction of identified cases by testing regime (H2 policy) and identified deaths")

text <- paste("minPop =",minPop,"\npolG1 =", paste(polG1, collapse = ", "),"\npolExcl =", paste(polExcl, collapse = ", "),
              "\nholidays =", paste(holidays, collapse = ", "),
              "\nmortMu =", mortMu, "\nmortSig =", mortSig, "\npOutl =", pOutl, "\nidgSig =", idgSig,"\n\n",
              paste(capture.output(print(fit, pars = c("pLagCase", "pLagDeath", "phiCase", "phiDeathRep","phiDeathTot", "idgLam1",
                                                       "idgLam2", "lmortality","fracCaseMu", "fracCaseSig", "fracDeathMu",
                                                       "fracDeathSig"))), collapse = "\n"))
fPars <- ggplot() + annotate("text", x = 0, y = 0, size=4, label = text, family = "mono") + theme_void()

pdf(paste0("output/charts-sup-", time.now, ".pdf"), width=10, height=10, onefile=T)
  print(fNew)
  print(fNewCase)
  print(fNewDeath)
  print(fPol)
  print(fPolCor)
  print(fGAll)
  print(fFrac)
  print(fPars)
dev.off()
