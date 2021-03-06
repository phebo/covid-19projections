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

nTPred <- 6

library(stats4)
library(xtable)
library(tidyverse)
library(rstan)
source("functions.R")

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
time.now <- format(Sys.time(), format='%y%m%d-%H%M%S')
print(time.now)
suppressWarnings(dir.create(file.path("output")))

dfJh <- read_csv("input/jh-database.csv")
dfEcon <- read_csv("input/econ-database.csv")
dfPop <- read_csv("input/econ-population.csv")
dfOx <- read_csv("input/oxford-policy.csv")
dfHol <- read_csv("input/holidays.csv")

lFull <- clean.data(dfJh, dfEcon, dfPop, dfOx, dfHol)
p <- lFull$p
dfP <- lFull$dfP
dfE <- lFull$dfE
vDate1 <- p$vDate[p$vDate <= max(dfE$date) - nTPred*7]
vDate2 <- p$vDate[p$vDate > max(dfE$date) - nTPred*7]
lData <- clean.data(dfJh, dfEcon, dfPop, dfOx, dfHol,
                    dates = c(as.Date("2020-02-01"), max(dfE$date) - nTPred*7), nTPred = nTPred)$lData
lData$mPol <- lFull$lData$mPol
lData$mPolChange <- lFull$lData$mPolChange
lData$mPolG1 <- lFull$lData$mPolG1
lData$mTest <- lFull$lData$mTest

m <- stan_model("model.stan")

fit <- sampling(m, data = lData, chains = 4, iter = 700, warmup = 500, thin = 2, control = list(adapt_delta = 0.9), seed = 99743)
#fit <- sampling(m, data = lData, chains = 2, iter = 300)
save(list = ls(), file = paste0("output/image-oos-", time.now, ".RData"))
print(fit, pars = c("deathAdj", "pLagCase", "pLagDeath", "phiCase", "phiDeathRep","phiDeathTot", "lmortality"))

sim <- rstan::extract(fit)
nIter <- length(sim$phiCase)

dfOutRaw <- bind_rows(
  expand_grid(date = p$vDate, geo = p$vGeo, iter = 1:nIter) %>%  mutate(name = "case", value = exp(as.vector(sim$logy))),
  expand_grid(date = p$vDate[-(1:lData$lagCaseMax)], geo = p$vGeo, iter = 1:nIter) %>%  mutate(name = "case", value = exp(as.vector(sim$lCaseEst))),
  expand_grid(date = p$vDate[-(1:lData$lagDeathMax)], geo = p$vGeo, iter = 1:nIter) %>%  mutate(name = "death", value = exp(as.vector(sim$lDeathEst))),
  expand_grid(date = vDate2, geo = p$vGeo, iter = 1:nIter) %>% mutate(name = "predCase", value = as.vector(sim$predCase)),
  expand_grid(date = vDate2, geo = p$vGeo, iter = 1:nIter) %>% mutate(name = "predDeath", value = as.vector(sim$predDeath)))
dfOut <- dfOutRaw %>% group_by(name, geo, date) %>%
  summarize(estimate = median(value), low = quantile(value, probs=0.025), high = quantile(value, probs=0.975)) %>% ungroup()

dfOutE <- dfOut %>%
  filter(name %in% c("predCase", "predDeath")) %>% mutate(name = ifelse(name == "predCase", "case", "death")) %>%
  left_join(dfE %>% select(geo:death) %>% pivot_longer(case:death) %>% rename(reported = value))

fPred <- dfOutE %>%
  ggplot(aes(x = date, color = name)) +
  geom_point(aes(y = reported)) + geom_line(aes(y=estimate)) +
  geom_line(aes(y=low), lty=2) + geom_line(aes(y=high), lty=2) +
  facet_wrap(~ geo, ncol = 5) + # ggtitle("New events per week (out of sample predictions; log scale)") + 
  # labs(subtitle = "Dot = reported; Line = model prediction; Dashed = 95% interval") +
  scale_color_manual(values = brewer.pal(4, "Set1")[c(2,1,4,3)]) +
  scale_x_date(breaks = as.Date(c(vDate2[1], vDate2[4])), date_labels = "%b %e", minor_breaks = vDate2) +
  scale_y_continuous(labels = scales::comma, trans="log10", limits = c(1, 1e7)) + 
  xlab(element_blank()) + ylab(element_blank()) +
  theme(legend.position="top", legend.title = element_blank())
ggsave(paste0("figures/fig-oos-full.png"), height = 9, width = 6.5)

fPerc <- dfOutE %>% group_by(name, date) %>% summarize(nRange = sum(reported >= low & reported <= high), n = n()) %>% ungroup() %>%
  mutate(fRange = nRange/n) %>%
  ggplot(aes(x = date, y = fRange)) + geom_col() + facet_wrap(~ name)  + 
  scale_y_continuous(labels = function(x) scales::percent(x, accuracy=1)) +
  ggtitle("Percentage of predictions within range") +
  xlab(element_blank()) + ylab(element_blank())

fPerc2 <- dfOutE %>% group_by(date) %>% summarize(nRange = sum(reported >= low & reported <= high), n = n()) %>% ungroup() %>%
  mutate(fRange = nRange/n) %>% # filter(name == "case") %>%
  ggplot(aes(x = date, y = fRange)) + geom_col() +
  scale_y_continuous(labels = function(x) scales::percent(x, accuracy=1)) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  geom_hline(yintercept = c(0,1)) +
  xlab(element_blank()) + ylab(element_blank())
ggsave(paste0("figures/fig-oos.png"), height = 3, width = 3)

pdf(paste0("output/charts-oos-", time.now, ".pdf"), width=10, height=10, onefile=T)
  print(fPred)
  print(fPerc)
dev.off()

tabOut <- dfOutE %>% group_by(name, date) %>% summarize(nRange = sum(reported >= low & reported <= high), n = n()) %>% ungroup() %>%
  mutate(week = format(date, "%b-%d"), fRange = nRange/n) %>% select(name, week, fRange) %>% xtable()
print(tabOut, type = "html", file = paste0("output/tab-oos-", time.now, ".html"))
