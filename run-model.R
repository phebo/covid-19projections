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

library(xtable)
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
suppressWarnings(dir.create(file.path("figures")))

saveAppData <- F
writeFigures <- T
selSc <- c("C1 - 1", "C2 - 1", "C3 - 2", "C4 - 3", "C5 - 0", "C6 - 1", "C7 - 1", "C8 - 3", "H1 - 2", "H2 - 2", "H3 - 2")

#### Read and clean data ####

dfJh <- read_csv("input/jh-database.csv")
dfEcon <- read_csv("input/econ-database.csv")
dfPop <- read_csv("input/econ-population.csv")
dfOx <- read_csv("input/oxford-policy.csv")
dfHol <- read_csv("input/holidays.csv") %>% select(-source)

l <- clean.data(dfJh, dfEcon, dfPop, dfOx, dfHol)
dfP <- l$dfP; dfE <- l$dfE; lData <- l$lData; p <- l$p

print(l$dfPCor %>% filter(cor < 1 & cor > 0.85) %>% arrange(-cor))
print(dfP %>% group_by(polCode, polName, level) %>% summarize(frac = sum(value) / n()) %>% group_by(polCode) %>%
        mutate(frac = frac - c(frac[-1], 0)) %>% filter(frac < 0.05))

#### Fit model ####

m <- stan_model("model.stan")

fit <- sampling(m, data = lData, chains = 4, iter = 700, warmup = 500, thin = 2, control = list(adapt_delta = 0.9, max_treedepth = 12), seed = 99743)
#fit <- sampling(m, data = lData, chains = 2, iter = 300)
save(list = ls(), file = paste0("output/image-", time.now, ".RData"))
print(fit, pars = c("deathAdj", "pLagCase", "pLagDeath", "phiCase", "phiDeathRep","phiDeathTot", "idgLam1", "idgLam2", 
                    "lmortality"))
if(writeFigures) print(xtable(
  summary(fit, pars=c("pLagCase", "pLagDeath", "phiCase", "phiDeathRep","phiDeathTot", "idgLam1", "idgLam2"))$summary[,c('50%','2.5%','97.5%','n_eff','Rhat')],
  digits=c(0,3,3,3,0,2)), type="html", file=paste0("figures/table.html"))

#### Process model output ####
sim <- rstan::extract(fit)
nIter <- length(sim$phiCase)

dfOutRaw <- bind_rows(
  expand_grid(date = p$vDate, geo = p$vGeo, iter = 1:nIter) %>%  mutate(name = "infection", value = exp(as.vector(sim$logy))),
  expand_grid(date = p$vDate[-(1:lData$lagCaseMax)], geo = p$vGeo, iter = 1:nIter) %>%  mutate(name = "case", value = exp(as.vector(sim$lCaseEst))),
  expand_grid(date = p$vDate[-(1:lData$lagDeathMax)], geo = p$vGeo, iter = 1:nIter) %>%  mutate(name = "death", value = exp(as.vector(sim$lDeathEst))),
  expand_grid(date = p$vDate[-(1:lData$lagDeathMax)], geo = p$vGeo, iter = 1:nIter) %>%  mutate(name = "deathTot", value = exp(as.vector(sim$lDeathTotEst))),
  expand_grid(pol = p$vPol, iter = 1:nIter) %>% mutate(name = "dg", value = as.vector(sim$dg)) %>%
    separate(pol, c("polCode", "polName", "level"), sep = " - ") %>%
    mutate(pol = paste(polCode, polName), level = as.integer(level)) %>% select(-c(polCode, polName)),
  expand_grid(level = 1:max(lData$mTest), geo = p$vGeo, iter = 1:nIter) %>% mutate(name = "fracCase", value = as.vector(sim$fracCase)),
  expand_grid(level = 1:max(lData$mTest), geo = p$vGeo, iter = 1:nIter) %>% mutate(name = "fracDeath", value = as.vector(sim$fracDeath)),
  expand_grid(geo = p$vGeo, iter = 1:nIter) %>% mutate(name = "g0", value = as.vector(sim$g0)),
  expand_grid(geo = p$vGeo, iter = 1:nIter) %>% mutate(name = "g1", value = as.vector(sim$g1)),
  expand_grid(date = p$vDate[-1], geo = p$vGeo, iter = 1:nIter) %>% mutate(name = "idg", value = as.vector(sim$idg)),
  expand_grid(date = p$vDate[-1], geo = p$vGeo, iter = 1:nIter) %>% mutate(name = "g", value = as.vector(sim$g)))
dfOutRaw <- dfOutRaw %>% bind_rows(
  dfOutRaw %>% filter(name == "dg") %>% mutate(name = "dgCum") %>%
    group_by(iter, pol) %>% mutate(value = cumsum(value)) %>% ungroup(),
  dfOutRaw %>% filter(name == "g") %>% mutate(name = "expg", value = expm1(value) * 100),
  dfOutRaw %>% filter(name == "idg") %>% left_join(dfOutRaw %>% filter(name == "g1") %>% select(geo, iter, g1 = value)) %>%
    mutate(name = "g1+idg", value = value + g1) %>% select(-g1) %>% filter(date >= as.Date("2020-04-01")))
dfOut <- dfOutRaw %>% group_by(name, geo, date, pol, level) %>%
  summarize(estimate = median(value), low = quantile(value, probs=0.025), high = quantile(value, probs=0.975)) %>% ungroup()

dfOutE <- dfOut %>% filter(name %in% c("infection", "case", "death", "deathTot"), date %in% p$vDate) %>%
  full_join(dfE %>% select(c(geo:deathTot)) %>% pivot_longer(case:deathTot, values_to = "reported")) %>%
  filter(name != "deathTot" | reported > 0) %>% select(-c(pol, level))
dfOutEPop <- left_join(dfOutE, dfPop) %>% mutate_at(vars(estimate:reported), ~ . / population * 1e4)
print(dfOutEPop %>% filter(name == "infection") %>% group_by(geo) %>%
        summarise_at(vars(estimate:high), ~ sum(.) / 1e4) %>% arrange(-estimate), n=100)

print(dfOutRaw %>% filter(name %in% c("g0", "g1", "idg")) %>% group_by(iter, name) %>% summarize(value = mean(value)) %>%
        group_by(name) %>% summarize(estimate = median(value), low = quantile(value, probs=0.025), high = quantile(value, probs=0.975))) %>%
  mutate(daily = expm1(estimate / 7))

maxLag <- 10
eps <- t(cbind(matrix(sim$eps, ncol = length(p$vDate) - 3), matrix(nrow = nIter * length(p$vGeo), ncol = maxLag)))
stopifnot(eps[1:(length(p$vDate)-3)] == sim$eps[1,1,])

if(saveAppData) {
  write_csv(dfOutRaw %>% filter(name %in% c("dg")) %>% select(iter, pol, level, value), "pol-app/dg-data.csv")
  write_csv(dfOutRaw %>% filter(name %in% c("g1+idg"), date == max(p$vDate)) %>% group_by(iter) %>%
              summarize(low = quantile(value, probs = 0.1), medium = quantile(value, probs = 0.5), high = quantile(value, probs = 0.9)) %>%
              pivot_longer(low:high),
            "pol-app/gbase-data.csv")
}

dfOutSc <- left_join(
  dfOutRaw %>% filter(name %in% c("g1+idg"), date == max(p$vDate)) %>% group_by(iter) %>%
    summarize(`Best 10%` = quantile(value, probs = 0.1), `Median` = quantile(value, probs = 0.5), `Worst 10%` = quantile(value, probs = 0.9)) %>%
    pivot_longer(`Best 10%`:`Worst 10%`) %>% rename(gBase = value),
  dfOutRaw %>% filter(name == "dgCum") %>% mutate(pol = paste(substr(pol, 1, 2), "-", level)) %>%
    filter(pol %in% selSc) %>% group_by(iter) %>% summarize(dgCum = sum(value)) %>% ungroup()
) %>% mutate(g = gBase - dgCum) %>%
  select(iter, group = name, gBase, g) %>% pivot_longer(gBase:g) %>%
  group_by(group, name) %>%
  summarize(estimate = median(value), low = quantile(value, probs=0.025), high = quantile(value, probs=0.975)) %>%
  ungroup() #%>% mutate(name = factor(name, levels = c("g","gBase"), labels = c(expression(italic(g)), expression(italic(g)^(base)))))

#### Make charts ####

ggplot(dfOutSc, aes(x = group, y = estimate, color = name)) + geom_point() +
  geom_errorbar(aes(x = group, ymin = low, ymax = high), width = 0.6) +
  geom_hline(yintercept = 0) +
  xlab("Jurisdiction quantile") + ylab(element_blank()) +
  scale_color_brewer(labels = c("With core\npolices", "Base"), palette = "Dark2",
                     guide = guide_legend(reverse = TRUE), name = "Growth rate") +
  theme(axis.title.x = element_text(vjust = -0.5))
if(writeFigures) ggsave(paste0("figures/fig-scenarios.png"), height = 4, width = 4)

fGPol <- dfOut %>% filter(name == "dgCum") %>%
  group_by(pol) %>% mutate(estimate = estimate - c(0, estimate[-length(estimate)])) %>%
  left_join(dfOut %>% filter(name == "dgCum") %>% group_by(pol) %>% summarize(levelMax = max(level))) %>%
  mutate(low = ifelse(level == levelMax, low, NA), high = ifelse(level == levelMax, high, NA)) %>%
  ggplot(aes(x = fct_rev(pol), y = estimate, ymin = low, ymax = high, fill = fct_rev(factor(level)))) +
  geom_col(width = 0.8) + geom_errorbar(width = 0.8) + coord_flip()  +
  scale_fill_brewer(palette = "YlOrRd", name="Policy\nlevel", direction = -1, guide = guide_legend(reverse = TRUE)) +
  xlab(element_blank()) + ylab(element_blank()) + theme(axis.text.y = element_text(hjust=0))
if(writeFigures) ggsave(paste0("figures/fig-gpol.png"), height = 4, width = 5)

fPolSum <- dfP %>% group_by(polCode, polName, level, date) %>% summarize(frac = sum(value) / n()) %>%
  mutate(pol = paste(polCode, polName, sep = " - ")) %>%
  group_by(pol, date) %>% mutate(frac = frac - c(frac[-1], 0)) %>%
  ggplot(aes(x = date, y = frac, fill = level)) + geom_area() + facet_wrap(~ pol, ncol = 2)  + xlab(element_blank()) + ylab(element_blank()) +
  scale_fill_brewer(palette = "YlOrRd", name="Policy\nlevel") +
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
  #labs(title = paste0("Current status (weekly rates as of ", format(max(dfOut$date, na.rm = T), "%d %b"),")"),
  #     subtitle = "Dot = median estimate; Line = 95% interval") +
  theme(axis.text.y = element_text(hjust=0))
if(writeFigures) ggsave(paste0("figures/fig-status.png"), height = 9, width = 6.5)

fGBase <- dfOut %>% filter(name == "g1+idg") %>%
  ggplot(aes(x = date, y = estimate, ymin = low, ymax = high)) + geom_ribbon(fill="grey70") + geom_line() +
  facet_wrap( ~ geo, ncol = 5) + xlab(element_blank()) + ylab(element_blank()) +
  #labs(title = "Weekly base growth rate: variation in policy effectiveness",
  #     subtitle = "Effect of level 1 policies for C1-C4") +
  scale_x_date(breaks = as.Date(c("2020-05-01", "2020-08-01")), date_labels = "%b%e", date_minor_breaks = "1 month") +
  theme(strip.text = element_text(size = 8.5))
if(writeFigures) ggsave(paste0("figures/fig-gbase.png"), height = 9, width = 6.5)

fNew <- dfOutE %>% mutate(reported = ifelse(reported == 0, NA, reported)) %>%
  ggplot(aes(x = date, color = name)) +
  geom_point(aes(y = reported), size = 0.4) + geom_line(aes(y=estimate)) +
  geom_line(aes(y=low), lty=2) + geom_line(aes(y=high), lty=2) +
  facet_wrap(~ geo, ncol = 5) +
  #ggtitle("New events per week (log scale)") + 
  #labs(subtitle = "Dot = reported; Line = model estimate; Dashed = 95% interval") +
  xlab(element_blank()) + ylab(element_blank()) + 
  scale_y_continuous(labels = scales::comma, trans="log10", limits = c(1, 1e6)) +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position="top", legend.title = element_blank())
if(writeFigures) ggsave(paste0("figures/fig-new.png"), height = 9, width = 6.5)

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

fPol <- dfP %>% mutate(pol = paste(polCode, polName)) %>%
  ggplot(aes(x = date, y = value, fill = pol)) + geom_col(width = 7) + facet_wrap(~ geo) +
  scale_x_date(date_breaks = "1 month", labels = function(x) format(x, format="%b")) +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  scale_fill_manual(name = element_blank(), values = brewer.pal(12,"Set3"), guide = guide_legend(nrow = 3)) +
  theme(legend.title=element_blank(), legend.position="top") +
  xlab(element_blank()) + ggtitle("Policy levels recorded in Oxford database")

dfPCorSum <- l$dfPCor %>% group_by(pol1) %>% summarise(value = sum(cor)) %>% arrange(-value) %>% mutate(pol = pol1) %>%
  separate(pol, c("polCode", "polName", "level"), sep = " - ") %>% mutate(pol2 = paste(polCode, level, sep = " - "))
fPolCor <- l$dfPCor %>% separate(pol2, c("polCode", "polName", "level"), sep = " - ") %>%
  mutate(pol2 = paste(polCode, level, sep = " - ")) %>%
  ggplot(aes(x = factor(pol2, levels = dfPCorSum$pol2), y = fct_rev(factor(pol1, levels = dfPCorSum$pol1)), fill = cor)) + geom_raster() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab(element_blank()) + ylab(element_blank()) +
  ggtitle("Correlation coefficient between policies")

fGAll <- dfOut %>% filter(name == "g") %>%
  ggplot(aes(x = date, y = estimate, ymin = low, ymax = high)) + geom_hline(yintercept = 0) +
  geom_ribbon(fill="grey70") + geom_line() + facet_wrap( ~ geo)  + xlab(element_blank()) + ylab(element_blank()) +
  ggtitle("Weekly growth rate estimates")

fFrac <- dfOut %>% filter(substr(name, 1, 4) == "frac") %>% select(-c(date, pol)) %>% 
  left_join(l$dfTest %>% filter(date == max(date)) %>% select(geo, level) %>% mutate(current = T)) %>%
  mutate(name = paste(name, level), current = ifelse(is.na(current), F, current)) %>%
  ggplot(aes(x = fct_rev(geo), y = estimate, ymin = low, ymax = high, color = fct_rev(factor(current)))) + geom_pointrange(size = 0.3) +
  facet_grid(~name) + coord_flip()  + xlab(element_blank()) + ylab(element_blank()) +
  scale_color_brewer(type = "qual", name = "Current testing\npolicy") +
  ggtitle("Fraction of identified cases by testing regime (H2 policy) and identified deaths")


options(width = 100)
text <- paste("minPop =",p$minPop,"\npolG1 =", paste(p$polG1, collapse = ", "),"\npolExcl =", paste(p$polExcl, collapse = ", "),
              "\nmortMu =", p$mortMu, "\nmortSig =", p$mortSig, "\npOutl =", p$pOutl, "\nidgSig =", p$idgSig,"\n\n",
              paste(capture.output(print(fit, pars = c("pLagCase", "pLagDeath", "phiCase", "phiDeathRep","phiDeathTot", "idgLam1",
                                                       "idgLam2", "lmortality"))), collapse = "\n"))
fPars <- ggplot() + annotate("text", x = 0, y = 0, size=4, label = text, family = "mono") + theme_void()

fPol2 <- dfP %>% mutate(pol = paste(polCode, polName, sep = " - ")) %>% filter(value == 1) %>% group_by(geo, pol, date) %>%
  summarize(level = max(level)) %>%
  ggplot(aes(x = date, y = fct_rev(geo), fill = level)) + geom_tile(height=0.8) +
  facet_grid(~ pol, labeller = label_wrap_gen(8)) +
  scale_fill_brewer(palette = "YlOrRd", name="Policy\nlevel") + ggtitle("With holidays") +
  ylab(element_blank()) + xlab(element_blank())
fPol2nh <- dfP %>% mutate(pol = paste(polCode, polName, sep = " - ")) %>% filter(valueOx == 1) %>% group_by(geo, pol, date) %>%
  summarize(level = max(level)) %>%
  ggplot(aes(x = date, y = fct_rev(geo), fill = level)) + geom_tile(height=0.8) +
  facet_grid(~ pol, labeller = label_wrap_gen(8)) +
  scale_fill_brewer(palette = "YlOrRd", name="Policy\nlevel") + ggtitle("No holidays") +
  ylab(element_blank()) + xlab(element_blank())

fPol2a <- dfP %>% mutate(pol = paste(polCode, polName, sep = " - ")) %>% filter(value == 1, polCode %in% c("C1","C2","C3","C4","C5","C6")) %>%
  group_by(geo, pol, date) %>% summarize(level = max(level)) %>%
  ggplot(aes(x = date, y = fct_rev(geo), fill = level)) + geom_tile(height=0.8) +
  facet_grid(~ pol, labeller = label_wrap_gen(8)) +
  scale_x_date(breaks = as.Date(c("2020-03-01", "2020-08-01")), date_labels = "%b%e", date_minor_breaks = "1 month") +
  scale_fill_brewer(palette = "YlOrRd", name="Policy\nlevel") + ylab(element_blank()) + xlab(element_blank()) +
  theme(legend.position = "none")
if(writeFigures) ggsave(paste0("figures/fig-pol-a.png"), height = 9, width = 6.5)

fPol2b <- dfP %>% mutate(pol = paste(polCode, polName, sep = " - ")) %>% filter(value == 1, !polCode %in% c("C1","C2","C3","C4","C5","C6")) %>%
  group_by(geo, pol, date) %>% summarize(level = max(level)) %>%
  ggplot(aes(x = date, y = fct_rev(geo), fill = level)) + geom_tile(height=0.8) +
  scale_x_date(breaks = as.Date(c("2020-03-01", "2020-08-01")), date_labels = "%b%e", date_minor_breaks = "1 month") +
  facet_grid(~ pol, labeller = label_wrap_gen(8)) +
  scale_fill_brewer(palette = "YlOrRd", name="Policy\nlevel") + ylab(element_blank()) + xlab(element_blank())
if(writeFigures) ggsave(paste0("figures/fig-pol-b.png"), height = 9, width = 6.5)

pdf(paste0("output/charts-main-", time.now, ".pdf"), width=5.5, height=7, onefile=T)
  print(fGPol)
  print(fPolSum)
  print(fDash)
  print(fGBase)
dev.off()

pdf(paste0("output/charts-sup-", time.now, ".pdf"), width=10, height=10, onefile=T)
  print(fNew)
  print(fNewCase)
  print(fNewDeath)
  print(fPol)
  print(fPolCor)
  print(fGAll)
  print(fFrac)
  print(fPars)
  par(mfrow=c(2,1))
  acf(as.vector(eps), na.action = na.pass, lag.max = maxLag)
  pacf(as.vector(eps), na.action = na.pass, lag.max = maxLag)
dev.off()

pdf(paste0("output/chart-pol.pdf"), width=12, height=8, onefile=T)
  print(fPol2)
  print(fPol2nh)
dev.off()


