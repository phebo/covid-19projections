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


library(stats4)
library(xtable)
library(tidyverse)
library(furrr)
library(rstan)
source("functions.R")

plan(multiprocess) # Allow parallel computing
time.now <- format(Sys.time(), format='%y%m%d-%H%M%S')
print(time.now)
suppressWarnings(dir.create(file.path("output")))

dfJh <- read_csv("input/jh-database.csv")
dfEcon <- read_csv("input/econ-database.csv")
dfPop <- read_csv("input/econ-population.csv")
dfOx <- read_csv("input/oxford-policy.csv")
dfHol <- read_csv("input/holidays.csv") %>% select(-source)

l <- list(
  clean.data(dfJh, dfEcon, dfPop, dfOx, dfHol),
  clean.data(dfJh, dfEcon, dfPop, dfOx, dfHol, minPop = 1e7),
  clean.data(dfJh, dfEcon, dfPop, dfOx, dfHol, minPop = 3e6),
  clean.data(dfJh, dfEcon, dfPop, dfOx, dfHol, pOutl = 1e-2),
  clean.data(dfJh, dfEcon, dfPop, dfOx, dfHol, pOutl = 1e-4),
  clean.data(dfJh, dfEcon, dfPop, dfOx, dfHol, idgSig = 0.01),
  clean.data(dfJh, dfEcon, dfPop, dfOx, dfHol, idgSig = 0.05),
  clean.data(dfJh, dfEcon, dfPop, dfOx, dfHol, geoExclude = c("Chile", "South Africa")),
  clean.data(dfJh, dfEcon, dfPop, dfOx, dfHol, dgSig = c(1,1e-4)),
  clean.data(dfJh, dfEcon, dfPop, dfOx, dfHol, dgMin = -2)
)
specNames <- c("Base", "Population > 10M", "Population > 3M", "P(outlier) = 0.01", "P(outlier) = 0.0001",
               "Stdev = 0.01", "Stdev = 0.05", "Excl Chile/S-Africa", "Spike & slab", "Negative Delta-g")

m <- stan_model("model.stan")
specs <- map(l, ~ list(m = m, data = .$lData, pars = "dg", iter = 700, warmup = 500, chains = 2, thin = 2))
chains <- make.chains(specs)
fits.chain <- future_map(chains, ~ do.chain(.))
fits <- cons.fits(fits.chain, chains)
save(list = ls(), file = paste0("output/image-sens-", time.now, ".RData"))
dgs <- map(fits, ~ extract(.)$dg)
names(dgs) <- specNames
dfRaw <- map2_dfr(dgs, l, ~ expand_grid(pol = .y$p$vPol, iter = 1:nrow(.x)) %>% mutate(value = as.vector(.x)), .id = "spec")
df <- dfRaw %>% group_by(spec, pol) %>%
  summarize(estimate = median(value), low = quantile(value, probs=0.025), high = quantile(value, probs=0.975)) %>% ungroup() %>%
  mutate(spec = factor(spec, levels = specNames))
fSens <- df %>% ggplot(aes(x = fct_rev(spec), y = estimate, ymin = low, ymax = high)) + geom_pointrange() +
  facet_wrap(~ pol, ncol = 4, labeller = label_wrap_gen(25)) + coord_flip() +
  xlab(element_blank()) + ylab(element_blank()) +
  theme(axis.text.y = element_text(hjust=0), strip.text = element_text(size = 8))
ggsave(paste0("figures/fig-sens.png"), height = 9, width = 6.5)
ggsave(paste0("output/chart-sens-", time.now, ".pdf"), width=12, height=8)
