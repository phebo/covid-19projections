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
  clean.data(dfJh, dfEcon, dfPop, dfOx, dfHol, minPop = 1e7),
  clean.data(dfJh, dfEcon, dfPop, dfOx, dfHol, minPop = 3e6),
  clean.data(dfJh, dfEcon, dfPop, dfOx, dfHol, pOutl = 1e-2),
  clean.data(dfJh, dfEcon, dfPop, dfOx, dfHol, pOutl = 1e-4),
  clean.data(dfJh, dfEcon, dfPop, dfOx, dfHol, idgSig = 0.01),
  clean.data(dfJh, dfEcon, dfPop, dfOx, dfHol, idgSig = 0.05),
  clean.data(dfJh, dfEcon, dfPop, dfOx, dfHol, geoExclude = c("Chile", "South Africa"))
)

m <- stan_model("model.stan")
specs <- map(l, ~ list(m = m, data = .$lData, pars = "dg", iter = 700, warmup = 500, chains = 2, thin = 2))
chains <- make.chains(specs)
fits.chain <- future_map(chains, ~ do.chain(.))
fits <- cons.fits(fits.chain, chains)
save(fits, file = paste0("output/image-sens-", time.now, ".RData"))
dgs <- map(fits, ~ extract(.)$dg)
dfRaw <- map2_dfr(dgs, l, ~ expand_grid(pol = .y$p$vPol, iter = 1:nrow(.x)) %>% mutate(value = as.vector(.x)), .id = "spec")
df <- dfRaw %>% group_by(spec, pol) %>%
  summarize(estimate = median(value), low = quantile(value, probs=0.025), high = quantile(value, probs=0.975)) %>% ungroup()
xt <- df %>% select(-estimate) %>% pivot_longer(c(low, high)) %>% pivot_wider(names_from = spec) %>% xtable()

print(xt, type = "html", include.rownames = F, include.colnames = F, file = paste0("output/table-sens-", time.now, ".html"))
