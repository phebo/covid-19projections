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


clean.data <- function(
  dfJh, dfEcon, dfPop, dfOx, dfHol,
  minPop = 5e6, geoExclude = NULL,
  dates = c(as.Date("2020-02-01"), Inf), nTPred = 0,
  polG1 = c("C1 - 1", "C2 - 1", "C3 - 1", "C4 - 1"),
  polExcl = c("C1 - 1", "C2 - 1", "C3 - 1", "C4 - 1", "C5 - 2", "C6 - 3", "C8 - 1"),
  lagCaseMax = 2, lagDeathMax = 4,
  mortMu = 0.01, mortSig = 0.5, # Parameters for log-normal distribution; mortSig = 0.5 means 95% interval of */exp(1.96*0.5)=2.6
  pOutl = 1e-3, # Probability of outlier (lower probability attaches more weight to extreme data points)
  idgSig = 0.02, # s.d. of change in idiosyncratic growth rate (AR(2) process)
  idgLam = c(0.9, 0.9),
  dgSig = c(0, 0), # If non-zero: use spike-and-slab-like prior standard deviations on dg (uninformative, "zero")
  dgMin = 0
) {

  vGeo <- sort(dfPop %>% filter(population >= minPop, geo != "US - New York City") %>% pull(geo))
  vGeo <- vGeo[!vGeo %in% geoExclude]
  
  vDate <- sort(unique(dfJh$date))
  vDate <- vDate[as.numeric(vDate) %% 7 == as.numeric(min(dfEcon$date)) %% 7]
  vDate <- vDate[vDate >= dates[1] & vDate <= dates[2]]
  stopifnot(as.numeric(vDate - lag(vDate))[-1] == 7)
  
  # Epidemiology data
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
  dfPHol <-
    bind_rows(
      dfHol %>% mutate(valueHol = ifelse(level >= 2, 1, 0), level = 1),
      dfHol %>% mutate(valueHol = ifelse(level >= 2, 1, 0), level = 2),
      dfHol %>% mutate(valueHol = ifelse(level == 3, 1, 0), level = 3)) %>%
    full_join(expand_grid(geo = vGeo, level = 1:3, date = vDate)) %>%
    arrange(geo, level, date) %>% group_by(geo, level) %>% fill(valueHol) %>% ungroup() %>%
    mutate(polCode = "C1", level = as.character(level))
  dfP <- left_join(dfP, dfPHol) %>%
    mutate(valueOx = value, value = ifelse(is.na(valueHol), valueOx, pmax(valueHol, valueOx)))
  stopifnot(dfP %>% group_by(geo, polCode, date) %>% mutate(test = value <= lag(value, default = 1)) %>% pull(test))
  
  dfPCor <- expand_grid(pol1 = vPolAll, pol2 = vPolAll) %>%
    mutate(cor = as.vector(cor(as.matrix(dfP %>% select(geo, pol, date, value) %>% pivot_wider(names_from = pol) %>% select(-c(geo, date))))))
  stopifnot(nrow(dfP) == length(vGeo) * length(vDate) * length(vPolAll), !is.na(dfP$value))
  
  mPol <- dfP %>% filter(!polS %in% c(polExcl)) %>%
    mutate(date = factor(date, levels = as.character(vDate))) %>%
    xtabs(value ~ geo + date + pol, .)
  vPol <- dimnames(mPol)$pol
  stopifnot(dimnames(mPol)$date == vDate, dimnames(mPol)$geo == vGeo,
            length(vPolAll) == length(vPol) + length(polExcl))
  mPolChange <- matrix(as.numeric(cbind(rep(TRUE, length(vGeo)), apply(mPol[,-length(vDate),] != mPol[,-1,], c(1,2), any))), nrow=length(vGeo))
  
  mPolG1 <- dfP %>% filter(polS %in% polG1) %>%
    group_by(geo, date) %>% summarize(value = sum(value)) %>% ungroup() %>%
    mutate(value = ifelse(value >= 1, 1, NA)) %>%
    group_by(geo) %>% fill(value) %>% ungroup() %>%
    mutate(value = ifelse(is.na(value), 0, 1)) %>%
    xtabs(value ~ geo + date, .)
  stopifnot(rownames(mPolG1) == vGeo, colnames(mPolG1) == vDate, apply(mPolG1[,-1] - mPolG1[,-ncol(mPolG1)], 1, FUN = min) == 0)

  dfTest <- dfOxSub %>% filter(polCode == "H2") %>%
    full_join(expand_grid(geo = vGeo, date = vDate)) %>% arrange(geo, date) %>%
    group_by(geo) %>% fill(level) %>% ungroup() %>% filter(date %in% vDate)
  mTest <- xtabs(level ~ geo + date, dfTest)
  mTest[mTest == 0] <- 1
  stopifnot(rownames(mTest) == vGeo, colnames(mTest) == vDate)
  
  # Estimate outlier distributions:
  nll <- function(lLam, lPhi) -sum(dnbinom(y, mu = exp(lLam), size = exp(lPhi), log=T))
  
  y <- as.vector(mCase)
  fit <- mle(nll, start = list(lLam = 0, lPhi = 0))
  outlCase <- exp(coef(fit))
  
  y <- as.vector(mDeathRep)
  fit <- mle(nll, start = list(lLam = 0, lPhi = 0))
  outlDeath <- exp(coef(fit))
  
  lData = list(mortMu = mortMu, mortSig = mortSig, lagDeathMax = lagDeathMax, lagCaseMax = lagCaseMax,
               nGeo = length(vGeo), nT = length(vDate), nTPred = nTPred, nPol = length(vPol), nTest = max(mTest),
               mPol = mPol, mPolChange = mPolChange, mPolG1 = mPolG1, mTest = mTest,
               mCase = mCase, mDeathRep = mDeathRep, mDeathTot = mDeathTot, mDeathExp = mDeathExp, outlCase = outlCase, 
               outlDeath = outlDeath, pOutl = pOutl, idgSig = idgSig, idgLam = idgLam, dgSig = dgSig, dgMin = dgMin)

  list(dfE = dfE, dfP = dfP, dfPCor = dfPCor, dfTest = dfTest, lData = lData,
       p = list(vDate = vDate, vGeo = vGeo, vPol = vPol, minPop = minPop, polG1 = polG1, polExcl = polExcl,
                mortMu = mortMu, mortSig = mortSig, pOutl = pOutl, idgSig = idgSig, idgLam = idgLam, dgSig = dgSig, dgMin = dgMin))
}

make.chains <- function(models){
  # Makes a list of chains from a list of model specifications
  # Requires model$chains with #of chains for that model
  chains <- list()
  n <- 1
  for(model in models){
    for(i in 1:model$chains){
      new.chain <- model
      new.chain$chain.id <- i
      new.chain$model.id <- n
      chains <- c(chains, list(new.chain))
    }
    n <- n + 1
  }
  chains
}

do.chain <- function(chain){
  require(rstan)
  sampling(chain$m, pars=chain$pars, data=chain$data, 
       chains=1, chain_id = chain$chain.id, seed=99743, iter=chain$iter, warmup=chain$warmup, thin=chain$thin,
       control = list(adapt_delta = 0.9, max_treedepth = 12))
}

cons.fits <- function(fits.chain, chains){
  # Consolidates multiple chains by fit into model fits
  model.ids <- sapply(chains, function(chain) chain$model.id)
  fits <- list()
  for(i in unique(model.ids)){
    fits[[i]] <- sflist2stanfit(fits.chain[model.ids == i])
  }
  fits
}

