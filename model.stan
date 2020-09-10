// Copyright (C) 2020, Phebo Wibbens, Wesley Koo, and Anita McGahan
// 
//   This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

functions {
  real fLag(real[] lInf, vector lpLag, int t0) { // t0 in <lagMax + 1, nT + 1>
    int lagMax = num_elements(lpLag);
    vector[lagMax] out;
    for(t in 1:lagMax) out[t] = lInf[t0 - t] + lpLag[t];
    return log_sum_exp(out);
  }
}

data {
  int<lower=1> lagCaseMax;
  int<lower=1> lagDeathMax;
  int<lower=1> nGeo;
  int<lower=1> nT;
  int<lower=0> nTPred;
  int<lower=1> nPol;
  int<lower=1> nTest;
  int<lower=0, upper=1> mPol[nGeo, nT, nPol];
  int<lower=0, upper=1> mPolChange[nGeo, nT]; // did policy change vs. previous period?
  int<lower=0, upper=1> mPolG1[nGeo, nT];
  int<lower=1, upper=nTest> mTest[nGeo, nT];
  int<lower=0> mCase[nGeo,nT];
  int<lower=0> mDeathRep[nGeo,nT];
  int<lower=-1> mDeathTot[nGeo,nT];
  real<lower=-1> mDeathExp[nGeo,nT];
  real<lower=0> outlCase[2];
  real<lower=0> outlDeath[2];
  real<lower=0,upper=1> mortMu;
  real<lower=0> mortSig;
  real<lower=0, upper=1> pOutl; // probability of outlier for cases
  real<lower=0> idgSig;
}

parameters {
  real logy0[nGeo]; // log infection rate at t=1
  real<lower=0,upper=3> g0[nGeo];  // Base weekly growth rate by geography (with continuous compounding)
  real<lower=-1,upper=2> g1[nGeo];  
  real<lower=0, upper=2> dg[nPol];
  simplex[nTest + 1] caseP[nGeo];
  simplex[nTest + 1] deathP[nGeo];
  real<upper=0> lmortality;
  real deathAdj[nGeo]; // Baseline death rate adjustment vs. historical trends
  real<lower=0> phiCase;
  real<lower=0> phiDeathRep;
  real<lower=0> phiDeathTot;
  simplex[lagCaseMax] pLagCase;
  simplex[lagDeathMax] pLagDeath;
  real<lower=-1,upper=1> idg[nGeo, nT+nTPred-1];
  real<lower=0,upper=1> idgLam1;
  real<lower=0,upper=1> idgLam2f; // Defining idgLam2 as fraction of idgLam1, ensuring idgLam2 <= idgLam1
}

transformed parameters{
  real logy[nGeo, nT+nTPred];  // log infection rate
  real<lower=0,upper=1> fracCase[nGeo, nTest]; // Fraction of infections reported 
  real<lower=0,upper=1> fracDeath[nGeo, nTest]; // Fraction of deaths reported 
  real g[nGeo, nT+nTPred-1];
  real lCaseEst[nGeo, nT + nTPred - lagCaseMax];
  real lDeathEst[nGeo, nT + nTPred - lagDeathMax];
  real lDeathTotEst[nGeo, nT - lagDeathMax];
  vector[lagCaseMax] lpLagCase;
  vector[lagDeathMax] lpLagDeath;
  real<lower=0,upper=1> idgLam2;
  real idgPhi[2];
  real eps[nGeo, nT+nTPred-3];

  lpLagCase = log(pLagCase);
  lpLagDeath = log(pLagDeath);
  idgLam2 = idgLam2f * idgLam1;
  idgPhi[1] = idgLam1 + idgLam2;
  idgPhi[2] = -idgLam1 * idgLam2;

  for(i in 1:nGeo) {
    real dgTot;
    logy[i, 1] = logy0[i];
    for(t in 1:(nT-1)){
      if(mPolChange[i, t]  == 1) {
        dgTot = 0;
        for(p in 1:nPol) if(mPol[i, t, p] == 1) dgTot -= dg[p];
      }
      if(mPolG1[i, t] == 0) g[i, t] = g0[i];
        else g[i, t] = g1[i];
      g[i, t] += dgTot + idg[i, t];
      logy[i, t + 1] = logy[i, t] + g[i, t];
    }
    for(t in nT:(nT + nTPred - 1)) {
      g[i, t] = g1[i] + dgTot + idg[i, t];
      logy[i, t + 1] = logy[i, t] + g[i, t];
    }

    // Fractions of infections reported from unit simplex (ensuring that fracCase[i, l-1] < fracCase[i, l])
    for(l in 1:nTest) {
      fracCase[i, l] = 0;
      fracDeath[i, l] = 0;
      for(l2 in 1:l) {
        fracCase[i, l] += caseP[i, l2];
        fracDeath[i, l] += deathP[i, l2];
      }
    }

    for(t in 1:(nT+nTPred-lagCaseMax)) {
      lCaseEst[i, t] = log(fracCase[i, mTest[i, min(nT, t + lagCaseMax)]]) + fLag(logy[i], lpLagCase, t + lagCaseMax);
    }
    for(t in 1:(nT+nTPred-lagDeathMax)) {
      real lDeath = fLag(logy[i], lpLagDeath, t + lagDeathMax) + lmortality;
      lDeathEst[i, t] = log(fracDeath[i, mTest[i, min(nT, t + lagDeathMax)]]) + lDeath;
      if(t + lagDeathMax <= nT) {
        if(mDeathExp[i, t + lagDeathMax] != -1) lDeathTotEst[i, t] = log(mDeathExp[i, t + lagDeathMax] + deathAdj[i] + exp(lDeath));
          else lDeathTotEst[i,t] = 0;
      }
    }
    for(t in 3:(nT+nTPred-1)) eps[i, t - 2] = idg[i, t] - idgPhi[1] * idg[i, t-1] - idgPhi[2] * idg[i, t-2];
  }
}

model {
  // Hyperprior for mortality
  lmortality ~ normal(log(mortMu), mortSig);
  
  // AR(2) model for idiosyncratic growth rate
  for(i in 1:nGeo) for(t in 1:(nT+nTPred-3)) eps[i, t] ~ normal(0, idgSig);

  // Likelihood for observations
  for(i in 1:nGeo) {
    for(t in (lagCaseMax + 1):nT) {
      target += log_sum_exp(log(1-pOutl) + neg_binomial_2_log_lpmf(mCase[i, t] | lCaseEst[i,t - lagCaseMax], phiCase),
                            log(pOutl) + neg_binomial_2_lpmf(mCase[i, t] | outlCase[1], outlCase[2]));
    }
    for(t in (lagDeathMax + 1):nT) {
      target += log_sum_exp(log(1-pOutl) + neg_binomial_2_log_lpmf(mDeathRep[i, t] | lDeathEst[i,t - lagDeathMax], phiDeathRep),
                            log(pOutl) + neg_binomial_2_lpmf(mDeathRep[i, t] | outlDeath[1], outlDeath[2]));
      if(mDeathTot[i, t] != -1 && mDeathExp[i, t] != -1) {
        mDeathTot[i, t] ~ neg_binomial_2_log(lDeathTotEst[i, t - lagDeathMax], phiDeathTot);
      }
    }
  }
}

generated quantities{
  int<lower=0> predCase[nGeo, nTPred];
  int<lower=0> predDeath[nGeo, nTPred];
  
  for(i in 1:nGeo) {
    for(t in 1:nTPred) {
      predCase[i, t] = neg_binomial_2_log_rng(lCaseEst[i, t + nT - lagCaseMax], phiCase);
      predDeath[i, t] = neg_binomial_2_log_rng(lDeathEst[i, t + nT - lagDeathMax], phiDeathRep);
    }
  }
}
