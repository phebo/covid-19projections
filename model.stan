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
  real to_real(int x) { return x; }
  real total_lag_weight(int lag_min, int lag_max, real mu, real sig) {
    real out = 0;
    for(i in lag_min:lag_max) {
      out += exp(normal_lpdf(to_real(i) | mu, sig));
    }
    return out;
  }
}

data {
  int<lower=1> nGeo;
  int<lower=1> nT;
  int<lower=1> nTPred;
  int<lower=1> nPol;
  int<lower=1> nK; // Number of reported event types; k=1 is Death, k=2 is Case
  int<lower=1> vTmax[nGeo];
  int<lower=0> lagMin[nK];
  int<lower=0> lagMax[nK];
  int<lower=0> mRep[nK, nGeo, nT-max(lagMax)];
  int<lower=1,upper=nT+1> mPol[nGeo, nPol];
  real lmaxDeath;
  real<lower=0,upper=1> mortality;
  real<lower=0> mortSig;
  real lr2;
  real<lower=0> phi2;
}

parameters {
  real lir1[nGeo]; // lir[i] at t=1
  real g0[nGeo];  // Base growth rate by geography
  real g0Mu;
  real<lower=0> g0Sig;
  real gMu;
  real<lower=0> gSig;
  vector<upper=0>[nPol] dg[nGeo];  // delta-growth rate due to various policies
  vector<upper=0>[nPol] dgMu;
  cov_matrix[nPol] dgSig;
  real<lower=0> phi;
  real<lower=0,upper=1> p; // Probability of outlier (with rate r2, and dispersion phi2)
  real<lower=0,upper=1> eventFrac[nK, nGeo]; // death & test fraction of total infections
  real<lower=min(lagMin),upper=max(lagMax)> muLag[nK];
  real<lower=0, upper=to_real(max(lagMax) - min(lagMin))/4> sigLag[nK];
  real g0Draw;
  vector[nPol] dgDraw;
}

transformed parameters {
  real lir[nGeo, nT];  // log infection rate
  real lrr[nK, nGeo, nT-max(lagMax)]; // log reported event rate
  real lLagWeight[nK, max(lagMax)];
  real gDraw[nPol + 1];
  
  for(k in 1:nK) {
    real x = log(total_lag_weight(lagMin[k], lagMax[k], muLag[k], sigLag[k]));
    for(t in 1:max(lagMax)) {
      lLagWeight[k,t] = normal_lpdf(to_real(t) | muLag[k], sigLag[k]) - x;
    }
  }
  
  for(i in 1:nGeo){
    lir[i, 1] = lir1[i];
    for(t in 1:(nT-1)){
      lir[i,t+1] = lir[i,t] + g0[i];
      for(l in 1:nPol) {
        if(t+1 >= mPol[i, l]) lir[i, t+1] += dg[i, l];
      }
    }
    for(t in 1:(nT-max(lagMax))) {
      for(k in 1:nK){
        real y[lagMax[k] - lagMin[k] + 1];
        for(dt in lagMin[k]:lagMax[k]){
          y[dt - lagMin[k] + 1] = lir[i, t - dt + max(lagMax)] + lLagWeight[k,dt];
        }
        lrr[k,i,t] = log_sum_exp(y) + log(eventFrac[k,i]);
      }
    }
  }
  for(l in 1:(nPol+1)) {
    gDraw[l] = g0Draw;
    for(l2 in 1:(l-1)) gDraw[l] += dgDraw[l2];
  }
}

model {
  eventFrac[1] ~ lognormal(log(mortality), mortSig);
  g0 ~ normal(g0Mu, g0Sig);

  // Diffuse priors to aid convergence
  phi ~ gamma(1,1);
  p ~ beta(0.001,1);

  for(i in 1:nGeo){
    dg[i] ~ multi_normal(dgMu, dgSig);
    g0[i] + sum(dg[i]) ~ normal(gMu, gSig);

    // Fit to data
    // Mixture model, see https://mc-stan.org/docs/2_22/stan-users-guide/summing-out-the-responsibility-parameter.html
    for(t in 1:(vTmax[i] - max(lagMax))) {
      for(k in 1:nK){
        target += log_sum_exp(log(1-p) + neg_binomial_2_log_lpmf(mRep[k, i, t] | lrr[k, i, t], phi),
                              log(p) + neg_binomial_2_log_lpmf(mRep[k, i, t] | lr2, phi2));
      }
    }
  }
  g0Draw ~ normal(g0Mu, g0Sig);
  dgDraw ~ multi_normal(dgMu, dgSig);
  g0Draw + sum(dgDraw) ~ normal(gMu, gSig);

}

generated quantities{
  real g[nGeo];
  real lirPred[nGeo, nTPred];  // log infection rate
  real lrrPred[nK, nGeo, nTPred]; // log reported event rate
  
  for(i in 1:nGeo) {
    g[i] = lir[i, nT] - lir[i, nT - 1];
    lirPred[i, 1] = lir[i, vTmax[i]] + g[i];
    for (t in 1:(nTPred-1)) {
      lirPred[i, t + 1] = lirPred[i, t] + g[i];
    }
    for(t in 1:nTPred) {
      for(k in 1:nK){
        real y[lagMax[k] - lagMin[k] + 1];
        for(dt in lagMin[k]:lagMax[k]){
          if(t - dt <= 0)
            y[dt - lagMin[k] + 1] = lir[i, t - dt + vTmax[i]] + lLagWeight[k,dt];
          else
            y[dt - lagMin[k] + 1] = lirPred[i, t - dt] + lLagWeight[k,dt];
        }
        lrrPred[k, i, t] = log_sum_exp(y) + log(eventFrac[k,i]);
      }
    }
  }
}
