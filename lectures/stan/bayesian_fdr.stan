data {
  int<lower=1> m;
  vector[m] z;
}
parameters {
  // Mixture weight: prior mass near 1 reflects GWAS sparsity (lecture: pi0 >> pi1)
  real<lower=0,upper=1> pi0;
  // Symmetric alternative shift; constrain to be nonnegative (half-normal prior)
  real<lower=0> mu;
  real<lower=0> sigma;
  real<lower=0> tau;
}
model {
  // Priors
  // Encourage sparsity: prior mean ~ 0.95 with moderate concentration
  pi0   ~ beta(20, 1);
  // Half-normal on mu via lower=0 constraint
  mu    ~ normal(0, 3);
  // Flexible tail under the alternative
  sigma ~ lognormal(0, 0.5);
  // Null spread tightly around 1 to reduce confounding (absorbs mild inflation)
  tau   ~ lognormal(0, 0.05);

  // Likelihood: mixture of null and symmetric alternative
  for (i in 1:m) {
    real lp0 = log(pi0)     + normal_lpdf(z[i] | 0, tau);
    real lpa = log1m(pi0) + (log_sum_exp(normal_lpdf(z[i] | mu,  sigma),
                                         normal_lpdf(z[i] | -mu, sigma)) - log(2));
    target += log_sum_exp(lp0, lpa);
  }
}
generated quantities {
  // Per-draw local false discovery rates: Pr(H=0 | z_i, parameters)
  vector[m] lfdr;
  for (i in 1:m) {
    real log_num   = log(pi0)     + normal_lpdf(z[i] | 0, tau);
    real log_alt   = log1m(pi0) + (log_sum_exp(normal_lpdf(z[i] | mu,  sigma),
                                               normal_lpdf(z[i] | -mu, sigma)) - log(2));
    real log_denom = log_sum_exp(log_num, log_alt);
    lfdr[i] = exp(log_num - log_denom);
  }
}
