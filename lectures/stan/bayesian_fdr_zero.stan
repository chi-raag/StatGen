data {
  int<lower=1> m;
  vector[m] z;
}
parameters {
  real<lower=0,upper=1> pi0;   // sparsity (proportion null)
  real<lower=0> sigma;         // alt SD (zero-centered)
  real<lower=0> tau;           // null SD (near 1)
}
model {
  // Priors: align with lecture (pi0 >> 1 - pi0; tau ~ 1)
  pi0   ~ beta(20, 1);
  sigma ~ lognormal(0, 0.5);
  tau   ~ lognormal(0, 0.05);

  // Likelihood: two-groups with zero-mean alternative
  for (i in 1:m) {
    real lp0 = log(pi0)     + normal_lpdf(z[i] | 0, tau);
    real lpa = log1m(pi0)   + normal_lpdf(z[i] | 0, sigma);
    target += log_sum_exp(lp0, lpa);
  }
}
generated quantities {
  // Local false discovery rates per SNP: Pr(H=0 | z_i, parameters)
  vector[m] lfdr;
  for (i in 1:m) {
    real log_num   = log(pi0)   + normal_lpdf(z[i] | 0, tau);
    real log_alt   = log1m(pi0) + normal_lpdf(z[i] | 0, sigma);
    real log_denom = log_sum_exp(log_num, log_alt);
    lfdr[i] = exp(log_num - log_denom);
  }
}

