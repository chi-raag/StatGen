data {
  int<lower=1> n;                 // samples
  int<lower=1> p;                 // predictors (SNPs)
  matrix[n, p] X;                 // standardized genotypes
  array[n] int<lower=0, upper=1> y; // binary outcome

  // Hyperparameters
  real<lower=0> sigma_spike;      // sd of spike (near 0)
  real<lower=0> sigma_slab;       // sd of slab (diffuse)
  real<lower=0> a_pi;             // Beta prior for inclusion probability
  real<lower=0> b_pi;
}

parameters {
  real alpha;                     // intercept
  vector[p] beta;                 // effects
  real<lower=0, upper=1> pi;      // prior inclusion probability
}

transformed parameters {
  vector[n] eta = alpha + X * beta;  // linear predictor
}

model {
  // Priors
  alpha ~ normal(0, 1);
  pi ~ beta(a_pi, b_pi);

  // Spike-and-slab mixture prior on coefficients (marginalized)
  for (j in 1:p) {
    target += log_mix(
      pi,
      normal_lpdf(beta[j] | 0, sigma_slab),
      normal_lpdf(beta[j] | 0, sigma_spike)
    );
  }

  // Probit likelihood via stable normal CDF/CDFc logs
  for (i in 1:n) {
    if (y[i] == 1)
      target += normal_lcdf(eta[i] | 0, 1);   // log Phi(eta)
    else
      target += normal_lccdf(eta[i] | 0, 1);  // log (1 - Phi(eta))
  }
}

generated quantities {
  // Posterior responsibility for each coefficient belonging to the slab
  vector[p] pip; // inclusion probabilities (per draw)
  for (j in 1:p) {
    real log_w1 = log(pi) + normal_lpdf(beta[j] | 0, sigma_slab);
    real log_w0 = log1m(pi) + normal_lpdf(beta[j] | 0, sigma_spike);
    pip[j] = exp(log_w1 - log_sum_exp(log_w1, log_w0));
  }
}
