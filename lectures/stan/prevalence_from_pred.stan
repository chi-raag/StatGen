data {
  int<lower=0> n;                       // sample size
  int<lower=0, upper=n> T;              // count predicted positive (\hat Y = 1)

  // Beta priors encoded as pseudo-counts (mean m, strength k => a=m*k, b=(1-m)*k)
  real<lower=0> a_pi;                   // prevalence prior alpha
  real<lower=0> b_pi;                   // prevalence prior beta
  real<lower=0> a_tpr;                  // TPR prior alpha
  real<lower=0> b_tpr;                  // TPR prior beta
  real<lower=0> a_fpr;                  // FPR prior alpha
  real<lower=0> b_fpr;                  // FPR prior beta
}

parameters {
  real<lower=0, upper=1> pi;            // prevalence Pr(Y=1)
  real<lower=0, upper=1> tpr;           // True Positive Rate Pr(\hat Y=1 | Y=1)
  real<lower=0, upper=1> fpr;           // False Positive Rate Pr(\hat Y=1 | Y=0)
}

transformed parameters {
  real<lower=0, upper=1> fnr = 1 - tpr; // False Negative Rate
  real<lower=0, upper=1> tnr = 1 - fpr; // True Negative Rate
  // Mixture probability for predicted positives in a representative sample
  real<lower=0, upper=1> p_pos = pi * tpr + (1 - pi) * fpr; // Pr(\hat Y = 1)
}

model {
  // Priors
  pi  ~ beta(a_pi, b_pi);
  tpr ~ beta(a_tpr, b_tpr);
  fpr ~ beta(a_fpr, b_fpr);

  // Likelihood: T predicted positives out of n
  T ~ binomial(n, p_pos);
}

generated quantities {
  // Naive rate and unbiased plug-in estimator for comparison/teaching
  real p_naive = T * 1.0 / n; 
  real pi_hat_unbiased = (p_naive - fpr) / ((tpr - fpr) + 1e-9);
  // Youden's J for stability intuition
  real J = tpr - fpr;
}
