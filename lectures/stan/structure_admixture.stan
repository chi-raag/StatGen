data {
  int<lower=1> N;                                // individuals
  int<lower=1> L;                                // loci
  int<lower=1> K;                                // ancestral pops
  array[N, L] int<lower=0, upper=2> y;           // genotypes (0,1,2)
  int<lower=0> L_soft;
  array[L_soft] int<lower=1, upper=L> soft_idx;  // e.g., {2}
  array[K] real<lower=0> a_theta;                // Beta 'a' for non-anchor loci
  array[K] real<lower=0> b_theta;                // Beta 'b' for non-anchor loci
  int<lower=1, upper=L> l_star;                  // anchor locus index
  real<lower=0> conc_pi;                         // shared Dirichlet conc. for pi
}

parameters {
  array[N] simplex[K] pi;                        // ancestry proportions
  ordered[K] eta;                                // ordered logits at anchor locus
  matrix<lower=0, upper=1>[K, L-1] theta_rest;   // allele freqs for non-anchor loci
}

transformed parameters {
  matrix<lower=0, upper=1>[K, L] theta;          // full allele-frequency matrix
  matrix[N, L] p_mix;                            // mixed allele frequency

  // anchor column (ordered)
  for (k in 1:K) theta[k, l_star] = inv_logit(eta[k]);

  // fill remaining columns
  {
    int c = 1;
    for (l in 1:L) {
      if (l == l_star) continue;
      for (k in 1:K) theta[k, l] = theta_rest[k, c];
      c += 1;
    }
  }

  // mixture expectations
  for (n in 1:N)
    for (l in 1:L)
      p_mix[n, l] = dot_product(pi[n], col(theta, l));
}

model {
  // priors
  eta ~ normal(0, 2.5);                          // weak prior; ordering gives ID
  for (k in 1:K)
    for (c in 1:(L - 1))
      theta_rest[k, c] ~ beta(a_theta[k], b_theta[k]);

  for (n in 1:N)
    pi[n] ~ dirichlet(rep_vector(conc_pi, K));

  for (c in 1:L_soft) {
  int l = soft_idx[c];
  if (l != l_star) {
    // comp 1 LOW, comp 2 HIGH at these loci (gentle)
    target += beta_lpdf(theta[1, l] | 2, 8);
    target += beta_lpdf(theta[2, l] | 8, 2);
  }
}

  // likelihood
  for (n in 1:N)
    for (l in 1:L)
      y[n, l] ~ binomial(2, p_mix[n, l]);
}

generated quantities {
  matrix[N, K] logit_pi;
  array[N, L] int y_rep;

  for (n in 1:N)
    for (k in 1:K)
      logit_pi[n, k] = logit(pi[n, k]);

  for (n in 1:N)
    for (l in 1:L)
      y_rep[n, l] = binomial_rng(2, p_mix[n, l]);
}
