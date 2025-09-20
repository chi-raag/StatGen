// STRUCTURE-style admixture model (fast classroom version)
// Bi-allelic markers; y[n, l] counts reference alleles (0, 1, 2)

data {
    int<lower=1> N;                               // individuals
    int<lower=1> L;                               // loci (bi-allelic)
    int<lower=1> K;                               // ancestral populations
    array[N, L] int<lower=0, upper=2> y;          // genotype matrix
    vector<lower=0>[K] alpha;                     // Dirichlet concentration for ancestries
    real<lower=0> a_theta;                        // Beta prior shape for theta
    real<lower=0> b_theta;                        // Beta prior shape for theta
}

parameters {
    array[N] simplex[K] pi;                       // individual ancestry proportions
    matrix<lower=0, upper=1>[K, L] theta;         // population-specific allele frequencies
}

transformed parameters {
    matrix[N, L] p_mix;                           // ancestry-weighted allele frequencies per individual
    for (n in 1:N) {
        for (l in 1:L) {
            real p = dot_product(pi[n], theta[, l]);
            p_mix[n, l] = fmin(1 - 1e-9, fmax(1e-9, p)); // guard in (0,1)
        }
    }
}

model {
    // Priors
    for (n in 1:N)
        pi[n] ~ dirichlet(alpha);
    for (k in 1:K)
        for (l in 1:L)
            theta[k, l] ~ beta(a_theta, b_theta);

    // Likelihood
    for (n in 1:N)
        for (l in 1:L)
            y[n, l] ~ binomial(2, p_mix[n, l]);
}

generated quantities {
    matrix[N, K] logit_pi;
    matrix[N, L] log_lik;                        // pointwise log-likelihood

    for (n in 1:N)
        for (k in 1:K)
            logit_pi[n, k] = logit(pi[n, k]);

    for (n in 1:N)
        for (l in 1:L)
            log_lik[n, l] = binomial_lpmf(y[n, l] | 2, p_mix[n, l]);
}
