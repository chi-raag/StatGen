// ABO phenotype model with latent genotype aggregation
// We model allele frequencies p = (pA, pB, pO) with a Dirichlet prior.
// Phenotype counts (nA, nAB, nB, nO) arise from genotype frequencies aggregated to phenotype.
// Genotype frequencies under HWE:
//  AA: pA^2
//  AO: 2 pA pO
//  BB: pB^2
//  BO: 2 pB pO
//  AB: 2 pA pB
//  OO: pO^2
// Phenotype probabilities:
//  A  : pA^2 + 2 pA pO
//  B  : pB^2 + 2 pB pO
//  AB : 2 pA pB
//  O  : pO^2
// We use a multinomial likelihood on phenotype counts.

functions {
  vector abo_pheno_probs(vector p) {
    // p[1]=pA, p[2]=pB, p[3]=pO
    real pA = p[1];
    real pB = p[2];
    real pO = p[3];
    vector[4] q;
    q[1] = pA * pA + 2 * pA * pO; // A
    q[2] = 2 * pA * pB;            // AB
    q[3] = pB * pB + 2 * pB * pO; // B
    q[4] = pO * pO;               // O
    return q;
  }
}

data {
  int<lower=0> nA;
  int<lower=0> nAB;
  int<lower=0> nB;
  int<lower=0> nO;
  vector<lower=0>[3] alpha; // Dirichlet prior hyperparameters
}

transformed data {
  int<lower=0> N = nA + nAB + nB + nO;
  array[4] int y = { nA, nAB, nB, nO };
}

parameters {
  simplex[3] p; // allele frequencies (pA, pB, pO)
}

transformed parameters {
  vector[4] q = abo_pheno_probs(p);
}

model {
  // Prior
  p ~ dirichlet(alpha);
  // Likelihood
  y ~ multinomial(q);
}

generated quantities {
  // Log likelihood for WAIC/LOO
  real log_lik = multinomial_lpmf(y | q);
  // Posterior predictive phenotype probabilities (already q)
  // Genotype frequencies if needed
  vector[6] geno_freq;
  real pA = p[1];
  real pB = p[2];
  real pO = p[3];
  geno_freq[1] = pA * pA;         // AA
  geno_freq[2] = 2 * pA * pO;     // AO
  geno_freq[3] = pB * pB;         // BB
  geno_freq[4] = 2 * pB * pO;     // BO
  geno_freq[5] = 2 * pA * pB;     // AB
  geno_freq[6] = pO * pO;         // OO
}
