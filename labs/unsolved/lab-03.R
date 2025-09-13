# Lab 03
# PUBH 8878 — Statistical Genetics
# ---------------------------------------------------------------
# Goals (aligned to Lecture 03):
# 1) Implement EM for ABO gene frequencies from phenotype counts
# 2) Two-point linkage: LOD and EM when EM is ACTUALLY necessary (partially informative mate)
# 3) Two-SNP haplotype EM under HWE + LD (D, r^2)
# 4) Mini single-marker association with basic QC
# ---------------------------------------------------------------

set.seed(8878)


library(dplyr)
library(ggplot2)
library(tibble)

# ---------------------------------------------------------------
# 1) EM for ABO gene frequencies
# ---------------------------------------------------------------

# Observed phenotype counts (from slides)
counts <- list(nA = 725, nAB = 72, nB = 258, nO = 1073)
N <- sum(unlist(counts))

# TODO 1.1: Write a function loglik_pheno(pA, pB) returning the log-likelihood
# under the phenotype model given global 'counts'.
# Hints:
#   pO = 1 - pA - pB (require pA>0, pB>0, pO>0)
#   P(A)  = pA^2 + 2 pA pO
#   P(AB) = 2 pA pB
#   P(B)  = pB^2 + 2 pB pO
#   P(O)  = pO^2
loglik_pheno <- function(pA, pB) {
  # TODO: implement and return a numeric scalar
}

# TODO 1.2: Implement EM to estimate (pA, pB, pO) given counts.
# Return a tibble with columns: iter, pA, pB, pO, logLik
# Tips:
#  - Track the observed-data log-likelihood each iteration (monotone ↑).
#  - Stop when max parameter change < tol. Also guard against invalid domains.
#  - Consider returning the whole trajectory for plotting.
em_abo <- function(pA0, pB0, max_iter = 200, tol = 1e-6) {
  # TODO: implement E- and M-steps per lecture and return trajectory
}

# TODO 1.3: Run EM from three starting values and compare convergence:
#  - Equal:     (1/3, 1/3)
#  - Unbalanced:(0.01, 0.98)
#  - Smart:     pO0 = sqrt(nO/N); use S=1-pO0 and c=nAB/(2N) to solve for (pA0, pB0)
# Produce:
#  - A small table of converged (pA, pB, pO) and iteration counts


# ---------------------------------------------------------------
# 2) Two-point linkage: LOD + EM (Transmitting parent Aa/Bb unknown phase; mate Aa/bb)
# ---------------------------------------------------------------

# We observe per child the two-locus genotypes (A-genotype, B-genotype).
# The mate is Aa/bb: the B allele is fully informative (comes from the transmitting parent),
# but when the child is Aa at A, the transmitting parent's A contribution is ambiguous.
#
# Use the six observable categories and their counts:
#   definitive: (AA,Bb), (AA,bb), (aa,Bb), (aa,bb)
#   ambiguous : (Aa,Bb), (Aa,bb)

counts2 <- c(
  n_AA_Bb = 12,
  n_AA_bb = 3,
  n_aa_Bb = 3,
  n_aa_bb = 12,
  n_Aa_Bb = 15,
  n_Aa_bb = 15
)
stopifnot(sum(counts2) == 60L)

# TODO 2.1: Write a function probs_partial(theta, w) that returns the
# category probabilities (named vector of length 6) under:
#  - theta in (0, 0.5)
#  - w = Pr(coupling phase) in [0,1]
#
# Hints:
#  Under coupling (AB/ab): hap probs from the transmitting parent are
#     P(AB)=(1-theta)/2, P(Ab)=theta/2, P(aB)=theta/2, P(ab)=(1-theta)/2.
#  Under repulsion (Ab/aB): swap AB<->Ab and ab<->aB in the above.
#  The mate transmits A vs a with prob 0.5 independently, which yields:
#    P(AA,Bb) = 0.5 * P(AB)
#    P(aa,Bb) = 0.5 * P(aB)
#    P(AA,bb) = 0.5 * P(Ab)
#    P(aa,bb) = 0.5 * P(ab)
#    P(Aa,Bb) = 0.5 * (P(AB) + P(aB))
#    P(Aa,bb) = 0.5 * (P(Ab) + P(ab))
probs_partial <- function(theta, w) {
  # TODO: implement and return a named vector with the six probabilities
}

# TODO 2.2: Write a function loglik_partial(theta, w, counts2) that computes the
# observed-data log-likelihood for the six categories using probs_partial.
loglik_partial <- function(theta, w, counts2) {
  # TODO: implement using safelog() and return a scalar
}

# TODO 2.3: Implement EM for (theta, w) with ambiguous categories.
# Follow the lecture: split ambiguous counts in the E-step using
#    q_AB = w*(1-theta) + (1-w)*theta  (prob an Aa,Bb child came from AB vs aB)
#    q_Ab = w*theta + (1-w)*(1-theta)  (prob an Aa,bb child came from Ab vs ab)
# Then form expected haplotype counts, compute R/NR under each phase, update w and theta.
em_twopoint_partial <- function(counts2, theta0 = 0.4, w0 = 0.5, max_iter = 200, tol = 1e-8) {
  # TODO: implement; return a tibble with iter, theta, w, logLik
}

# TODO 2.4: Run EM, report (theta_hat, w_hat), and compute the LOD vs theta=0.5:
#   LOD = log10{ L(theta_hat,w_hat) / L(theta=0.5) }.
# (At theta=0.5 the phase mixture is irrelevant.) Also, make a small convergence plot.


# ---------------------------------------------------------------
# 3) Two-SNP haplotype EM + LD (D, r^2)
# ---------------------------------------------------------------

# Simulate genotypes for two SNPs from true hap frequencies
p_true <- c(ab = .40, aB = .10, Ab = .25, AB = .25)
n <- 1000

draw_hap <- function(p) {
  sample(names(p), 2, replace = TRUE, prob = unname(p))
}

hap_to_vec <- function(h) {
  map1 <- list(ab = c(0, 0), aB = c(0, 1), Ab = c(1, 0), AB = c(1, 1))
  v1 <- map1[[h[1]]]
  v2 <- map1[[h[2]]]
  c(g1 = v1[1] + v2[1], g2 = v1[2] + v2[2])
}

haps <- replicate(n, draw_hap(p_true))
G <- t(apply(haps, 2, hap_to_vec))
G <- as_tibble(G)

# TODO 3.1: Implement EM for two-SNP haplotypes under HWE
# p0 is a named vector c(ab=..., aB=..., Ab=..., AB=...)
# Hints:
#  - Only (g1=1,g2=1) is ambiguous between {AB+ab} vs {Ab+aB}.
#  - E-step: weight the two resolutions by p_AB*p_ab vs p_Ab*p_aB (the factor of 2 cancels).
#  - M-step: normalize expected hap counts to sum to 1.
em_hap2 <- function(G, p0 = rep(0.25, 4), tol = 1e-8, maxit = 500) {
  # TODO: implement E-step expected hap counts per person and M-step normalize
}

# TODO 3.2: Estimate p via EM and compare to p_true (print and plot).
# TODO 3.3: Compute LD measures from p_hat: D and r^2
#   Hints: pA = p(Ab)+p(AB); pB = p(aB)+p(AB); p11 = p(AB); D = p11 - pA*pB;
#          r2 = D^2 / (pA*(1-pA)*pB*(1-pB)).


# ---------------------------------------------------------------
# 4) Mini association demo + QC
# ---------------------------------------------------------------

# Simulate one SNP under HWE and a binary trait with covariate
set.seed(1203)
n <- 2000
maf <- 0.30
geno <- rbinom(n, 2, maf) # 0/1/2 copies
age <- rnorm(n, 50, 10)
linpred <- -3.0 + 0.40 * scale(age)[, 1] + 0.35 * geno # log-odds
pr <- plogis(linpred)
case <- rbinom(n, 1, pr)
dat <- tibble(geno, age, case)

# TODO 4.1: HWE test on controls only (subset case==0)
#  - compute counts (n_AA, n_Aa, n_aa) and a chi-squared HWE p-value
#  - function hwe_chisq(nAA, nAa, naa): return list(stat, pval)

# TODO 4.2: Fit logistic regression case ~ geno + age and report
#  - OR for geno (per allele) and 95% CI
#  - Hint: glm(case ~ geno + scale(age), family=binomial)

# TODO 4.3: Apply QC rule: if HWE p < 1e-6 (in controls), flag the SNP.

# End of Lab 03
