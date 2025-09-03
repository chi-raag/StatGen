############################################################
# PUBH 8878 — Lab 02: Linear modeling for heritability concepts
# Topic: Heritability (A/C/E), liability-threshold, ascertainment.
#
# Theme: Use familiar modeling tools (lm, glm with probit, optional REML via lme4)
# to answer questions. Avoid custom MLE implementations from scratch.
#
# You will:
#   1) Build a simple pedigree and A (2\phi) matrix
#   2) Simulate quantitative and binary traits on the pedigree
#   3) Estimate h^2 via:
#       - Parent–offspring regression (lm)
#       - Twins (MZ/DZ) via regression slopes
#       - Optional: REML LMM using a package (no from-scratch code)
#   4) For binary traits, fit probit GLMs and apply scale conversion
#   5) Examine ascertainment bias using probit GLMs
############################################################

## ---- 0) Setup ---------------------------------------------------------------

library(kinship2) # pedigree and kinship matrix
library(mvtnorm) # multivariate normal probabilities
library(dplyr) # data manipulation
library(lme4) # optional LMM (if installed)
library(sandwich) # robust SEs

set.seed(8878)
options(stringsAsFactors = FALSE)

## ---- 1) Build a pedigree and the additive relationship matrix A ------------

# We simulate many unrelated nuclear families (founders unrelated across families).
# kinship2::pedigree() expects: id, dadid, momid, sex (1=male, 2=female).

make_nuclear_families <- function(n_fam = 100, kids_min = 2, kids_max = 4) {
  stopifnot(kids_min >= 1, kids_max >= kids_min)
  ids <- character()
  dadid <- character()
  momid <- character()
  sex <- integer()

  famid <- character()
  role <- character()

  for (f in seq_len(n_fam)) {
    nk <- sample(kids_min:kids_max, size = 1)
    dad <- sprintf("F%03d", f)
    mom <- sprintf("M%03d", f)
    kids <- sprintf("K%03d_%02d", f, seq_len(nk))

    id_vec <- c(dad, mom, kids)
    dad_vec <- c(NA, NA, rep(dad, nk))
    mom_vec <- c(NA, NA, rep(mom, nk))
    sex_vec <- c(1L, 2L, sample(c(1L, 2L), nk, replace = TRUE))

    ids <- c(ids, id_vec)
    dadid <- c(dadid, dad_vec)
    momid <- c(momid, mom_vec)
    sex <- c(sex, sex_vec)

    famid <- c(famid, rep(sprintf("Fam%03d", f), length(id_vec)))
    role <- c(role, c("father", "mother", rep("child", nk)))
  }

  ped_df <- data.frame(
    id = ids, dadid = dadid, momid = momid, sex = sex,
    fam = famid, role = role, stringsAsFactors = FALSE
  )
  ped_df
}

ped_df <- make_nuclear_families(n_fam = 120, kids_min = 2, kids_max = 4)

# Build a kinship2 pedigree object and compute the kinship matrix Φ.
ped_obj <- with(ped_df, pedigree(id = id, dadid = dadid, momid = momid, sex = sex))
Phi <- kinship(ped_obj) # kinship coefficients
A <- 2 * as.matrix(Phi) # additive relationship
dim(A) # sanity check

# Index helpers
id <- ped_df$id
n <- length(id)

# Add a continuous covariate (e.g., age) and a binary covariate (sex)
set.seed(42)
ped_df$age <- scale(rnorm(n, mean = 50, sd = 10))
ped_df$sex01 <- as.integer(ped_df$sex == 1L) # male=1, female=0

## ---- 2) Simulate a quantitative trait from the pedigree ---------------------

true_beta <- c(`(Intercept)` = 2.0, age = 0.10, sex01 = -0.25)
true_h2 <- 0.60
true_VA <- true_h2
true_VE <- 1 - true_h2

X <- cbind(Intercept = 1, age = ped_df$age, sex01 = ped_df$sex01)
V <- true_VA * A + true_VE * diag(n)

# Simulate via Cholesky
R <- chol(V) # upper triangular s.t. V = R'R
z <- rnorm(n)
y <- as.vector(X %*% true_beta + t(R) %*% z)

ped_df$y <- y

## ---- 3) Linear-model estimators (no custom MLE) -----------------------------

# A) Parent–offspring regression (lm)
#  - Build a parent–offspring dataset from ped_df
#  - Regress offspring y on parent y; h2_hat ≈ 2 * slope
# TODO: construct a data.frame with columns: fam, parent_y, child_y (one child per family or all children)
# TODO: fit: lm(child ~ midparent)

# TODO: compute h2_hat as coef(fit)[2]
# TODO: fit mixed model with lmer4, compare confidence intervals using confint

# B) Twins (regression on standardized phenotypes)
#  - Create synthetic MZ/DZ pairs from the simulated pedigree or simulate separately
#  - Standardize y within zygosity; fit lm(y1_z ~ y2_z) in MZ and DZ subsets to estimate r_MZ, r_DZ
#  - Compute h2_hat = 2*(r_MZ - r_DZ); c2_hat = 2*r_DZ - r_MZ; e2_hat = 1 - h2_hat - c2_hat

# Parameters for twin simulation (feel free to change)
VA <- 0.4 # additive genetic variance
VC <- 0.3 # shared environment variance
VE <- 0.3 # individual environment variance
VP <- VA + VC + VE
n_MZ <- 1500
n_DZ <- 1500

set.seed(1202)

# Helper to simulate one cohort of twin pairs given genetic correlation rA
simulate_twins <- function(n_pairs, rA, VA, VC, VE, label) {
  # Shared C per pair
  C <- rnorm(n_pairs, mean = 0, sd = sqrt(VC))
  # Additive genetic A1, A2 with Corr = rA
  if (rA == 1) {
    A_shared <- rnorm(n_pairs, mean = 0, sd = sqrt(VA))
    A1 <- A_shared
    A2 <- A_shared
  } else {
    SigmaA <- matrix(c(VA, rA * VA, rA * VA, VA), nrow = 2)
    A <- mvtnorm::rmvnorm(n_pairs, sigma = SigmaA)
    A1 <- A[, 1]
    A2 <- A[, 2]
  }
  # Individual-specific E
  E1 <- rnorm(n_pairs, mean = 0, sd = sqrt(VE))
  E2 <- rnorm(n_pairs, mean = 0, sd = sqrt(VE))
  # Phenotypes
  Y1 <- A1 + C + E1
  Y2 <- A2 + C + E2
  data.frame(
    zygosity = label,
    pair = seq_len(n_pairs),
    y1 = Y1,
    y2 = Y2
  )
}

tw_MZ <- simulate_twins(n_MZ, rA = 1.0, VA = VA, VC = VC, VE = VE, label = "MZ")
tw_DZ <- simulate_twins(n_DZ, rA = 0.5, VA = VA, VC = VC, VE = VE, label = "DZ")
tw <- bind_rows(tw_MZ, tw_DZ)

# Standardize within zygosity
tw <- tw |>
  group_by(zygosity) |>
  mutate(y1_z = scale(y1)[, 1], y2_z = scale(y2)[, 1]) |>
  ungroup()

# TODO: get correlation for r_MZ

# TODO: get correlation for r_DZ

# Compute h^2 = 2*(r_MZ - r_DZ)

# compute bootstrap CIs for h2_hat



# Optional REML using a package
#  - If available, demonstrate lme4::lmer with a family random intercept (captures C)
#  - For kinship-based A random effects, consider coxme::lmekin (if installed)
# if (requireNamespace("lme4", quietly = TRUE)) {
#   library(lme4)
#   # Example: random intercept for family (shared environment):
#   # lmer(y ~ age + sex01 + (1|fam), data = ped_df, REML = TRUE)
# }


## ---- 4) Binary trait: probit GLM + recurrence risk -------------------------

# Simulate liabilities with correlation structure h_ell^2 * A + (1 - h_ell^2) * I,
# threshold at T = qnorm(1-K) to get prevalence K.

simulate_binary_from_liability <- function(A, h2_liab = 0.5, K = 0.10) {
  n <- nrow(A)
  V <- h2_liab * A + (1 - h2_liab) * diag(n) # variance normalized to 1
  R <- chol(V)
  L <- as.vector(t(R) %*% rnorm(n)) # liabilities
  T <- qnorm(1 - K)
  ybin <- as.integer(L > T)
  list(L = L, y = ybin, T = T, K_emp = mean(ybin))
}

set.seed(91)
h2_liab_true <- 0.50
K_true <- 0.10
bin_dat <- simulate_binary_from_liability(A, h2_liab = h2_liab_true, K = K_true)
ped_df$case <- bin_dat$y

# Probit GLM sanity check (intercept-only): prevalence via probit link
# TODO: fit <- glm(case ~ 1, family = binomial(link = "probit"), data = ped_df)
# TODO: recover prevalence_hat = 1 - pnorm(-coef(fit)[1]) and compare to mean(case)

# compare to empirical prevalence
