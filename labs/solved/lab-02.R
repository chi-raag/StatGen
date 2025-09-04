library(kinship2) # pedigree and kinship matrix
library(mvtnorm) # multivariate normal probabilities
library(dplyr) # data manipulation
library(lme4qtl) # LMM with custom relationship matrices
library(lme4)
set.seed(8878)

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
        # Use NA for unknown parents (kinship2 expects 0 or NA; avoid "0" as a string)
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
Phi <- kinship(ped_obj) # kinship coefficients φ_ij
A <- 2 * as.matrix(Phi) # additive relationship A = 2Φ
dim(A) # sanity check

# Index helpers
id <- ped_df$id
n <- length(id)

# Add a continuous covariate (e.g., age) and a binary covariate (sex)
set.seed(42)
ped_df$age <- scale(rnorm(n, mean = 50, sd = 10))
ped_df$sex01 <- as.integer(ped_df$sex == 1L) # male=1, female=0

## ---- 2) Simulate a quantitative trait from the pedigree ---------------------

# We simulate Y ~ N(Xβ, V) with V = σ_A^2 * A + σ_E^2 * I.
# Choose true parameters and compute h^2 = σ_A^2 / (σ_A^2 + σ_E^2).

true_beta <- c(`(Intercept)` = 2.0, age = 0.10, sex01 = -0.25)
true_h2 <- 0.60
true_VA <- true_h2
true_VE <- 1 - true_h2

X <- cbind(Intercept = 1, age = ped_df$age, sex01 = ped_df$sex01)
V <- true_VA * A + true_VE * diag(n)

# Simulate via Cholesky
R <- chol(V) # upper triangular s.t. V = R'R
z <- rnorm(n)
y <- as.vector(X %*% true_beta + t(R) %*% z) # note: R' * z is N(0, V)

ped_df$y <- y

## ---- 3) Linear-model estimators -----------------------------

# A) Parent–offspring regression (lm)
#  - Build a parent–offspring dataset from ped_df
#  - Regress offspring y on parent y; h2_hat ≈ 2 * slope
# TODO: construct a data.frame with columns: fam, parent_y, child_y (one child per family or all children)

parents <- ped_df %>%
    filter(role %in% c("father", "mother")) %>%
    group_by(fam) %>%
    filter(n() == 2) %>% # require both parents
    summarise(midparent = mean(y, .groups = "drop"))

children <- ped_df %>%
    filter(role == "child") %>%
    transmute(fam, child = y)

df_child <- children %>%
    left_join(parents, by = "fam") %>%
    filter(!is.na(midparent))

# TODO: fit: lm(child ~ midparent)
fit <- lm(child ~ midparent, data = df_child)

# TODO: compute confidence intervals for h2_hat
confint(fit, "midparent", level = 0.95)

# TODO: fit mixed model with lme4qtl (retmatLmer)

# The ID for the random effect must be a factor
ped_df$id_factor <- factor(ped_df$id)

# The (1|id_factor) term is linked to the relationship matrix A
fit_lmer <- relmatLmer(y ~ age + sex01 + (1 | id_factor),
    data = ped_df,
    relmat = list(id_factor = A)
)

# Extract variance components
vars <- as.data.frame(VarCorr(fit_lmer))
sigma2_A <- vars$vcov[1] # Additive genetic variance
sigma2_E <- vars$vcov[2] # Residual variance

# Calculate heritability
h2_est_lmer <- sigma2_A / (sigma2_A + sigma2_E)
h2_est_lmer

# TODO: get confidence intervals for h2_hat

# Profile likelihood for variance components
pr <- profile(fit_lmer, which = "theta_", prof.scale = "varcov")

# Convert the profile to “proportion of variance” scale
pr_prop <- lme4qtl::varpropProf(pr)

# Get profile-likelihood CIs; pick the row for the additive component
ci_all <- confint(pr_prop, level = 0.95)

# In a one‑random‑effect model, the additive proportion is `.sigprop01`
ci_h2 <- ci_all[setdiff(
    grep("sigprop", rownames(ci_all), value = TRUE),
    ".sigmaprop"
), , drop = FALSE]

ci_h2


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
r_MZ <- cor(tw$y1_z[tw$zygosity == "MZ"], tw$y2_z[tw$zygosity == "MZ"])

# TODO: get correlation for r_DZ
r_DZ <- cor(tw$y1_z[tw$zygosity == "DZ"], tw$y2_z[tw$zygosity == "DZ"])

# Compute h^2 = 2*(r_MZ - r_DZ)
h2_hat <- 2 * (r_MZ - r_DZ)
h2_hat
VA

# compute bootstrap CIs for h2_hat
bootstrap_h2 <- c()
for (i in 1:1000) {
    tw_boot <- tw[sample(nrow(tw), replace = TRUE), ]
    r_MZ_boot <- cor(tw_boot$y1_z[tw_boot$zygosity == "MZ"], tw_boot$y2_z[tw_boot$zygosity == "MZ"])
    r_DZ_boot <- cor(tw_boot$y1_z[tw_boot$zygosity == "DZ"], tw_boot$y2_z[tw_boot$zygosity == "DZ"])
    h2_boot <- 2 * (r_MZ_boot - r_DZ_boot)
    bootstrap_h2 <- c(bootstrap_h2, h2_boot)
}
ci_twins <- quantile(bootstrap_h2, probs = c(0.025, 0.975))
ci_twins

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
h2_liab_true <- 0.2
K_true <- 0.1
bin_dat <- simulate_binary_from_liability(A, h2_liab = h2_liab_true, K = K_true)
ped_df$case <- bin_dat$y
# TODO: Estimate h^2 from offspring-parent regression of case on midparent case
#  - Construct a data.frame with columns: fam, parent_case, child_case (all children)
#  - Fit glm(child_case ~ parent_case, family=binomial(link="probit"))
#  - Estimate h2_liab_hat

# TODO: construct `parents_bin`, a data.frame with columns: fam, midparent_case
parents_bin <- ped_df %>%
    dplyr::filter(role %in% c("father", "mother")) %>%
    dplyr::group_by(fam) %>%
    dplyr::filter(dplyr::n() == 2) %>% # require both parents
    dplyr::summarise(midparent_case = mean(case), .groups = "drop")

# TODO: construct `children_bin`, a data.frame with columns: fam, child_case (all children)
children_bin <- ped_df %>%
    dplyr::filter(role == "child") %>%
    dplyr::transmute(fam, child = case) %>%
    dplyr::ungroup()

# TODO: merge to get `df_child_bin` with columns: fam, midparent_case, child_case
df_child_bin <- children_bin %>%
    dplyr::left_join(parents_bin, by = "fam") %>%
    dplyr::filter(!is.na(midparent_case))

# TODO: K <- prevalence (mean) of child_case
K <- mean(df_child_bin$child)

# TODO: get T, which is the threshold on the liability scale
T <- qnorm(1 - K)

# TODO: get phiT = density of N(0,1) at T
phiT <- dnorm(T)

# TODO: Build midparent liability score s_mp = phi(T)*(m - K) / [K(1-K)]
# assign to df_child_bin$mp_score
df_child_bin$mp_score <- phiT * (df_child_bin$midparent_case - K) / (K * (1 - K))

# TODO: fit a probit GLM of child_case on mp_score
fit_mp <- glm(child ~ mp_score, family = binomial(link = "probit"), data = df_child_bin)

# Retrieve h2_hat
h2_hat <- coef(fit_mp)[["mp_score"]] # this is on the liability scale

h2_hat

# TODO: get standard error for h2_hat and compute 95% CI using sandwich::vcovCL (cluster by familiy)
vcl <- vcovCL(fit_mp, cluster = df_child_bin$fam)
se_b <- sqrt(diag(vcl))[["mp_score"]]
se_h2 <- as.numeric(se_b * K * (1 - K) / phiT) # delta method
ci_h2_bin <- h2_hat + c(-1.96, 1.96) * se_h2
ci_h2_bin
