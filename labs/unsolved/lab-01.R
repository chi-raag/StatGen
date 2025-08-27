# Lab 01
# PUBH 8878 — Statistical Genetics (Foundations)
# ------------------------------------------------------------------
# Goal: Hands-on practice with HWE testing and allele-frequency CIs.
# ------------------------------------------------------------------

library(HardyWeinberg)
library(ggplot2)
library(dplyr)
library(fastR2)
# --- lets practice some simulations ---

# TODO: simulate a normal distribution with mean 0 and sd 1, n=1000

# TODO: plot histogram of simulated values

allele_freq_hat <- function(n_AA, n_Aa, n_aa) {
    # TODO: returns \hat p = freq(A)
}

hwe_expected_counts <- function(n_AA, n_Aa, n_aa) {
    # TODO: return expected counts under HWE
}

chisq_hwe <- function(n_AA, n_Aa, n_aa) {
    # TODO: return list with chisq statistic, df, p.value, expected counts
}

allele_freq_wald_se <- function(p, n) {
    # TODO: return SE of allele frequency estimate
}

# SNP1: common variant
snp1 <- list(id = "SNP1", n_AA = 175, n_Aa = 33, n_aa = 4)

# SNP2: rarer variant; expect small expected count in one cell
snp2 <- list(id = "SNP2", n_AA = 198, n_Aa = 12, n_aa = 0)

# TODO: TEST HWE for each SNP

# TODO: 95% CI for allele frequency for each SNP

# TODO: Parametric bootstrap under HWE null ------------------------
