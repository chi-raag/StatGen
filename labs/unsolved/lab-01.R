# Lab 01 — HWE: χ² vs Exact, and Small-Sample Pitfalls
# PUBH 8878 — Statistical Genetics (Foundations)
# ------------------------------------------------------------------
# Goal: Hands-on practice with HWE testing and allele-frequency CIs.
# ------------------------------------------------------------------


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

snps <- list(snp1, snp2)

# TODO: TEST HWE for each SNP
