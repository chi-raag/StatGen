# Lab 03 — SOLVED
# PUBH 8878 — Statistical Genetics
# ---------------------------------------------------------------
# Solutions aligned to Lecture 03
# ---------------------------------------------------------------

set.seed(8878)

library(dplyr)
library(ggplot2)
library(tibble)


# ---------------------------------------------------------------
# 1) EM for ABO gene frequencies
# ---------------------------------------------------------------

counts <- list(nA = 725, nAB = 72, nB = 258, nO = 1073)
N <- sum(unlist(counts))

loglik_pheno <- function(pA, pB) {
  pO <- 1 - pA - pB
  if (pA <= 0 || pB <= 0 || pO <= 0) {
    return(-Inf)
  }
  PA <- pA^2 + 2 * pA * pO
  PAB <- 2 * pA * pB
  PB <- pB^2 + 2 * pB * pO
  PO <- pO^2
  if (min(PA, PAB, PB, PO) <= 0) {
    return(-Inf)
  }
  with(counts, nA * log(PA) + nAB * log(PAB) + nB * log(PB) + nO * log(PO))
}

em_abo <- function(pA0, pB0, max_iter = 200, tol = 1e-6) {
  pA <- pA0
  pB <- pB0
  pO <- 1 - pA - pB
  out <- tibble(iter = 0, pA = pA, pB = pB, pO = pO, logLik = loglik_pheno(pA, pB))
  for (t in 1:max_iter) {
    if (pA <= 0 || pB <= 0 || (1 - pA - pB) <= 0) break
    pO <- 1 - pA - pB
    PA <- pA^2 + 2 * pA * pO
    PBp <- pB^2 + 2 * pB * pO
    # E-step weights
    wA_AA <- (pA^2) / PA
    wA_AO <- (2 * pA * pO) / PA
    wB_BB <- (pB^2) / PBp
    wB_BO <- (2 * pB * pO) / PBp
    nAA <- counts$nA * wA_AA
    nAO <- counts$nA * wA_AO
    nBB <- counts$nB * wB_BB
    nBO <- counts$nB * wB_BO
    # M-step
    pA_new <- (2 * nAA + nAO + counts$nAB) / (2 * N)
    pB_new <- (2 * nBB + nBO + counts$nAB) / (2 * N)
    pO_new <- 1 - pA_new - pB_new
    ll_new <- loglik_pheno(pA_new, pB_new)
    out <- add_row(out, iter = t, pA = pA_new, pB = pB_new, pO = pO_new, logLik = ll_new)
    if (max(abs(c(pA_new - pA, pB_new - pB, pO_new - pO))) < tol) break
    pA <- pA_new
    pB <- pB_new
  }
  out
}

# Starts
pO0 <- sqrt(counts$nO / N)
S <- 1 - pO0
cprod <- counts$nAB / (2 * N)
disc <- max(S^2 - 4 * cprod, 0)
root1 <- (S + sqrt(disc)) / 2
root2 <- (S - sqrt(disc)) / 2
pA_smart <- max(root1, root2)
pB_smart <- S - pA_smart

res_equal <- em_abo(1 / 3, 1 / 3)
res_unbal <- em_abo(0.01, 0.98)
res_smart <- em_abo(pA_smart, pB_smart)

all_res <- bind_rows(
  mutate(res_equal, strategy = "Equal"),
  mutate(res_unbal, strategy = "Unbalanced"),
  mutate(res_smart, strategy = "Smart")
) %>%
  group_by(strategy) %>%
  mutate(max_ll = max(logLik), gap = max_ll - logLik)

conv_tbl <- all_res %>%
  group_by(strategy) %>%
  summarise(
    Iterations = max(iter),
    pA = last(pA), pB = last(pB), pO = last(pO),
    logLik = last(logLik), .groups = "drop"
  )
print(conv_tbl)



# ---------------------------------------------------------------
# 2) Two-point linkage: LOD + EM where EM is necessary
#    (Transmitting parent Aa/Bb unknown phase; mate Aa/bb)
# ---------------------------------------------------------------

counts2 <- c(
  n_AA_Bb = 12,
  n_AA_bb = 3,
  n_aa_Bb = 3,
  n_aa_bb = 12,
  n_Aa_Bb = 15,
  n_Aa_bb = 15
)

probs_partial <- function(theta, w) {
  stopifnot(theta > 0, theta < 0.5, w >= 0, w <= 1)
  # Hap probs under coupling
  AB_c <- (1 - theta) / 2
  Ab_c <- theta / 2
  aB_c <- theta / 2
  ab_c <- (1 - theta) / 2
  # Hap probs under repulsion
  AB_r <- theta / 2
  Ab_r <- (1 - theta) / 2
  aB_r <- (1 - theta) / 2
  ab_r <- theta / 2
  # Mixture over phase
  AB <- w * AB_c + (1 - w) * AB_r
  Ab <- w * Ab_c + (1 - w) * Ab_r
  aB <- w * aB_c + (1 - w) * aB_r
  ab <- w * ab_c + (1 - w) * ab_r
  # Category probabilities (mate transmits A/a with prob 0.5)
  out <- c(
    n_AA_Bb = 0.5 * AB,
    n_AA_bb = 0.5 * Ab,
    n_aa_Bb = 0.5 * aB,
    n_aa_bb = 0.5 * ab,
    n_Aa_Bb = 0.5 * (AB + aB),
    n_Aa_bb = 0.5 * (Ab + ab)
  )
  stopifnot(abs(sum(out) - 1) < 1e-12)
  out
}

loglik_partial <- function(theta, w, counts2) {
  p <- probs_partial(theta, w)
  sum(counts2 * log(p))
}

em_twopoint_partial <- function(counts2, theta0 = 0.4, w0 = 0.5, max_iter = 200, tol = 1e-8) {
  th <- theta0
  w <- w0
  N <- sum(counts2)
  traj <- tibble(iter = 0, theta = th, w = w, logLik = loglik_partial(th, w, counts2))
  for (t in 1:max_iter) {
    # E-step: split ambiguous categories
    qAB <- w * (1 - th) + (1 - w) * th # P(AB | Aa,Bb; th, w)
    qAb <- w * th + (1 - w) * (1 - th) # P(Ab | Aa,bb; th, w)
    # Expected hap counts
    nAB <- counts2["n_AA_Bb"] + qAB * counts2["n_Aa_Bb"]
    naB <- counts2["n_aa_Bb"] + (1 - qAB) * counts2["n_Aa_Bb"]
    nAb <- counts2["n_AA_bb"] + qAb * counts2["n_Aa_bb"]
    nab <- counts2["n_aa_bb"] + (1 - qAb) * counts2["n_Aa_bb"]
    # Recombinant vs nonrecombinant under each phase
    R_c <- nAb + naB
    NR_c <- nAB + nab
    R_r <- nAB + nab
    NR_r <- nAb + naB
    # Update w using stable log-likelihood ratio
    logLc <- NR_c * log(1 - th) + R_c * log(th)
    logLr <- NR_r * log(1 - th) + R_r * log(th)
    m <- max(logLc, logLr)
    w_new <- exp(logLc - m) / (exp(logLc - m) + exp(logLr - m))
    # Update theta
    th_new <- (w_new * R_c + (1 - w_new) * R_r) / N
    ll_new <- loglik_partial(th_new, w_new, counts2)
    traj <- add_row(traj, iter = t, theta = th_new, w = w_new, logLik = ll_new)
    if (max(abs(c(th_new - th, w_new - w))) < tol) break
    th <- th_new
    w <- w_new
  }
  traj
}

fit2 <- em_twopoint_partial(counts2, theta0 = 0.1, w0 = 0.5)
print(tail(fit2, 3))

ggplot(fit2, aes(iter, theta)) +
  geom_line() +
  geom_point() +
  labs(x = "Iteration", y = expression(hat(theta))) +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank())

theta_hat <- tail(fit2$theta, 1)
w_hat <- tail(fit2$w, 1)
ll_hat <- loglik_partial(theta_hat, w_hat, counts2)
ll_null <- loglik_partial(0.5 - 1e-8, 0.5, counts2) # theta=0.5 (phase irrelevant)
LOD2 <- (ll_hat - ll_null) / log(10)
cat(sprintf(
  "Two-point partial-info EM: theta_hat=%.4f, w_hat=%.3f, LOD=%.3f\n",
  theta_hat, w_hat, LOD2
))

# ---------------------------------------------------------------
# 3) Two-SNP haplotype EM + LD
# ---------------------------------------------------------------

p_true <- c(ab = .40, aB = .10, Ab = .25, AB = .25)
n <- 1000

draw_hap <- function(p) sample(names(p), 2, replace = TRUE, prob = unname(p))
hap_to_vec <- function(h) {
  map1 <- list(ab = c(0, 0), aB = c(0, 1), Ab = c(1, 0), AB = c(1, 1))
  v1 <- map1[[h[1]]]
  v2 <- map1[[h[2]]]
  c(g1 = v1[1] + v2[1], g2 = v1[2] + v2[2])
}

haps <- replicate(n, draw_hap(p_true))
G <- t(apply(haps, 2, hap_to_vec))
G <- as_tibble(G)

em_hap2 <- function(G, p0 = rep(0.25, 4), tol = 1e-8, maxit = 500) {
  names(p0) <- c("ab", "aB", "Ab", "AB")
  p <- p0 / sum(p0)
  n <- nrow(G)
  traj <- tibble(iter = 0, p_ab = p["ab"], p_aB = p["aB"], p_Ab = p["Ab"], p_AB = p["AB"])
  for (t in 1:maxit) {
    exp_counts <- c(ab = 0, aB = 0, Ab = 0, AB = 0)
    for (i in 1:n) {
      g1 <- G$g1[i]
      g2 <- G$g2[i]
      if (g1 == 0 && g2 == 0) {
        exp_counts["ab"] <- exp_counts["ab"] + 2
      } else if (g1 == 2 && g2 == 2) {
        exp_counts["AB"] <- exp_counts["AB"] + 2
      } else if (g1 == 0 && g2 == 2) {
        exp_counts["aB"] <- exp_counts["aB"] + 2
      } else if (g1 == 2 && g2 == 0) {
        exp_counts["Ab"] <- exp_counts["Ab"] + 2
      } else if (g1 == 1 && g2 == 0) {
        exp_counts["Ab"] <- exp_counts["Ab"] + 1
        exp_counts["ab"] <- exp_counts["ab"] + 1
      } else if (g1 == 1 && g2 == 2) {
        exp_counts["AB"] <- exp_counts["AB"] + 1
        exp_counts["aB"] <- exp_counts["aB"] + 1
      } else if (g1 == 2 && g2 == 1) {
        exp_counts["AB"] <- exp_counts["AB"] + 1
        exp_counts["Ab"] <- exp_counts["Ab"] + 1
      } else if (g1 == 0 && g2 == 1) {
        exp_counts["aB"] <- exp_counts["aB"] + 1
        exp_counts["ab"] <- exp_counts["ab"] + 1
      } else if (g1 == 1 && g2 == 1) {
        # ambiguous double het: {AB+ab} vs {Ab+aB}
        w1 <- p["AB"] * p["ab"]
        w2 <- p["Ab"] * p["aB"]
        s <- w1 + w2
        if (s == 0) {
          w1 <- w2 <- 0.5
        } else {
          w1 <- w1 / s
          w2 <- 1 - w1
        }
        exp_counts["AB"] <- exp_counts["AB"] + w1
        exp_counts["ab"] <- exp_counts["ab"] + w1
        exp_counts["Ab"] <- exp_counts["Ab"] + w2
        exp_counts["aB"] <- exp_counts["aB"] + w2
      } else {
        stop("Unexpected genotype combination.")
      }
    }
    p_new <- exp_counts / (2 * n)
    if (max(abs(p_new - p)) < tol) {
      p <- p_new
      traj <- add_row(traj, iter = t, p_ab = p["ab"], p_aB = p["aB"], p_Ab = p["Ab"], p_AB = p["AB"])
      break
    }
    p <- p_new
    traj <- add_row(traj, iter = t, p_ab = p["ab"], p_aB = p["aB"], p_Ab = p["Ab"], p_AB = p["AB"])
  }
  list(p = p, traj = traj)
}

fit_h <- em_hap2(G)
p_hat <- fit_h$p
print(rbind(true = p_true, est = p_hat))

# LD
pA <- p_hat["Ab"] + p_hat["AB"]
pB <- p_hat["aB"] + p_hat["AB"]
p11 <- p_hat["AB"]
D <- as.numeric(p11 - pA * pB)
r2 <- as.numeric(D^2 / (pA * (1 - pA) * pB * (1 - pB)))
cat(sprintf("LD: D=%.4f, r^2=%.4f\n", D, r2))

# ---------------------------------------------------------------
# 4) Mini association demo + QC
# ---------------------------------------------------------------

set.seed(1203)
n <- 2000
maf <- 0.30
geno <- rbinom(n, 2, maf) # 0/1/2 copies
age <- rnorm(n, 50, 10)
linpred <- -3.0 + 0.40 * scale(age)[, 1] + 0.35 * geno # log-odds
pr <- plogis(linpred)
case <- rbinom(n, 1, pr)
dat <- tibble(geno, age, case)

hwe_chisq <- function(nAA, nAa, naa) {
  n <- nAA + nAa + naa
  p <- (2 * nAA + nAa) / (2 * n)
  expAA <- n * p^2
  expAa <- 2 * n * p * (1 - p)
  expaa <- n * (1 - p)^2
  stat <- (nAA - expAA)^2 / expAA + (nAa - expAa)^2 / expAa + (naa - expaa)^2 / expaa
  pval <- pchisq(stat, df = 1, lower.tail = FALSE)
  list(stat = stat, pval = pval)
}

ctrl <- dat %>% filter(case == 0)
nAA <- sum(ctrl$geno == 2)
nAa <- sum(ctrl$geno == 1)
naa <- sum(ctrl$geno == 0)
hwe <- hwe_chisq(nAA, nAa, naa)
cat(sprintf("Controls HWE: X2=%.3f, p=%.3g\n", hwe$stat, hwe$pval))

fit <- glm(case ~ geno + scale(age), family = binomial(), data = dat)
or <- exp(coef(fit)["geno"])
se <- sqrt(vcov(fit)["geno", "geno"])
ci <- exp(coef(fit)["geno"] + c(-1, 1) * 1.96 * se)
cat(sprintf("Per-allele OR: %.3f (95%% CI %.3f–%.3f)\n", or, ci[1], ci[2]))

if (hwe$pval < 1e-6) {
  message("QC FLAG: SNP fails HWE in controls (p < 1e-6).")
} else {
  message("QC PASS: SNP HWE in controls is acceptable.")
}

# End of Lab 03 (solved)
