# Shared functions. Genotypes are coded -1 / 0 / 1 (rrBLUP convention);
# inbred maize lines are homozygous, so simulated genotypes are -1 or 1.

suppressPackageStartupMessages({
  library(rrBLUP)
  library(ggplot2)
})

# ---------------------------------------------------------------- simulation
simulate_genotypes <- function(n_lines, n_snps, n_chr, n_pops, fst, missing_rate, ld_rho) {
  # Balding-Nichols model: subpopulation allele frequencies drift from an
  # ancestral frequency (population structure). Linkage disequilibrium is
  # introduced with a Gaussian copula: along each chromosome, a latent AR(1)
  # process with neighbour correlation `ld_rho` is thresholded at each SNP's
  # subpopulation allele frequency. This is a simplified LD model, not a
  # coalescent simulation (see README limitations).
  p_anc <- pmin(pmax(rbeta(n_snps, 0.7, 0.7), 0.02), 0.98)
  a <- p_anc * (1 - fst) / fst
  b <- (1 - p_anc) * (1 - fst) / fst
  p_pop <- sapply(seq_len(n_pops), function(k) rbeta(n_snps, a, b))

  pop <- sample(rep_len(seq_len(n_pops), n_lines))
  chr <- sort(rep_len(seq_len(n_chr), n_snps))

  Z <- matrix(rnorm(n_lines * n_snps), n_lines, n_snps)
  for (j in seq_len(n_snps)[-1]) {
    if (chr[j] == chr[j - 1]) Z[, j] <- ld_rho * Z[, j - 1] + sqrt(1 - ld_rho^2) * Z[, j]
  }
  thr <- qnorm(p_pop[, pop, drop = FALSE])             # n_snps x n_lines
  X <- ifelse(Z < t(thr), 1, -1)                       # P(+1) = subpop frequency

  line_ids <- sprintf("L%03d", seq_len(n_lines))
  snp_ids  <- sprintf("SNP%04d", seq_len(n_snps))
  dimnames(X) <- list(line_ids, snp_ids)

  pos <- ave(seq_len(n_snps), chr, FUN = function(i) as.integer(round(seq(1, 3e8 / n_chr, length.out = length(i)))))
  map <- data.frame(snp = snp_ids, chr = chr, pos = pos)

  X_true <- X
  X[sample(length(X), round(missing_rate * length(X)))] <- NA
  list(X = X, X_true = X_true, pop = setNames(pop, line_ids), map = map)
}

simulate_phenotype <- function(X_true, n_qtl, h2, mean_y) {
  p   <- (colMeans(X_true) + 1) / 2
  maf <- pmin(p, 1 - p)
  qtl <- sample(which(maf > 0.10), n_qtl)
  eff <- rnorm(n_qtl)
  g   <- as.numeric(X_true[, qtl, drop = FALSE] %*% eff)
  g   <- (g - mean(g)) / sd(g)               # true genetic values, variance 1
  e   <- rnorm(length(g), 0, sqrt((1 - h2) / h2))
  list(y = mean_y + g + e, tbv = g, qtl = qtl, qtl_effect = eff)
}

# ------------------------------------------------------------------------ QC
qc_genotypes <- function(X, maf_min, snp_miss_max, ind_miss_max) {
  n0 <- nrow(X); m0 <- ncol(X)
  keep_i <- rowMeans(is.na(X)) <= ind_miss_max
  X <- X[keep_i, , drop = FALSE]

  snp_miss <- colMeans(is.na(X))
  p   <- (colMeans(X, na.rm = TRUE) + 1) / 2
  maf <- pmin(p, 1 - p)
  stats <- data.frame(snp = colnames(X), maf = maf, missing = snp_miss,
                      kept = snp_miss <= snp_miss_max & maf >= maf_min)
  X <- X[, stats$kept, drop = FALSE]

  # Mean imputation (simple; see README limitations).
  n_imputed <- sum(is.na(X))
  cm <- colMeans(X, na.rm = TRUE)
  for (j in seq_len(ncol(X))) X[is.na(X[, j]), j] <- cm[j]

  summary <- data.frame(
    step = c("lines_before", "lines_after", "snps_before", "snps_after", "calls_imputed"),
    value = c(n0, nrow(X), m0, ncol(X), n_imputed))
  list(X = X, stats = stats, summary = summary)
}

# ------------------------------------------------ genomic relationship matrix
# VanRaden (2008) method 1: G = ZZ' / (2 * sum p(1-p)), Z = (W - 2p), W = X + 1
vanraden_G <- function(X) {
  p <- (colMeans(X) + 1) / 2
  Z <- sweep(X, 2, 2 * p - 1)
  G <- tcrossprod(Z) / (2 * sum(p * (1 - p)))
  dimnames(G) <- list(rownames(X), rownames(X))
  G
}

# ------------------------------------------------------------------- models
# rrBLUP: marker effects as random effects (ridge regression BLUP)
fit_predict_rrblup <- function(X, y, train, test) {
  fit <- mixed.solve(y[train], Z = X[train, , drop = FALSE])
  as.numeric(X[test, , drop = FALSE] %*% fit$u + as.numeric(fit$beta))
}

# GBLUP: genotype effects as random effects with covariance G.
# Test phenotypes are masked (NA); mixed.solve returns effects for all lines.
fit_predict_gblup <- function(G, y, train, test) {
  y_masked <- y
  y_masked[test] <- NA
  fit <- mixed.solve(y_masked, K = G)
  as.numeric(fit$u[test] + as.numeric(fit$beta))
}

# ------------------------------------------------------------------ metrics
score_prediction <- function(obs, pred, tbv = NULL) {
  data.frame(
    r      = suppressWarnings(cor(obs, pred)),            # predictive ability
    r_tbv  = if (is.null(tbv)) NA_real_ else suppressWarnings(cor(tbv, pred)),
    rmse   = sqrt(mean((obs - pred)^2)),
    slope  = if (var(pred) < 1e-8) NA_real_ else cov(obs, pred) / var(pred))  # ~1 = unbiased
}

run_split <- function(X, G, y, tbv, train, test, meta) {
  preds <- list(rrblup = fit_predict_rrblup(X, y, train, test),
                gblup  = fit_predict_gblup(G, y, train, test))
  m <- do.call(rbind, lapply(names(preds), function(nm) {
    cbind(meta, model = nm, n_train = length(train), n_test = length(test),
          score_prediction(y[test], preds[[nm]], tbv[test]))
  }))
  list(metrics = m,
       preds = data.frame(id = test, obs = y[test], rrblup = preds$rrblup, gblup = preds$gblup))
}

make_folds <- function(n, k) sample(rep_len(seq_len(k), n))

new_meta <- function(scenario, rep = NA, fold = NA, group = NA, train_size = NA) {
  data.frame(scenario = scenario, rep = rep, fold = fold, group = group, train_size = train_size)
}

# --------------------------------------------------------------------- plots
okabe_ito <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9")

theme_paper <- function() {
  theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey92"))
}

save_fig <- function(plot, name, width = 6, height = 4.5) {
  ggsave(file.path("results", "figures", name), plot, width = width, height = height, dpi = 300)
}

check_repo_root <- function() {
  if (!file.exists("config.R") || !dir.exists("R"))
    stop("Run scripts from the repository root (the folder containing config.R).", call. = FALSE)
}
