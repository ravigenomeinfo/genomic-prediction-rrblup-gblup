make_small <- function(seed = 1) {
  set.seed(seed)
  g <- simulate_genotypes(n_lines = 80, n_snps = 400, n_chr = 4, n_pops = 2, fst = 0.1, missing_rate = 0.02, ld_rho = 0.95)
  p <- simulate_phenotype(g$X_true, n_qtl = 20, h2 = 0.6, mean_y = 8)
  list(g = g, p = p)
}

test_that("simulation is reproducible for a fixed seed", {
  a <- make_small(7); b <- make_small(7)
  expect_identical(a$g$X, b$g$X)
  expect_identical(a$p$y, b$p$y)
})

test_that("QC removes rare and missing-heavy SNPs and leaves no NA", {
  X <- make_small()$g$X
  X[, 1] <- 1                         # monomorphic
  X[1:30, 2] <- NA                    # >10% missing
  qc <- qc_genotypes(X, maf_min = 0.05, snp_miss_max = 0.10, ind_miss_max = 0.20)
  expect_false(any(c("SNP0001", "SNP0002") %in% colnames(qc$X)))
  expect_false(anyNA(qc$X))
  expect_true(all(qc$stats$maf[qc$stats$kept] >= 0.05))
})

test_that("VanRaden G matches rrBLUP::A.mat", {
  X <- qc_genotypes(make_small()$g$X, 0.05, 0.1, 0.2)$X
  G <- vanraden_G(X)
  expect_equal(unname(G), unname(A.mat(X)), tolerance = 1e-6)
  expect_equal(G, t(G))
})

test_that("rrBLUP and GBLUP give near-identical predictions with the same G", {
  s <- make_small(3)
  X <- qc_genotypes(s$g$X, 0.05, 0.1, 0.2)$X
  G <- vanraden_G(X); y <- s$p$y
  test <- 1:16; train <- 17:80
  expect_gt(cor(fit_predict_rrblup(X, y, train, test), fit_predict_gblup(G, y, train, test)), 0.99)
})

test_that("genomic prediction beats chance on simulated data", {
  s <- make_small(5)
  X <- qc_genotypes(s$g$X, 0.05, 0.1, 0.2)$X
  G <- vanraden_G(X)
  set.seed(1)
  r <- replicate(5, { f <- make_folds(80, 5); t <- which(f == 1)
    cor(s$p$tbv[t], fit_predict_gblup(G, s$p$y, which(f != 1), t)) })
  expect_gt(mean(r), 0.2)
})

test_that("folds are balanced and metrics are sane", {
  f <- make_folds(103, 5)
  expect_equal(length(f), 103); expect_lte(diff(range(table(f))), 1)
  m <- score_prediction(1:10, 1:10)
  expect_equal(m$r, 1); expect_equal(m$rmse, 0); expect_equal(m$slope, 1)
})
