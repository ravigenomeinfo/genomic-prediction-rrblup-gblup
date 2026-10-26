# Step 4: fit rrBLUP and GBLUP on ALL lines (descriptive fit, not validation).
if (!file.exists("config.R")) stop("Run scripts from the repository root (the folder containing config.R).", call. = FALSE)
source("config.R"); source("R/functions.R")

X <- readRDS("data/processed/X_qc.rds")
G <- readRDS("data/processed/G.rds")
y <- readRDS("data/processed/pheno.rds")$grain_yield_sim

rr <- mixed.solve(y, Z = X)
gb <- mixed.solve(y, K = G)

h2_hat <- gb$Vu / (gb$Vu + gb$Ve)
vc <- data.frame(model = c("GBLUP", "rrBLUP"),
                 Vu = c(gb$Vu, rr$Vu), Ve = c(gb$Ve, rr$Ve),
                 h2_genomic_estimate = c(h2_hat, NA))
write.csv(vc, "results/tables/variance_components.csv", row.names = FALSE)

gebv <- data.frame(line_id = rownames(X), phenotype = y,
                   gebv_gblup = gb$u, gebv_rrblup = as.numeric(X %*% rr$u))
gebv$rank_gebv <- rank(-gebv$gebv_gblup)
write.csv(gebv[order(gebv$rank_gebv), ], "results/tables/gebv_all_lines.csv", row.names = FALSE)

me <- data.frame(snp = colnames(X), effect = rr$u)
write.csv(me, "results/tables/rrblup_marker_effects.csv", row.names = FALSE)

cat("GBLUP genomic h2 estimate:", round(h2_hat, 3), "(simulated truth:", config$h2, ")\n")
cat("Correlation of GEBV (GBLUP vs rrBLUP):", round(cor(gebv$gebv_gblup, gebv$gebv_rrblup), 4), "\n")
cat("NOTE: these GEBVs are in-sample; use the cross-validation results for accuracy.\n")
