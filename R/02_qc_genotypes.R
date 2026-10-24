# Step 2: genotype QC (line/SNP missingness, MAF) and mean imputation.
# QC uses genotypes only (no phenotypes), so it cannot leak phenotype information.
check_repo_root()
source("config.R"); source("R/functions.R")

X_raw <- as.matrix(read.csv("data/simulated/genotypes_raw.csv", row.names = 1, check.names = FALSE))
pheno <- read.csv("data/simulated/phenotypes.csv")
stopifnot(identical(rownames(X_raw), pheno$line_id))

qc <- qc_genotypes(X_raw, config$maf_min, config$snp_miss_max, config$ind_miss_max)
pheno <- pheno[match(rownames(qc$X), pheno$line_id), ]
stopifnot(identical(rownames(qc$X), pheno$line_id), !anyNA(qc$X))

dir.create("data/processed", showWarnings = FALSE)
saveRDS(qc$X, "data/processed/X_qc.rds")
saveRDS(pheno, "data/processed/pheno.rds")
write.csv(qc$stats, "results/tables/qc_marker_stats.csv", row.names = FALSE)
write.csv(qc$summary, "results/tables/qc_summary.csv", row.names = FALSE)
print(qc$summary)
