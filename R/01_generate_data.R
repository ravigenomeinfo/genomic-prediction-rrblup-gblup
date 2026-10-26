# Step 1: simulate a maize-like SNP panel and a quantitative trait.
# ALL DATA ARE SIMULATED. Ground truth is saved separately in truth.csv.
if (!file.exists("config.R")) stop("Run scripts from the repository root (the folder containing config.R).", call. = FALSE)
source("config.R"); source("R/functions.R")
set.seed(config$seed)

geno <- simulate_genotypes(config$n_lines, config$n_snps, config$n_chr,
                           config$n_pops, config$fst, config$missing_rate, config$ld_rho)
ph <- simulate_phenotype(geno$X_true, config$n_qtl, config$h2, config$trait_mean)

write.csv(geno$X, "data/simulated/genotypes_raw.csv", quote = FALSE)
write.csv(data.frame(line_id = rownames(geno$X), grain_yield_sim = round(ph$y, 4)),
          "data/simulated/phenotypes.csv", row.names = FALSE, quote = FALSE)
write.csv(geno$map, "data/simulated/marker_map.csv", row.names = FALSE, quote = FALSE)
write.csv(data.frame(line_id = rownames(geno$X), true_pop = geno$pop,
                     true_breeding_value = round(ph$tbv, 4)),
          "data/simulated/truth_lines.csv", row.names = FALSE, quote = FALSE)
write.csv(data.frame(snp = colnames(geno$X)[ph$qtl], effect = ph$qtl_effect),
          "data/simulated/truth_qtl.csv", row.names = FALSE, quote = FALSE)

message("Simulated ", nrow(geno$X), " lines x ", ncol(geno$X), " SNPs; ",
        "mean missing rate = ", round(mean(is.na(geno$X)), 4))
