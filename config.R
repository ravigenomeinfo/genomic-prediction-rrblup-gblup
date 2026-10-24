# Central configuration for the whole workflow.
# Every script sources this file, so all settings live in one place.

config <- list(
  # ---- reproducibility ----
  seed = 2026,

  # ---- simulated maize-like panel ----
  n_lines      = 400,    # inbred lines
  n_snps       = 2000,   # SNP markers
  n_chr        = 10,     # maize has 10 chromosomes
  n_pops       = 3,      # simulated subpopulations (population structure)
  fst          = 0.10,   # divergence among subpopulations (Balding-Nichols model)
  ld_rho       = 0.97,   # latent neighbour correlation controlling LD along chromosomes
  missing_rate = 0.02,   # fraction of genotype calls set to NA
  n_qtl        = 50,     # causal loci (drawn from the marker set)
  h2           = 0.50,   # true narrow-sense heritability (entry-mean basis)
  trait_mean   = 8,      # simulated mean grain yield (t/ha); purely illustrative

  # ---- genotype quality control ----
  maf_min      = 0.05,   # drop SNPs with minor allele frequency below this
  snp_miss_max = 0.10,   # drop SNPs missing in more than 10% of lines
  ind_miss_max = 0.20,   # drop lines missing more than 20% of SNP calls

  # ---- validation ----
  k_folds      = 5,
  n_reps       = 10,
  train_sizes  = c(50, 100, 150, 200, 240)
)
