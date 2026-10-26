# Step 3: VanRaden genomic relationship matrix, PCA and inferred clusters.
if (!file.exists("config.R")) stop("Run scripts from the repository root (the folder containing config.R).", call. = FALSE)
source("config.R"); source("R/functions.R")
set.seed(config$seed)

X <- readRDS("data/processed/X_qc.rds")
G <- vanraden_G(X)
saveRDS(G, "data/processed/G.rds")

pca <- prcomp(X, center = TRUE, scale. = FALSE)
pcs <- pca$x[, 1:5]
var_expl <- pca$sdev^2 / sum(pca$sdev^2)

# Clusters are inferred from genotypes only (as in a real study).
cl <- kmeans(pcs[, 1:3], centers = config$n_pops, nstart = 25)$cluster
saveRDS(list(pcs = pcs, var_expl = var_expl[1:5], cluster = cl), "data/processed/pca.rds")

truth <- read.csv("data/simulated/truth_lines.csv")
print(table(inferred_cluster = cl, simulated_pop = truth$true_pop[match(rownames(X), truth$line_id)]))
cat("Mean diagonal of G (expected ~2 for fully inbred lines):", round(mean(diag(G)), 3), "\n")
cat("Variance explained by PC1-PC2:", round(sum(var_expl[1:2]) * 100, 1), "%\n")
