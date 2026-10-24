# Step 5: validation under three designs.
#   A) repeated random k-fold CV
#   B) accuracy vs training-set size (random training subsets)
#   C) leave-one-cluster-out CV (predict a genetically distinct group)
# r_tbv (correlation with true breeding values) exists only because data are simulated.
check_repo_root()
source("config.R"); source("R/functions.R")

X <- readRDS("data/processed/X_qc.rds")
G <- readRDS("data/processed/G.rds")
y <- readRDS("data/processed/pheno.rds")$grain_yield_sim
cl <- readRDS("data/processed/pca.rds")$cluster
truth <- read.csv("data/simulated/truth_lines.csv")
tbv <- truth$true_breeding_value[match(rownames(X), truth$line_id)]
n <- nrow(X)

all_metrics <- list(); oof_rep1 <- list()

# ---- A) repeated random k-fold ----
set.seed(config$seed + 1)
for (rep in seq_len(config$n_reps)) {
  folds <- make_folds(n, config$k_folds)
  for (f in seq_len(config$k_folds)) {
    test <- which(folds == f); train <- which(folds != f)
    out <- run_split(X, G, y, tbv, train, test, new_meta("random_kfold", rep, f))
    all_metrics[[length(all_metrics) + 1]] <- out$metrics
    if (rep == 1) oof_rep1[[length(oof_rep1) + 1]] <- cbind(out$preds, fold = f)
  }
}
oof <- do.call(rbind, oof_rep1)
oof$line_id <- rownames(X)[oof$id]
write.csv(oof, "results/tables/oof_predictions_rep1.csv", row.names = FALSE)

# ---- B) training-set size ----
set.seed(config$seed + 2)
for (s in config$train_sizes) for (rep in seq_len(config$n_reps)) {
  train <- sample(n, s); test <- setdiff(seq_len(n), train)
  out <- run_split(X, G, y, tbv, train, test, new_meta("train_size", rep, NA, NA, s))
  all_metrics[[length(all_metrics) + 1]] <- out$metrics
}

# ---- C) leave-one-cluster-out ----
for (k in sort(unique(cl))) {
  test <- which(cl == k); train <- which(cl != k)
  out <- run_split(X, G, y, tbv, train, test, new_meta("leave_cluster_out", NA, NA, k))
  all_metrics[[length(all_metrics) + 1]] <- out$metrics
}

metrics <- do.call(rbind, all_metrics)
write.csv(metrics, "results/tables/cv_metrics_all.csv", row.names = FALSE)

# ---- summary table ----
agg <- aggregate(cbind(r, r_tbv, rmse, slope) ~ scenario + model + train_size,
                 data = transform(metrics, train_size = ifelse(is.na(train_size), -1, train_size)),
                 FUN = function(v) mean(v, na.rm = TRUE))
sds <- aggregate(r ~ scenario + model + train_size,
                 data = transform(metrics, train_size = ifelse(is.na(train_size), -1, train_size)),
                 FUN = function(v) sd(v, na.rm = TRUE))
agg$r_sd <- sds$r[match(paste(agg$scenario, agg$model, agg$train_size),
                        paste(sds$scenario, sds$model, sds$train_size))]
# Slope is unstable when prediction variance is tiny, so summarise it by the median.
slp <- aggregate(slope ~ scenario + model + train_size,
                 data = transform(metrics, train_size = ifelse(is.na(train_size), -1, train_size)),
                 FUN = function(v) median(v, na.rm = TRUE))
agg$slope <- slp$slope[match(paste(agg$scenario, agg$model, agg$train_size),
                             paste(slp$scenario, slp$model, slp$train_size))]
agg$train_size[agg$train_size == -1] <- NA
agg[, c("r", "r_tbv", "rmse", "slope", "r_sd")] <- round(agg[, c("r", "r_tbv", "rmse", "slope", "r_sd")], 3)
write.csv(agg, "results/tables/cv_summary.csv", row.names = FALSE)
print(agg)
