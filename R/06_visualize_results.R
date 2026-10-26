# Step 6: publication-quality figures (300 dpi PNG).
if (!file.exists("config.R")) stop("Run scripts from the repository root (the folder containing config.R).", call. = FALSE)
source("config.R"); source("R/functions.R")

stats   <- read.csv("results/tables/qc_marker_stats.csv")
G       <- readRDS("data/processed/G.rds")
pca     <- readRDS("data/processed/pca.rds")
metrics <- read.csv("results/tables/cv_metrics_all.csv")
oof     <- read.csv("results/tables/oof_predictions_rep1.csv")
me      <- read.csv("results/tables/rrblup_marker_effects.csv")
map     <- read.csv("data/simulated/marker_map.csv")
qtl     <- read.csv("data/simulated/truth_qtl.csv")
truth   <- read.csv("data/simulated/truth_lines.csv")
model_lab <- c(gblup = "GBLUP", rrblup = "rrBLUP")
metrics$model <- factor(model_lab[metrics$model], levels = model_lab)

# Fig 1: MAF distribution with QC threshold
f1 <- ggplot(stats, aes(maf, fill = kept)) +
  geom_histogram(bins = 40, colour = "white") +
  geom_vline(xintercept = config$maf_min, linetype = 2) +
  scale_fill_manual(values = c(`TRUE` = okabe_ito[1], `FALSE` = "grey60"),
                    labels = c(`TRUE` = "retained", `FALSE` = "removed"), name = NULL) +
  labs(x = "Minor allele frequency", y = "Number of SNPs",
       title = "SNP minor allele frequency (simulated data)") + theme_paper()
save_fig(f1, "fig01_maf_distribution.png")

# Fig 2: PCA
d <- data.frame(pca$pcs[, 1:2], line_id = rownames(pca$pcs))
d$pop <- factor(truth$true_pop[match(d$line_id, truth$line_id)])
f2 <- ggplot(d, aes(PC1, PC2, colour = pop)) + geom_point(size = 1.8, alpha = 0.8) +
  scale_colour_manual(values = okabe_ito, name = "Simulated\nsubpopulation") +
  labs(x = sprintf("PC1 (%.1f%%)", 100 * pca$var_expl[1]),
       y = sprintf("PC2 (%.1f%%)", 100 * pca$var_expl[2]),
       title = "Population structure") + theme_paper()
save_fig(f2, "fig02_pca.png")

# Fig 3: genomic relationship heatmap (lines ordered by hierarchical clustering)
ord <- hclust(dist(G))$order
gl <- expand.grid(i = seq_along(ord), j = seq_along(ord))
gl$value <- G[cbind(ord[gl$i], ord[gl$j])]
f3 <- ggplot(gl, aes(i, j, fill = value)) + geom_raster() +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-1, 2), oob = scales::squish, name = "G") +
  coord_fixed() + labs(x = "Lines", y = "Lines", title = "Genomic relationship matrix (VanRaden)") +
  theme_paper() + theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank())
save_fig(f3, "fig03_grm_heatmap.png", 5.5, 4.8)

# Fig 4: out-of-fold predicted vs observed (repetition 1)
long <- rbind(data.frame(oof[, c("line_id", "obs")], pred = oof$gblup, model = "GBLUP"),
              data.frame(oof[, c("line_id", "obs")], pred = oof$rrblup, model = "rrBLUP"))
ann <- do.call(rbind, lapply(split(long, long$model), function(x)
  data.frame(model = x$model[1], label = sprintf("r = %.2f", cor(x$obs, x$pred)))))
f4 <- ggplot(long, aes(pred, obs)) + geom_point(alpha = 0.6, colour = okabe_ito[1]) +
  geom_smooth(method = "lm", se = FALSE, colour = okabe_ito[2], linewidth = 0.8, formula = y ~ x) +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = label), hjust = -0.2, vjust = 1.5, inherit.aes = FALSE) +
  facet_wrap(~model) +
  labs(x = "Predicted value (out-of-fold)", y = "Observed phenotype (simulated t/ha)",
       title = "Five-fold cross-validation, repetition 1") + theme_paper()
save_fig(f4, "fig04_predicted_vs_observed.png", 7, 3.8)

# Fig 5: accuracy distribution, random k-fold
a <- subset(metrics, scenario == "random_kfold")
f5 <- ggplot(a, aes(model, r, fill = model)) + geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.12, size = 0.8, alpha = 0.5) +
  scale_fill_manual(values = okabe_ito[1:2], guide = "none") +
  labs(x = NULL, y = "Predictive ability (r)",
       title = sprintf("Random %d-fold CV, %d repetitions", config$k_folds, config$n_reps)) + theme_paper()
save_fig(f5, "fig05_cv_accuracy.png", 4, 4)

# Fig 6: training-set size
ts <- aggregate(r ~ train_size + model, data = subset(metrics, scenario == "train_size"),
                FUN = function(v) c(m = mean(v), s = sd(v)))
ts <- data.frame(ts[, 1:2], m = ts$r[, "m"], s = ts$r[, "s"])
f6 <- ggplot(ts, aes(train_size, m, colour = model)) +
  geom_line(position = position_dodge(4)) + geom_point(position = position_dodge(4)) +
  geom_errorbar(aes(ymin = m - s, ymax = m + s), width = 5, position = position_dodge(4)) +
  scale_colour_manual(values = okabe_ito[1:2], name = NULL) +
  labs(x = "Training population size", y = "Predictive ability (r, mean \u00b1 SD)",
       title = "Effect of training-set size") + theme_paper()
save_fig(f6, "fig06_training_size.png", 5.5, 4)

# Fig 7: validation design matters
b <- rbind(
  transform(aggregate(r ~ rep + model, data = subset(metrics, scenario == "random_kfold"), FUN = mean),
            design = "Random 5-fold CV")[, c("model", "r", "design")],
  transform(subset(metrics, scenario == "leave_cluster_out"),
            design = "Leave-one-cluster-out")[, c("model", "r", "design")])
b$design <- factor(b$design, levels = c("Random 5-fold CV", "Leave-one-cluster-out"))
f7 <- ggplot(b, aes(design, r, colour = model)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.5), size = 2) +
  scale_colour_manual(values = okabe_ito[1:2], name = NULL) +
  labs(x = NULL, y = "Predictive ability (r)",
       title = "Random vs leave-one-cluster-out validation",
       subtitle = "Points: CV repetitions (random) or held-out clusters") +
  coord_cartesian(ylim = c(0, 0.6)) + theme_paper()
save_fig(f7, "fig07_validation_design.png", 6, 4)

# Fig 8: rrBLUP marker effects along the genome, true QTL highlighted
m <- merge(me, map, by = "snp"); m <- m[order(m$chr, m$pos), ]
m$x <- seq_len(nrow(m)); m$is_qtl <- m$snp %in% qtl$snp
f8 <- ggplot(m, aes(x, abs(effect), colour = factor(chr %% 2))) + geom_point(size = 0.7, alpha = 0.7) +
  geom_point(data = subset(m, is_qtl), colour = okabe_ito[2], size = 1.8, shape = 17) +
  scale_colour_manual(values = c("grey40", "grey70"), guide = "none") +
  labs(x = "Marker (ordered by chromosome and position)", y = "|rrBLUP marker effect|",
       title = "Estimated marker effects", subtitle = "Orange triangles: simulated causal loci that passed QC") +
  theme_paper()
save_fig(f8, "fig08_marker_effects.png", 7, 3.8)
message("Saved 8 figures to results/figures/")
