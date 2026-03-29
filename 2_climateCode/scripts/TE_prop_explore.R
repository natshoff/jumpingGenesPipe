# TE Analysis Pipeline
# Script: TE_prop_explore.R
#
# Geographic partitioning analyses for TE family profiles.
# Asks: is there significant ecoregion / geographic structure in TE composition?
#
# Analyses:
#   1. PERMANOVA (adonis2)           — does ecoregion explain TE variation?
#   2. Homogeneity of dispersion     — validates PERMANOVA assumption of equal spread
#   3. ANOSIM                        — non-parametric alternative group test
#   4. Mantel test                   — isolation-by-distance   ≥≥≥≥≥≥≥≥÷
#   5. Population-averaged heatmap   — TE family fingerprints per population
#   6. Geographic tSNE map           — k-means clusters mapped to coordinates
#   7. Variance partitioning         — climate vs geography vs overlap
#
# Input:  2_climateCode/data/MR_fam.csv
#         2_climateCode/data/MR_sup.csv
# Output: 2_climateCode/results/TE_prop_explore/  (created if absent)


###################
# Library loading #
###################

library(tidyverse)
library(vegan)      # adonis2, anosim, mantel, varpart, vegdist, decostand
library(pheatmap)   # clustered heatmap
library(Rtsne)      # tSNE
library(maps)       # US state outlines for geographic map
library(patchwork)  # combining plots


############################
# Parameters — edit here   #
############################

setwd("~/Documents/GitHub/jumpingGenesPipe")

fam_path <- "2_climateCode/data/MR_fam.csv"
sup_path <- "2_climateCode/data/MR_sup.csv"
out_dir  <- "2_climateCode/results/TE_prop_explore"

# Number of k-means clusters for geographic tSNE map
n_clusters <- 5

# Number of top-variable families to display in heatmap (keeps it readable)
n_heatmap_fams <- 60

# tSNE perplexity — try 15, 30, 50 to check stability
tsne_perplexity <- 30

# tSNE random seed for reproducibility
tsne_seed <- 142


####################
# Setup & load data #
####################

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

MR_fam <- read.csv(fam_path, header = TRUE) %>%
  mutate(Ecoregion = as.factor(Ecoregion))

MR_sup <- read.csv(sup_path, header = TRUE) %>%
  mutate(Ecoregion = as.factor(Ecoregion))

# Isolate TE family proportion matrix (individuals × families)
te_mat <- MR_fam %>%
  dplyr::select(ends_with("_mean")) %>%
  as.matrix()

# Hellinger transformation: square root of relative proportions
# Reduces the influence of very abundant families and makes Euclidean
# distance ecologically meaningful (Legendre & Gallagher 2001)
te_hel <- decostand(te_mat, method = "hellinger")

# Environmental matrices for variance partitioning
climate_mat <- MR_fam %>%
  dplyr::select(bio1:alt) %>%
  as.data.frame()

geo_mat <- MR_fam %>%
  dplyr::select(Latitude, Longitude) %>%
  as.data.frame()

message("Data loaded. Individuals: ", nrow(MR_fam),
        " | TE families: ", ncol(te_mat),
        " | Ecoregions: ", nlevels(MR_fam$Ecoregion))


#####################################################
# 1. PERMANOVA — does ecoregion predict TE profiles? #
#####################################################
# adonis2 partitions variance in the Bray-Curtis dissimilarity matrix
# by ecoregion membership.  R² = fraction of total variance explained.
# Significance assessed by 999 permutations.

message("\n========== 1. PERMANOVA ==========")
perm_result <- adonis2(
  te_mat ~ Ecoregion,
  data         = MR_fam,
  permutations = 999,
  method       = "bray"
)
print(perm_result)
write.csv(as.data.frame(perm_result),
          file.path(out_dir, "permanova_result.csv"))


########################################################
# 2. Homogeneity of dispersion — PERMANOVA assumption   #
########################################################
# If ecoregions differ in *spread* (not just centroid position) PERMANOVA
# can give false positives. betadisper tests this. A non-significant result
# here strengthens the PERMANOVA conclusion.

message("\n========== 2. Homogeneity of dispersion ==========")
bc_dist    <- vegdist(te_mat, method = "bray")
disp_result <- betadisper(bc_dist, MR_fam$Ecoregion)
disp_anova  <- anova(disp_result)
print(disp_anova)

# Visualise dispersion per ecoregion
disp_df <- data.frame(
  Ecoregion  = MR_fam$Ecoregion,
  Distance   = disp_result$distances
)
p_disp <- ggplot(disp_df, aes(x = Ecoregion, y = Distance, fill = Ecoregion)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 4) +
  theme_bw() +
  labs(title = "Multivariate dispersion per ecoregion",
       subtitle = paste0("betadisper ANOVA p = ",
                         round(disp_anova$`Pr(>F)`[1], 4)),
       y = "Distance to centroid") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 20, hjust = 1))
p_disp
ggsave(file.path(out_dir, "dispersion_boxplot.png"),
       p_disp, width = 8, height = 5, dpi = 300)


###################################
# 3. ANOSIM — non-parametric test #
###################################
# Compares within-group vs between-group dissimilarities.
# R statistic: 0 = no separation, 1 = perfect separation.
# Less sensitive to unequal group sizes than PERMANOVA.

message("\n========== 3. ANOSIM ==========")
anosim_result <- anosim(bc_dist, MR_fam$Ecoregion, permutations = 999)
print(summary(anosim_result))


###########################
# 4. Mantel test          #
###########################
# Tests whether TE profile dissimilarity (Bray-Curtis) correlates with
# geographic distance (Euclidean on lat/lon).
# A significant positive correlation = isolation by distance.
# NOTE: lat/lon Euclidean distance is an approximation; fine for relative
# comparisons across a continental-US range.

message("\n========== 4. Mantel test ==========")
geo_dist    <- dist(geo_mat)
mantel_result <- mantel(bc_dist, geo_dist,
                        method = "spearman", permutations = 999)
print(mantel_result)

# Scatter: pairwise geographic distance vs Bray-Curtis dissimilarity
mantel_df <- data.frame(
  geo  = as.vector(geo_dist),
  bc   = as.vector(bc_dist)
)
p_mantel <- ggplot(mantel_df, aes(x = geo, y = bc)) +
  geom_point(alpha = 0.05, size = 0.8) +
  geom_smooth(method = "lm", color = "#1B4332", se = TRUE) +
  theme_bw() +
  labs(title = "Mantel: geographic distance vs TE dissimilarity",
       subtitle = paste0("Mantel r = ", round(mantel_result$statistic, 3),
                         "  p = ", round(mantel_result$signif, 4)),
       x = "Geographic distance (Euclidean lat/lon)",
       y = "Bray-Curtis dissimilarity")
p_mantel
ggsave(file.path(out_dir, "mantel_scatter.png"),
       p_mantel, width = 7, height = 5, dpi = 300)


##################################################
# 5. Population-averaged heatmap                  #
##################################################
# Average TE family proportions per population, then cluster.
# Only the top N most variable families are shown for readability.
# Row annotation = ecoregion; column annotation = top contributing families.

message("\n========== 5. Heatmap ==========")

pop_avg <- MR_fam %>%
  group_by(Population, Ecoregion) %>%
  summarise(across(ends_with("_mean"), mean, .names = "{.col}"),
            .groups = "drop")

pop_mat <- pop_avg %>%
  dplyr::select(ends_with("_mean")) %>%
  as.matrix()
rownames(pop_mat) <- pop_avg$Population

# Keep only top N families by variance across populations
fam_var   <- apply(pop_mat, 2, var)
top_fams  <- names(sort(fam_var, decreasing = TRUE))[1:min(n_heatmap_fams, ncol(pop_mat))]
pop_mat_sub <- pop_mat[, top_fams]
# Clean up column names: remove "_mean" suffix
colnames(pop_mat_sub) <- str_remove(colnames(pop_mat_sub), "_mean$")

# Annotation: ecoregion per population
anno_row <- data.frame(
  Ecoregion = as.character(pop_avg$Ecoregion),
  row.names = pop_avg$Population
)

# Colour palette per ecoregion (matches plotting scripts)
eco_levels <- sort(unique(anno_row$Ecoregion))
eco_colors <- setNames(
  c("#54286f","#a989b9","#8b0000","#9ab87a","#304c00",
    "#cdb79e","#85683c","#1066bc","#8abff5","#BDDBF9")[seq_along(eco_levels)],
  eco_levels
)
anno_colors <- list(Ecoregion = eco_colors)

png(file.path(out_dir, "heatmap_top_families.png"),
    width = 14, height = 10, units = "in", res = 300)
pheatmap(
  pop_mat_sub,
  annotation_row  = anno_row,
  annotation_colors = anno_colors,
  scale           = "column",
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  fontsize_row    = 7,
  fontsize_col    = 6,
  main            = paste0("Population TE fingerprints (top ", n_heatmap_fams,
                           " most variable families, z-scored)")
)
dev.off()
message("Heatmap saved.")


##################################################
# 6. Geographic tSNE map                          #
##################################################
# Run tSNE on shared TE families, k-means cluster the embedding,
# then plot population centroids on a US map coloured by cluster.
# This reveals whether tSNE clusters have geographic coherence.

message("\n========== 6. Geographic tSNE map ==========")

# --- prepare family data: keep only families present in ALL individuals ---
te_tsne <- MR_fam %>%
  dplyr::select(ends_with("_mean")) %>%
  as.matrix()
te_tsne[te_tsne == 0] <- NA
shared_cols <- colSums(is.na(te_tsne)) == 0
te_tsne_shared <- te_tsne[, shared_cols]
te_tsne_shared[is.na(te_tsne_shared)] <- 0  # safety: shouldn't trigger
message(paste("Shared TE families for tSNE:", ncol(te_tsne_shared)))

# --- run tSNE ---
set.seed(tsne_seed)
tsne_fit <- Rtsne(
  scale(te_tsne_shared),
  perplexity  = tsne_perplexity,
  check_duplicates = FALSE
)

tsne_df <- as.data.frame(tsne_fit$Y) %>%
  rename(tSNE1 = V1, tSNE2 = V2) %>%
  bind_cols(MR_fam %>% dplyr::select(ID, Population, Ecoregion,
                                      Latitude, Longitude, State))

# --- k-means cluster tSNE embedding ---
set.seed(tsne_seed)
km <- kmeans(tsne_df[, c("tSNE1", "tSNE2")], centers = n_clusters, nstart = 25)
tsne_df$Cluster <- as.factor(km$cluster)

# tSNE coloured by cluster
p_tsne_cluster <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2,
                                       color = Cluster)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text(aes(label = Population), size = 2.5, nudge_x = 0.4,
            check_overlap = TRUE) +
  theme_bw() +
  labs(title = paste0("tSNE coloured by k-means cluster (k = ", n_clusters, ")"),
       subtitle = paste0("Perplexity = ", tsne_perplexity))

# tSNE coloured by ecoregion (for comparison)
p_tsne_eco <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2,
                                   color = Ecoregion)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text(aes(label = Population), size = 2.5, nudge_x = 0.4,
            check_overlap = TRUE) +
  scale_color_manual(values = c(
    "#54286f","#a989b9","#8b0000","#9ab87a","#304c00",
    "#cdb79e","#85683c","#1066bc","#8abff5","#BDDBF9")) +
  theme_bw() +
  labs(title = "tSNE coloured by ecoregion")

p_tsne_cluster + p_tsne_eco

ggsave(file.path(out_dir, "tsne_cluster_vs_ecoregion.png"),
       p_tsne_cluster + p_tsne_eco,
       width = 16, height = 6, dpi = 300)

# --- map: population centroids coloured by tSNE cluster ---
pop_geo <- tsne_df %>%
  group_by(Population, Cluster, Ecoregion) %>%
  summarise(Latitude  = mean(Latitude),
            Longitude = mean(Longitude),
            .groups = "drop")

us_map <- map_data("state")

p_geo_map <- ggplot() +
  geom_polygon(data = us_map, aes(x = long, y = lat, group = group),
               fill = "grey92", color = "white", linewidth = 0.3) +
  geom_point(data = pop_geo,
             aes(x = Longitude, y = Latitude, color = Cluster),
             size = 4, alpha = 0.85) +
  geom_text(data = pop_geo,
            aes(x = Longitude, y = Latitude, label = Population),
            size = 2, vjust = -0.8, check_overlap = TRUE) +
  coord_fixed(ratio = 1.3,
              xlim = range(pop_geo$Longitude, na.rm = TRUE) + c(-2, 2),
              ylim = range(pop_geo$Latitude,  na.rm = TRUE) + c(-2, 2)) +
  theme_bw() +
  labs(title = paste0("Population TE cluster membership (k = ", n_clusters, ")"),
       subtitle = "Clusters derived from tSNE embedding of shared TE family proportions",
       x = "Longitude", y = "Latitude")

p_geo_map

ggsave(file.path(out_dir, "geographic_cluster_map.png"),
       p_geo_map, width = 12, height = 7, dpi = 300)

# Also save cluster assignments
write.csv(tsne_df %>% dplyr::select(ID, Population, Ecoregion, State, Cluster),
          file.path(out_dir, "tsne_cluster_assignments.csv"),
          row.names = FALSE)
message("Geographic map saved.")


##################################################
# 7. Variance partitioning                        #
##################################################
# Partitions variance in Hellinger-transformed TE family data into:
#   [a]     = climate alone (bio1–alt)
#   [b]     = shared climate + geography
#   [c]     = geography alone (lat/lon)
#   [d]     = unexplained
# Fractions [a] and [c] are the "pure" effects of each driver.

message("\n========== 7. Variance partitioning ==========")

# Drop rows with any NA in climate or geo (needed for RDA)
complete_idx <- complete.cases(climate_mat, geo_mat)
message(paste("Complete cases for varpart:", sum(complete_idx),
              "of", nrow(MR_fam)))

vp <- varpart(
  te_hel[complete_idx, ],
  climate_mat[complete_idx, ],
  geo_mat[complete_idx, ]
)
print(vp)

# Plot Venn diagram of variance fractions
png(file.path(out_dir, "varpart_venn.png"),
    width = 7, height = 6, units = "in", res = 300)
plot(vp,
     bg      = c("#AED9E0", "#FFA69E"),
     Xnames  = c("Climate\n(bio1–alt)", "Geography\n(lat/lon)"),
     main    = "Variance partitioning of TE family profiles")
dev.off()

# Test significance of each pure fraction via RDA + permutation test
message("\n-- RDA significance: climate (controlling for geography) --")
rda_climate <- rda(te_hel[complete_idx, ] ~
                     . + Condition(Latitude + Longitude),
                   data = cbind(climate_mat, geo_mat)[complete_idx, ])
print(anova(rda_climate, permutations = 999))

message("\n-- RDA significance: geography (controlling for climate) --")
rda_geo <- rda(te_hel[complete_idx, ] ~
                 Latitude + Longitude +
                 Condition(bio1 + bio2 + bio3 + bio4 + bio5 + bio6 +
                             bio7 + bio8 + bio9 + bio10 + bio11 +
                             bio12 + bio13 + bio14 + bio15 + bio16 +
                             bio17 + bio18 + bio19 + alt),
               data = cbind(climate_mat, geo_mat)[complete_idx, ])
print(anova(rda_geo, permutations = 999))


####################
# Summary printout #
####################

message("\n========== SUMMARY ==========")
message(paste("PERMANOVA  R² =", round(perm_result$R2[1], 3),
              " p =", perm_result$`Pr(>F)`[1]))
message(paste("ANOSIM     R  =", round(anosim_result$statistic, 3),
              " p =", anosim_result$signif))
message(paste("Mantel     r  =", round(mantel_result$statistic, 3),
              " p =", round(mantel_result$signif, 4)))
message(paste("Varpart    Climate R² =",
              round(vp$part$indfract$Adj.R.squared[1], 3),
              " | Geography R² =",
              round(vp$part$indfract$Adj.R.squared[3], 3),
              " | Shared =",
              round(vp$part$indfract$Adj.R.squared[2], 3)))
message(paste("Output written to:", out_dir))
