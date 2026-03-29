# TE Analysis Pipeline
# Script: TE_pop_explore.R
#
# Population-level controls on TE dissimilarity.
# Asks: do individuals from the same population share more similar TE profiles
# than individuals from different populations?
#
# Analyses:
#   1. Sample size filter         — retain populations with >= min_n individuals
#   2. PERMANOVA (population)     — does population membership explain TE variation?
#   3. Within vs between          — intra vs inter-population Bray-Curtis comparison
#   4. Consistency scores         — mean intra-population dissimilarity per population
#   5. tSNE small multiples       — per-population panels with convex hulls
#   6. tSNE overview with hulls   — single plot, all populations, colored hulls
#   7. Ecoregion map              — CEC Level II polygons + population points by k-means cluster
#   8. SIMPER                     — which TE families drive between-population differences?
#   9. Pairwise population heatmap — population-averaged Bray-Curtis distance matrix
#
# Input:  2_climateCode/data/MR_fam.csv
# Output: 2_climateCode/results/TE_pop_explore/  (created if absent)


###################
# Library loading #
###################

library(tidyverse)
library(vegan)      # adonis2, betadisper, simper, vegdist, decostand
library(pheatmap)   # pairwise distance heatmap
library(Rtsne)      # tSNE
library(sf)         # read ecoregion shapefile
library(patchwork)  # combining plots


############################
# Parameters — edit here   #
############################

setwd("~/Documents/GitHub/jumpingGenesPipe")

fam_path <- "2_climateCode/data/MR_fam.csv"
out_dir  <- "2_climateCode/results/TE_pop_explore"

# Minimum individuals per population to be included
min_n <- 3

# Ecoregion color palette — used consistently across all plots
eco_palette <- c(
  "8.1"  = "#54286f",
  "8.2"  = "#a989b9",
  "8.3"  = "#8b0000",
  "8.4"  = "#9ab87a",
  "8.5"  = "#304c00",
  "9.2"  = "#cdb79e",
  "9.4"  = "#85683c",
  "10.2" = "#1066bc",
  "11.1" = "#8abff5",
  "13.1" = "#BDDBF9"
)

# tSNE settings
tsne_perplexity <- 30
tsne_seed       <- 142

# k-means clusters
n_clusters <- 5

# Number of TE families to show in SIMPER output
n_simper_fams <- 15

# Number of top-variable families for pairwise heatmap
n_heatmap_fams <- 60


####################
# Setup & load data #
####################

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

MR_fam_all <- read.csv(fam_path, header = TRUE) %>%
  mutate(Ecoregion  = as.factor(Ecoregion),
         Population = as.factor(Population))


######################################################
# 1. Filter to populations with >= min_n individuals #
######################################################

pop_counts <- MR_fam_all %>%
  count(Population, name = "n_individuals")

pops_keep <- pop_counts %>%
  filter(n_individuals >= min_n) %>%
  pull(Population)

MR_fam <- MR_fam_all %>%
  filter(Population %in% pops_keep) %>%
  mutate(Population = droplevels(Population))

excluded <- pop_counts %>%
  filter(n_individuals < min_n)

message("========== Sample size filter ==========")
message(paste("Populations retained (>=", min_n, "individuals):", nlevels(MR_fam$Population)))
message(paste("Populations excluded (<", min_n, "individuals):", nrow(excluded)))
if (nrow(excluded) > 0) {
  message(paste("  Excluded:", paste(excluded$Population, collapse = ", ")))
}
message(paste("Individuals retained:", nrow(MR_fam)))

write.csv(pop_counts, file.path(out_dir, "population_sample_sizes.csv"), row.names = FALSE)

# TE matrix and Bray-Curtis distance for filtered dataset
te_mat <- MR_fam %>%
  dplyr::select(ends_with("_mean")) %>%
  as.matrix()

bc_dist <- vegdist(te_mat, method = "bray")
te_hel  <- decostand(te_mat, method = "hellinger")


##############################################################
# 2. PERMANOVA — does population membership explain TE variation?
##############################################################
# Compare to ecoregion PERMANOVA from TE_prop_explore.R:
# Population R² > Ecoregion R² suggests fine-scale structure dominates.

message("\n========== 2. PERMANOVA (Population) ==========")
perm_pop <- adonis2(
  te_mat ~ Population,
  data         = MR_fam,
  permutations = 999,
  method       = "bray"
)
print(perm_pop)
write.csv(as.data.frame(perm_pop),
          file.path(out_dir, "permanova_population.csv"))

# Homogeneity of dispersion check
disp_pop    <- betadisper(bc_dist, MR_fam$Population)
disp_anova  <- anova(disp_pop)
message(paste("betadisper ANOVA p =", round(disp_anova$`Pr(>F)`[1], 4),
              " (non-significant = equal dispersion = PERMANOVA is trustworthy)"))


####################################################
# 3. Within vs between population dissimilarity    #
####################################################
# For each pair of individuals, label whether they are from the
# same population (within) or different populations (between).
# A clear separation confirms population-level structure.

message("\n========== 3. Within vs between population dissimilarity ==========")

bc_mat   <- as.matrix(bc_dist)
pop_vec  <- as.character(MR_fam$Population)
n        <- nrow(bc_mat)

within_vals  <- c()
between_vals <- c()

for (i in 1:(n - 1)) {
  for (j in (i + 1):n) {
    val <- bc_mat[i, j]
    if (pop_vec[i] == pop_vec[j]) {
      within_vals  <- c(within_vals,  val)
    } else {
      between_vals <- c(between_vals, val)
    }
  }
}

wb_df <- bind_rows(
  data.frame(type = "Within population",  BC = within_vals),
  data.frame(type = "Between population", BC = between_vals)
)

# Wilcoxon test: are within-population distances smaller than between?
wb_test <- wilcox.test(within_vals, between_vals, alternative = "less")
message(paste("Wilcoxon test (within < between): W =",
              round(wb_test$statistic, 0),
              " p =", format(wb_test$p.value, digits = 3)))

p_wb <- ggplot(wb_df, aes(x = type, y = BC, fill = type)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  scale_fill_manual(values = c("#1B4332", "#AED9E0")) +
  theme_bw() +
  labs(title = "Within vs between population TE dissimilarity",
       subtitle = paste0("Wilcoxon p = ", format(wb_test$p.value, digits = 3),
                         "  (within < between, one-sided)"),
       x = NULL, y = "Bray-Curtis dissimilarity") +
  theme(legend.position = "none")
p_wb
ggsave(file.path(out_dir, "within_between_dissimilarity.png"),
       p_wb, width = 7, height = 5, dpi = 300)


####################################################
# 4. Population consistency scores                 #
####################################################
# Mean intra-population Bray-Curtis per population.
# Low score = individuals in that population have similar TE profiles.
# Colored by ecoregion; sized by sample count.
# A correlation with sample size would suggest a sampling artefact.

message("\n========== 4. Population consistency scores ==========")

consistency <- lapply(levels(MR_fam$Population), function(pop) {
  idx    <- which(pop_vec == pop)
  if (length(idx) < 2) return(NULL)
  sub_bc <- bc_mat[idx, idx]
  # Upper triangle only
  mean_bc <- mean(sub_bc[upper.tri(sub_bc)])
  data.frame(
    Population = pop,
    mean_intra_BC = mean_bc,
    n = length(idx)
  )
}) %>%
  bind_rows() %>%
  left_join(MR_fam %>% distinct(Population, Ecoregion), by = "Population")

write.csv(consistency,
          file.path(out_dir, "population_consistency_scores.csv"),
          row.names = FALSE)

# Correlation between consistency score and sample size
cor_test <- cor.test(consistency$n, consistency$mean_intra_BC,
                     method = "spearman")
message(paste("Spearman correlation (n vs mean intra-pop BC): r =",
              round(cor_test$estimate, 3), " p =", round(cor_test$p.value, 4)))

# Scatter: consistency vs sample size
p_consist_n <- ggplot(consistency,
                      aes(x = n, y = mean_intra_BC, color = Ecoregion)) +
  geom_point(size = 3, alpha = 0.85) +
  geom_smooth(method = "lm", color = "black", se = TRUE, linewidth = 0.7) +
  scale_color_manual(values = eco_palette) +
  theme_bw() +
  labs(title = "Population consistency vs sample size",
       subtitle = paste0("Spearman r = ", round(cor_test$estimate, 3),
                         "  p = ", round(cor_test$p.value, 4)),
       x = "Individuals per population",
       y = "Mean intra-population Bray-Curtis")
p_consist_n
ggsave(file.path(out_dir, "consistency_vs_samplesize.png"),
       p_consist_n, width = 8, height = 5, dpi = 300)

# Barplot: consistency scores ranked
p_consist_bar <- consistency %>%
  arrange(mean_intra_BC) %>%
  mutate(Population = factor(Population, levels = Population)) %>%
  ggplot(aes(x = Population, y = mean_intra_BC, fill = Ecoregion)) +
  geom_col() +
  scale_fill_manual(values = eco_palette) +
  theme_bw() +
  labs(title = "Population consistency scores (lower = more similar within population)",
       x = "Population", y = "Mean intra-population Bray-Curtis") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
p_consist_bar
ggsave(file.path(out_dir, "consistency_barplot.png"),
       p_consist_bar, width = 12, height = 5, dpi = 300)


####################################################
# 5. tSNE small multiples with convex hulls        #
####################################################
# Each panel shows one population. All individuals are shown as grey
# background points for context. The focal population is highlighted
# in color with a convex hull. Visually reveals whether individuals
# from the same population cluster together or scatter across tSNE space.

message("\n========== 5. tSNE small multiples ==========")

# Prepare shared TE family matrix for tSNE
te_tsne <- te_mat
te_tsne[te_tsne == 0] <- NA
shared_cols    <- colSums(is.na(te_tsne)) == 0
te_tsne_shared <- te_tsne[, shared_cols]
message(paste("Shared TE families for tSNE:", ncol(te_tsne_shared)))

set.seed(tsne_seed)
tsne_fit <- Rtsne(
  scale(te_tsne_shared),
  perplexity       = tsne_perplexity,
  check_duplicates = FALSE
)

tsne_df <- as.data.frame(tsne_fit$Y) %>%
  rename(tSNE1 = V1, tSNE2 = V2) %>%
  bind_cols(MR_fam %>% dplyr::select(ID, Population, Ecoregion, State,
                                      Latitude, Longitude))

# --- k-means cluster tSNE embedding ---
set.seed(tsne_seed)
km <- kmeans(tsne_df[, c("tSNE1", "tSNE2")], centers = n_clusters, nstart = 25)
tsne_df$Cluster <- as.factor(km$cluster)

# tSNE coloured by k-means cluster
p_tsne_cluster <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = Cluster)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text(aes(label = Population), size = 2.5, nudge_x = 0.4,
            check_overlap = TRUE) +
  theme_bw() +
  labs(title = paste0("tSNE coloured by k-means cluster (k = ", n_clusters, ")"),
       subtitle = paste0("Filtered to populations >= ", min_n,
                         " individuals | Perplexity = ", tsne_perplexity))

# tSNE coloured by ecoregion (for comparison)
p_tsne_eco <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = Ecoregion)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text(aes(label = Population), size = 2.5, nudge_x = 0.4,
            check_overlap = TRUE) +
  scale_color_manual(values = eco_palette) +
  theme_bw() +
  labs(title = "tSNE coloured by ecoregion")

# tSNE coloured by population (each population a distinct color)
p_tsne_pop <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = Population)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(title = "tSNE coloured by population (legend suppressed — too many levels)",
       subtitle = "Use small multiples below to examine individual populations")

p_tsne_cluster + p_tsne_eco
ggsave(file.path(out_dir, "tsne_cluster_vs_ecoregion.png"),
       p_tsne_cluster + p_tsne_eco,
       width = 16, height = 6, dpi = 300)

p_tsne_pop
ggsave(file.path(out_dir, "tsne_by_population.png"),
       p_tsne_pop, width = 10, height = 7, dpi = 300)

# Save cluster assignments alongside population/ecoregion info
write.csv(
  tsne_df %>% dplyr::select(ID, Population, Ecoregion, State, Cluster),
  file.path(out_dir, "tsne_cluster_assignments.csv"),
  row.names = FALSE
)

# Compute convex hulls per population
pop_hulls <- tsne_df %>%
  group_by(Population, Ecoregion) %>%
  filter(n() >= 3) %>%                   # hull needs >= 3 points
  slice(chull(tSNE1, tSNE2)) %>%
  ungroup()

# Background points (all individuals, used in every panel)
bg_all <- tsne_df %>% dplyr::select(tSNE1, tSNE2)

# Build small multiples: grey background + coloured focal population + hull
p_small_mult <- ggplot() +
  # Grey background — all individuals in every panel
  geom_point(data = bg_all,
             aes(x = tSNE1, y = tSNE2),
             color = "grey80", size = 0.8, alpha = 0.5) +
  # Focal population points
  geom_point(data = tsne_df,
             aes(x = tSNE1, y = tSNE2, color = Ecoregion),
             size = 2, alpha = 0.9) +
  # Convex hull
  geom_polygon(data = pop_hulls,
               aes(x = tSNE1, y = tSNE2, fill = Ecoregion, group = Population),
               alpha = 0.15, linewidth = 0.5, color = "black") +
  scale_color_manual(values = eco_palette) +
  scale_fill_manual(values  = eco_palette) +
  facet_wrap(~ Population, ncol = 8) +
  theme_bw(base_size = 8) +
  theme(legend.position  = "bottom",
        strip.text       = element_text(size = 7, face = "bold"),
        axis.text        = element_blank(),
        axis.ticks       = element_blank(),
        panel.grid       = element_blank()) +
  labs(title = paste0("tSNE small multiples: within-population clustering (n >= ", min_n, ")"),
       subtitle = "Grey = all individuals. Color + hull = focal population.",
       x = "tSNE1", y = "tSNE2")

n_pops   <- nlevels(MR_fam$Population)
n_cols   <- 8
n_rows   <- ceiling(n_pops / n_cols)
fig_h    <- max(8, n_rows * 2.2)

ggsave(file.path(out_dir, "tsne_small_multiples.png"),
       p_small_mult,
       width = 20, height = fig_h, dpi = 200, limitsize = FALSE)
message(paste("Small multiples saved:", n_pops, "panels,", n_cols, "columns."))


####################################################
# 6. tSNE overview — all hulls on one plot         #
####################################################
# Complements the small multiples. Shows all population hulls together,
# colored by ecoregion. Hull overlap indicates populations with similar
# TE profiles; isolated hulls indicate distinctive populations.

p_hull_overview <- ggplot() +
  geom_polygon(data = pop_hulls,
               aes(x = tSNE1, y = tSNE2,
                   fill = Ecoregion, group = Population),
               alpha = 0.12, color = "grey40", linewidth = 0.3) +
  geom_point(data = tsne_df,
             aes(x = tSNE1, y = tSNE2, color = Ecoregion),
             size = 1.5, alpha = 0.7) +
  scale_color_manual(values = eco_palette) +
  scale_fill_manual(values  = eco_palette) +
  theme_bw() +
  labs(title = "tSNE overview: per-population convex hulls colored by ecoregion",
       subtitle = "Overlapping hulls = similar TE profiles; isolated hulls = distinctive populations")
p_hull_overview
ggsave(file.path(out_dir, "tsne_hull_overview.png"),
       p_hull_overview, width = 12, height = 8, dpi = 300)


####################################################
# 7. Ecoregion map with population k-means clusters #
####################################################
# EPA Level III Ecoregion polygons as the base layer, filled by ecoregion.
# Population centroids overlaid as points colored by tSNE k-means cluster.
# Directly links the TE-based clustering back to geography and environment.

message("\n========== 7. Ecoregion map ==========")

eco_shp_path  <- "2_climateCode/data/ecoregions/na_cec_eco_l2/NA_CEC_Eco_Level2.shp"
eco_code_field <- "NA_L2CODE"   # update if shapefile uses a different field name

if (!file.exists(eco_shp_path)) {
  warning("Ecoregion shapefile not found at: ", eco_shp_path,
          "\nSkipping map.")
} else {

  # Load shapefiles and ensure both are in WGS84
  ecoregions_sf <- st_read(eco_shp_path, quiet = TRUE) %>%
    st_transform(4326)

  states_sf <- st_read("2_climateCode/data/s_16ap26/conus.shp", quiet = TRUE) %>%
    st_transform(4326)

  message("Shapefile columns: ", paste(names(ecoregions_sf), collapse = ", "))

  eco_levels <- sort(unique(as.character(MR_fam$Ecoregion)))

  # All codes in shapefile: data ecoregions get named colors, extras get grey
  all_l2_codes    <- unique(as.character(ecoregions_sf[[eco_code_field]]))
  extra_codes     <- setdiff(all_l2_codes, eco_levels)
  eco_fill_values <- c(eco_palette,
                       setNames(rep("grey85", length(extra_codes)), extra_codes))

  # Population centroids: one row per (Population × Cluster) for coloring points.
  # Populations split across clusters will appear as multiple points.
  pop_centroids <- tsne_df %>%
    group_by(Population, Cluster, Ecoregion) %>%
    summarise(Longitude = mean(Longitude, na.rm = TRUE),
              Latitude  = mean(Latitude,  na.rm = TRUE),
              n         = n(),
              .groups   = "drop") %>%
    filter(!is.na(Longitude) & !is.na(Latitude))

  # One label per population regardless of cluster split — placed at the
  # overall geographic centroid of all individuals in that population
  pop_labels <- tsne_df %>%
    group_by(Population) %>%
    summarise(Longitude = mean(Longitude, na.rm = TRUE),
              Latitude  = mean(Latitude,  na.rm = TRUE),
              .groups   = "drop") %>%
    filter(!is.na(Longitude) & !is.na(Latitude)) %>%
    mutate(label = str_remove(Population, "^p_"))

  p_eco_map <- ggplot() +
    geom_sf(data = ecoregions_sf,
            aes(fill = .data[[eco_code_field]]),
            color = "white", linewidth = 0.15, alpha = 0.55) +
    geom_sf(data = states_sf,
            fill = NA, color = "grey30", linewidth = 0.35) +
    geom_point(data = pop_centroids,
               aes(x = Longitude, y = Latitude, color = Cluster, size = n),
               alpha = 0.9, shape = 21, stroke = 1.2) +
    geom_text(data = pop_labels,
              aes(x = Longitude, y = Latitude, label = label),
              size = 2.2, vjust = -0.8, color = "black") +
    scale_fill_manual(values   = eco_fill_values,
                      name     = "CEC L2 Ecoregion",
                      na.value = "grey85") +
    scale_size_continuous(range = c(3, 7), name = "Individuals") +
    coord_sf(xlim = range(pop_centroids$Longitude) + c(-3, 3),
             ylim = range(pop_centroids$Latitude)  + c(-3, 3)) +
    theme_bw() +
    labs(title    = "Population TE cluster membership on CEC Level II Ecoregions",
         subtitle = paste0("Point color = tSNE k-means cluster (k = ", n_clusters,
                           ")  |  Fill = Level II ecoregion  |  Size = n individuals"),
         x = "Longitude", y = "Latitude",
         color = paste0("k-means\ncluster (k=", n_clusters, ")")) +
    theme(legend.position = "right")

  p_eco_map
  ggsave(file.path(out_dir, "ecoregion_cluster_map.png"),
         p_eco_map, width = 14, height = 8, dpi = 300)
  message("Ecoregion map saved.")
}


####################################################
# 8. SIMPER — which families drive differences?    #
####################################################
# Identifies TE families that contribute most to Bray-Curtis dissimilarity
# between populations. Useful for knowing *what* the biological signal is.
# Run on population-averaged data to keep it tractable.

message("\n========== 8. SIMPER ==========")
# SIMPER identifies which TE families contribute most to Bray-Curtis
# dissimilarity *between populations* (not ecoregions).
# Uses individual-level data grouped by Population.
# With many populations this produces many pairwise comparisons;
# we aggregate across all pairs to find universally important families.

# Join ecoregion for downstream use
pop_eco <- MR_fam %>%
  distinct(Population, Ecoregion) %>%
  arrange(Population)

eco_levels <- sort(unique(as.character(pop_eco$Ecoregion)))

# Run SIMPER at the population level using individual observations
# permutations = 0 skips significance testing (fast); increase for p-values
simper_result <- simper(te_mat,
                        group        = MR_fam$Population,
                        permutations = 0)

simper_summary <- summary(simper_result, ordered = TRUE)
n_pairs <- length(simper_summary)
message(paste("SIMPER computed across", n_pairs, "population pairs"))

# --- Extract and aggregate contributions across all pairwise comparisons ---
all_contribs <- lapply(names(simper_summary), function(pair) {
  df         <- as.data.frame(simper_summary[[pair]])
  df$Family  <- rownames(df)
  df$Pair    <- pair
  # Split pair name into the two populations being compared
  pops       <- str_split(pair, "_", n = 2)[[1]]
  df$Pop_A   <- pops[1]
  df$Pop_B   <- pops[2]
  df
}) %>% bind_rows()

# Overall top families: mean contribution across ALL pairwise comparisons
top_fams_simper <- all_contribs %>%
  group_by(Family) %>%
  summarise(
    mean_contrib  = mean(average,  na.rm = TRUE),
    sd_contrib    = sd(average,    na.rm = TRUE),
    n_pairs_present = sum(!is.na(average)),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_contrib)) %>%
  slice_head(n = n_simper_fams) %>%
  mutate(Family = str_remove(Family, "_mean$"))

write.csv(top_fams_simper,
          file.path(out_dir, "simper_top_families_populations.csv"),
          row.names = FALSE)

# Bar chart: top families overall
p_simper <- ggplot(top_fams_simper,
                   aes(x = reorder(Family, mean_contrib), y = mean_contrib)) +
  geom_col(fill = "#1B4332") +
  geom_errorbar(aes(ymin = mean_contrib - sd_contrib,
                    ymax = mean_contrib + sd_contrib),
                width = 0.35, color = "grey40") +
  coord_flip() +
  theme_bw() +
  labs(title    = paste0("Top ", n_simper_fams,
                         " TE families driving between-population dissimilarity (SIMPER)"),
       subtitle = paste0("Mean ± SD contribution to Bray-Curtis across ",
                         n_pairs, " population pairs"),
       x = "TE Family", y = "Mean contribution to dissimilarity")
p_simper
ggsave(file.path(out_dir, "simper_top_families.png"),
       p_simper, width = 9, height = 6, dpi = 300)


# --- Heatmap: proportion of top SIMPER families per population ---
# Shows WHICH populations have elevated/depleted levels of the key families.
# Rows = populations sorted by ecoregion; columns = top SIMPER families.
# Cell color = z-scored mean proportion (red = elevated, blue = depleted).

top_fam_cols <- paste0(top_fams_simper$Family, "_mean")
top_fam_cols <- top_fam_cols[top_fam_cols %in% colnames(te_mat)]

pop_top_mat <- MR_fam %>%
  group_by(Population) %>%
  summarise(across(all_of(top_fam_cols), mean), .groups = "drop") %>%
  left_join(pop_eco, by = "Population") %>%
  arrange(Ecoregion, Population) %>%
  column_to_rownames("Population")

# Separate annotation and matrix
anno_simper <- data.frame(
  Ecoregion = as.character(pop_top_mat$Ecoregion),
  row.names = rownames(pop_top_mat)
)
pop_top_mat <- pop_top_mat %>%
  dplyr::select(all_of(top_fam_cols)) %>%
  as.matrix()
colnames(pop_top_mat) <- str_remove(colnames(pop_top_mat), "_mean$")

png(file.path(out_dir, "simper_families_by_population.png"),
    width = 14, height = 10, units = "in", res = 300)
pheatmap(
  pop_top_mat,
  annotation_row    = anno_simper,
  annotation_colors = list(Ecoregion = eco_palette[eco_levels]),
  scale             = "column",            # z-score within each family
  cluster_rows      = TRUE,
  cluster_cols      = FALSE,               # keep families in SIMPER rank order
  clustering_distance_rows = "correlation",
  color             = colorRampPalette(c("#1066bc","white","#8b0000"))(100),
  fontsize_row      = 7,
  fontsize_col      = 8,
  main              = paste0("Top ", n_simper_fams,
                             " SIMPER families across populations (z-scored proportion)",
                             "\nColumns ranked by contribution to between-population dissimilarity")
)
dev.off()
message("SIMPER heatmap saved.")


####################################################
# 9. Pairwise population distance heatmap          #
####################################################
# Bray-Curtis dissimilarity between population-averaged TE profiles.
# Clustering reveals which populations have the most similar fingerprints.
# Row/column annotation = ecoregion.

message("\n========== 9. Pairwise population heatmap ==========")

# Top N most variable families for readability
fam_var   <- apply(pop_avg_mat, 2, var)
top_fams  <- names(sort(fam_var, decreasing = TRUE))[1:min(n_heatmap_fams, ncol(pop_avg_mat))]
pop_sub   <- pop_avg_mat[, top_fams]
colnames(pop_sub) <- str_remove(colnames(pop_sub), "_mean$")

# Bray-Curtis distance between populations
pop_bc    <- as.matrix(vegdist(pop_sub, method = "bray"))

# Annotation
eco_levels <- sort(unique(as.character(pop_eco$Ecoregion)))
anno_pop <- data.frame(
  Ecoregion  = as.character(pop_eco$Ecoregion),
  row.names  = pop_eco$Population
)
# Align annotation to matrix row order
anno_pop <- anno_pop[rownames(pop_bc), , drop = FALSE]

png(file.path(out_dir, "pairwise_population_heatmap.png"),
    width = 14, height = 12, units = "in", res = 300)
pheatmap(
  pop_bc,
  annotation_row    = anno_pop,
  annotation_col    = anno_pop,
  annotation_colors = list(Ecoregion = eco_palette[eco_levels]),
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  color             = colorRampPalette(c("#f7fbff","#2171b5","#08306b"))(100),
  fontsize_row      = 6,
  fontsize_col      = 6,
  main              = "Pairwise Bray-Curtis dissimilarity between populations (top variable families)"
)
dev.off()
message("Pairwise heatmap saved.")


####################
# Summary printout #
####################

message("\n========== SUMMARY ==========")
message(paste("Populations analyzed:", nlevels(MR_fam$Population),
              " | Individuals:", nrow(MR_fam)))
message(paste("PERMANOVA (Population) R² =",
              round(perm_pop$R2[1], 3),
              " p =", perm_pop$`Pr(>F)`[1]))
message(paste("betadisper p =",
              round(disp_anova$`Pr(>F)`[1], 4)))
message(paste("Wilcoxon within < between: p =",
              format(wb_test$p.value, digits = 3)))
message(paste("Mean within-pop BC:",
              round(mean(within_vals), 4),
              " | Mean between-pop BC:",
              round(mean(between_vals), 4)))
message(paste("Most consistent population:",
              consistency$Population[which.min(consistency$mean_intra_BC)],
              "(mean BC =", round(min(consistency$mean_intra_BC), 4), ")"))
message(paste("Least consistent population:",
              consistency$Population[which.max(consistency$mean_intra_BC)],
              "(mean BC =", round(max(consistency$mean_intra_BC), 4), ")"))
message(paste("Top SIMPER family:", top_fams_simper$Family[1]))
message(paste("Output written to:", out_dir))
