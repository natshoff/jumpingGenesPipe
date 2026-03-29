# TE Analysis Pipeline
# Script 2: TE Proportion Plots
#
# Description: Creates plots from master data files prepared in 1_TEpipe_MRfileprep.R
# Input:  2_climateCode/data/MR_sup.csv
#         2_climateCode/data/MR_fam.csv
# Output: Various PCA, tSNE, boxplot, and variance summary CSVs


###################
# Library loading #
###################

library(tidyverse)
library(PerformanceAnalytics)
library(recipes)
library(Rtsne)
library(tidytext)


##############################
# Load in MR files from csvs #
##############################

setwd("~/Documents/GitHub/jumpingGenesPipe")

# SUPERFAMILY
MR_sup <- read.csv("2_climateCode/data/MR_sup.csv", header = TRUE)
# NOTE: optional filtering of outlier (add filter pipe to MR_sup before building MR_sup_add)
# MR_sup <- MR_sup %>% filter(ID != "JG_CO_054_F")

# Convert ecoregion to factor
MR_sup$Ecoregion <- as.factor(MR_sup$Ecoregion)

# Add (1) LTR proportion and (2) total TE content to master file
# Ecoregion_num converts the ecoregion factor to a numeric for ordered axis plotting
MR_sup_add <- MR_sup %>%
  mutate(
    LTR_mean      = Copia_mean + Gypsy_mean,
    TE_tot        = rowSums(across(ends_with("_mean"))),
    Ecoregion_num = as.numeric(Ecoregion)
  )


# FAMILY
MR_fam <- read.csv("2_climateCode/data/MR_fam.csv", header = TRUE)
# NOTE: optional filtering of outlier
# MR_fam <- MR_fam %>% filter(ID != "JG_CO_054_F")

# Convert ecoregion to factor and add numeric version
MR_fam <- MR_fam %>%
  mutate(
    Ecoregion     = as.factor(Ecoregion),
    Ecoregion_num = as.numeric(Ecoregion)
  )


######################################
# Inter / Intra population variation #
######################################

###########
# BOXPLOT #
###########

# LTR proportion across populations, ordered and colored by ecoregion
ggplot() +
  geom_boxplot(data = MR_sup_add,
               mapping = aes(x = reorder(Population, Ecoregion_num),
                             y = LTR_mean, color = Ecoregion),
               outlier.shape = 4) +
  geom_jitter(data = MR_sup_add,
              mapping = aes(x = reorder(Population, Ecoregion_num),
                            y = LTR_mean, color = Ecoregion),
              width = 0.35, size = 3, alpha = 0.5) +
  scale_color_manual(values = c(
    "#54286f",  # 8.1
    "#a989b9",  # 8.2
    "#8b0000",  # 8.3
    "#9ab87a",  # 8.4
    "#304c00",  # 8.5
    "#cdb79e",  # 9.2
    "#85683c",  # 9.4
    "#1066bc",  # 10.2
    "#8abff5",  # 11.1
    "#BDDBF9")) +
  theme_bw() +
  ggtitle("LTR proportions across populations") +
  theme(axis.text.x = element_text(angle = 20))

# Total TE proportion across populations
ggplot() +
  geom_boxplot(data = MR_sup_add,
               mapping = aes(x = reorder(Population, Ecoregion_num),
                             y = TE_tot, color = Ecoregion),
               outlier.shape = 4) +
  geom_jitter(data = MR_sup_add,
              mapping = aes(x = reorder(Population, Ecoregion_num),
                            y = TE_tot, color = Ecoregion),
              width = 0.25, size = 3, alpha = 0.5) +
  scale_color_manual(values = c(
    "#54286f",
    "#a989b9",
    "#8b0000",
    "#9ab87a",
    "#304c00",
    "#cdb79e",
    "#85683c",
    "#1066bc",
    "#8abff5",
    "#BDDBF9")) +
  theme_bw() +
  ggtitle("Total TE proportions across populations") +
  xlab("Population Number") +
  ylab("Total TE Proportion") +
  theme(axis.text.x = element_text(angle = 20))

# Total TE proportion — transparent background version
ggplot() +
  geom_boxplot(data = MR_sup_add,
               mapping = aes(x = reorder(Population, Ecoregion_num),
                             y = TE_tot, color = Ecoregion),
               outlier.shape = 4, fill = "white", alpha = 0.5, linewidth = 0.75) +
  scale_color_manual(values = c(
    "#54286f",
    "#a989b9",
    "#8b0000",
    "#9ab87a",
    "#304c00",
    "#cdb79e",
    "#85683c",
    "#1066bc",
    "#8abff5",
    "#BDDBF9")) +
  geom_jitter(data = MR_sup_add,
              mapping = aes(x = reorder(Population, Ecoregion_num),
                            y = TE_tot, fill = Ecoregion),
              width = 0.25, size = 3, color = "black", shape = 21) +
  scale_fill_manual(values = alpha(c(
    "#54286f",
    "#a989b9",
    "#8b0000",
    "#9ab87a",
    "#304c00",
    "#cdb79e",
    "#85683c",
    "#1066bc",
    "#8abff5",
    "#BDDBF9"), 0.6)) +
  theme_bw() +
  ggtitle("Total TE proportions across populations") +
  xlab("Population Number") +
  ylab("Total TE Proportion") +
  theme(axis.text.x = element_text(angle = 20),
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA))

ggsave("TEprop_pops.png", height = 8, width = 15, dpi = 300, bg = "transparent")


# Just populations 045 and 229
MR_sup_add_reduced <- MR_sup_add %>%
  filter(Population %in% c("p_045", "p_229"))

ggplot() +
  geom_boxplot(data = MR_sup_add_reduced,
               mapping = aes(x = reorder(Population, Ecoregion_num),
                             y = TE_tot, color = Ecoregion),
               outlier.shape = 4, fill = "white", alpha = 0.5, linewidth = 0.75) +
  scale_color_manual(values = c(
    "#54286f",
    "#a989b9",
    "#8b0000",
    "#9ab87a",
    "#304c00",
    "#cdb79e",
    "#85683c",
    "#1066bc",
    "#8abff5",
    "#BDDBF9")) +
  geom_jitter(data = MR_sup_add_reduced,
              mapping = aes(x = reorder(Population, Ecoregion_num),
                            y = TE_tot, fill = Ecoregion),
              width = 0.25, size = 3, color = "black", shape = 21) +
  scale_fill_manual(values = alpha(c(
    "#54286f",
    "#a989b9",
    "#8b0000",
    "#9ab87a",
    "#304c00",
    "#cdb79e",
    "#85683c",
    "#1066bc",
    "#8abff5",
    "#BDDBF9"), 0.6)) +
  theme_bw() +
  ggtitle("Total TE proportions: populations 045 and 229") +
  xlab("Population Number") +
  ylab("Total TE Proportion") +
  theme(axis.text.x = element_text(angle = 20),
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA))

ggsave("TEprop_pops_reduced.png", height = 8, width = 15, dpi = 300, bg = "transparent")


############
# Variance #
############

# Average proportion of each superfamily across all individuals
# NOTE: update column names here if your superfamily set differs
Sup_mean <- MR_sup_add %>%
  summarise(Copia       = mean(Copia_mean),
            Gypsy       = mean(Gypsy_mean),
            Harbinger   = mean(Harbinger_mean),
            Helitron    = mean(Helitron_mean),
            L1          = mean(L1_mean),
            MSAT        = mean(MSAT_mean),
            MuDR        = mean(MuDR_mean),
            RTE         = mean(RTE_mean),
            SINE2.tRNA  = mean(SINE2.tRNA_mean),
            hAT         = mean(hAT_mean),
            unclassified = mean(unclassified_mean))

Sup_sd <- MR_sup_add %>%
  summarise(Copia       = sd(Copia_mean),
            Gypsy       = sd(Gypsy_mean),
            Harbinger   = sd(Harbinger_mean),
            Helitron    = sd(Helitron_mean),
            L1          = sd(L1_mean),
            MSAT        = sd(MSAT_mean),
            MuDR        = sd(MuDR_mean),
            RTE         = sd(RTE_mean),
            SINE2.tRNA  = sd(SINE2.tRNA_mean),
            hAT         = sd(hAT_mean),
            unclassified = sd(unclassified_mean))
write.csv(Sup_sd, "Superfamily_sd.csv", row.names = FALSE)


# Individuals
Var_ind <- MR_sup_add %>%
  summarise(TE_tot_mean  = mean(TE_tot),
            TE_tot_sd    = sd(TE_tot),
            TE_tot_min   = min(TE_tot),
            TE_tot_max   = max(TE_tot),
            TE_tot_range = TE_tot_max - TE_tot_min,
            TE_tot_var   = var(TE_tot))

# Within populations
Var_intra <- MR_sup_add %>%
  group_by(Population) %>%
  summarise(TE_tot_mean  = mean(TE_tot),
            TE_tot_sd    = sd(TE_tot),
            TE_tot_min   = min(TE_tot),
            TE_tot_max   = max(TE_tot),
            TE_tot_range = TE_tot_max - TE_tot_min,
            TE_tot_var   = var(TE_tot),
            Sample_per   = n())
Var_mean_intra <- Var_intra %>%
  summarise(TE_tot_mean_var   = mean(TE_tot_var, na.rm = TRUE),
            TE_tot_mean_range = mean(TE_tot_range, na.rm = TRUE))

# Between populations
Var_inter <- Var_intra %>%
  summarise(TE_pop_mean  = mean(TE_tot_mean),
            TE_pop_sd    = sd(TE_tot_mean),
            TE_pop_min   = min(TE_tot_mean),
            TE_pop_max   = max(TE_tot_mean),
            TE_pop_range = TE_pop_max - TE_pop_min,
            TE_pop_var   = var(TE_tot_mean, na.rm = TRUE))

write.csv(Var_ind,   "AI_individual_var.csv",  row.names = FALSE)
write.csv(Var_intra, "AI_intrapop_var.csv",    row.names = FALSE)
write.csv(Var_inter, "AI_interpop_var.csv",    row.names = FALSE)

# Intrapopulation variance
Pop_var <- MR_sup_add %>%
  group_by(Population) %>%
  summarise(TEtot_var = var(TE_tot),
            LTR_var   = var(LTR_mean))

Mean_intra_var <- Pop_var %>%
  summarise(TEtot_var_mean = mean(TEtot_var, na.rm = TRUE),
            TEtot_var_sd   = sd(TEtot_var, na.rm = TRUE),
            LTR_var_mean   = mean(LTR_var, na.rm = TRUE),
            LTR_var_sd     = sd(LTR_var, na.rm = TRUE))

# Interpopulation variance
Ind_var <- MR_sup_add %>%
  summarise(TEtot_var = var(TE_tot),
            LTR_var   = var(LTR_mean))

Pop_means <- MR_sup_add %>%
  group_by(Population) %>%
  summarise(TEtot_pop_mean = mean(TE_tot),
            TEtot_pop_min  = min(TE_tot),
            TEtot_pop_max  = max(TE_tot),
            LTR_pop_mean   = mean(LTR_mean),
            LTR_pop_min    = min(LTR_mean),
            LTR_pop_max    = max(LTR_mean)) %>%
  mutate(TEtot_range   = TEtot_pop_max - TEtot_pop_min,
         LTR_pop_range = LTR_pop_max   - LTR_pop_min)

Mean_inter_var <- Pop_means %>%
  summarise(TEtot_var = var(TEtot_pop_mean),
            LTR_var   = var(LTR_pop_mean))


##################################
# PCA of Superfamily differences #
##################################

######################
# Correlation matrix #
######################

# Spearman correlation on all superfamily proportions (non-parametric)
MR_sup %>%
  dplyr::select(ends_with("_mean")) %>%
  chart.Correlation(histogram = TRUE, pch = 19, method = "spearman")


###################
# Superfamily PCA #
###################

# Subset to metadata + superfamily mean proportions
sup_pca_data <- MR_sup %>%
  dplyr::select(ID, Population, Individual, State, Locality,
                Latitude, Longitude, Ecoregion, ends_with("_mean"))

sup_recipe <-
  recipe(~., data = sup_pca_data) %>%
  update_role(ID, Population, Individual, State, Locality,
              Latitude, Longitude, Ecoregion, new_role = "id") %>%
  step_naomit(all_predictors()) %>%
  step_normalize(all_predictors()) %>%
  step_pca(all_predictors(), id = "pca") %>%
  prep()

# Execute PCA
sup_pca <- sup_recipe %>% tidy(id = "pca")
sup_pca

# Visualize variance explained by each component
sup_recipe %>%
  tidy(id = "pca", type = "variance") %>%
  dplyr::filter(terms == "percent variance") %>%
  ggplot(aes(x = component, y = value)) +
  geom_col(fill = "#1B4332") +
  xlim(c(0, 15)) +
  ylab("% of total variance")

# Plot PCA loadings to identify TE types driving each PC
sup_pca %>%
  mutate(terms = tidytext::reorder_within(terms, abs(value), component)) %>%
  ggplot(aes(abs(value), terms, fill = value > 0)) +
  geom_col() +
  facet_wrap(~component, scales = "free_y") +
  tidytext::scale_y_reordered() +
  scale_fill_manual(values = c("#52B788", "#1B4332")) +
  labs(x = "Absolute value of contribution", y = NULL, fill = "Positive?")

# Widen PCA loadings for biplot arrows
sup_wider <- sup_pca %>%
  tidyr::pivot_wider(names_from = component, id_cols = terms)

# PC1 vs PC2 biplot
arrow_style <- arrow(length = unit(0.05, "inches"), type = "closed")

sup_plot12 <-
  bake(sup_recipe, new_data = NULL) %>%
  ggplot(aes(PC1, PC2, label = Population, color = Ecoregion)) +
  geom_point(aes(color = Ecoregion), size = 5) +
  geom_text(check_overlap = TRUE, size = 4, hjust = "left", nudge_x = 0.07) +
  labs(color = NULL) +
  geom_segment(inherit.aes = FALSE, data = sup_wider,
               aes(xend = PC1 * 5, yend = PC2 * 5),
               x = 0, y = 0, arrow = arrow_style) +
  geom_text(inherit.aes = FALSE, data = sup_wider,
            aes(x = PC1 * 5, y = PC2 * 5, label = terms),
            hjust = "outward", nudge_x = 0.00,
            vjust = 0, nudge_y = 0.00, size = 5, color = "black")

sup_plot12 +
  ggtitle("PCA of TE Proportions") +
  theme_bw() +
  xlab("PC1") +
  ylab("PC2")


#############################
# PCA of Family differences #
#############################

fam_pca_data <- MR_fam %>%
  dplyr::select(ID, Population, Individual, State, Locality,
                Latitude, Longitude, Ecoregion, ends_with("_mean"))

fam_recipe <-
  recipe(~., data = fam_pca_data) %>%
  update_role(ID, Population, Individual, State, Locality,
              Latitude, Longitude, Ecoregion, new_role = "id") %>%
  step_naomit(all_predictors()) %>%
  step_normalize(all_predictors()) %>%
  step_pca(all_predictors(), id = "pca") %>%
  prep()

fam_pca <- fam_recipe %>% tidy(id = "pca")
fam_pca

# Variance explained
fam_recipe %>%
  tidy(id = "pca", type = "variance") %>%
  dplyr::filter(terms == "percent variance") %>%
  ggplot(aes(x = component, y = value)) +
  geom_col(fill = "#1B4332") +
  xlim(c(0, 350)) +
  ylab("% of total variance")

# Widen loadings for arrows
fam_wider <- fam_pca %>%
  tidyr::pivot_wider(names_from = component, id_cols = terms)

# PC1 vs PC2 biplot
fam_plot12 <-
  bake(fam_recipe, new_data = NULL) %>%
  ggplot(aes(PC1, PC2, label = Population, color = Ecoregion)) +
  geom_point(aes(color = Ecoregion), size = 5) +
  geom_text(check_overlap = TRUE, size = 4, hjust = "left", nudge_x = 0.07) +
  labs(color = NULL)
  #geom_segment(inherit.aes = FALSE, data = fam_wider,
  #             aes(xend = PC1 * 500, yend = PC2 * 500),
  #             x = 0, y = 0, arrow = arrow_style) +
  #geom_text(inherit.aes = FALSE, data = fam_wider,
  #          aes(x = PC1 * 500, y = PC2 * 500, label = terms),
  #          hjust = "outward", nudge_x = 0.00,
  #          vjust = 0, nudge_y = 0.00, size = 5, color = "black")

fam_plot12 +
  ggtitle("PCA of Family TE Proportions") +
  theme_bw() +
  xlab("PC1") +
  ylab("PC2")

# Zoom plot
fam_plot12 +
  ggtitle("PCA of Family TE Proportions") +
  theme_bw() +
  xlab("PC1") +
  ylab("PC2") +
  coord_cartesian(ylim = c(-30, -10), xlim = c(-40, -5))


##############
# tSNE plots #
##############

# ---------- FAMILY tSNE ----------

# Select metadata + family mean proportions, add row ID for rejoining after tSNE
tSNE_MR_fam <- MR_fam %>%
  dplyr::select(ID:Ecoregion, ends_with("_mean")) %>%
  mutate(ID2 = row_number())

# Store metadata for rejoining
tSNE_MR_fam_meta <- tSNE_MR_fam %>%
  dplyr::select(ID:Ecoregion, ID2)

# Work only with TE proportion columns for filtering and tSNE
te_cols_fam <- tSNE_MR_fam %>%
  dplyr::select(ends_with("_mean"), ID2)

# Replace 0s with NA — families absent from a sample are treated as missing
te_cols_fam[te_cols_fam == 0] <- NA

# Keep only columns (families) present in ALL individuals
test_rep <- te_cols_fam[, colSums(is.na(te_cols_fam)) == 0]
message(paste("TE families shared across all samples:", ncol(test_rep) - 1))

# Perform tSNE on shared families
set.seed(142)
tSNE_rep_fit <- test_rep %>%
  column_to_rownames("ID2") %>%
  scale() %>%
  Rtsne()

# Extract components and rejoin metadata
tSNE_rep_df <- tSNE_rep_fit$Y %>%
  as.data.frame() %>%
  rename(tSNE1 = "V1", tSNE2 = "V2") %>%
  mutate(ID2 = row_number()) %>%
  inner_join(tSNE_MR_fam_meta, by = "ID2")

# Colored by ecoregion
tSNE_eco_rep <- tSNE_rep_df %>%
  ggplot(aes(x = tSNE1, y = tSNE2, color = Ecoregion)) +
  geom_point(size = 5) +
  scale_color_manual(values = c(
    "#54286f", "#a989b9", "#8b0000", "#9ab87a", "#304c00",
    "#cdb79e", "#85683c", "#1066bc", "#8abff5", "#BDDBF9")) +
  geom_text(mapping = aes(label = State), nudge_x = 0.5) +
  ggtitle("tSNE clustering of TE families (shared), colored by ecoregion") +
  theme_bw()
tSNE_eco_rep

# Transparent ecoregion tSNE
ggplot() +
  geom_point(data = tSNE_rep_df,
             mapping = aes(x = tSNE1, y = tSNE2, fill = Ecoregion),
             size = 5, shape = 21, color = "black") +
  scale_fill_manual(values = alpha(c(
    "#54286f", "#a989b9", "#8b0000", "#9ab87a", "#304c00",
    "#cdb79e", "#85683c", "#1066bc", "#8abff5", "#BDDBF9"), 0.8)) +
  geom_text(data = tSNE_rep_df,
            mapping = aes(x = tSNE1, y = tSNE2, label = State, color = Ecoregion),
            nudge_x = 0.5, show.legend = FALSE) +
  scale_color_manual(values = c(
    "#54286f", "#a989b9", "#8b0000", "#9ab87a", "#304c00",
    "#cdb79e", "#85683c", "#1066bc", "#8abff5", "#BDDBF9")) +
  ggtitle("tSNE clustering of TE families (shared), colored by ecoregion") +
  theme_bw() +
  theme(legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA))
ggsave("TEfam_tSNE.png", height = 8, width = 15, dpi = 300, bg = "transparent")

# Colored by population
tSNE_pop_rep <- tSNE_rep_df %>%
  ggplot(aes(x = tSNE1, y = tSNE2, color = Population)) +
  geom_point(size = 5) +
  geom_text(mapping = aes(label = State), nudge_x = 0.5) +
  ggtitle("tSNE clustering of TE families (shared), colored by population") +
  theme_bw()
tSNE_pop_rep

# Merge family tSNE output with superfamily proportions to assess clustering drivers
tSNE_rep_df_merge <- tSNE_rep_df %>%
  left_join(MR_sup_add)
tSNE_rep_df_merge <- tSNE_rep_df_merge %>%
  mutate(TE_tot_log = log(TE_tot))

# Colored by log-transformed TE total
tSNE_sup_rep <- tSNE_rep_df_merge %>%
  ggplot(aes(x = tSNE1, y = tSNE2, color = TE_tot_log)) +
  geom_point(size = 5) +
  geom_text(mapping = aes(label = State), nudge_x = 0.5) +
  ggtitle("tSNE clustering of TE families (shared), colored by TE total") +
  theme_bw()
tSNE_sup_rep

# Outlier assessment via IQR
boxplot.stats(tSNE_rep_df_merge$TE_tot)$out
tSNE_rep_df_merge %>%
  ggplot(aes(x = tSNE1, y = tSNE2, color = TE_tot)) +
  geom_point(size = 5) +
  geom_text(mapping = aes(label = State), nudge_x = 0.5) +
  ggtitle("tSNE clustering of TE families (shared), colored by TE total") +
  theme_bw()


# ---------- SUPERFAMILY tSNE ----------

tSNE_MR_sup <- MR_sup %>%
  dplyr::select(ID:Ecoregion, ends_with("_mean")) %>%
  mutate(ID2 = row_number())

tSNE_MR_sup_meta <- tSNE_MR_sup %>%
  dplyr::select(ID:Ecoregion, ID2)

set.seed(142)
tSNE_sup_fit <- tSNE_MR_sup %>%
  dplyr::select(ends_with("_mean"), ID2) %>%
  column_to_rownames("ID2") %>%
  scale() %>%
  Rtsne()

tSNE_sup_df <- tSNE_sup_fit$Y %>%
  as.data.frame() %>%
  rename(tSNE1 = "V1", tSNE2 = "V2") %>%
  mutate(ID2 = row_number()) %>%
  inner_join(tSNE_MR_sup_meta, by = "ID2")

# Colored by ecoregion
tSNE_eco_sup <- tSNE_sup_df %>%
  ggplot(aes(x = tSNE1, y = tSNE2, color = Ecoregion)) +
  geom_point(size = 5) +
  scale_color_manual(values = c(
    "#54286f", "#a989b9", "#8b0000", "#9ab87a", "#304c00",
    "#cdb79e", "#85683c", "#1066bc", "#8abff5", "#BDDBF9")) +
  geom_text(mapping = aes(label = State), nudge_x = 0.5) +
  ggtitle("tSNE clustering of TE superfamilies, colored by ecoregion") +
  theme_bw()
tSNE_eco_sup
