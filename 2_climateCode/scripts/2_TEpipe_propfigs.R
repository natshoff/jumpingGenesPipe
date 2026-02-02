# TE Analysis Pipeline

# Script 2: TE Proportion Plots

# Description: This script creates plots from the master data file prepared
# in 1_TEpipe_MRfileprep.R


#########################
# Set working directory #
#########################

# NOTE: change this depending on the directory you are using
setwd("~/Research/McKain Lab/CPING/TE/R/TE Pipeline/April 2024 Analysis")



###################
# Library loading #
###################

# Data manipulation library
library(tidyverse)
# Correlation matrix library
library(PerformanceAnalytics)
# Recipes for PCA
library(recipes)
# tSNE analysis
library(Rtsne)


##############################
# Load in MR files from csvs #
##############################
# Population metadata
Pop_meta <- read.csv("Pop_meta.csv", header = T)


# SUPERFAMILY
MR_sup <- read.csv("MR_sup.csv", header = T) 
  # NOTE: optional filtering of outlier (add pipe above)
 # filter(ID != "JG_CO_054_F")
# Convert ecoregion to factor
MR_sup$Ecoregion <- as.factor(MR_sup$Ecoregion)

# Add (1) LTR proprtion and (2) total TE content to master file
# Added Ecoregion_num for graphing assistance
MR_sup_add <- MR_sup %>%
  # add LTR and total column
  mutate(LTR_mean = Copia_mean + Gypsy_mean,
         TE_tot = rowSums(across(Copia_mean:unclassified_mean)),
         Ecoregion_num = Ecoregion)
# Convert Ecoregion_num to number from factor
MR_sup_add$Ecoregion_num <- as.numeric(MR_sup_add$Ecoregion_num)
# Reorder Seq_run levels
MR_sup_add$Seq_run <- factor(MR_sup_add$Seq_run,
                              c("10.08.21","01.24.22","03.23.22"))


# FAMILY
MR_fam <- read.csv("MR_fam.csv", header = T) %>% 
  # Add numeric ecoregion column for easier plotting later
  mutate(Ecoregion_num = Ecoregion)
# NOTE: optional filtering of outlier (add pipe above)
# filter(ID != "JG_CO_054_F")
# Convert one ecoregion to factor
MR_fam$Ecoregion <- as.factor(MR_fam$Ecoregion)
# Reorder Seq_run levels
MR_fam$Seq_run <- factor(MR_fam$Seq_run,
                             c("10.08.21","01.24.22","03.23.22"))




#####################
# TE Pipeline Stats #
#####################

################################
# DF Manipulation and creation #
################################
# We need four new working dataframes fro these plots
# (1) A pivot longer form of the MR_sup df, lengthened across pipe stage (TE_pipe_long)
# (2) A summary df of the mean read pairs of each stage (TE_pipe_sum)
# (3) A summary df of the mean read pairs of each stage, grouped by sequence run (TE_pipe_sum_seq)
# (4) A pivot longer form of that summary file, lengthened across pipe stage (TE_sum_seq_long)

# Make long version of MR_sup file
TE_pipe_long <- MR_sup_add %>% 
  pivot_longer(Trimmed:Deduplicated, names_to = "Stage", values_to = "Read_pairs")
# Reorder Stage and Seq_run levels for plotting in order
TE_pipe_long$Stage <- factor(TE_pipe_long$Stage,
                             c("Trimmed","Organella","Fungal","Bacterial","Deduplicated"))

# Make summary file of pipeline stage means
TE_pipe_sum <- MR_sup_add %>% 
  summarise(Trimmed_sd = sd(Trimmed),
            Organella_sd = sd(Organella),
            Fungal_sd = sd(Fungal),
            Bacterial_sd = sd(Bacterial),
            Deduplicated_sd = sd(Deduplicated),
            Trimmed = mean(Trimmed),
            Organella = mean(Organella),
            Fungal = mean(Fungal),
            Bacterial = mean(Bacterial),
            Deduplicated = mean(Deduplicated))

# Make summary file of pipeline stage means
TE_pipe_sum_seq <- MR_sup_add %>% 
  group_by(Seq_run) %>% 
  summarise(Trimmed = mean(Trimmed),
            Organella = mean(Organella),
            Fungal = mean(Fungal),
            Bacterial = mean(Bacterial),
            Deduplicated = mean(Deduplicated))

# Make long version of stage summary file
TE_sum_seq_long <- TE_pipe_sum_seq %>% 
  pivot_longer(Trimmed:Deduplicated, names_to = "Stage", values_to = "Mean_rp")


#################
# Scatter plots #
#################

# Scatter plot with total TE content and deduplicated reads
ggplot(MR_sup_add) +
  geom_point(mapping = aes(x = Deduplicated, y = TE_tot,
                           color = Seq_run)) +
  geom_smooth(mapping = aes(x = Deduplicated, y = TE_tot,
                            color = Seq_run),
              method = lm) +
  geom_smooth(mapping = aes(x = Deduplicated, y = TE_tot),
              method = lm, color = "black")
# Scatter plot with LTR content and deduplicated reads
ggplot(MR_sup_add) +
  geom_point(mapping = aes(x = Deduplicated, y = LTR_mean,
                           color = Seq_run)) +
  geom_smooth(mapping = aes(x = Deduplicated, y = LTR_mean,
                            color = Seq_run),
              method = lm) +
  geom_smooth(mapping = aes(x = Deduplicated, y = LTR_mean),
              method = lm, color = "black")

#################
# Density plots #
#################

# Density plot of TE_total reads 
ggplot(MR_sup_add) +
  geom_density(mapping = aes(x = TE_tot, fill = Seq_run),
               alpha = 0.5)
# Density plot of deduplicated across sequencing runs 
ggplot(MR_sup_add) +
  geom_density(mapping = aes(x = Deduplicated, fill = Seq_run),
               alpha = 0.5)

#############
# Box plots #
#############

# Box plots across filtering levels all
ggplot(TE_pipe_long) +
  geom_boxplot(mapping = aes(x = Stage, y = Read_pairs)) +
  geom_point(data = TE_sum_seq_long, 
             mapping = aes(x = Stage, y = Mean_rp, color = Seq_run),
             shape=4, size=3,
             position = position_dodge2(width = 0.75,
                                        preserve = "single")) +
  geom_line(data = TE_sum_seq_long, 
            mapping = aes(x = Stage, y = Mean_rp, color = Seq_run, group = Seq_run),
            size = 1,
            position = position_dodge2(width = 0.75,
                                       preserve = "single"))

# Box plots across filtering levels (subset by sequencing run)
ggplot(TE_pipe_long) +
  geom_boxplot(mapping = aes(x = Stage, y = Read_pairs, color = Seq_run)) +
  geom_point(data = TE_sum_seq_long, 
             mapping = aes(x = Stage, y = Mean_rp, color = Seq_run),
             shape=4, size=3,
             position = position_dodge2(width = 0.75,
                                        preserve = "single")) +
  geom_line(data = TE_sum_seq_long, 
            mapping = aes(x = Stage, y = Mean_rp, color = Seq_run, group = Seq_run),
            size = 1,
            position = position_dodge2(width = 0.75,
                                       preserve = "single"))

# Violin plot with total TE proportion for each sequencing run
ggplot(TE_pipe_long) +
  geom_violin(mapping = aes(x = Seq_run, y = TE_tot, color = Seq_run)) +
  geom_jitter(mapping = aes(x = Seq_run, y = TE_tot, color = Seq_run),
              alpha = 0.5)

######################################
# Inter / Intra population variation #
######################################


###########
# BOXPLOT #
###########

# LTR proportion
ggplot() +
  geom_boxplot(data = MR_sup_add, mapping = aes(x = reorder(Population, Ecoregion_num), 
                                                y = LTR_mean,
                                                color = Ecoregion),
               outlier.shape = 4) +
  geom_jitter(data = MR_sup_add, mapping = aes(x = reorder(Population, Ecoregion_num), 
                                               y = LTR_mean,
                                               color = Ecoregion),
              width = 0.35, size = 3, alpha = 0.5) +
  scale_color_manual(values = c(
    "#54286f", #8.1
    "#a989b9", #8.2
    "#8b0000", #8.3
    "#9ab87a", #8.4
    "#304c00", #8.5
    "#cdb79e", #9.2
    "#85683c", #9.4
    "#1066bc", #10.2
    "#8abff5", #11.1
    "#BDDBF9")) +
  theme_bw() +
  ggtitle("LTR proportions across populations") +
  theme(axis.text.x = element_text(angle = 20))

# Total TE proportion
ggplot() +
  geom_boxplot(data = MR_sup_add, mapping = aes(x = reorder(Population, Ecoregion_num), 
                                                y = TE_tot,
                                                color = Ecoregion),
               outlier.shape = 4) +
  geom_jitter(data = MR_sup_add, mapping = aes(x = reorder(Population, Ecoregion_num), 
                                               y = TE_tot,
                                               color = Ecoregion),
              width = 0.25, size = 3, alpha = 0.5) +
  scale_color_manual(values = c(
    "#54286f", #8.1
    "#a989b9", #8.2
    "#8b0000", #8.3
    "#9ab87a", #8.4
    "#304c00", #8.5
    "#cdb79e", #9.2
    "#85683c", #9.4
    "#1066bc", #10.2
    "#8abff5", #11.1
    "#BDDBF9")) +
  theme_bw() +
  ggtitle("Total TE proportions across populations") +
  xlab("Population Number") +
  ylab("Total TE Proportion") +
  theme(axis.text.x = element_text(angle = 20))

# Total TE proportion but transparent background
ggplot() +
  geom_boxplot(data = MR_sup_add, mapping = aes(x = reorder(Population, Ecoregion_num), 
                                                y = TE_tot,
                                                color = Ecoregion),
               outlier.shape = 4, fill = "white", alpha = 0.5, lwd = 0.75) +
  scale_color_manual(values = c(
    "#54286f", #8.1
    "#a989b9", #8.2
    "#8b0000", #8.3
    "#9ab87a", #8.4
    "#304c00", #8.5
    "#cdb79e", #9.2
    "#85683c", #9.4
    "#1066bc", #10.2
    "#8abff5", #11.1
    "#BDDBF9")) +
  geom_jitter(data = MR_sup_add, mapping = aes(x = reorder(Population, Ecoregion_num), 
                                               y = TE_tot,
                                               fill = Ecoregion),
              width = 0.25, size = 3, color = "black", shape = 21) +
  scale_fill_manual(values = alpha(c(
    "#54286f", #8.1
    "#a989b9", #8.2
    "#8b0000", #8.3
    "#9ab87a", #8.4
    "#304c00", #8.5
    "#cdb79e", #9.2
    "#85683c", #9.4
    "#1066bc", #10.2
    "#8abff5", #11.1
    "#BDDBF9"), 0.6)) +
  theme_bw() +
  ggtitle("Total TE proportions across populations") +
  xlab("Population Number") +
  ylab("Total TE Proportion") +
  theme(axis.text.x = element_text(angle = 20),
        legend.background = element_rect(fill = 'transparent'),
        legend.key = element_rect(fill = 'transparent'),
        panel.background = element_rect(fill = 'transparent'),
        plot.background = element_rect(fill = 'transparent', color=NA))

ggsave("TEprop_pops.png", height = 8, width = 15, dpi = 300, bg='transparent')

# Just population 045 and 229
# Filter just pop 229 and 045
MR_sup_add_reduced <- MR_sup_add %>% 
  filter(Population == "p_045"|"p_229")

ggplot() +
  geom_boxplot(data = MR_sup_add, mapping = aes(x = reorder(Population, Ecoregion_num), 
                                                y = TE_tot,
                                                color = Ecoregion),
               outlier.shape = 4, fill = "white", alpha = 0.5, lwd = 0.75) +
  scale_color_manual(values = c(
    "#54286f", #8.1
    "#a989b9", #8.2
    "#8b0000", #8.3
    "#9ab87a", #8.4
    "#304c00", #8.5
    "#cdb79e", #9.2
    "#85683c", #9.4
    "#1066bc", #10.2
    "#8abff5", #11.1
    "#BDDBF9")) +
  geom_jitter(data = MR_sup_add, mapping = aes(x = reorder(Population, Ecoregion_num), 
                                               y = TE_tot,
                                               fill = Ecoregion),
              width = 0.25, size = 3, color = "black", shape = 21) +
  scale_fill_manual(values = alpha(c(
    "#54286f", #8.1
    "#a989b9", #8.2
    "#8b0000", #8.3
    "#9ab87a", #8.4
    "#304c00", #8.5
    "#cdb79e", #9.2
    "#85683c", #9.4
    "#1066bc", #10.2
    "#8abff5", #11.1
    "#BDDBF9"), 0.6)) +
  theme_bw() +
  ggtitle("Total TE proportions across populations") +
  xlab("Population Number") +
  ylab("Total TE Proportion") +
  theme(axis.text.x = element_text(angle = 20),
        legend.background = element_rect(fill = 'transparent'),
        legend.key = element_rect(fill = 'transparent'),
        panel.background = element_rect(fill = 'transparent'),
        plot.background = element_rect(fill = 'transparent', color=NA))

ggsave("TEprop_pops.png", height = 8, width = 15, dpi = 300, bg='transparent')



# Total TE proportion, colored by sequencing run
ggplot() +
  geom_boxplot(data = MR_sup_add, mapping = aes(x = reorder(Population, Ecoregion_num), 
                                                y = TE_tot,
                                                color = Seq_run),
               outlier.shape = 4) +
  geom_jitter(data = MR_sup_add, mapping = aes(x = reorder(Population, Ecoregion_num), 
                                               y = TE_tot,
                                               color = Seq_run),
              width = 0.25, size = 3, alpha = 0.5) +
  theme_bw() +
  ggtitle("Total TE proportions across populations") +
  theme(axis.text.x = element_text(angle = 20))

############
# Variance #
############

# Calculate average proportion of each superfamily across all individuals
Sup_mean <- MR_sup_add %>% 
  summarise(Copia = mean(Copia_mean),
            Gypsy = mean(Gypsy_mean),
            Harbinger = mean(Harbinger_mean),
            Helitron = mean(Helitron_mean),
            L1 = mean(L1_mean),
            MSAT = mean(MSAT_mean),
            MuDR = mean(MuDR_mean),
            RTE = mean(RTE_mean),
            SINE2.tRNA = mean(SINE2.tRNA_mean),
            hAT = mean(hAT_mean),
            unclassified = mean(unclassified_mean))
Sup_sd <- MR_sup_add %>% 
  summarise(Copia = sd(Copia_mean),
            Gypsy = sd(Gypsy_mean),
            Harbinger = sd(Harbinger_mean),
            Helitron = sd(Helitron_mean),
            L1 = sd(L1_mean),
            MSAT = sd(MSAT_mean),
            MuDR = sd(MuDR_mean),
            RTE = sd(RTE_mean),
            SINE2.tRNA = sd(SINE2.tRNA_mean),
            hAT = sd(hAT_mean),
            unclassified = sd(unclassified_mean))
write.csv(Sup_sd, "Superfamily_sd.csv", row.names = F)            


# Individuals
Var_ind <- MR_sup_add %>% 
  summarise(TE_tot_mean = mean(TE_tot),
            TE_tot_sd = sd(TE_tot),
            TE_tot_min = min(TE_tot),
            TE_tot_max = max(TE_tot),
            TE_tot_range = TE_tot_max - TE_tot_min,
            TE_tot_var = var(TE_tot))

# Within populations
Var_intra <- MR_sup_add %>% 
  group_by(Population) %>% 
  summarise(TE_tot_mean = mean(TE_tot),
            TE_tot_sd = sd(TE_tot),
            TE_tot_min = min(TE_tot),
            TE_tot_max = max(TE_tot),
            TE_tot_range = TE_tot_max - TE_tot_min,
            TE_tot_var = var(TE_tot),
            Sample_per = n())
Var_mean_intra <- Var_intra %>% 
  summarise(TE_tot_mean_var = mean(TE_tot_var, na.rm = T),
            TE_tot_mean_range = mean(TE_tot_range, na.rm = T))

# Between populations
Var_inter <- Var_intra %>% 
  summarise(TE_pop_mean = mean(TE_tot_mean),
            TE_pop_sd = sd(TE_tot_mean),
            TE_pop_min = min(TE_tot_mean),
            TE_pop_max = max(TE_tot_mean),
            TE_pop_range = TE_pop_max - TE_pop_min,
            TE_pop_var = var(TE_tot_mean, na.rm = T))

# CSVs of output data
write.csv(Var_ind, "AI_individual_var.csv", row.names = F)
write.csv(Var_intra, "AI_intrapop_var.csv", row.names = F)
write.csv(Var_inter, "AI_interpop_var.csv", row.names = F)

# Intrapopulation metrics
# Calculate variance within populations
Pop_var <- MR_sup_add %>%
  group_by(Population) %>% 
  summarise(TEtot_var = var(TE_tot),
            LTR_var = var(LTR_mean))
# Calculate average intrapopulation variance
Mean_intra_var <- Pop_var %>% 
  summarise(TEtot_var_mean = mean(TEtot_var, na.rm = T),
            TEtot_var_sd = sd(TEtot_var, na.rm = T),
            LTR_var_mean = mean(LTR_var, na.rm = T),
            LTR_var_sd = sd(LTR_var, na.rm = T))

# Interpopulation metrics
# Variance across all samples
Ind_var <- MR_sup_add %>% 
  summarise(TEtot_var = var(TE_tot),
            LTR_var = var(LTR_mean))
# Variance between population means
# Calculate population averages
Pop_means <- MR_sup_add %>% 
  group_by(Population) %>% 
  summarise(TEtot_pop_mean = mean(TE_tot),
            TEtot_pop_min = min(TE_tot),
            TEtot_pop_max = max(TE_tot),
            LTR_pop_mean = mean(LTR_mean),
            LTR_pop_min = min(LTR_mean),
            LTR_pop_max = max(LTR_mean)) %>% 
  mutate(TEtot_range = TEtot_pop_max - TEtot_pop_min,
         LTR_pop_range = LTR_pop_max - LTR_pop_min)
# Calculate interpopulation variance
Mean_inter_var <- Pop_means %>% 
  summarise(TEtot_var = var(TEtot_pop_mean),
            LTR_var = var(LTR_pop_mean))


##################################
# PCA of Superfamily differences #
##################################

######################
# Correlation matrix #
######################

# Correlation matrix on all superfamily proportions
# Use spearman because of non-parametric data
MR_sup %>% 
  dplyr::select(c(Copia_mean:unclassified_mean)) %>% 
  chart.Correlation(histogram=T, pch=19, method = "spearman")

###################
# Superfamily PCA #
###################

# PCA of all superfamily proportions
sup_recipe <-
  # Subset data to only include individual metadata(1:6) and TE proportion means(37:48)
  recipe(~., data = MR_sup[c(1:6,37:48)]) %>%
  # converts categorical vars to dummy vars (we'll use this at some point)
  update_role(c(ID:Ecoregion),new_role = "id") %>% 
  # Omits NAs
  step_naomit(all_predictors()) %>%
  # Normalize for SD = 1, mean = 0
  step_normalize(all_predictors()) %>%
  # Creates PCA out of predictors
  step_pca(all_predictors(), id = "pca") %>% 
  prep()

# Execute PCA with recipe
sup_pca <- 
  sup_recipe %>% 
  tidy(id = "pca") 
sup_pca

# Visualize variation explained by each component
# NOTE: run the first three lines of this to determine the % variance explained
# by each PC
sup_recipe %>% 
  tidy(id = "pca", type = "variance") %>% 
  dplyr::filter(terms == "percent variance") %>% 
  ggplot(aes(x = component, y = value)) + 
  geom_col(fill = "#1B4332") + 
  xlim(c(0, 15)) + 
  ylab("% of total variance")

# Plot PCA loadings
# This helps identify the TE types that contribute to each PC
sup_pca %>%
  mutate(terms = tidytext::reorder_within(terms, 
                                          abs(value), 
                                          component)) %>%
  ggplot(aes(abs(value), terms, fill = value > 0)) +
  geom_col() +
  facet_wrap(~component, scales = "free_y") +
  tidytext::scale_y_reordered() +
  scale_fill_manual(values = c("#52B788", "#1B4332")) +
  labs(
    x = "Absolute value of contribution",
    y = NULL, fill = "Positive?"
  ) 

# Widen PCA - converts pca for better graphing
sup_wider <- sup_pca %>% 
  tidyr::pivot_wider(names_from = component, id_cols = terms)

# Plot PC1 vs PC2
# define arrow style
arrow_style <- arrow(length=unit(0.05, "inches"),type = "closed")
# Plots points of PCA
# Specifies colors and all that and pipes in recipe information
sup_plot12 <-
  juice(sup_recipe) %>%
  ggplot(aes(PC1, PC2, label = Population, color = Ecoregion)) +
  geom_point(aes(color = Ecoregion), size = 5) +
  geom_text(check_overlap = TRUE, size = 4, hjust = "left", nudge_x = 0.07) +
  labs(color = NULL) + 
  geom_segment(inherit.aes = F, data = sup_wider,
               aes(xend = PC1*5, yend = PC2*5), 
               x = 0, 
               y = 0, 
               arrow = arrow_style) + 
  geom_text(inherit.aes = F, data = sup_wider,
            aes(x = PC1*5, y = PC2*5, label = terms), 
            hjust = "outward",
            nudge_x = 0.00,
            vjust = 0,
            nudge_y = 0.00,
            size = 5, 
            color = 'black')
# Whole plot
sup_plot12 + 
  ggtitle("PCA of TE Proportions") +
  theme_bw() +
  xlab("PC1 (42.6%)") +
  ylab("PC2 (23.5%)")


#############################
# PCA of Family differences #
#############################

##############
# Family PCA #
##############

# PCA with all families
fam_recipe <-
  # Subset data to only include individual metadata(1:6) and TE proportion means(37:1633)
  recipe(~., data = MR_fam[c(1:6,37:1633)]) %>%
  # converts categorical vars to dummy vars (we'll use this at some point)
  update_role(c(ID:Ecoregion),new_role = "id") %>% 
  # Omits NAs
  step_naomit(all_predictors()) %>%
  # Normalize for SD = 1, mean = 0
  step_normalize(all_predictors()) %>%
  # Creates PCA out of predictors
  step_pca(all_predictors(), id = "pca") %>% 
  prep()

# Execute PCA with recipe
fam_pca <- 
  fam_recipe %>% 
  tidy(id = "pca") 
fam_pca

# Visualize variation explained by each component
# NOTE: run the first three lines of this to determine the % variance explained
# by each PC
fam_recipe %>% 
  tidy(id = "pca", type = "variance") %>% 
  dplyr::filter(terms == "percent variance") %>% 
  ggplot(aes(x = component, y = value)) + 
  geom_col(fill = "#1B4332") + 
  xlim(c(0, 350)) + 
  ylab("% of total variance")

# Widen PCA - converts pca for better graphing
fam_wider <- fam_pca %>% 
  tidyr::pivot_wider(names_from = component, id_cols = terms)

# Plot PC1 vs PC2
# define arrow style
arrow_style <- arrow(length=unit(0.05, "inches"),type = "closed")
# Plots points of PCA
# Specifies colors and all that and pipes in recipe information
fam_plot12 <-
  juice(fam_recipe) %>%
  ggplot(aes(PC1, PC2, label = Population, color = Ecoregion)) +
  geom_point(aes(color = Ecoregion), size = 5) +
  geom_text(check_overlap = TRUE, size = 4, hjust = "left", nudge_x = 0.07) +
  labs(color = NULL) + 
  geom_segment(inherit.aes = F, data = fam_wider,
               aes(xend = PC1*500, yend = PC2*500), 
               x = 0, 
               y = 0, 
               arrow = arrow_style) + 
  geom_text(inherit.aes = F, data = fam_wider,
            aes(x = PC1*500, y = PC2*500, label = terms), 
            hjust = "outward",
            nudge_x = 0.00,
            vjust = 0,
            nudge_y = 0.00,
            size = 5, 
            color = 'black')
# Whole plot
fam_plot12 + 
  ggtitle("PCA of Family TE Proportions") +
  theme_bw() +
  xlab("PC1 (9.83%)") +
  ylab("PC2 (6.24%)")

# Zoom plot
fam_plot12 + 
  ggtitle("PCA of Family TE Proportions") +
  theme_bw() +
  xlab("PC1 (9.83%)") +
  ylab("PC2 (6.24%)") +
  coord_cartesian(ylim = c(-30, -10), xlim = c(-40,-5))



##############
# tSNE plots #
##############

# Family
# Add ID column for joining metadata after tSNE
tSNE_MR_fam <- MR_fam %>%
  # Select only metadata and mean proportions
  dplyr::select(c(1:1633)) %>% 
  mutate(ID2 = row_number())

# Reduced family dataset test
test_fam <- tSNE_MR_fam

# Replace all 0s with NAs
test_fam[test_fam == 0] <- NA

# Remove columns that contain NAs
test_rep <- test_fam[, colSums(is.na(test_fam)) == 0]
# this dataset now contains only families shared by all individuals (406)

# Reduced tSNE with 406 shared TE families
# Store metadata
tSNE_MR_fam_meta <- tSNE_MR_fam %>%
  dplyr::select(ID:Deduplicated,ID2)

# Perform tSNE
set.seed(142)
# run tSNE
tSNE_rep_fit <- test_rep %>%
  # exclude bioclims and lat/long as predictors
  dplyr::select(-c(ID:Deduplicated)) %>% 
  # only assess numeric columns
  dplyr::select(where(is.numeric)) %>%
  column_to_rownames("ID2") %>% 
  # standardize data
  scale() %>% 
  # perform tSNE
  Rtsne()

# Extract tSNE components
tSNE_rep_df <- tSNE_rep_fit$Y %>% 
  as.data.frame() %>%
  rename(tSNE1="V1",
         tSNE2="V2") %>%
  mutate(ID2=row_number())

# Join tSNE components to metadata
tSNE_rep_df <- tSNE_rep_df %>%
  inner_join(tSNE_MR_fam_meta, by="ID2")

# Colored by ecoregion
tSNE_eco_rep <- tSNE_rep_df %>%
  ggplot(aes(x = tSNE1, 
             y = tSNE2,
             color = Ecoregion)) +
  geom_point(size = 5) +
  scale_color_manual(values = c(
    "#54286f", #8.1
    "#a989b9", #8.2
    "#8b0000", #8.3
    "#9ab87a", #8.4
    "#304c00", #8.5
    "#cdb79e", #9.2
    "#85683c", #9.4
    "#1066bc", #10.2
    "#8abff5", #11.1
    "#BDDBF9")) +
  geom_text(tSNE_rep_df, mapping = aes(label = State), nudge_x = 0.5) +
  ggtitle("tSNE clustering of TE families (406), colored by level II ecoregion") +
  theme_bw()
tSNE_eco_rep
# transparent tSNE ecoregion plot
ggplot() +
  geom_point(data = tSNE_rep_df, mapping = aes(x = tSNE1,
                                               y = tSNE2,
                                               fill = Ecoregion), 
             size = 5, shape = 21, color = "black") +
  scale_fill_manual(values = alpha(c(
    "#54286f", #8.1
    "#a989b9", #8.2
    "#8b0000", #8.3
    "#9ab87a", #8.4
    "#304c00", #8.5
    "#cdb79e", #9.2
    "#85683c", #9.4
    "#1066bc", #10.2
    "#8abff5", #11.1
    "#BDDBF9"),0.8)) +
  geom_text(tSNE_rep_df, mapping = aes(x = tSNE1, y = tSNE2, 
                                       label = State, color = Ecoregion), 
            nudge_x = 0.5, show.legend = F) +
  scale_color_manual(values = c(
    "#54286f", #8.1
    "#a989b9", #8.2
    "#8b0000", #8.3
    "#9ab87a", #8.4
    "#304c00", #8.5
    "#cdb79e", #9.2
    "#85683c", #9.4
    "#1066bc", #10.2
    "#8abff5", #11.1
    "#BDDBF9")) +
  ggtitle("tSNE clustering of TE families (406), colored by level II ecoregion") +
  theme_bw() +
  theme(legend.background = element_rect(fill = 'transparent'),
        legend.key = element_rect(fill = 'transparent'),
        panel.background = element_rect(fill = 'transparent'),
        plot.background = element_rect(fill = 'transparent', color=NA))
# save file
ggsave("TEfam_tSNE.png", height = 8, width = 15, dpi = 300, bg='transparent')




# Colored by sequencing run
tSNE_seq_rep <- tSNE_rep_df %>%
  ggplot(aes(x = tSNE1, 
             y = tSNE2,
             color = Seq_run)) +
  geom_point(size = 5) +
  geom_text(tSNE_rep_df, mapping = aes(label = State), nudge_x = 0.5) +
  ggtitle("tSNE clustering of TE families (406), colored by sequencing run") +
  theme_bw()
tSNE_seq_rep

# Colored by population
tSNE_pop_rep <- tSNE_rep_df %>%
  ggplot(aes(x = tSNE1, 
             y = tSNE2,
             color = Population)) +
  geom_point(size = 5) +
  geom_text(tSNE_rep_df, mapping = aes(label = State), nudge_x = 0.5) +
  ggtitle("tSNE clustering of TE families (406), colored by population") +
  theme_bw()
tSNE_pop_rep

# Merge TE family tSNE output with superfamily proportions to assess potential
# drivers of clustering
tSNE_rep_df_merge <- tSNE_rep_df %>% 
  left_join(MR_sup_add)
# Log transform TE_total to reduce skew for better visualization
tSNE_rep_df_merge <- tSNE_rep_df_merge %>% 
  mutate(TE_tot_log = log(TE_tot))

# colored by log transformed TE total
tSNE_sup_rep <- tSNE_rep_df_merge %>%
  ggplot(aes(x = tSNE1, 
             y = tSNE2,
             color = TE_tot_log)) +
  geom_point(size = 5) +
  geom_text(tSNE_rep_df_merge, mapping = aes(label = State), nudge_x = 0.5) +
  ggtitle("tSNE clustering of TE families (406), colored by TE total") +
  theme_bw()
tSNE_sup_rep

# Filter few samples with reallly high TE proportions to assess
# Outlier test based on IQR criterion
boxplot.stats(tSNE_rep_df_merge$TE_tot)$out
# Identifies > 0.733 as outliers
# Plot and color by TE total
tSNE_rep_df_merge %>% 
 # filter(TE_tot < 0.733) %>% 
  ggplot(aes(x = tSNE1, 
             y = tSNE2,
             color = TE_tot)) +
  geom_point(size = 5) +
  geom_text(mapping = aes(label = State), nudge_x = 0.5) +
  ggtitle("tSNE clustering of TE families (406), colored by TE total") +
  theme_bw()


# Superfamily
# Add ID column for joining metadata after tSNE
tSNE_MR_sup <- MR_sup %>%
  # Select only metadata and mean proportions
  dplyr::select(c(1:48)) %>% 
  mutate(ID2 = row_number())

# Store metadata
tSNE_MR_sup_meta <- tSNE_MR_fam %>%
  dplyr::select(ID:Deduplicated,ID2)

# Perform tSNE
set.seed(142)
# run tSNE
tSNE_sup_fit <- tSNE_MR_sup %>%
  # exclude bioclims and lat/long as predictors
  dplyr::select(-c(ID:Deduplicated)) %>% 
  # only assess numeric columns
  dplyr::select(where(is.numeric)) %>%
  column_to_rownames("ID2") %>% 
  # standardize data
  scale() %>% 
  # perform tSNE
  Rtsne()

# Extract tSNE components
tSNE_sup_df <- tSNE_sup_fit$Y %>% 
  as.data.frame() %>%
  rename(tSNE1="V1",
         tSNE2="V2") %>%
  mutate(ID2=row_number())

# Join tSNE components to metadata
tSNE_sup_df <- tSNE_sup_df %>%
  inner_join(tSNE_MR_sup_meta, by="ID2")

# Colored by ecoregion
tSNE_eco_sup <- tSNE_sup_df %>%
  ggplot(aes(x = tSNE1, 
             y = tSNE2,
             color = Ecoregion)) +
  geom_point(size = 5) +
  scale_color_manual(values = c(
    "#54286f", #8.1
    "#a989b9", #8.2
    "#8b0000", #8.3
    "#9ab87a", #8.4
    "#304c00", #8.5
    "#cdb79e", #9.2
    "#85683c", #9.4
    "#1066bc", #10.2
    "#8abff5", #11.1
    "#BDDBF9")) +
  geom_text(tSNE_sup_df, mapping = aes(label = State), nudge_x = 0.5) +
  ggtitle("tSNE clustering of TE superfamilies, colored by level II ecoregion") +
  theme_bw()
tSNE_eco_sup

# Colored by sequencing run
tSNE_seq_sup <- tSNE_sup_df %>%
  ggplot(aes(x = tSNE1, 
             y = tSNE2,
             color = Seq_run)) +
  geom_point(size = 5) +
  geom_text(tSNE_sup_df, mapping = aes(label = State), nudge_x = 0.5) +
  ggtitle("tSNE clustering of TE families, colored by sequencing run") +
  theme_bw()
tSNE_seq_sup



