#########################
# Set working directory #
#########################

# NOTE: change this depending on the directory you are using
setwd("~/Research/McKain Lab/CPING/TE/R/TE Pipeline")



###################
# Library loading #
###################

# Data manipulation library
library(tidyverse)
# Correlation matrix library
library(PerformanceAnalytics)
# Recipes for PCA
library(recipes)
# Patchworks
library(patchwork)
# Mapping libraries
library(rgdal)
library(raster)
# Partial linear regression
library(MASS)
library(ppcor)
library(ggpubr)


##############################
# Load in MR files from csvs #
##############################
# Population metadata
MR_bio <- read.csv("Pop_meta.csv", header = T)
# Convert ecoregion to factor
MR_bio$Ecoregion <- as.factor(MR_bio$Ecoregion)
# Add in envirem data too
MR_envm <- read.csv("Pop_envirem.csv", header = T)
# Convert ecoregion to factor
MR_envm$Ecoregion <- as.factor(MR_envm$Ecoregion)


##############################
# PCA of Population BioClims #
##############################

######################
# Correlation matrix #
######################

# Correlation matrix on all bioclims + alt
# Use spearman because of non-parametric data
MR_envm %>% 
  dplyr::select(-c(Population:Longitude,Ecoregion)) %>% 
  chart.Correlation(histogram=F, pch=19, method = "spearman")

# Reduce bioclims to less correlated (<0.80) variables
#MR_map %>% 
#  dplyr::select(-c(Population:Longitude,Ecoregion)) %>%
#  dplyr::select(c(bio1,bio2,bio3,bio12,bio15)) %>%
#  chart.Correlation(histogram=F, pch=19, method = "spearman")


###################
# 0.5 arc sec PCA #
###################

# Reduced PCA with just a couple bioclim vars
red_recipe <-
  recipe(~., data = MR_envm[c(1:24)]) %>%
  # converts categorical vars to dummy vars (we'll use this at some point)
  update_role(c(Population:Ecoregion),new_role = "id") %>% 
  # Omits NAs
  step_naomit(all_predictors()) %>%
  # Normalize for SD = 1, mean = 0
  step_normalize(all_predictors()) %>%
  # Creates PCA out of predictors
  step_pca(all_predictors(), id = "pca") %>% 
  prep()

# Execute PCA with recipe
red_pca <- 
  red_recipe %>% 
  tidy(id = "pca") 
red_pca

# Visualize variation explained by each component
red_recipe %>% 
  tidy(id = "pca", type = "variance") %>% 
  dplyr::filter(terms == "percent variance") %>% 
  ggplot(aes(x = component, y = value)) + 
  geom_col(fill = "#1B4332") + 
  xlim(c(0, 19)) + 
  ylab("% of total variance")

# Plot PCA loadings
# This helps identify the bioclim vars that contribute to each PC
red_pca %>%
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
red_wider <- red_pca %>% 
  tidyr::pivot_wider(names_from = component, id_cols = terms)

# Plot PC1 vs PC2
# define arrow style
arrow_style <- arrow(length=unit(0.05, "inches"),type = "closed")
# Plots points of PCA
# Specifies colors and all that and pipes in recipe information
red_plot12 <-
  juice(red_recipe) %>%
  ggplot(aes(PC1, PC2, label = State, color = Ecoregion)) +
  geom_point(aes(color = Ecoregion), size = 5) +
  geom_text(check_overlap = TRUE, size = 4, hjust = "left", nudge_x = 0.07) +
  labs(color = NULL) + 
  geom_segment(inherit.aes = F, data = red_wider,
               aes(xend = PC1*5, yend = PC2*5), 
               x = 0, 
               y = 0, 
               arrow = arrow_style) + 
  geom_text(inherit.aes = F, data = red_wider,
            aes(x = PC1*5, y = PC2*5, label = terms), 
            hjust = "outward",
            nudge_x = 0.00,
            vjust = 0,
            nudge_y = 0.00,
            size = 5, 
            color = 'black') +
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
    "#BDDBF9"))
# Whole plot
red_plot12 + 
  ggtitle("PCA of BioClim variables across TE populations") +
  theme_bw() +
  xlab("PC1 (53.8%)") +
  ylab("PC2 (18.5%)")

# Juice recipe to extract PC axes
pop_PC <- juice(red_recipe)



#############################
# Add in TE proportion data #
#############################

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

# Add PC axes to MR file
MR_sup_add_pc <- MR_sup_add %>% 
  left_join(pop_PC)


###########################################################
# TE Proportion correlation with envioronmental variables #
###########################################################

#########################
# LTR pop summary plots #
#########################

# Summarise Population Means
TE_pop_sum <- MR_sup_add_pc %>%
  # Select columns for grouping/summarizing
  group_by(Population) %>% 
  summarise(LTR_pop_mean = mean(LTR_mean),
            LTR_pop_sd = sd(LTR_mean),
            TEtot_pop_mean = mean(TE_tot),
            TEtot_pop_sd = sd(TE_tot)) %>% 
  # rejoin to pop environmental data
  left_join(MR_envm) %>%
  # join to PCs
  left_join(pop_PC) %>% 
  relocate(Population,State:Ecoregion,PC1:PC5)
# Convert ecoregion to factor again
TE_pop_sum$Ecoregion <- as.factor(TE_pop_sum$Ecoregion)


# Convert to long form for faceting
TE_pop_sum_long <- TE_pop_sum %>%
  dplyr::select(c(Population:Ecoregion,
                  PC1:PC5,
                  LTR_pop_mean:TEtot_pop_sd,
                  annualPET:topoWet)) %>% 
  pivot_longer(c(Latitude:Longitude,PC1:PC5,annualPET:topoWet),
               names_to = "EnvChar",
               values_to = "EnvVals")



#################
# Scatter plots #
#################

# LTR proportions
ggplot() +
  geom_point(data = TE_pop_sum_long, mapping = aes(x = EnvVals, y = LTR_pop_mean,
                                                   color = Ecoregion), size = 3) +
  geom_errorbar(data = TE_pop_sum_long, mapping = aes(x = EnvVals, y = LTR_pop_mean, 
                                                      ymin=LTR_pop_mean-LTR_pop_sd, 
                                                      ymax=LTR_pop_mean+LTR_pop_sd,
                                                      color = Ecoregion), width = 0) +
  geom_smooth(data = TE_pop_sum_long, mapping = aes(x = EnvVals, y = LTR_pop_mean),
              color = "black", method = lm) +
  facet_wrap(~EnvChar, scales = "free") +
  ggtitle("Mean population LTR read proportions across Envirem variables") +
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
  xlab("Environmental Variable") +
  ylab("LTR Proportion of Reads")

# TE total proportions
ggplot() +
  geom_point(data = TE_pop_sum_long, mapping = aes(x = EnvVals, y = TEtot_pop_mean,
                                                   color = Ecoregion), size = 3) +
  geom_errorbar(data = TE_pop_sum_long, mapping = aes(x = EnvVals, y = TEtot_pop_mean, 
                                                      ymin=TEtot_pop_mean-TEtot_pop_sd, 
                                                      ymax=TEtot_pop_mean+TEtot_pop_sd,
                                                      color = Ecoregion), width = 0) +
  geom_smooth(data = TE_pop_sum_long, mapping = aes(x = EnvVals, y = TEtot_pop_mean),
              color = "black", method = lm) +
  facet_wrap(~EnvChar, scales = "free") +
  ggtitle("Mean population TE total read proportions across BiolClim variables") +
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
  xlab("Environmental Variable") +
  ylab("Total TE Proportion of Reads")

#  Just TWI figure
ggplot() +
  geom_point(data = TE_pop_sum, mapping = aes(x = topoWet, y = TEtot_pop_mean,
                                                   color = Ecoregion), size = 3) +
  geom_errorbar(data = TE_pop_sum, mapping = aes(x = topoWet, y = TEtot_pop_mean, 
                                                      ymin=TEtot_pop_mean-TEtot_pop_sd, 
                                                      ymax=TEtot_pop_mean+TEtot_pop_sd,
                                                      color = Ecoregion), width = 0) +
  geom_smooth(data = TE_pop_sum, mapping = aes(x = topoWet, y = TEtot_pop_mean),
              color = "black", method = lm) +
  ggtitle("Mean population TE total read proportions across Topographic Wetness Index") +
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
  xlab("TWI") +
  ylab("Total TE Proportion of Reads")


