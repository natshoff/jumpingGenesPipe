# TE Analysis Pipeline

# Script 3: Population Plots

# Description: This script creates plots from the master data file prepared
# in 1_TEpipe_MRfileprep.R



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
MR_map <- read.csv("Pop_meta.csv", header = T)
# Convert ecoregion to factor
MR_map$Ecoregion <- as.factor(MR_map$Ecoregion)



########
# MAPS #
########
# Download shapefile of US
us <- getData("GADM", country = "USA", level = 1)
# Remove Alaska and Hawaii
us_sub <- us %>% 
  subset(NAME_1 != "Alaska") %>% 
  subset(NAME_1 != "Hawaii")

# Download WorldClim data on 2.5 arc min scale
bio2.5 <- getData("worldclim", var = 'bio', res = 2.5)

# Bound bioclims to US
# Bio1
bio1 <- crop(bio2.5[[1]], bbox(us_sub))

# Bio12
bio12 <- crop(bio2.5[[12]], bbox(us_sub))

# Bio3
bio3 <- crop(bio2.5[[3]], bbox(us_sub))

# Bio15
bio15 <- crop(bio2.5[[15]], bbox(us_sub))

# Convert raster to dataframe
# Bio1 
# raster to points
bio1_point <- rasterToPoints(bio1)
# points to df
bio1_df <- data.frame(bio1_point)
# rescale data
bio1_df$bio1 <- (bio1_df$bio1)*0.1

# Bio12 
# raster to points
bio12_point <- rasterToPoints(bio12)
# points to df
bio12_df <- data.frame(bio12_point)

# Bio3 
# raster to points
bio3_point <- rasterToPoints(bio3)
# points to df
bio3_df <- data.frame(bio3_point)

# Bio15 
# raster to points
bio15_point <- rasterToPoints(bio15)
# points to df
bio15_df <- data.frame(bio15_point)

# Maps of important bioclim variables
# Bio1
b1_g <- ggplot(data=bio1_df, aes(x=x, y=y)) +
  # Plots bioclim data
  geom_raster(aes(fill=bio1)) +
  scale_fill_gradient(low = "#8ABFF5", high = "#8B0000") +
  # Plots US border over raster data
  borders("state", colour = "white", size = 0.5) +
  geom_point(data = MR_map, aes(x = Longitude, y = Latitude, color = Ecoregion), size = 2)  +
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
  ggtitle("Average Annual Temperature (C) - Bio1") +
  theme_bw() +
  theme(legend.position="none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.background = element_rect(fill = 'transparent'),
        legend.key = element_rect(fill = 'transparent'),
        panel.background = element_rect(fill = 'transparent'),
        plot.background = element_rect(fill = 'transparent', color=NA))
# Bio3
b3_g <- ggplot(data=bio3_df, aes(x=x, y=y)) +
  # Plots bioclim data
  geom_raster(aes(fill=bio3)) +
  scale_fill_gradient(low = "#D4B4C8", high = "#3D0029") +
  # Plots US border over raster data
  borders("state", colour = "white", size = 0.5) +
  geom_point(data = MR_map, aes(x = Longitude, y = Latitude, color = Ecoregion), size = 2) +
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
  ggtitle("Bio3 (Isothermality)") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())

# Bio12
b12_g <- ggplot(data=bio12_df, aes(x=x, y=y)) +
  # Plots bioclim data
  geom_raster(aes(fill=bio12)) +
  scale_fill_gradient2(low = "#CDB79E", mid = "#5D782E", high = "#3A4919", midpoint = 1800) +
  # Plots US border over raster data
  borders("state", colour = "white", size = 0.5) +
  geom_point(data = MR_map, aes(x = Longitude, y = Latitude, color = Ecoregion), size = 2) +
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
  ggtitle("Average Annual Precipitation (mm) - Bio12") +
  theme_bw() +
  theme(legend.position="none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.background = element_rect(fill = 'transparent'),
        legend.key = element_rect(fill = 'transparent'),
        panel.background = element_rect(fill = 'transparent'),
        plot.background = element_rect(fill = 'transparent', color=NA))

# Bio15
b15_g <- ggplot(data=bio15_df, aes(x=x, y=y)) +
  # Plots bioclim data
  geom_raster(aes(fill=bio15)) +
  scale_fill_gradient(low = "#EDE5E5", high = "#520000") +
  # Plots US border over raster data
  borders("state", colour = "white", size = 0.5) +
  geom_point(data = MR_map, aes(x = Longitude, y = Latitude, color = Ecoregion), size = 2) +
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
  ggtitle("Precipitation Seasonality - Bio15") +
  theme_bw() +
  theme(legend.position="none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.background = element_rect(fill = 'transparent'),
        legend.key = element_rect(fill = 'transparent'),
        panel.background = element_rect(fill = 'transparent'),
        plot.background = element_rect(fill = 'transparent', color=NA))
# Bio1 Map
b1_g
# save file
ggsave("bio1_map.png", height = 3.5, width = 5.86, dpi = 300, bg='transparent',
       units = "in")

# Bio12 Map
b12_g
# save file
ggsave("bio12_map.png", height = 3.5, width = 5.86, dpi = 300, bg='transparent',
       units = "in")

# Bio15 Map
b15_g
# save file
ggsave("bio15_map.png", height = 3.5, width = 5.86, dpi = 300, bg='transparent',
       units = "in")



##############################
# PCA of Population BioClims #
##############################

######################
# Correlation matrix #
######################

# Correlation matrix on all bioclims + alt
# Use spearman because of non-parametric data
MR_map %>% 
  dplyr::select(-c(Population:Longitude,Ecoregion)) %>% 
  chart.Correlation(histogram=F, pch=19, method = "spearman")

# Reduce bioclims to less correlated (<0.80) variables
MR_map %>% 
  dplyr::select(-c(Population:Longitude,Ecoregion)) %>%
  dplyr::select(c(bio1,bio2,bio3,bio12,bio15)) %>%
  chart.Correlation(histogram=F, pch=19, method = "spearman")


###################
# 0.5 arc sec PCA #
###################

# Reduced PCA with just a couple bioclim vars
red_recipe <-
  recipe(~., data = MR_map[c(1:5,6:8,17,20,26)]) %>%
  # converts categorical vars to dummy vars (we'll use this at some point)
  update_role(c(Population:Longitude,Ecoregion),new_role = "id") %>% 
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
  xlim(c(0, 7)) + 
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
  xlab("PC1 (68.5%)") +
  ylab("PC2 (19.1%)")

# Transparent version
juice(red_recipe) %>%
  ggplot(aes(PC1, PC2, label = State)) +
  geom_point(aes(fill = Ecoregion), shape = 21, size = 5) +
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
    "#BDDBF9"), 0.8)) +
  geom_text(check_overlap = TRUE, size = 4, hjust = "left", nudge_x = 0.07) +
  labs(color = NULL) + 
  geom_segment(inherit.aes = F, data = red_wider,
               aes(xend = PC1*5, yend = PC2*5), 
               x = 0, 
               y = 0, 
               arrow = arrow_style) + 
#  geom_text(inherit.aes = F, data = red_wider,
#            aes(x = PC1*5, y = PC2*5, label = terms), 
#            hjust = "outward",
#            nudge_x = 0.00,
#            vjust = 0,
#            nudge_y = 0.00,
#            size = 5, 
#            color = 'black') + 
  ggtitle("PCA of BioClim variables across TE populations") +
  xlab("PC1 (68.5%)") +
  ylab("PC2 (19.1%)") +
  theme_bw() +
  theme(legend.background = element_rect(fill = 'transparent'),
        legend.key = element_rect(fill = 'transparent'),
        panel.background = element_rect(fill = 'transparent'),
        plot.background = element_rect(fill = 'transparent', color=NA))
# save file
ggsave("TEpop_pca.png", height = 7.25, width = 7.5, dpi = 300, bg='transparent',
       units = "in")


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
  left_join(MR_map) %>%
  # join to PCs
  left_join(pop_PC) %>% 
  relocate(Population,State:Ecoregion,PC1:PC5)
# Convert ecoregion to factor again
TE_pop_sum$Ecoregion <- as.factor(TE_pop_sum$Ecoregion)
# Add ecoregion number column
TE_pop_sum <- TE_pop_sum %>% 
  mutate(Ecoregion_num2 = as.numeric(Ecoregion))


# Convert to long form for faceting
TE_pop_sum_long <- TE_pop_sum %>%
  dplyr::select(c(Population:Longitude,Ecoregion,
           LTR_pop_mean:TEtot_pop_sd,
           bio1,bio2,bio3,bio12,bio15,
           PC1:PC2)) %>% 
  pivot_longer(c(Latitude:Longitude,bio1:bio15,PC1:PC2),
               names_to = "EnvChar",
               values_to = "EnvVals")



#################
# Scatter plots #
#################

# All pops mean
ggplot() +
  geom_point(data = TE_pop_sum, mapping = aes(x = reorder(Population, Ecoregion_num2), 
                                              y = TEtot_pop_mean,
                                              fill = Ecoregion), 
             size = 5, shape = 21, color = "black", alpha = 0.8) +
  geom_errorbar(data = TE_pop_sum, mapping = aes(x = reorder(Population, Ecoregion_num2), 
                                                      y = TEtot_pop_mean, 
                                                      ymin=TEtot_pop_mean-TEtot_pop_sd, 
                                                      ymax=TEtot_pop_mean+TEtot_pop_sd,
                                                      color = Ecoregion), width = 0) +
  ggtitle("Mean population TE total read proportions") +
  scale_fill_manual(values = c(
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
  xlab("Population Number") +
  ylab("Total TE Proportion of Reads") +
  theme(axis.text.x = element_text(angle = 20),
        legend.background = element_rect(fill = 'transparent'),
        legend.key = element_rect(fill = 'transparent'),
        panel.background = element_rect(fill = 'transparent'),
        plot.background = element_rect(fill = 'transparent', color=NA))

# save
ggsave("TEprop_pop_mean.png", height = 8, width = 15, dpi = 300, bg='transparent')




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
  ggtitle("Mean population LTR read proportions across BiolClim variables") +
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

# Add in ag info and color by ag
pop_ag <- read.csv("TEpops_ag.csv", header = T) %>% 
  rename(Population = population)
# Join to TE_pop_sum_long
TE_pop_sum_long_ag <- TE_pop_sum_long %>% 
  left_join(pop_ag)
# TE total proportions
ggplot() +
  geom_point(data = TE_pop_sum_long_ag, mapping = aes(x = EnvVals, y = TEtot_pop_mean,
                                                   color = ag), size = 3) +
  geom_errorbar(data = TE_pop_sum_long_ag, mapping = aes(x = EnvVals, y = TEtot_pop_mean, 
                                                      ymin=TEtot_pop_mean-TEtot_pop_sd, 
                                                      ymax=TEtot_pop_mean+TEtot_pop_sd,
                                                      color = ag), width = 0) +
  geom_smooth(data = TE_pop_sum_long_ag, mapping = aes(x = EnvVals, y = TEtot_pop_mean,
                                                       color = ag), method = lm) +
  facet_wrap(~EnvChar, scales = "free") +
  ggtitle("Mean population TE total read proportions across BiolClim variables") +
  theme_bw() +
  xlab("Environmental Variable") +
  ylab("Total TE Proportion of Reads")



# Individual sequencing runs related to env vars
ggplot() +
  geom_point(data = TE_pop_sum, mapping = aes(x = bio1, y = TEtot_pop_mean,
                                              color = Ecoregion), size = 3) +
  geom_smooth(data = TE_pop_sum, mapping = aes(x = bio1, y = TEtot_pop_mean),
              method = lm) +
  theme_bw()



#############################
# Partial linear regression #
#############################

# POPULATION LEVEL #
# Test correlation between PC1, PC2, bioclim vars and TE total proportion
# Use Latitude and Longitude as partial predictors
pcor_pop_full <- TE_pop_sum %>% 
  dplyr::select(c(Population:Longitude,Ecoregion,
                  LTR_pop_mean:TEtot_pop_sd,
                  bio1,bio2,bio3,bio12,bio15,
                  PC1:PC2))
# reduce data to only predictors
pcor_pop_red <- pcor_pop_full %>% 
  dplyr::select(c(Latitude:Longitude,
                  TEtot_pop_mean,
            #      bio1:bio15,
                  PC1:PC2))
# Partial linear regression
pcor_pop <- pcor(pcor_pop_red, method = "pearson")
pcor_pop

# Population level with bioclims as single predictors
pcor_pop_red_bio <- pcor_pop_full %>% 
  dplyr::select(c(Latitude:Longitude,
                  TEtot_pop_mean,
                  bio1:bio15))
# Partial linear regression
pcor_pop_bio <- pcor(pcor_pop_red_bio, method = "pearson")
pcor_pop_bio


# INDIVIDUAL LEVEL #
pcor_ind_full <- MR_sup_add_pc %>% 
  dplyr::select(c(ID:Ecoregion,
                  Latitude:Longitude,
                  LTR_mean:TE_tot,
                  bio1,bio2,bio3,bio12,bio15,
                  PC1:PC2))
# reduce data to only predictors
pcor_ind_red <- pcor_ind_full %>% 
  dplyr::select(c(Latitude:Longitude,
                  TE_tot,
              #    bio1:bio15,
                  PC1:PC2))
# Partial linear regression
pcor_ind <- pcor(pcor_ind_red, method = "pearson")
pcor_ind
