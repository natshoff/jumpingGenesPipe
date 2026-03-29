# TE Analysis Pipeline

# Script 1: File Preparation

# Description: This script is used to prep the TE proportion and metadata
# files for downstream analyses.


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



#############################
# Data loading and cleaning #
#############################


######################################
# TE SUPERFAMILY mean/sd proportions #
######################################

# 306 samples from April 2024 analysis
TEsup_avg_raw <- read.delim("Raw Data/apr21_allsamples_TE_annotation_results_average.txt", 
                             header = T, sep = "\t")
TEsup_sd_raw <- read.delim("Raw Data/apr21_allsamples_TE_annotation_results_stdev.txt", 
                            header = T, sep = "\t")

# We have several issues here
# 1) The ID column has extra information "4_Transposome/"
# 2) There is a non-TE proportion column to be removed (X) that tranposome just spits out
# 3) The column headers are shifted over to the right one space

# Issue 1 and 2 (remove path + extra column)
# Average
TEsup_avg_cor <- TEsup_avg_raw %>%
  # Remove path information from ID column
  mutate(Species = str_remove_all(Species, "4_Transposome/Sorghum_halepense-")) %>%
  # Drop useless column that contains non-TE proportions
  # Drop unlcassified column because we'll shift it over
  dplyr::select(-c(X, unclassified))
# Standard deviation
TEsup_sd_cor <- TEsup_sd_raw %>%
  # Remove path information from ID column
  mutate(Species = str_remove_all(Species, "4_Transposome/Sorghum_halepense-")) %>%
  # Drop useless column that contains non-TE proportions
  # Drop unlcassified column because we'll shift it over
  dplyr::select(-c(X, unclassified))

# Issue 3 (fix column headers)
# Extract raw column names
TEsup_head_raw <- colnames(TEsup_avg_raw)
# Remove blank column names (X and X.)
TEsup_head_correct <- TEsup_head_raw[-c(2:3)]
# Rename columns with correct headers
# Mean
names(TEsup_avg_cor) <- c(TEsup_head_correct) %>% 
  paste0("_mean")
# Rename Species_mean column header as ID
TEsup_avg_cor <- TEsup_avg_cor %>% 
  as_tibble() %>% 
  dplyr::rename(ID = Species_mean)
# Standard dev
names(TEsup_sd_cor) <- c(TEsup_head_correct) %>% 
  paste0("_sd")
# Rename Species_sd column header as ID
TEsup_sd_cor <- TEsup_sd_cor %>% 
  as_tibble() %>% 
  dplyr::rename(ID = Species_sd)

# Join all TE superfamily mean/sd data
TE_sup <- TEsup_avg_cor %>% 
  left_join(TEsup_sd_cor) %>% 
  # Split ID column into pop and individual
  separate(ID, c("Species","Collection","Population","Individual"),remove = F) %>% 
  dplyr::select(-c(Species, Collection))
# Add "p_" to all population columns to preserve 0 in front of some pop IDs
# when saving to csv
TE_sup$Population <- sub("^","p_",TE_sup$Population)

# Optional save cleaned file to csv for easier import into other scripts
write.csv(TE_sup, "TE_superfamily_317.csv", row.names = F)


#################################
# TE FAMILY mean/sd proportions #
#################################

# Individuals from 2024 re-analysis
TEfam_avg_raw <- read.delim("apr21_allsamples_TE_annotation_results_average.txt", 
                             header = T, sep = "\t")
TEfam_sd_raw <- read.delim("apr21_allsamples_TE_annotation_results_stdev.txt", 
                            header = T, sep = "\t")

# There are the same issues with these files headers
# The fix is the same EXCEPT you need to change the names of the columns to drop
# when shifting over SEE BELOW

# Issue 1 and 2 (remove path + extra column)
# Average
TEfam_avg_cor <- TEfam_avg_raw %>%
  # Remove path information from ID column
  mutate(Species = str_remove_all(Species, "4_Transposome/")) %>%
  # Drop useless column that contains non-TE proportions (X)
  # Drop tRNASAT.1_ZM column because we'll shift it over
  # THIS IS THE CHANGE
  dplyr::select(-c(X, tRNASAT.1_ZM))
# Standard deviation
TEfam_sd_cor <- TEfam_sd_raw %>%
  # Remove path information from ID column
  mutate(Species = str_remove_all(Species, "4_Transposome/")) %>%
  # Drop useless column that contains non-TE proportions
  # Drop unlcassified column because we'll shift it over
  # THIS IS THE CHANGE
  dplyr::select(-c(X, tRNASAT.1_ZM))

# Issue 3 (fix column headers)
# Extract raw column names
TEfam_head_raw <- colnames(TEfam_avg_raw)
# Remove blank column names (X and X.)
TEfam_head_correct <- TEfam_head_raw[-c(2:3)]
# Rename columns with correct headers
# Mean
names(TEfam_avg_cor) <- c(TEfam_head_correct) %>% 
  paste0("_mean")
# Rename Species_mean column header as ID
TEfam_avg_cor <- TEfam_avg_cor %>% 
  as_tibble() %>% 
  dplyr::rename(ID = Species_mean)
# Standard dev
names(TEfam_sd_cor) <- c(TEfam_head_correct) %>% 
  paste0("_sd")
# Rename Species_sd column header as ID
TEfam_sd_cor <- TEfam_sd_cor %>% 
  as_tibble() %>% 
  dplyr::rename(ID = Species_sd)

# Join all TE superfamily mean/sd data
TE_fam <- TEfam_avg_cor %>% 
  left_join(TEfam_sd_cor) %>% 
  # Split ID column into pop and individual
  separate(ID, c("Species","Collection","Population","Individual"),remove = F) %>% 
  dplyr::select(-c(Species, Collection))
# Add "p_" to all population columns to preserve 0 in front of some pop IDs
# when saving to csv
TE_fam$Population <- sub("^","p_",TE_fam$Population)

# DON'T IGNORE THIS I AM REPLACING NAs WITH 0s HERE
# Because these data weren't analyzed all together there are differneces in
# what TE families were recovered for the samples in the 05.03 and 03.04 transponsome runs.
# I am replacing the NAs (not recovered TE family) with 0s here so the tSNE analysis can run
TE_fam[is.na(TE_fam)] <- 0

# Optional save cleaned file to csv for easier import into other scripts
write.csv(TE_fam, "TE_family_317.csv", row.names = F)



#######################
# Population Metadata #
#######################

###############################
# 0.5 arc minute BioClim Data #
###############################
# Read in BioClim data
# NOTE: this data was extracted for each pop using the BioClim_Blitz.R script
Pop_bio0.5 <- read.csv("Prepped Files/joinbioclim0.5_alt.csv", header = T) %>%  
  # Had to get rid of some dummy column in this
  dplyr::select(-c(X)) %>% 
  as_tibble()
# Rescale BioClim data to actual units
# BioClim scales
# Apply scale factor to raw bioclim data
Pop_bio0.5$bio1 <- (Pop_bio0.5$bio1)*0.1
Pop_bio0.5$bio2 <- (Pop_bio0.5$bio2)*0.1
Pop_bio0.5$bio4 <- (Pop_bio0.5$bio4)*0.01
Pop_bio0.5$bio5 <- (Pop_bio0.5$bio5)*0.1
Pop_bio0.5$bio6 <- (Pop_bio0.5$bio6)*0.1
Pop_bio0.5$bio7 <- (Pop_bio0.5$bio7)*0.1
Pop_bio0.5$bio8 <- (Pop_bio0.5$bio8)*0.1
Pop_bio0.5$bio9 <- (Pop_bio0.5$bio9)*0.1
Pop_bio0.5$bio10 <- (Pop_bio0.5$bio10)*0.1
Pop_bio0.5$bio11 <- (Pop_bio0.5$bio11)*0.1


##############################
# Lvl III Ecoregion Metadata #
##############################
# Add in metadata with ecoregions
Pop_ecoIII <- read.csv("JGI_grant_map_meta.csv", header = T) %>% 
  select(c(population_number,state:longitude,ecoregion)) %>%  
  as_tibble() %>% 
  rename(Population = population_number,
         State = state,
         Locality = locality,
         Latitude = latitude,
         Longitude = longitude,
         Ecoregion = ecoregion)
# Set ecoregion to be factor
Pop_ecoIII$Ecoregion <- as.factor(Pop_ecoIII$Ecoregion)


############################
# Master Pop Metadata file #
############################
# Join ecoregion and bioclims for master file
Pop_meta <- Pop_bio0.5 %>% 
  left_join(Pop_ecoIII)
# Write to csv for easy loading downstream
write.csv(Pop_meta, "Pop_meta.csv", row.names = F)


######################
# Molecular Metadata #
######################
# This is the isolation, library, and sequencing metadate from all individuals
# This is useful for checking batch effects
TE_mol <- read.csv("TE_molecular_meta.csv", header = T)


#####################
# TE Pipeline Stats #
#####################
# These are the number of read pairs at each stage of the TE pipeline
# Trimmed = post-trimmomatic read pairs
# Organella = post-organellar (mitochondria/chloroplast) Bowtie filtering
# Fungal = post-fungal contamination bowtie filtering
# Bacterial = post-bacterial contamination bowtie filtering
# Deduplicated = post-nubeamdedup (PCR duplicates) filtering

# Read in files from different TE pipeline runs
#TE_pipe1 <- read.delim("TE_prep_stats_02.22.22_212.txt", 
#                             header = T, sep = "\t")
#TE_pipe2 <- read.delim("TE_prep_stats_05.03.22_107.txt", 
#                       header = T, sep = "\t")
# Rowbind data together
#TE_pipe <- bind_rows(TE_pipe1, TE_pipe2)



####################################
# Bind TE proportions and metadata #
####################################
# Bind to SUPERFAMILY
MR_sup <- TE_sup %>% 
  left_join(Pop_bio0.5) %>% 
  left_join(Pop_ecoIII) %>% 
  left_join(TE_mol) %>% 
  left_join(TE_pipe) %>% 
  relocate(ID:Individual,State:Locality,Ecoregion,bio1:alt,
           Iso_ID,Lib_ID,Primer_501,Primer_701,Seq_run,
           Trimmed:Deduplicated)
# Bind to FAMILY
MR_fam <- TE_fam %>% 
  left_join(Pop_bio0.5) %>% 
  left_join(Pop_ecoIII) %>% 
  left_join(TE_mol) %>% 
  left_join(TE_pipe) %>% 
  relocate(ID:Individual,State:Locality,Ecoregion,bio1:alt,
           Iso_ID,Lib_ID,Primer_501,Primer_701,Seq_run,
           Trimmed:Deduplicated)



###################################################
# Write master CSV files for downstream analayses #
###################################################
# SUPERFAMILY
write.csv(MR_sup, "MR_sup.csv", row.names = F)
# FAMILY
write.csv(MR_fam, "MR_fam.csv", row.names = F)


