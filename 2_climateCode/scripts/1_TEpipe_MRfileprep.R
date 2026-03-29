# TE Analysis Pipeline
# Script 1: File Preparation
#
# Preps TE proportion files (superfamily and family) and joins to
# population metadata from allEnv.csv to create master analysis files.
#
# Input files:
#   - data/output/*_superfamily_TE_annotation_results_average.txt
#   - data/output/*_superfamily_TE_annotation_results_stdev.txt
#   - data/output/*_family_TE_annotation_results_average.txt
#   - data/output/*_family_TE_annotation_results_stdev.txt
#   - 2_climateCode/data/allEnv.csv
#
# Output files:
#   - MR_sup.csv  (superfamily proportions + metadata)
#   - MR_fam.csv  (family proportions + metadata)


###################
# Library loading #
###################

library(tidyverse)


#############################
# File paths - update here  #
#############################

setwd("~/Documents/GitHub/jumpingGenesPipe")

sup_avg_path <- "data/output/feb1_286_superfamily_TE_annotation_results_average.txt"
sup_sd_path  <- "data/output/feb1_286_superfamily_TE_annotation_results_stdev.txt"
fam_avg_path <- "data/output/feb1_286_family_TE_annotation_results_average.txt"
fam_sd_path  <- "data/output/feb1_286_family_TE_annotation_results_stdev.txt"
env_path     <- "2_climateCode/data/allEnv.csv"


#' Clean raw TE annotation output file
#'
#' Handles the two issues in the pipeline output:
#'   1. Species column contains path prefix "4_Transposome/"
#'   2. The "-" TE type column is non-informative and dropped
#'
#' @param path     Path to the .txt file
#' @param suffix   Column suffix to append: "_mean" or "_sd"
#' @return Tibble with ID, Population, Individual columns + TE proportion columns
clean_te_file <- function(path, suffix) {
  
  raw <- read.delim(path, header = TRUE, sep = "\t")
  
  raw %>%
    # Remove path prefix from ID column
    mutate(Species = str_remove(Species, "^4_Transposome/")) %>%
    # Drop the "-" TE type column (R reads it as "X.")
    dplyr::select(-any_of("X.")) %>%
    # Rename TE columns with suffix
    rename_with(~ paste0(.x, suffix), -Species) %>%
    # Rename Species as ID
    rename(ID = Species) %>%
    # Parse population and individual from ID (format: JG_CO_011_L)
    separate(ID, into = c("sp", "col", "Population", "Individual"),
             sep = "_", remove = FALSE) %>%
    dplyr::select(-sp, -col) %>%
    # Add "p_" prefix to population to preserve leading zeros
    mutate(Population = paste0("p_", Population)) %>%
    as_tibble()
}


###############################
# STEP 1: Load population env #
###############################

Pop_env <- read.csv(env_path, header = TRUE) %>%
  dplyr::select(-Lat, -Long) %>%
  as_tibble()


######################################
# STEP 2: TE SUPERFAMILY proportions #
######################################

TEsup_avg <- clean_te_file(sup_avg_path, "_mean")
TEsup_sd  <- clean_te_file(sup_sd_path,  "_sd")

# Join mean and sd, drop duplicate metadata columns from sd
TE_sup <- TEsup_avg %>%
  left_join(TEsup_sd %>% dplyr::select(-Population, -Individual),
            by = "ID")


#################################
# STEP 3: TE FAMILY proportions #
#################################

TEfam_avg <- clean_te_file(fam_avg_path, "_mean")
TEfam_sd  <- clean_te_file(fam_sd_path,  "_sd")

# Join mean and sd
TE_fam <- TEfam_avg %>%
  left_join(TEfam_sd %>% dplyr::select(-Population, -Individual),
            by = "ID")

# DON'T IGNORE THIS I AM REPLACING NAs WITH 0s HERE
# Because these data weren't analyzed all together there are differneces in
# what TE families were recovered for the samples in the 05.03 and 03.04 transponsome runs.
# I am replacing the NAs (not recovered TE family) with 0s here so the tSNE analysis can run
TE_fam[is.na(TE_fam)] <- 0


##############################################
# STEP 4: Join TE data to population metadata #
##############################################

MR_sup <- TE_sup %>%
  left_join(Pop_env, by = "Population") %>%
  # Put metadata columns first, then TE proportions
  relocate(ID, Population, Individual,
           State, Locality, Latitude, Longitude,
           Ecoregion, bio1:alt)

MR_fam <- TE_fam %>%
  left_join(Pop_env, by = "Population") %>%
  relocate(ID, Population, Individual,
           State, Locality, Latitude, Longitude,
           Ecoregion, bio1:alt)


####################
# STEP 5: QC check #
####################

message("========== FILE PREP SUMMARY ==========")
message(paste("Superfamily samples:", nrow(MR_sup)))
message(paste("Family samples:     ", nrow(MR_fam)))
message(paste("Populations in env: ", nrow(Pop_env)))
message(paste("Matched in MR_sup:  ", sum(!is.na(MR_sup$Latitude))))
message(paste("Matched in MR_fam:  ", sum(!is.na(MR_fam$Latitude))))

unmatched <- MR_sup %>%
  filter(is.na(Latitude)) %>%
  distinct(Population) %>%
  pull(Population)
if (length(unmatched) > 0) {
  message(paste("WARNING - populations with no env match:", 
                paste(unmatched, collapse = ", ")))
}
message("========================================")


#################################################
# STEP 6: Write master CSVs for downstream work #
#################################################

write.csv(MR_sup, "2_climateCode/data/MR_sup.csv", row.names = FALSE)
write.csv(MR_fam, "2_climateCode/data/MR_fam.csv", row.names = FALSE)

message("Output written:")
message("  2_climateCode/data/MR_sup.csv")
message("  2_climateCode/data/MR_fam.csv")
