# Bioclim blitz #
# Description: code to extract bioclimatic raster data from WorldClim database
# and elevation data using modern R spatial packages.
#
# Note: this is still a work in progress to extract all neceassry climate
# data in one single script. Proceeding with allEnv.csv for now
# since it has soils + other data incorporated already
#
# Input: population file with Latitude in column 4 and Longitude in column 5
# Author: Nate Hofford
# Updated: 2026 - migrated from rgdal/raster to sf/terra/geodata


#setwd("~/Research/McKain Lab/CPING/TE/R/TE Pipeline/April 2024 Analysis")


#####################
# STEP 0: Libraries #
#####################

library(tidyverse)
library(sf)       # Vector data handling (replaces sp/rgdal)
library(terra)    # Raster data handling (replaces raster)
library(geodata)  # WorldClim and SRTM elevation data downloads



########################
# STEP 1: Load in data #
########################

TE_pops <- read.csv("data/TE_pops_all.csv", header = T)



#####################################################
# STEP 2: Create functions for pulling bioclim data #
#####################################################

#' Pull WorldClim bioclimatic variables
#' Loads from cached CSV if already downloaded, otherwise downloads and saves.
#'
#' @param geotable   Data frame with Latitude in column 4, Longitude in column 5
#' @param res        Resolution: 0.5, 2.5, 5, or 10 (arc-minutes)
#' @param cache_csv  Path to cache CSV (default: data/rawbioclim_<res>.csv)
#' @param data_path  Directory to cache raster tiles (default: tempdir())
#' @return Data frame with Longitude, Latitude, and bio1-bio19
pull_bioclim <- function(geotable, res = 0.5,
                         cache_csv = paste0("data/rawbioclim_", res, ".csv"),
                         data_path = tempdir()) {

  if (file.exists(cache_csv)) {
    message(paste("Loading cached bioclim data from:", cache_csv))
    return(read.csv(cache_csv))
  }

  message(paste("Downloading WorldClim data at", res, "arc-minute resolution..."))
  wc_data <- geodata::worldclim_global(var = "bio", res = res, path = data_path)

  geotable_sf <- st_as_sf(geotable, coords = c(5, 4), crs = 4326)

  message("Extracting values at sample locations...")
  extracted_values <- terra::extract(wc_data, vect(geotable_sf))

  coords <- st_coordinates(geotable_sf)
  final_df <- data.frame(
    Longitude = coords[, 1],
    Latitude  = coords[, 2],
    extracted_values[, -1]  # Drop the ID column from extract()
  )
  colnames(final_df) <- c("Longitude", "Latitude", paste0("bio", 1:19))
  row.names(final_df) <- row.names(geotable)

  message(paste("Saving bioclim data to:", cache_csv))
  write.csv(final_df, file = cache_csv, row.names = FALSE)

  message(paste("Extracted bioclim data for", nrow(final_df), "locations"))
  return(final_df)
}


#' Pull elevation data using geodata (SRTM 30 arc-second ~1km resolution)
#' Loads from cached CSV if already downloaded, otherwise downloads and saves.
#'
#' @param geotable   Data frame with Latitude in column 4, Longitude in column 5
#' @param cache_csv  Path to cache CSV (default: data/srtm.csv)
#' @param data_path  Directory to cache raster tiles (default: tempdir())
#' @return Data frame with Longitude, Latitude, and alt (meters)
pull_elevation <- function(geotable,
                           cache_csv = "data/srtm.csv",
                           data_path = tempdir()) {

  if (file.exists(cache_csv)) {
    message(paste("Loading cached elevation data from:", cache_csv))
    return(read.csv(cache_csv))
  }

  message("Downloading SRTM elevation data (30 arc-second resolution, ~1km)...")
  elev_raster <- geodata::elevation_global(res = 0.5, path = data_path)

  geotable_sf <- st_as_sf(geotable, coords = c(5, 4), crs = 4326)

  extracted_values <- terra::extract(elev_raster, vect(geotable_sf))

  coords <- st_coordinates(geotable_sf)
  final_df <- data.frame(
    Longitude = coords[, 1],
    Latitude  = coords[, 2],
    alt       = extracted_values[, 2]
  )
  row.names(final_df) <- row.names(geotable)

  message(paste("Saving elevation data to:", cache_csv))
  write.csv(final_df, file = cache_csv, row.names = FALSE)

  message(paste("Extracted elevation for", nrow(final_df), "locations"))
  return(final_df)
}



#############################
# STEP 3: Execute functions #
#############################

# Pull bioclim - loads from data/rawbioclim_0.5.csv if it exists
TE_bioclim0.5 <- pull_bioclim(TE_pops, res = 0.5)

# Optional: other resolutions (uncomment as needed)
# TE_bioclim2.5 <- pull_bioclim(TE_pops, res = 2.5)
# TE_bioclim5.0 <- pull_bioclim(TE_pops, res = 5)

# Pull SRTM elevation - loads from data/srtm.csv if it exists
TE_srtm <- pull_elevation(TE_pops)

# Assign EPA Level 3 Ecoregions via spatial join with local shapefile
message("Loading EPA Level III Ecoregions...")
ecoregions <- st_read("data/ecoregions/us_eco_l3/us_eco_l3.shp", quiet = TRUE) %>%
  st_transform(4326) %>%
  st_make_valid()

geotable_sf <- st_as_sf(TE_pops, coords = c(5, 4), crs = 4326)
eco_joined  <- st_join(geotable_sf, ecoregions, join = st_within)
eco_coords  <- st_coordinates(geotable_sf)

TE_ecoregions <- data.frame(
  Longitude = eco_coords[, 1],
  Latitude  = eco_coords[, 2],
  L3_KEY    = eco_joined$US_L3CODE,
  L3_NAME   = eco_joined$US_L3NAME
)

n_assigned <- sum(!is.na(TE_ecoregions$L3_NAME))
n_missing  <- sum(is.na(TE_ecoregions$L3_NAME))
message(paste("Assigned L3 ecoregions to", n_assigned, "locations"))
if (n_missing > 0) {
  message(paste("Warning:", n_missing, "locations outside US ecoregion boundaries"))
}



################################################
# STEP 4: Join population metadata and bioclim #
################################################

TE_joined0.5 <- TE_pops %>%
  rename(Latitude = Lat, Longitude = Long) %>%
  inner_join(TE_bioclim0.5, by = c("Longitude", "Latitude")) %>%
  tibble()

TE_joinedalt <- TE_pops %>%
  rename(Latitude = Lat, Longitude = Long) %>%
  inner_join(TE_srtm, by = c("Longitude", "Latitude")) %>%
  tibble()

# Join all data together: bioclim + elevation + L3 ecoregions
TE_bio0.5_alt <- TE_joined0.5 %>%
  left_join(TE_joinedalt,    by = c("Longitude", "Latitude")) %>%
  left_join(TE_ecoregions,   by = c("Longitude", "Latitude"))



############################
# STEP 5: Write final CSVs #
############################

# Final joined file: population metadata + bio1-bio19 + elevation + L3 ecoregion
write.csv(TE_bio0.5_alt,
          file = "data/joinbioclim0.5_alt_ecoregion.csv",
          row.names = FALSE)

message("\n========== FINAL DATASET SUMMARY ==========")
message(paste("Total samples:       ", nrow(TE_bio0.5_alt)))
message(paste("Unique L3 ecoregions:", length(unique(na.omit(TE_bio0.5_alt$L3_NAME)))))
message(paste("Total columns:       ", ncol(TE_bio0.5_alt)))
message("Variables: bio1-bio19, alt, L3_KEY, L3_NAME")
message("============================================\n")

# YOU NOW HAVE YOUR BIOCLIM, ALTITUDE, AND ECOREGION DATA
