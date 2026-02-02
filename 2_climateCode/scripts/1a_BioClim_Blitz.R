# Bioclim blitz #
# Description: code to extract bioclimatic raster data from WorldClim database
# and elevation data using modern R spatial packages.
#
# Input: population file with Latitude in column 4 and Longitude in column 5
# Author: Nate Hofford
# Updated: 2026 - migrated from rgdal/raster to sf/terra/geodata

#setwd("~/Research/McKain Lab/CPING/TE/R/TE Pipeline/April 2024 Analysis")


#####################
# STEP 0: Libraries #
#####################

# Data management
library(tidyverse)

# Modern spatial packages (replacing deprecated rgdal/raster/sp)
library(sf)        # Vector data handling (replaces sp/rgdal)
library(terra)     # Raster data handling (replaces raster)
library(geodata)   # WorldClim and elevation data downloads
library(elevatr)   # Alternative elevation data source



########################
# STEP 1: Load in data #
########################

# Read in csv of populations
TE_pops <- read.csv("Raw Data/TE_pops_all.csv", header = T)



#####################################################
# STEP 2: Create functions for pulling bioclim data #
#####################################################

#' Pull WorldClim bioclimatic variables
#' 
#' @param geotable Data frame with Latitude in column 4, Longitude in column 5
#' @param res Resolution: 0.5, 2.5, 5, or 10 (arc-minutes)
#' @param data_path Directory to cache downloaded data (default: tempdir())
#' @return Data frame with coordinates and 19 bioclim variables
pull_bioclim <- function(geotable, res = 0.5, data_path = tempdir()) {
  
  message(paste("Downloading WorldClim data at", res, "arc-minute resolution..."))
  
  # Download WorldClim bioclimatic variables
  wc_data <- geodata::worldclim_global(var = "bio", res = res, path = data_path)
  
  # Convert geotable to sf object (column 5 = Longitude, column 4 = Latitude)
  geotable_sf <- st_as_sf(geotable, coords = c(5, 4), crs = 4326)
  
  # Extract raster values at point locations
  message("Extracting values at sample locations...")
  extracted_values <- terra::extract(wc_data, vect(geotable_sf))
  
  # Combine coordinates and extracted values
  coords <- st_coordinates(geotable_sf)
  final_df <- data.frame(
    Longitude = coords[, 1],
    Latitude = coords[, 2],
    extracted_values[, -1]  # Drop the ID column from extract()
  )
  
  # Rename columns for clarity
  colnames(final_df) <- c(
    "Longitude", "Latitude", 
    paste0("bio", 1:19)
  )
  
  row.names(final_df) <- row.names(geotable)
  
  message(paste("Extracted bioclim data for", nrow(final_df), "locations"))
  return(final_df)
}


#' Pull elevation data using elevatr (AWS Terrain Tiles)
#' 
#' @param geotable Data frame with Latitude in column 4, Longitude in column 5
#' @param z Zoom level for resolution (1-14, higher = finer resolution)
#' @return Data frame with coordinates and elevation
pull_elevation <- function(geotable, z = 9) {
  
  message("Downloading elevation data...")
  
  # Convert geotable to sf object
  geotable_sf <- st_as_sf(geotable, coords = c(5, 4), crs = 4326)
  
  # Get elevation using elevatr (uses AWS terrain tiles)
  elev_data <- get_elev_point(geotable_sf, src = "aws", z = z)
  
  # Extract coordinates and elevation
  coords <- st_coordinates(geotable_sf)
  final_df <- data.frame(
    Longitude = coords[, 1],
    Latitude = coords[, 2],
    alt = elev_data$elevation
  )
  
  row.names(final_df) <- row.names(geotable)
  
  message(paste("Extracted elevation for", nrow(final_df), "locations"))
  return(final_df)
}


#' Alternative: Pull elevation using geodata (SRTM 30s resolution)
#' 
#' @param geotable Data frame with Latitude in column 4, Longitude in column 5
#' @param data_path Directory to cache downloaded data
#' @return Data frame with coordinates and elevation
pull_elevation_srtm <- function(geotable, data_path = tempdir()) {
  
  message("Downloading SRTM elevation data...")
  
  # Download SRTM elevation data (30 arc-second resolution)
  elev_raster <- geodata::elevation_global(res = 0.5, path = data_path)
  
  # Convert geotable to sf object
  geotable_sf <- st_as_sf(geotable, coords = c(5, 4), crs = 4326)
  
  # Extract elevation values
  extracted_values <- terra::extract(elev_raster, vect(geotable_sf))
  
  # Combine coordinates and elevation
  coords <- st_coordinates(geotable_sf)
  final_df <- data.frame(
    Longitude = coords[, 1],
    Latitude = coords[, 2],
    alt = extracted_values[, 2]
  )
  
  row.names(final_df) <- row.names(geotable)
  
  message(paste("Extracted elevation for", nrow(final_df), "locations"))
  return(final_df)
}


# Legacy wrapper for backwards compatibility
pull_bioclim0.5 <- function(geotable) {
  pull_bioclim(geotable, res = 0.5)
}

pull_alt0.5 <- function(geotable) {
  pull_elevation(geotable)
}



# Additional resolution wrappers (if needed)
pull_bioclim2.5 <- function(geotable) {
  pull_bioclim(geotable, res = 2.5)
}

pull_bioclim5.0 <- function(geotable) {
  pull_bioclim(geotable, res = 5)
}

pull_bioclim10 <- function(geotable) {
  pull_bioclim(geotable, res = 10)
}



#############################
# STEP 3: Execute functions #
#############################

# Pull bioclim data using functions
# Note: First run will download data (~1GB for 0.5 resolution)
message("Starting bioclim extraction...")
TE_bioclim0.5 <- pull_bioclim(TE_pops, res = 0.5)

# Optional: other resolutions (uncomment as needed)
# TE_bioclim2.5 <- pull_bioclim(TE_pops, res = 2.5)
# TE_bioclim5.0 <- pull_bioclim(TE_pops, res = 5)

# Pull elevation data
message("Starting elevation extraction...")
TE_srtm <- pull_elevation(TE_pops)
# Alternative using SRTM: TE_srtm <- pull_elevation_srtm(TE_pops)

# Save raw 0.5 bioclim data as csv so this doesn't have to happen again
write.csv(TE_bioclim0.5, file = "Raw Data/rawbioclim0.5.csv", row.names = FALSE)

# OPTIONAL: read in previously saved bioclim data
# TE_bioclim0.5 <- read.csv("Raw Data/rawbioclim0.5.csv", header = TRUE)


################################################
# STEP 4: Join population metadata and bioclim #
################################################

# Join bioclim and pop location info by lat/long
TE_joined0.5 <- TE_pops %>% 
  # Rename lat/long columns so variables are standardized
  rename(Latitude = Lat, Longitude = Long) %>% 
  # Inner join to filter non-TE pops

  inner_join(TE_bioclim0.5, by = c("Longitude", "Latitude")) %>% 
  # Convert to tibble
  tibble()

# Join altitude data
TE_joinedalt <- TE_pops %>% 
  rename(Latitude = Lat, Longitude = Long) %>% 
  inner_join(TE_srtm, by = c("Longitude", "Latitude")) %>% 
  tibble()

# Join bio0.5 and alt together to simplify export
TE_bio0.5_alt <- TE_joined0.5 %>% 
  left_join(TE_joinedalt, by = c("Longitude", "Latitude"))

# Optional: join other resolutions
# TE_joined2.5 <- TE_pops %>% 
#   rename(Latitude = Lat, Longitude = Long) %>% 
#   inner_join(TE_bioclim2.5, by = c("Longitude", "Latitude")) %>% 
#   tibble()


############################
# STEP 5: Write final CSVs #
############################

# Write joined files to a csv for saving
#write.csv(TE_joined0.5, file = "joinbioclim0.5.csv")
#write.csv(TE_joined2.5, file = "joinbioclim2.5.csv")
#write.csv(TE_joined5.0, file = "joinbioclim5.0.csv")
#write.csv(TE_joinedalt, file = "joinalt.csv")
write.csv(TE_bio0.5_alt, file = "Prepped Files/joinbioclim0.5_alt.csv")

# YOU NOW HAVE YOUR BIOCLIM AND ALTITUDE DATA
