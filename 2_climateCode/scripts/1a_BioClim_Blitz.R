# Bioclim blitz #
# Description: code to extract 0.5, 2.5 and 5.0 arc second raster data
# from WorldClim database
# Need to provide (1) input population file with lat in column 4 and long in
# column 5.
# Author: Nate Hofford

setwd("~/Research/McKain Lab/CPING/TE/R/TE Pipeline/April 2024 Analysis")


#####################
# STEP 0: Libraries #
#####################

# Data management
library(tidyverse)
# Packages for BioClim variable extraction
# Note: these can be finicky and rgdal is going to be deprecated in 2023
# Might need to reload r / check compatible version are loading
library(rgdal)



########################
# STEP 1: Load in data #
########################

# Read in csv of populations
TE_pops <- read.csv("Raw Data/TE_pops_all.csv", header = T)



#####################################################
# STEP 2: Create functions for pulling bioclim data #
#####################################################

# Create pull bioclim function for 0.5 scale
# Note that the 0.5 'getData()' command requires explicit lat/long
pull_bioclim0.5 <- function(geotable) {
  final_df <- data.frame(matrix(ncol=21))
  colnames(final_df)<-c("Longitude","Latitude","bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18", "bio19")
  for (i in 1:nrow(geotable)) {
    w<-getData('worldclim', var='bio', res=0.5, lon=geotable[i,5], lat=geotable[i,4])
    coords<-data.frame(x=geotable[i,5], y=geotable[i,4])
    points <- SpatialPoints(coords, proj4string = w@crs)
    values <- extract(w,points)
    df <- cbind.data.frame(coordinates(points),values)
    colnames(df) <- colnames(final_df)
    print (df)
    final_df <- rbind(final_df, df);
  }
  final_df <- final_df[-1,]
  row.names(final_df) <- row.names(geotable)
  row
  return(final_df)
}

# Pull altitude data
# hole-filled CGIAR-SRTM (90 m resolution)
# Elevation data is in meters (m)
pull_alt0.5 <- function(geotable) {
  final_df <- data.frame(matrix(ncol=3))
  colnames(final_df)<-c("Longitude","Latitude","alt")
  for (i in 1:nrow(geotable)) {
    w<-getData('SRTM',lon=geotable[i,5], lat=geotable[i,4])
    coords<-data.frame(x=geotable[i,5], y=geotable[i,4])
    points <- SpatialPoints(coords, proj4string = w@crs)
    values <- extract(w,points)
    df <- cbind.data.frame(coordinates(points),values)
    colnames(df) <- colnames(final_df)
    print (df)
    final_df <- rbind(final_df, df);
  }
  final_df <- final_df[-1,]
  row.names(final_df) <- row.names(geotable)
  row
  return(final_df)
}


library(climateR)
library(terra)
library(sf)

library(geodata)
library(terra)
library(sf)

pull_bioclim0.5 <- function(geotable) {
  # Download WorldClim bioclimatic variables at 0.5-degree resolution using geodata
  wc_data <- geodata::worldclim_global(var = "bio", res = 0.5, path = tempdir())
  
  # Convert geotable into an sf object (assuming columns 4 and 5 are Latitude and Longitude)
  geotable_sf <- st_as_sf(geotable, coords = c(5, 4), crs = 4326)
  
  # Extract raster values at point locations
  extracted_values <- terra::extract(wc_data, vect(geotable_sf), cells = FALSE)
  
  # Combine coordinates and extracted values into a final data frame
  coords <- st_coordinates(geotable_sf)
  final_df <- data.frame(
    Longitude = coords[, 1],
    Latitude = coords[, 2],
    extracted_values[, -1]  # Drop the ID column
  )
  
  # Rename columns for clarity
  colnames(final_df) <- c(
    "Longitude", "Latitude", "bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7",
    "bio8", "bio9", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", 
    "bio16", "bio17", "bio18", "bio19"
  )
  
  row.names(final_df) <- row.names(geotable)
  return(final_df)
}



# Function for 2.5 arc sec scale
# Note that lat/long is not required for 'getData()'
#pull_bioclim2.5 <- function(geotable) {
#  final_df <- data.frame(matrix(ncol=21))
#  colnames(final_df)<-c("Longitude","Latitude","bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18", "bio19")
#  for (i in 1:nrow(geotable)) {
#    w<-getData('worldclim', var='bio', res=2.5)
#    coords<-data.frame(x=geotable[i,5], y=geotable[i,4])
#    points <- SpatialPoints(coords, proj4string = w@crs)
#    values <- extract(w,points)
#    df <- cbind.data.frame(coordinates(points),values)
#    colnames(df) <- colnames(final_df)
#    print (df)
#    final_df <- rbind(final_df, df);
#  }
#  final_df <- final_df[-1,]
#  row.names(final_df) <- row.names(geotable)
#  row
#  return(final_df)
#}

# Function for 5 arc sec scale
# Note that lat/long is not required for 'getData()'
#pull_bioclim5.0 <- function(geotable) {
#  final_df <- data.frame(matrix(ncol=21))
#  colnames(final_df)<-c("Longitude","Latitude","bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18", "bio19")
#  for (i in 1:nrow(geotable)) {
#    w<-getData('worldclim', var='bio', res=5.0)
#    coords<-data.frame(x=geotable[i,5], y=geotable[i,4])
#    points <- SpatialPoints(coords, proj4string = w@crs)
#    values <- extract(w,points)
#    df <- cbind.data.frame(coordinates(points),values)
#    colnames(df) <- colnames(final_df)
#    print (df)
#    final_df <- rbind(final_df, df);
#  }
#  final_df <- final_df[-1,]
#  row.names(final_df) <- row.names(geotable)
#  row
#  return(final_df)
#}



#############################
# STEP 3: Execute functions #
#############################

# Pull bioclim data using functions
# Note: this will take a minute
TE_bioclim0.5 <- pull_bioclim0.5(TE_pops)
#TE_bioclim2.5 <- pull_bioclim2.5(TE_pops)
#TE_bioclim5.0 <- pull_bioclim5.0(TE_pops)
TE_srtm <- pull_alt0.5(TE_pops)

# Save raw 0.5 bioclim data as csv so this doesn't have to happen again
write.csv(TE_bioclim0.5, file = "Raw Data/rawbioclim0.5.csv")
# OPTIONAL read in of raw 0.5 bioclim data
#TE_bioclim0.5 <- read.csv("rawbioclim0.5.csv", header = T) %>% 
 # dplyr::select(-c(1))


################################################
# STEP 4: Join population metadata and bioclim #
################################################

# Join bioclim and pop location info by lat/long
# I remove columns 6-11 here because they contain status information
# that we don't need for this
TE_joined0.5 <- TE_pops %>% 
  # Rename lat/long columns so variables are standardized
  rename(Latitude = Lat, Longitude = Long) %>% 
  # Inner join to filter non-TE pops
  inner_join(TE_bioclim0.5) %>% 
  # Convert to tibble
  tibble()
# Join 2.5 and 5.0 as well
#TE_joined2.5 <- TE_pops %>% 
  # Rename lat/long columns so variables are standardized
#  rename(Latitude = Lat, Longitude = Long) %>% 
  # Inner join to filter non-TE pops
#  inner_join(TE_bioclim2.5) %>% 
  # Convert to tibble
#  tibble()
#TE_joined5.0 <- TE_pops[-c(6:11)] %>% 
  # Rename lat/long columns so variables are standardized
#  rename(Latitude = Lat, Longitude = Long) %>% 
  # Inner join to filter non-TE pops
#  inner_join(TE_bioclim5.0) %>% 
  # Convert to tibble
#  tibble()
# Join altitude data
TE_joinedalt <- TE_pops %>% 
  # Rename lat/long columns so variables are standardized
  rename(Latitude = Lat, Longitude = Long) %>% 
  # Inner join to filter non-TE pops
  inner_join(TE_srtm) %>% 
  # Convert to tibble
  tibble()

# Join bio0.5 and alt together to simplify export
TE_bio0.5_alt <- TE_joined0.5 %>% 
  left_join(TE_joinedalt)


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
