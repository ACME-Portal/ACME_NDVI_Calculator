### Oak Bay Deer - NDVI

# Libraries 
library(mapedit)
library(tidyverse)
library(lme4)
library(lattice)
library(MuMIn)
library(AICcmodavg)
library(faraway)
library(knitr)
library(AICcmodavg)
library(DT)
library(boot)
library(dotwhisker)
library(stargazer)
library(cowplot)
library(broom)
library(sf)
library(tmap)
library(raster)
library(blockCV)
library(shinyjs)
library(tmaptools)
library(precrec)
library(mapview)
library(visreg)
library(MASS)
library(stars)
library(lwgeom)
library(dplyr)
library(GmAMisc)
library(BAMMtools)
library(Hmisc)


## Steps ---
# Get data for CRD, summer 2018
# Download data, load into R
# Stack datasets 
# Calculate NDVI
# Average NDVI for summer months

### Part 1 -- Load Landsat Layers ------ 
# Source - https://www.earthdatascience.org/courses/earth-analytics/multispectral-remote-sensing-data/landsat-data-in-r-geotiff/
options(stringsAsFactors = FALSE)  # Turn off factors

# Folder LC080470262018050501T1-SC20200204205109 -- Image: LC08_L1TP_047026_20180505_20180517
landsat_2 <- list.files("Landsat/Wilmore/LE070450222009072501T1-SC20200703175144.tar", 
                      pattern = glob2rx("*band*.tif$"),
                      full.names = TRUE)


### Stack bands in Landsat Images ---
landsat_2 <- stack(landsat_2) %>% brick()  

### View data ---
mapview::viewRGB(landsat_2, r = 4, g = 3, b = 2)
# Looks dark, but good. Veg is green, water is blue, urban is grey, snow is white, Frazer plume is silty.



### Part 2 -- NDVI ------

# Option 1) NDVI manually 
landsat_2_ndvi <- (landsat_2[[5]] - landsat_2[[4]]) / (landsat_2[[5]] + landsat_2[[4]])
x11(10,10)
plot(landsat_2_ndvi)  # Some error here. Bands?

backup_landsat_2_ndvi <- landsat_2_ndvi

## Option 2) NDVI with function
# Create a function for NDVI
ndvi_function <- function(b5, b4){
  diff <- (b5 - b4) / (b5 + b4)
  return(diff)
}
# Apply the function to the raster layers with overlay()
landsat_ndvi_ov <- overlay(landsat_2[[5]],
                           landsat_2[[4]],
                           fun = ndvi_function)

## Fix bands error
landsat_2[[5]]  # Band 5
landsat_2[[4]]  # Band 4

# Calculation looks right. What is wrong?
# NDVI = (NIR - Red) / (NIR + Red)
# source - https://www.usgs.gov/land-resources/nli/landsat/landsat-normalized-difference-vegetation-index?qt-science_support_page_related_con=0#qt-science_support_page_related_con
### RANGE of NDVI is between -1 and 1. This would solve the plotting problem.
# HIGH = 470, LOW = -Inf. This would explain it. Frequency plot?

# freq_landsat_2_ndvi <- freq(landsat_2_ndvi)
# hist_landsat_2_ndvi <- hist(freq(landsat_2_ndvi))   # SLOW STEP - probably not neccesary 

# Some Landsat 8 pixels are saturated. This is a known problem. 
# source - https://gis.stackexchange.com/questions/155178/landsat8-ndvi-values-not-in-the-1-to-1-range
# source - https://www.researchgate.net/post/Why_Landsat_8_NDVI_Values_are_out_of_Range_Not_in_between-1_to_1
# Call all pixels outside +-1 range NA. Will this cause problems? Yes. Call them -1. Or.. multiply them? Hmmm. Help!

# source - https://www.researchgate.net/post/How_might_one_scale_down_the_reflectance_values_in_the_range_0-1_from_landsat_8/1
# You can try different GIS tools to rescale the values from 0 to 1 based on the SR product guide of Landsat 8. Here is a way to rescale the values using Raster Calculator tool in ArcGIS
# SetNull("%Landsat8_Band%","%Landsat8_Band%", "VALUE < 0 OR VALUE > 10000" ) * 0.0001 
# This basically sets all values below 0 and above 10000 as null and converts the rest of the valid reflectance values in the 0 to 1 range.

# Get frequency tables for the two NIR and R bands
freq(landsat_2$LC08_L1TP_047026_20180505_20180517_01_T1_sr_band5)   # NOT from zero to 1. 
freq(landsat_2$LC08_L1TP_047026_20180505_20180517_01_T1_sr_band4)

# Set values below zero as NA 
landsat_2_ndvi <- reclassify(landsat_2_ndvi, cbind(-Inf, 0, NA), right=FALSE)    # right=FALSE avoids setting values of 0 to NA
freq(landsat_2_ndvi)

mapview(landsat_2_ndvi)  # Closer. Would adjusting the values above 1 help?

# Make all values above 1 NA
landsat_2_ndvi <- reclassify(landsat_2_ndvi, cbind(1, 470, NA), right=TRUE)
# test using +Inf instead of a real number here -- reclassify(landsat_2_ndvi, cbind(1, +Inf, NA), right=TRUE)
mapview(landsat_2_ndvi)

# # Make all values above 1 multiplied by 0.0001
# landsat_2_ndvi[landsat_2_ndvi > 1] * 0.0001
# mapview(landsat_2_ndvi)



### Part 3 - Load the rest of the 2018 images -----

### 2018_04_26 ---
landsat_1 <- list.files("Data/Landsat/2018/LC080480262018042601T1-SC20200204205123.tar", 
                        pattern = glob2rx("*band*.tif$"),
                        full.names = TRUE)
# Stack bands
landsat_1 <- stack(landsat_1) %>% brick()  

# Calculate NDVI 
# <Note> bands will be different if using Landsat_7 data rather than Landsat_8 
landsat_1_ndvi <- (landsat_1[[5]] - landsat_1[[4]]) / (landsat_1[[5]] + landsat_1[[4]])
hist(freq(landsat_1_ndvi))



# Reclassify below zero
landsat_1_ndvi <- reclassify(landsat_1_ndvi, cbind(-Inf, 0, NA), right=FALSE)

# Rescale above 1
# test <- landsat_1_ndvi[landsat_1_ndvi > 1, drop = FALSE]
## This DOESN'T work. DROP = FALSE is good though - it keeps the raster instead of converting to a matrix!
mapview(test)

# Convert values above 1 to NA
landsat_1_ndvi <- reclassify(landsat_1_ndvi, cbind(1, 961, NA), right=TRUE)
# test using +Inf instead of a real number here - max pixel value of 961 is arbitrary
mapview(landsat_1_ndvi)
mapview(landsat_2_ndvi)


### 2018_07_15  ---
landsat_3 <- list.files("Data/Landsat/2018/LC080480262018071501T1-SC20200204205116.tar", 
                        pattern = glob2rx("*band*.tif$"),
                        full.names = TRUE)
# Stack bands
landsat_3 <- stack(landsat_3) %>% brick()  

# Calculate NDVI 
landsat_3_ndvi <- (landsat_3[[5]] - landsat_3[[4]]) / (landsat_3[[5]] + landsat_3[[4]])
hist(freq(landsat_3_ndvi))

# Reclassify below zero
landsat_3_ndvi <- reclassify(landsat_3_ndvi, cbind(-Inf, 0, NA), right=FALSE)

# Reclassify above zero
landsat_3_ndvi <- reclassify(landsat_3_ndvi, cbind(1, 961, NA), right=TRUE)
mapview(landsat_3_ndvi)



### 2018_07_24  ---
landsat_4 <- list.files("Data/Landsat/2018/LC080470262018072401T1-SC20200204205121.tar", 
                        pattern = glob2rx("*band*.tif$"),
                        full.names = TRUE)
# Stack bands
landsat_4 <- stack(landsat_4) %>% brick()  

# Calculate NDVI 
landsat_4_ndvi <- (landsat_4[[5]] - landsat_4[[4]]) / (landsat_4[[5]] + landsat_4[[4]])
hist(freq(landsat_4_ndvi))

# Reclassify below zero
landsat_4_ndvi <- reclassify(landsat_4_ndvi, cbind(-Inf, 0, NA), right=FALSE)

# Convert values above 1 to NA
landsat_4_ndvi <- reclassify(landsat_4_ndvi, cbind(1, 961, NA), right=TRUE)
# test using +Inf instead of a real number here



### 2018_08_09  ---
landsat_5 <- list.files("Data/Landsat/2018/LC080470262018080901T1-SC20200204205117.tar", 
                        pattern = glob2rx("*band*.tif$"),
                        full.names = TRUE)
# Stack bands
landsat_5 <- stack(landsat_5) %>% brick()  

# Calculate NDVI 
landsat_5_ndvi <- (landsat_5[[5]] - landsat_5[[4]]) / (landsat_5[[5]] + landsat_5[[4]])
hist(freq(landsat_5_ndvi))

# Reclassify below zero
landsat_5_ndvi <- reclassify(landsat_5_ndvi, cbind(-Inf, 0, NA), right=FALSE)

# Convert values above 1 to NA
landsat_5_ndvi <- reclassify(landsat_5_ndvi, cbind(1, 961, NA), right=TRUE)



### 2018_09_26  ---
landsat_6 <- list.files("Data/Landsat/2018/LC080470262018092601T1-SC20200204205031.tar", 
                        pattern = glob2rx("*band*.tif$"),
                        full.names = TRUE)
# Stack bands
landsat_6 <- stack(landsat_6) %>% brick()  

# Calculate NDVI 
landsat_6_ndvi <- (landsat_6[[5]] - landsat_6[[4]]) / (landsat_6[[5]] + landsat_6[[4]])
hist(freq(landsat_6_ndvi))

# Reclassify below zero
landsat_6_ndvi <- reclassify(landsat_6_ndvi, cbind(-Inf, 0, NA), right=FALSE)

# Convert values above 1 to NA
landsat_6_ndvi <- reclassify(landsat_6_ndvi, cbind(1, 961, NA), right=TRUE)



### Part 4 - Crop Images, Stack images, Calculate Mean  ------

mapview(extent(landsat_1_ndvi) %>% as("SpatialPolygons")) + 
  mapview(extent(landsat_2_ndvi) %>% as("SpatialPolygons")) + 
  mapview(extent(landsat_3_ndvi) %>% as("SpatialPolygons")) + 
  mapview(extent(landsat_4_ndvi) %>% as("SpatialPolygons")) + 
  mapview(extent(landsat_5_ndvi) %>% as("SpatialPolygons")) + 
  mapview(extent(landsat_6_ndvi) %>% as("SpatialPolygons"))
# These extents do not line up - need to crop to study area first

# Create study area
study_area_small <- mapview(landsat_1_ndvi) %>% editMap()
# editMap lets you draw on an html map interactivey, and create a polygon for your study site
study_area_small$finished %>% mapview()

# Create a study_area from your mapedit polygon
study_area <- study_area_small$finished
mapview(study_area)

# Reproject study area to raster projection
crs_landsat = " +proj=utm +zone=10 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
study_area <- st_transform(study_area, crs = crs_landsat)

## Crop rasters --
study_area_sp <- study_area %>% as("Spatial")

landsat_1_ndvi_crop <- crop(landsat_1_ndvi, study_area_sp)
mapview(landsat_1_ndvi_crop)
### Note - have to crop 0 values from the coastline to make the focal statistics work
landsat_2_ndvi_crop <- crop(landsat_2_ndvi, study_area_sp)
landsat_3_ndvi_crop <- crop(landsat_3_ndvi, study_area_sp)
landsat_4_ndvi_crop <- crop(landsat_4_ndvi, study_area_sp)
landsat_5_ndvi_crop <- crop(landsat_5_ndvi, study_area_sp)
landsat_6_ndvi_crop <- crop(landsat_6_ndvi, study_area_sp)

## Stack rasters --
ndvi_stack <- stack(landsat_1_ndvi_crop, landsat_2_ndvi_crop, landsat_3_ndvi_crop, landsat_4_ndvi_crop, landsat_5_ndvi_crop, landsat_6_ndvi_crop)
mapview(ndvi_stack$layer.6)  # plot these and make sure they're different


## Calculate mean --
ndvi_mean <- calc(ndvi_stack, fun = mean, na.rm = T)
mapview(ndvi_mean)
# Substitude SUM for mean if you want, or use other raster math 

writeRaster(ndvi_mean, "Data/ndvi_mean.tif", "GTiff")
# Write final raster