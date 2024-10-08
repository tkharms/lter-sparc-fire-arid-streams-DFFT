## ------------------------------------------------------- ##
  # SPARC Fire & Aridlands - Extract Climate Covariate
## ------------------------------------------------------- ##
# Written by: Nick J Lyon

# Purpose:
## Using the provided shapefile(s), extract the following data:
## PALMER DROUGHT SEVERITY INDEX (PDSI)

# Data Source
## Dai Global Palmer Drought Severity Index (PDSI)
## https://rda.ucar.edu/datasets/ds299.0/

## -------------------------------- ##
          # Housekeeping ----
## -------------------------------- ##

# Read needed libraries
# install.packages("librarian")
librarian::shelf(tidyverse, sf, ncdf4, terra, exactextractr, geojsonio, scicomptools, googledrive)

# Silence `dplyr::summarize` preemptively
options(dplyr.summarise.inform = FALSE)

# Clear environment / collect garbage
rm(list = ls()); gc()

# Identify path to location of shared data
(path <- scicomptools::wd_loc(local = F, 
                              remote_path = file.path('/', "home", "shares",
                                                      "lter-sparc-fire-arid"),
                              local_path = getwd()))

## -------------------------------- ##
        # Extraction Prep ----
## -------------------------------- ##

# Load in the catchment delineations (stored as GeoJSON)
sf_file <- geojsonio::geojson_read(x = file.path(path, "catchment-geojsons",
                                                 "aridland_fire_catchments.geojson"),
                                   what = "sp") %>%
  ## Convert to simple features object
  sf::st_as_sf(x = .)

# Check result
dplyr::glimpse(sf_file)

# Exploratory plot
plot(sf_file["usgs_site"], axes = T)

# Identify the grouping columns
(group_cols <- c(setdiff(x = names(sf_file), y = c("geometry", "geom"))))

# Define filename as an object (easier to change if/as needed)
pdsi_src <- "pdsisc.monthly.maps.1850-2018.fawc-1.r2.5x2.5.ipe-2.nc"

# Read in the netCDF file and examine for context on units / etc.
pdsi_nc <- ncdf4::nc_open(filename = file.path(path, "raw-spatial-data", "ncar_pdsi", pdsi_src))

# Look at this
print(pdsi_nc)

# Read it as a raster too
## This format is more easily manipulable for our purposes
pdsi_rast <- terra::rast(x = file.path(path, "raw-spatial-data", "ncar_pdsi", pdsi_src))

# Check names
names(pdsi_rast)

# Check out just one of those
print(pdsi_rast$sc_PDSI_pm_1)

# Visual check for overlap
plot(pdsi_rast$sc_PDSI_pm_1, axes = T, reset = F)
plot(sf_file["usgs_site"], axes = T, add = T)

## -------------------------------- ##
            # Extract ----
## -------------------------------- ##

# Create an empty list for storing extracted information
out_list <- list()

# Identify how many layers are in this
(layer_ct <- length(names(pdsi_rast)))

# We'll need to strip each layer separately
for(k in 1:layer_ct){

  # Build name of layer
  focal_layer <- paste0("sc_PDSI_pm_", k)
  
  # Starting message
  message("Extraction begun for '", focal_layer, "'")
  
  # Identify time of this layer
  layer_time <- terra::time(x = pdsi_rast[[focal_layer]])
  
  # Strip out the relevant bit
  small_out_df <- exactextractr::exact_extract(x = pdsi_rast[[focal_layer]], 
                                               y = sf_file,
                                               include_cols = group_cols,
                                               progress = T) %>%
    # Above returns a list so switch it to a dataframe
    purrr::list_rbind(x = .) %>%
    # Filter out NAs
    dplyr::filter(!is.na(value)) %>%
    # Average within existing groups
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::summarize(value_avg = mean(value, na.rm = T)) %>%
    dplyr::ungroup() %>%
    # Add a column for what timestamp this is
    dplyr::mutate(time = layer_time, 
                  layer_num = k,
                  .before = value_avg)
  
  # Add it to the list
  out_list[[focal_layer]] <- small_out_df
  
  # Success message
  message("Processing complete for ", layer_time, " (number ", k, " of ", layer_ct, ")") }

## -------------------------------- ##
            # Wrangle ----
## -------------------------------- ##

# Unlist the output of that loop for easier wrangling
pdsi_v1 <- purrr::list_rbind(x = out_list) %>% 
  # Rename extracted column more intuitively
  dplyr::rename(ncar_pdsi = value_avg)

# Check structure
dplyr::glimpse(pdsi_v1)

# Wrangle time information too
pdsi_v2 <- pdsi_v1 %>% 
  # Time is just the year so rename column
  dplyr::rename(year = time) %>% 
  # Remove unwanted years (for faster/easier processing from here on out)
  dplyr::filter(year >= 1978) %>% 
  # Exactly 12 layers / year means that's how month is stored
  dplyr::group_by(dplyr::across(dplyr::all_of(group_cols)), year) %>% 
  dplyr::mutate(month = dplyr::row_number(),
                .after = year) %>% 
  dplyr::ungroup() %>% 
  # Get month names too in case they are wanted
  dplyr::mutate(month_name = dplyr::case_when(
    month == 1 ~ "jan", month == 2 ~ "feb", month == 3 ~ "mar",
    month == 4 ~ "apr", month == 5 ~ "may", month == 6 ~ "jun",
    month == 7 ~ "jul", month == 8 ~ "aug", month == 9 ~ "sep",
    month == 10 ~ "oct", month == 11 ~ "nov", month == 12 ~ "dec",
    T ~ 'x'), .after = month) %>% 
  # Drop layer number column
  dplyr::select(-layer_num)

# Make sure there are no non-month numbers
pdsi_v2 %>% 
  dplyr::filter(month_name == "x") %>% 
  dplyr::select(year, month) %>% 
  dplyr::distinct()

# Check structure again
dplyr::glimpse(pdsi_v2)

## -------------------------------- ##
              # Export ----
## -------------------------------- ##

# Pick final object name
final_pdsi <- pdsi_v2

# Create folder to export to
dir.create(path = file.path(path, "extracted-data"), showWarnings = F)

# Define file path for CSV
pdsi_path <- file.path(path, "extracted-data", "fire-arid_pdsi.csv")

# Export the summarized data
write.csv(x = final_pdsi, na = '', row.names = F, file = pdsi_path)

# Upload to GoogleDrive
googledrive::drive_upload(media = pdsi_path, overwrite = T,
                          path = googledrive::as_id("https://drive.google.com/drive/u/0/folders/1XxvY56h1cMmaYatF7WhVrbYbaOgdRBGC"))

# End ----

