# retrieve the downloaded climate values on specific lat/lon's from downloaded CHELSA files

# Script depends on three arguments from the command line: 
# 1) location of the downloaded files - The files function is recursive, so you can name the first directory 
#    where downstream there are only CHELSA files
# 2) the file with lat-lon information - formatting of the coordinates df (#1) would need changing for other use
# 3) output location for the climate csv

# MAKE SURE TO CHECK THE LOG FILES AT THE END. SCRIPT WILL PRINT FILES WHERE IT IS UNABLE TO EXTRACT DATA 
# Usually this is due to a corrupted download

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(raster))
suppressPackageStartupMessages(library(RCurl))

# arguments from command line
args     <- commandArgs(trailingOnly = TRUE)
download.location <- args[1]
lat.lon.file <- args[2]
output_dir <- args[3]

print(paste("location of downloads:", download.location))
print(paste("location of lat/lon file:", lat.lon.file))
print(paste("output_directory:", output_dir))

print(paste("work directory:", getwd()))

# 1. set up lat/lon data ------------------------------------------------------

# coord
coord_df <- read.csv(lat.lon.file,
                     stringsAsFactors = F) %>% 
  dplyr::select( Lat, Lon, SpeciesAuthor, MatrixPopulation ) %>% 
  rename( Longitude = Lon,
          Latitude  = Lat )

# site coord MATRIX (to feed "raster")
site_coord <-  matrix(c(coord_df$Longitude,
                        coord_df$Latitude),
                      dimnames = list(rep('value',nrow(coord_df)),
                                      c('Long', 'Lat')),
                      byrow = FALSE, nrow = nrow(coord_df) )


# 2. List downloaded CHELSA files using the text file used for downloading  --------------------------
files <- list.files(path = download.location, pattern = ".tif$", recursive = T, full.names = T)

# 3. Extract CHELSA DATA ------------------------------------------------------------------------

# create function to get climate information
get_clim <- function(file) {
  
  # get info on month & variable of this file
  variables <- stringr::str_split(regmatches(file, regexec("CHELSAcruts_(.+)_V.1.0.tif", file))[[1]][2], 
                                  "[[:punct:]]")
  
  # read raster file
  repP        <- raster( file,  )
  
  # extract info from raster file
  values_clim <- raster::extract(repP, site_coord, method = 'bilinear')
  
  # get all relevant data into dataframe
  clim_df     <- coord_df %>%
    mutate( variable = variables[[1]][1],
            year     = variables[[1]][3],
            month    = variables[[1]][2],
            value    = values_clim)
  
  return(clim_df)
}

climate <- lapply(as.list(files),function(x) tryCatch({get_clim(x)}, error = function(e) {message(paste("problem file:", x))}) ) %>% bind_rows()



# put it all out
write.csv(climate,
          file.path(output_dir, 'All_populations.csv'),
          row.names=F)

