# automatic download of CHELSA's climate data
# this need be run on the cluster
# 1. Set up lat/lon data
# 2. Set up variables to download chelsa data directly
# 3. Extract chelsa data
library(dplyr)
library(tidyr)
library(raster)
library(RCurl)

# arguments from command line
args     <- commandArgs(trailingOnly = TRUE)

# climate variable
job_n    <- 1
out_file <- args[1]

# pipe-able Reduce_rbind and grep functions
rbind_l     <- function(df_list){ Reduce(function(...) rbind(...), df_list) }


# 1. set up lat/lon data ------------------------------------------------------

# coord
coord_df <- read.csv('data/lagged/species_authors.csv',
                     stringsAsFactors = F) %>% 
              dplyr::select( Latitude, Longitude ) %>% 
              mutate( Longitude = as.numeric(Longitude),
                      Latitude  = as.numeric(Latitude) )

# site coord MATRIX (to feed "raster")
site_coord <-  matrix(c(coord_df$Longitude,
                        coord_df$Latitude),
                      dimnames = list(rep('value',nrow(coord_df)),
                                      c('Long', 'Lat')),
                      byrow = FALSE, nrow = nrow(coord_df) )


# 2. Set up variables to download chelsa data directly --------------------------

# what do I need from CHELSA?
chelsa_df <- expand.grid( variable = c('prec', 'tmax', 'tmin'),
                          year     = c(1901:2019),
                          month    = c(paste0(1:9),
                          # month    = c(paste0('0',1:9),
                                       10,11,12),
                          stringsAsFactors = F) %>%
                arrange(variable, year, month)

# set up reading
read_dir  <- 'https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V1/chelsa_cruts/'

# produce file name based on index 'ii'
produce_file_name <- function(ii){

  paste0(chelsa_df$variable[ii],
         '/CHELSAcruts_',
         chelsa_df$variable[ii],'_',
         chelsa_df$month[ii],'_',
         chelsa_df$year[ii],'_',
         'V.1.0.tif')

}


# create a directory to download file in
set.seed( job_n )
pp_v <- runif(1, 1, 200) %>% paste0()
dir.create( pp_v )

# get all file links (from file name)
file_names <- sapply( job_n, produce_file_name)
file_links <- paste0( read_dir, file_names)
file_dest  <- paste0( pp_v, gsub('^prec','',file_names ) )


# 3. Extract CHELSA DATA ------------------------------------------------------------------------

# download file
file_path <- file_links

# get the size of files BEFORE downloading them
get_file_size <- function( file_path ){
  
  getURL(file_path, nobody=1L, header=1L) %>% 
    strsplit(., "\r\n") %>% 
    .[[1]] %>% 
    grep( 'Content-Length: ',. ,value=T) %>% 
    gsub( 'Content-Length: ','',. ) %>% 
    as.numeric()
  
} 


# keep downloading until the file is actually downloaded
test_file <- NULL
file_s_t  <- TRUE


while( length(test_file) == 0 & file_s_t ){
    
  file_siz <- get_file_size(file_path)
  
  # tryCatch does not end the process upon error
  tryCatch(download.file(file_path, 
                         file_dest, 
                         mode = "wb"), 
           error = function(e) print('http error 500 I assume')
           )
           
  # test whether the file was downloaded or not
  test_file <- grepl('.tif', list.files( pp_v ) )

  # test whether downloaded file is large enough
  file_n    <- grep('.tif', list.files( pp_v ), value = TRUE )
  file_d_s  <- file.info( paste0( pp_v, '/', file_n) )$size
  file_s_t  <- file_d_s != file_siz
  
  # if file downloaded only partially, remove it and start over
  if( file_s_t ){
    file_r   <- grep('.tif$',list.files( pp_v ),value=T)
    file.remove( paste0(pp_v, '/', file_r) )
  }
  
}


# get climate information

# read raster file
print( 'read.raster.file' )
raster_file <- grep('.tif',list.files( pp_v ), value=T)
if( length(raster_file) > 1 ) warning( "More than one download in this dir" )
repP        <- raster( paste0( pp_v, '/', raster_file ) )

# extract info from raster file
print( 'extract.raster.file' )
values_clim <- raster::extract(repP, site_coord, method = 'bilinear')
clim_df     <- coord_df %>%
                  mutate( variable = chelsa_df$variable[job_n],
                          year     = chelsa_df$year[job_n],
                          month    = chelsa_df$month[job_n],
                          value    = values_clim)

# remove file
print( 'remove.download' )
file_r   <- grep('.tif$',list.files( pp_v ),value=T)
file.remove( paste0(pp_v, '/', file_r) )

# remove temp directory
file.remove( pp_v, recursive=T )


# print(clim_df)
# paste0(out_file, '_slice_',job_n,'.csv')

# put it all out
write.csv(clim_df,
          paste0(out_file, '_slice_',job_n,'.csv'),
          row.names=F)
