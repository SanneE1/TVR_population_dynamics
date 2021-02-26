## Create file with download links for all required CHELSA files 
## Will be used in the wget download command
## Customization:
## Change chelsa_df for what you need 
## Check that download adress (L19) is still correct
## Change location of text file output location in L36

library(dplyr)

# what do I need from CHELSA?
chelsa_df <- expand.grid( variable = c('prec', 'tmax', 'tmin'),
                          year     = c(1990:2010),
                          month    = c(paste0(1:9),
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

file_names <- sapply( c(1:nrow(chelsa_df)), produce_file_name)
file_links <- paste0( read_dir, file_names)

write.table(file_links, "analysis/compadre_studies/CHELSA_files.txt", row.names = F, col.names = F, quote = F)

