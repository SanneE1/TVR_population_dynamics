library(dplyr)

## Get info from command line

args     <- commandArgs(trailingOnly = TRUE)

other_job <- args[1]

read_dir <- paste0("mpms_climate-", other_job)


## Read all month/year seperate csv and create one csv
files <- list.files(read_dir, pattern = "All_populations_")
  
all <- lapply(files, read.csv) %>% bind_rows %>% mutate(SpeciesAuthor = factor(SpeciesAuthor),
                                                        MatrixPopulation = factor(MatrixPopulation))

## Split up main file into seperate per population
per_pop <- split(all, list(all$SpeciesAuthor, all$MatrixPopulation), drop = T)

## Save per population csv
lapply(per_pop, function(x) write.csv(x, file = file.path(read_dir, paste0("climate_", names(x), ".csv")), row.names = F))



