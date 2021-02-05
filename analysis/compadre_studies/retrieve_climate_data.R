# library(dplyr)
library(tidyverse)
library(raster)
library(RCurl)
library(prism)
library(pbapply)
library(lubridate)
library(popbio)

load('data/COMPADRE_v.X.X.X.4.RData')

# Select 2x2 mpm's  ---- I should probably get more than only the Dalgleish studies, but just taking all is not a good idea either

## Subsection from Dalgleish
id <- which(grepl("Dalgleish", compadre$metadata$Authors) & compadre$metadata$MatrixDimension == 2)

species <- unique(compadre$metadata$SpeciesAuthor[id])


## Start for loop through all species --------------------------

# for(i in species){
#   
#   print(i)
#   id2 <- id[which(compadre$metadata$SpeciesAuthor[id] == i)]
#   
#   # Population coordinates
#   pop <- distinct(data.frame(name = compadre$metadata$MatrixPopulation[id2],
#                              lat = compadre$metadata$Lat[id2],
#                              lon = compadre$metadata$Lon[id2],
#                              start = compadre$metadata$StudyStart[id2],
#                              end = compadre$metadata$StudyEnd[id2]))
#   
#   ## prism data
#   
#   # Build in a check to download the prism data if the years have not yet been downloaded?
#   
#   # get_prism_monthlys(type = "tmean", years = 1900:2019, mon = 1:12, keepZip = F)
#   # get_prism_monthlys(type = "ppt", years = 1900:2019, mon = 1:12, keepZip = F)
#   
#   ## years needed
#   if(length(c(pop$start:pop$end)) >= 30) {
#     yrs <-   c(pop$start:pop$end)
#   } else {
#     yrs <-   c((pop$end - 29):pop$end)
#   }
#   
#   ## appropriate file names
#   
#   files <- grep(paste(yrs, collapse = "|"), ls_prism_data()[,1], value = T)  
#   
#   climate <- prism_slice(c(pop$lon, pop$lat), files)
#   
#   clim <- climate$data %>% rownames_to_column(var = "variable") %>%
#     mutate(date = as.Date(date)) %>%
#     rename(value = data)%>% 
#     mutate(variable = stringr::str_match(variable, "PRISM_(.*?)_")[,2]) %>% 
#     pivot_wider(names_from = variable, values_from = value) %>%
#     arrange(date) %>%
#     group_by(month(date)) %>%
#     mutate(tmean = scale(tmean),
#            ppt = scale(ppt)) %>%
#     ungroup() %>%
#     dplyr::select(-'month(date)')
#   
#   write.csv(clim, file = paste0("data/climate_compadre/", i, "_climate.csv"))
#   
# }


### run regression for each cell/vitalrate with recent or lagged yearly climate anomalies
i = species[1]

clim <- read.csv(paste0("data/climate_compadre/", i, "_climate.csv")) %>%
  dplyr::select(-X) %>%
  mutate(date = as.Date(date, "%Y-%m-%d")) %>%
  group_by(Year = year(date)) %>%
  summarise(tmean = mean(tmean, na.rm = T),
            ppt = mean(ppt, na.rm = T))

clim_lag <- clim %>%
  transmute(Year = Year +1,
            tmean_lagged = tmean,
            ppt_lagged = ppt)

clim <- left_join(clim, clim_lag)


acf(clim$ppt, plot = F)
acf(clim$tmean, plot = F)


id2 <- id[which(compadre$metadata$SpeciesAuthor[id] == i & compadre$metadata$MatrixComposite[id] == "Individual")]

Amat_list <- lapply(as.list(id2), function(x) compadre$mat[x][[1]]$matA)  

# get growth, survival and fecundity values and correlate with climate      
lambdas <- data.frame(Year = compadre$metadata$MatrixEndYear[id2],
                      lambda = sapply(Amat_list, lambda))


lambdas <- left_join(lambdas, clim)

summary(lm(lambda ~ tmean, data = lambdas))
summary(lm(lambda ~ tmean_lagged, data = lambdas))

summary(lm(lambda ~ ppt, data = lambdas))
summary(lm(lambda ~ ppt_lagged, data = lambdas))



