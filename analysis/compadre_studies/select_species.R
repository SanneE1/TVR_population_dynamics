library(tidyverse)
library(popbio)
library(Rage)

load("data/COMPADRE_v.X.X.X.4.RData")


# With these criteria there are no annual species
multiple_Ind <- compadre$metadata %>%
  filter(MatrixComposite == "Individual" & 
           MatrixSplit == "Divided" &
           MatrixFec == "Yes") %>%
  select(SpeciesAuthor, MatrixPopulation, OrganismType, Lat, Lon) %>%
  group_by(SpeciesAuthor, MatrixPopulation) %>%
  summarise(OrganismType = OrganismType,
            Lat = Lat,
            Lon = Lon,
            n = n())  %>%
  distinct() %>%
  filter(n >= 10)

# Only more than 10 individual matrices for annual species if you do not group per population
# per population; 
#             Setaria_faberi - 3 pops - 8 yrs each
#             Hordeum_spontaneum - 1 pop - 4 yrs
# All others are only 2 yrs

annual <- compadre$metadata %>%
  filter(OrganismType == "Annual",
         MatrixSplit == "Divided",
         MatrixFec == "Yes",
         MatrixComposite == "Individual") %>%
  select(SpeciesAuthor, MatrixPopulation, OrganismType, Lat, Lon, MatrixComposite) %>%
  group_by(SpeciesAuthor) %>%
  summarise(OrganismType = OrganismType,
            Lat = Lat,
            Lon = Lon,
            n = n())  %>%
  distinct() %>%
  filter(n >= 10)

multiple_Ind <- rbind(multiple_Ind, annual)

#-----------------------------------------------------------
# Gather Life histories
#-----------------------------------------------------------
get_all <- function(SpeciesAuthor, MatrixPopulation) {
  if(any(
    compadre$metadata$OrganismType[which(compadre$metadata$SpeciesAuthor == SpeciesAuthor)] == "Annual")) {
    id <- which(compadre$metadata$SpeciesAuthor == SpeciesAuthor & 
                  compadre$metadata$MatrixComposite == "Individual")
  
    }else{
  id <- which(compadre$metadata$SpeciesAuthor == SpeciesAuthor & 
                compadre$metadata$MatrixPopulation == MatrixPopulation &
                compadre$metadata$MatrixComposite == "Individual")
    }
  
  # Retrieve all full matrices
  Amats <- suppressMessages(
    lapply(as.list(id), function(x) as.vector(compadre$mat[x][[1]]$matA)) %>% bind_cols)
  ## matrix cell position in vector -> c([1,1], [2,1], etc., [1,2], [2,2], etc.)
  Umats <- suppressMessages(
    lapply(as.list(id), function(x) as.vector(compadre$mat[x][[1]]$matU)) %>% bind_cols) 
  Fmats <- suppressMessages(
    lapply(as.list(id), function(x) as.vector(compadre$mat[x][[1]]$matF)) %>% bind_cols)
  
  dim <- unique(compadre$metadata$MatrixDimension[id])
  
  if(length(dim) == 1){
    
    Acell_values <- Amats %>%
      summarise(mean = apply(., 1, mean),
                sd = apply(., 1, sd))%>%
      mutate(across(everything(), ~ replace(., is.na(.), 0)))
    
    Ucell_values <- (Umats) %>%
      mutate(across(everything(), ~ case_when(. > 6 ~ 6,
                                              . < -6 ~ -6,
                                              between(., -6, 6) ~ .))) %>%
      summarise(mean = apply(., 1, mean),
                sd = apply(., 1, sd)) %>%
      mutate(across(everything(), ~ replace(., is.na(.), 0)))
    
    Fcell_values <- (Fmats) %>%
      mutate(across(everything(), ~ case_when(. < -12 ~ -12,
                                              . >= -12 ~ .))) %>%
      summarise(mean = apply(., 1, mean),
                sd = apply(., 1, sd))%>%
      mutate(across(everything(), ~ replace(., is.na(.), 0)))
    
    meanA <- matrix(Acell_values$mean, nrow = dim, ncol = dim)
    meanU <- matrix(Ucell_values$mean, nrow = dim, ncol = dim)
    meanF <- matrix(Fcell_values$mean, nrow = dim, ncol = dim)
    
    lifetable <- makeLifeTable(meanU, meanF, nSteps = 50) %>% drop_na()
    
    R0 <- R0(meanU, meanF)                                               ## Net reproductive value
    La <- lifeTimeRepEvents(meanU, meanF)$La                             ## mean age at maturity
    Gen_time <- sum(lifetable$lxmx * lifetable$x) / sum(lifetable$lxmx)  ## Generation time
    
  }else{
    
    R0 <- NA                                              
    La <- NA
    Gen_time <- NA
    
    print(paste("Different sized matrices in same study. Species:", SpeciesAuthor, "population:", MatrixPopulation))
  
    }
  
    return(data.frame(SpeciesAuthor = SpeciesAuthor,
                      MatrixPopulation = MatrixPopulation,
                      R0 = R0,
                      La = La,
                      Gen_time = Gen_time))
  
  
  
}

LH <- map2(multiple_Ind$SpeciesAuthor, multiple_Ind$MatrixPopulation, get_all) %>% bind_rows

species <- left_join(multiple_Ind, LH)

write.csv(species, "data/species_authors.csv", row.names = F)


