library(tidyverse)
library(popbio)
library(Rage)

load("data/COMPADRE_v.6.21.1.0.RData")


# With these criteria there are no annual species
multiple_Ind <- compadre$metadata %>%
  filter(MatrixComposite == "Individual" & 
           MatrixSplit == "Divided" &
           MatrixFec == "Yes",
         MatrixTreatment == "Unmanipulated") %>%
  select(SpeciesAuthor, MatrixPopulation, StudyStart, OrganismType, Lat, Lon) %>%
  group_by(SpeciesAuthor, MatrixPopulation) %>%
  summarise(OrganismType = OrganismType,
            Start = StudyStart,
            Lat = Lat,
            Lon = Lon,
            n = n())  %>%
  distinct() %>%
  filter(n >= 10)


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
    
    R0 <- R0(meanU, meanF)                                               ## Net reproductive value
    La <- lifeTimeRepEvents(meanU, meanF)$La                             ## mean age at maturity
    Gen_time <- generation.time(meanA)                                   ## Generation time
    
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


