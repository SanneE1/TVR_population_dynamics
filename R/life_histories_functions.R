## functions for the lifehistories dataframe

# create matrix or a list of matrices from vital rates
create_mpm_from_vr <- function(df) {
  
  if(nrow(df) == 1) {
    mat <- matrix(0,2,2)
    mat[1,1] <- df$Sj*(1-df$gamma)
    mat[2,1] <- df$Sj*df$gamma
    mat[1,2] <- df$phi
    mat[2,2] <- df$Sa
  } else {
    mat <- lapply(as.list(1:nrow(df)), function(x) {
      mat <- matrix(0,2,2)
      mat[1,1] <- df$Sj[x]*(1-df$gamma[x])
      mat[2,1] <- df$Sj[x]*df$gamma[x]
      mat[1,2] <- df$phi[x]
      mat[2,2] <- df$Sa[x]
      mat
    } )
    
  }
  return(mat)
}

# calculate elasticities
single_row_elas <- function(ii) {
  df <- mat_df[ii,] %>%
    add_column(elas_Sj = NA,
               elas_Sa = NA,
               elas_gamma = NA,
               elas_phi = NA)
  pert <- 0.01
  
  for(vr in c("Sj", "Sa", "phi", "gamma")) {
    
    df_h <- df_l <- df
    df_h[,vr] <- df_h[, vr] + pert
    df_l[,vr] <- df_l[, vr] - pert
    
    mean <- mat_df[ii,] %>% create_mpm_from_vr(.) %>% popbio::lambda(.)
    high <- create_mpm_from_vr(df_h) %>% popbio::lambda(.)
    low <- create_mpm_from_vr(df_l) %>% popbio::lambda(.)
    
    sens <- (high - low) / (2 * pert)
    elas <- sens * abs(df[,vr]) / mean
    
    df[grep(paste0("elas_", vr), colnames(df))] <- elas
    
  }
  
  return(df)
}

df_elas_per_vr <- function(vr_df) {
  
  lapply(as.list(1:nrow(vr_df)), single_row_elas) %>% bind_rows
  
}

# Calculate the maximum possible coefficient of variation for a probability (Morris & Doak 2004)
CVmax <- function(mean) {
  sqrt((1-mean)/mean)
}
