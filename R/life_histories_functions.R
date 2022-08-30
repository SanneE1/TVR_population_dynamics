## functions for the lifehistories dataframe

# create matrix or a list of matrices from vital rates
create_mpm_from_vr <- function(df) {
  
  if(nrow(df) == 1) {
    mat <- matrix(0,2,2)
    mat[1,1] <- df$Sj*(1-df$gamma)
    mat[2,1] <- df$Sj*df$gamma
    mat[1,2] <- df$rho
    mat[2,2] <- df$Sa
  } else {
    mat <- lapply(seq_len(nrow(df)), function(x) {
      mat <- matrix(0,2,2)
      mat[1,1] <- df$Sj[x]*(1-df$gamma[x])
      mat[2,1] <- df$Sj[x]*df$gamma[x]
      mat[1,2] <- df$rho[x]
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
               elas_rho = NA)

  matrix.elements <- expression(Sj*(1-gamma),  rho,
                                Sj*gamma,      Sa)
  vr_sens_elas <- popbio::vitalsens(matrix.elements, list(Sj = df$Sj,
                                                          Sa = df$Sa,
                                                          gamma = df$gamma,
                                                          rho = df$rho) )
  df[,c(5:8)] <- t(vr_sens_elas$elasticity)
  
  return(df)
}

df_elas_per_vr <- function(vr_df) {
  
  lapply(as.list(1:nrow(vr_df)), single_row_elas) %>% bind_rows
  
}

# Calculate the maximum possible coefficient of variation for a probability (Morris & Doak 2004)
CVmax <- function(mean) {
  sqrt((1-mean)/mean)
}

Vmax <- function(mean) {
  mean*(1-mean)
}
