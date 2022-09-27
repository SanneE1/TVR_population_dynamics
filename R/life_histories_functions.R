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

# calculate life history traits
single_row_lh_traits <- function(ii, vr_df) {
  
  df <- vr_df[ii,] %>%
    add_column()
  
  mat <- create_mpm_from_vr(df)
  matU <- mat
  matU[1,2] <- 0
  matR <- matrix(0,2,2)
  matR[1,2] <- mat[1,2]
  
  
  df <- df %>% add_column(
    damp.ratio = damping.ratio(mat),
    iteroparity = entropy_d(lx = mpm_to_lx(matU), mx = mpm_to_mx(matU, matR)),
    life.expect = life_expect_mean(matU)
  ) 
  
  df <- do.call(data.frame, lapply(df,
                                   function(x) replace(x, is.infinite(x), NA)))
  
    return(df)
  
}


df_lh_traits <- function(vr_df) {
  
  lapply(as.list(seq_len(nrow(vr_df))), function(x) single_row_lh_traits(x, vr_df = vr_df)) %>% bind_rows
  
}




# Calculate the maximum possible coefficient of variation for a probability (Morris & Doak 2004)
CVmax <- function(mean) {
  sqrt((1-mean)/mean)
}

