# Create dataframe for full range of lifehistories, according to Koons et al., 2016

library(tidyverse)

source("R/life_histories_functions.R")


# Create dataframe with sequences of juvenile/adult survival (Sj, Sa) & J->A transition probability (gamma)
# Then solve fecundity (rho) so that the 2x2 matrix has lambda ~ 1 
# Using the formula of determinant, solve for rho. Otto and Day, Eq. P2.6 at page 238
mat_df <- expand.grid(
  Sj = seq(0.05, 0.95, 0.15),
  Sa = seq(0.05, 0.95, 0.15),
  gamma = c(0.2, 0.5, 0.8)
) %>%    
  mutate( rho = -(Sj*(1-gamma) + Sa - 1 - (Sj*(1-gamma)*Sa)) / (Sj*(gamma)) )


# Calculate proportional measure of buffering (tau) as defined by Koons et al.
df <- df_elas_per_vr(mat_df) %>% 
  rowwise() %>%
  mutate(
    tau_Sj = (1-elas_Sj)/max(elas_Sj, elas_Sa, elas_rho),
    tau_Sa = (1-elas_Sa)/max(elas_Sj, elas_Sa, elas_rho),
    tau_gamma = (1-elas_gamma)/max(elas_Sj, elas_Sa, elas_rho),
    tau_rho = (1-elas_rho)/max(elas_Sj, elas_Sa, elas_rho)
      ) %>%
# Calculate variance for each vital rate
  mutate(
    Sj_sd = tau_Sj * 0.75 * CVmax(Sj) * Sj,
    Sa_sd = tau_Sa * 0.75 * CVmax(Sa) * Sa,
    rho_sd = tau_rho * 0.75 * 1 * rho
  ) %>%
  select(-c(contains("elas"), contains("tau"))) %>%
  rowid_to_column(var = "lh_id")

## save df
write.csv(df, file = "results/life_histories_df.csv", row.names = F)


