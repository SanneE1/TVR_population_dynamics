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
    tau_Sj = (1-elas_Sj)/max(elas_Sj, elas_Sa, elas_gamma, elas_rho),
    tau_Sa = (1-elas_Sa)/max(elas_Sj, elas_Sa, elas_gamma, elas_rho),
    tau_gamma = (1-elas_gamma)/max(elas_Sj, elas_Sa, elas_gamma, elas_rho),
    tau_rho = (1-elas_rho)/max(elas_Sj, elas_Sa, elas_gamma, elas_rho)
      ) %>%
  # mutate(
  #   tau_Sj = 1,
  #   tau_Sa = 1,
  #   tau_gamma = 1,
  #   tau_rho = 1
  # ) %>%
# Calculate variance for each vital rate
  mutate(
    Sj_var = tau_Sj * 0.5 * Vmax(Sj),
    Sa_var = tau_Sa * 0.5 * Vmax(Sa),
    gamma_var = tau_gamma * 0.5 * Vmax(gamma),
    rho_var = tau_rho * 0.5 * 1
  ) %>%
  mutate(
    sum_elas = sum(elas_Sj, elas_Sa, elas_gamma, elas_rho),
    Sj_alpha = dampack::beta_params(Sj, sqrt(Sj_var))$alpha,
    Sj_beta = dampack::beta_params(Sj, sqrt(Sj_var))$beta,
    Sa_alpha = dampack::beta_params(Sa, sqrt(Sa_var))$alpha,
    Sa_beta = dampack::beta_params(Sa, sqrt(Sa_var))$beta,
    gamma_alpha = dampack::beta_params(gamma, sqrt(gamma_var))$alpha,
    gamma_beta = dampack::beta_params(gamma, sqrt(gamma_var))$beta,
    rho_alpa = dampack::gamma_params(rho, sqrt(rho_var))$shape
    ) %>% 
  # select(-c(contains("elas"), contains("tau"))) %>%
  rowid_to_column(var = "lh_id")

## save df
write.csv(df, file = "results/life_histories_df.csv", row.names = F)


