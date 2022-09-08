
library(tidyverse)
library(popbio)

source("R/life_histories_functions.R")



# big_same <- lapply(as.list(c(0.05,0.25,0.5,1)), function(x) {
#   read.csv(file.path("results/01_Simulations_mpm_same_directions/", 
#                      paste("life_histories_lambda_same_dir_", x, ".csv", sep = ""))) %>%
#     mutate(sig.strength = x)
# }
# ) %>% bind_rows %>%
#   df_lh_traits
# 
# big_opp <- lapply(as.list(c(0.05,0.25,0.5,1)), function(x) {
#   read.csv(file.path("results/02_Simulations_mpm_opposing_directions",
#                      paste("life_histories_lambda_diff_dir_", x, ".csv", sep = ""))) %>%
#     mutate(sig.strength = x)
# }
# ) %>% bind_rows %>%
#   df_lh_traits

big_same <- read.csv("results/01_Simulations_mpm_same_directions/life_histories_all.csv") %>%
  mutate(auto_cat = cut(clim_auto, breaks = 3, labels = c("-0.9", "0", "0.9")),
         vr_cov = "positive") %>%
  group_by(lh_id, clim_sd, auto_cat, sig.strength, lag_type) %>% 
  mutate(rep = row_number()) %>%
  pivot_wider(names_from = lag_type, values_from = lambda) %>%
  mutate(l_diff = Umatrix - none,
         rel_dif = l_diff/abs(none))


big_opp <- read.csv("results/02_Simulations_mpm_opposing_directions/life_histories_all.csv") %>%
  mutate(auto_cat = cut(clim_auto, breaks = 3, labels = c("-0.9", "0", "0.9")),
         vr_cov = "negative") %>%
  group_by(lh_id, clim_sd, auto_cat, sig.strength, lag_type) %>% 
  mutate(rep = row_number()) %>%
  pivot_wider(names_from = lag_type, values_from = lambda) %>%
  mutate(l_diff = Umatrix - none,
         rel_dif = l_diff/abs(none))

# anti_join(big_same %>% select(lh_id, clim_sd, auto_cat, sig.strength, rep),
#           big_opp %>% select(lh_id, clim_sd, auto_cat, sig.strength, rep))

df <- rbind(big_same, big_opp) 

write.csv(df, "results/lambdas_life_histories.csv", row.names = F)

df <- df %>%
  mutate(auto = as.numeric(levels(auto_cat))[auto_cat],
         elas_Sj = scale(elas_Sj), 
         elas_Sa = scale(elas_Sa), 
         elas_gamma = scale(elas_gamma))

ggplot(df) + 
  geom_point(aes(x = clim_sd, y = rel_dif, colour = auto_cat), position = position_dodge(width = 0.1)) +
  facet_grid(rows = vars(vr_cov))



mod_diff <- lm(l_diff ~ gen.time + damp.ratio + 
                 elas_Sj + elas_Sa + elas_gamma + elas_rho +
                 clim_sd + auto + sig.strength + vr_cov + 
                 clim_sd:auto + clim_sd:vr_cov + auto:vr_cov,
               data = df)
summary(mod_diff)

drop1(mod_diff, test = "Chisq") %>% 
  rownames_to_column("dropped_covariate") %>%
  rename(p_value = `Pr(>Chi)`) %>%
  mutate(
    p_adj_BH = p.adjust(p_value,method="BH",n=96),
    p_adj_holm = p.adjust(p_value,method="holm",n=96),
    include = ifelse(p_value < 0.05 & p_adj_BH < 0.05 & p_adj_holm < 0.05, "yes", "no")
  ) %>% as_tibble()


mod_rell <- lm(rel_dif ~ gen.time + damp.ratio + 
            elas_Sj + elas_Sa + elas_gamma + elas_rho +
            clim_sd + auto + sig.strength + vr_cov + 
              clim_sd:auto + clim_sd:vr_cov + auto:vr_cov,
          data = df)
summary(mod_rell)

drop1(mod_rell, test = "Chisq") %>% 
  rownames_to_column("dropped_covariate") %>%
  rename(p_value = `Pr(>Chi)`) %>%
  mutate(
    p_adj_BH = p.adjust(p_value,method="BH",n=96),
    p_adj_holm = p.adjust(p_value,method="holm",n=96),
    include = ifelse(p_value < 0.05 & p_adj_BH < 0.05 & p_adj_holm < 0.05, "yes", "no")
  ) %>% as_tibble()

mod <- lm(l_diff ~ clim_sd + auto + sig.strength + vr_cov +
            clim_sd:auto + clim_sd:vr_cov + auto:vr_cov +
            clim_sd:auto:vr_cov,
          data = df)
summary(mod)

mod <- lm(l_diff ~ clim_sd*auto*sig.strength*vr_cov,
          data = df)


