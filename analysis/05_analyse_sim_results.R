
library(tidyverse)
library(popbio)
library(lme4)

source("R/life_histories_functions.R")

# # Combine all result files from the array jobs into one file
# s0.05 <- lapply(list.files("/work/evers/simulate_lag-20986126/", full.names = T),
#                 read.csv) %>% bind_rows() %>% mutate(sig.strength = 0.05, vr_cov = "positive")
# s0.25 <- lapply(list.files("/work/evers/simulate_lag-20986273/", full.names = T),
#                 read.csv) %>% bind_rows() %>% mutate(sig.strength = 0.25, vr_cov = "positive")
# s0.5 <- lapply(list.files("/work/evers/simulate_lag-20986350/", full.names = T),
#                read.csv) %>% bind_rows() %>% mutate(sig.strength = 0.5, vr_cov = "positive")
# s1 <- lapply(list.files("/work/evers/simulate_lag-20986426/", full.names = T),
#              read.csv) %>% bind_rows() %>% mutate(sig.strength = 1, vr_cov = "positive")
# 
# o0.05 <- lapply(list.files("/work/evers/simulate_oposing-20986717/", full.names = T),
#                 read.csv) %>% bind_rows() %>% mutate(sig.strength = 0.05, vr_cov = "negative")
# o0.25 <- lapply(list.files("/work/evers/simulate_oposing-20986794/", full.names = T),
#                 read.csv) %>% bind_rows() %>% mutate(sig.strength = 0.25, vr_cov = "negative")
# o0.5 <- lapply(list.files("/work/evers/simulate_oposing-20987011/", full.names = T),
#                read.csv) %>% bind_rows() %>% mutate(sig.strength = 0.5, vr_cov = "negative")
# o1 <- lapply(list.files("/work/evers/simulate_oposing-20987088/", full.names = T),
#              read.csv) %>% bind_rows() %>% mutate(sig.strength = 1, vr_cov = "negative")
# 
# df <- rbind(s0.05, s0.25, s0.5, s1, o0.05, o0.25, o0.5, o1)
# 
# write.csv(df, file = "results/lambdas_life_histories.csv", row.names = F)

df <- read.csv("results/lambdas_life_histories.csv") %>%
  left_join(., 
            read.csv("results/life_histories_df.csv")) %>%
  mutate(auto_cat = cut(clim_auto, breaks = 7, labels = c(-0.6, -0.3, -0.1, 0, 0.1, 0.3, 0.6)),
         life.expect = scale(log(life.expect)),
         iteroparity = scale(iteroparity)) 



# check that all simulations produced correct number of matrices
if(any(df$n_mats != 10002)) warning("some sequences do not have the correct number of matrices")

lamb_pos <- lmer(lambda ~ clim_sd + I(clim_sd^2) + clim_auto + sig.strength + 
                   lag_type + lag_type:clim_sd + lag_type:clim_auto + lag_type:sig.strength + 
            (I(clim_sd^2) - 1|lh_id), data = df %>% filter(vr_cov == "positive"))
lamb_neg <- lmer(lambda ~ clim_sd + I(clim_sd^2) + clim_auto + sig.strength + 
                   lag_type + lag_type:clim_sd + lag_type:clim_auto + lag_type:sig.strength +  
                   (I(clim_sd^2) - 1|lh_id), data = df %>% filter(vr_cov == "negative"))

summary(lamb_pos)
summary(lamb_neg)

saveRDS(object = data.frame(Estimate_pos = round(fixef(lamb_pos), 4), 
                            Estimate_neg = round(fixef(lamb_neg), 4)), 
        file = "results/table_lambda.rds")



df1 <- df %>%
  group_by(lh_id, clim_sd, auto_cat, sig.strength, lag_type, vr_cov) %>%
  mutate(rep = row_number()) %>%
  pivot_wider(names_from = lag_type, values_from = lambda) 

df2 <- df1 %>%
  select(-clim_auto) %>% 
  filter(clim_sd == 0.01 | clim_sd == 1) %>% 
  group_by(lh_id, rep, auto_cat) %>% 
  pivot_wider(names_from = "clim_sd", values_from = c("Umatrix", "none")) %>% 
  mutate(rel_decrease = (Umatrix_1 - Umatrix_0.01)/(none_1 - none_0.01),
         rel_diff = none_1 - none_0.01,
         auto_cat = as.numeric(levels(auto_cat))[auto_cat])


rel_dec_pos <- lmer(log(rel_decrease) ~ life.expect + iteroparity +
                      auto_cat + sig.strength + 
                      auto_cat:life.expect + auto_cat:iteroparity +
                      (1|lh_id), 
                    data = df2 %>% filter(vr_cov == "positive"))
rel_dec_neg <- lmer(log(rel_decrease) ~ life.expect + iteroparity +
                      auto_cat + sig.strength + 
                      auto_cat:life.expect + auto_cat:iteroparity +
                      (1|lh_id), 
                    data = df2 %>% filter(vr_cov == "negative"))

summary(rel_dec_pos)
summary(rel_dec_neg)

saveRDS(object = data.frame(Estimate_pos = round(fixef(rel_dec_pos), 4), 
                            Estimate_neg = round(fixef(rel_dec_neg), 4)), 
        file = "results/table_rel_decr.rds")

# 
# df3 <- df1 %>%
#   filter(clim_sd == 1) %>%
#   rowwise() %>%
#   mutate(diff = (abs(none) - abs(Umatrix))/abs(none),
#          auto_cat = as.numeric(levels(auto_cat))[auto_cat])
# 
# lower_bound_pos <- quantile(df3$diff[which(df3$vr_cov == "positive")], 0.01, na.rm = T)
# upper_bound_pos <- quantile(df3$diff[which(df3$vr_cov == "positive")], 0.99, na.rm = T)
# incl_id_pos <- which(!(df3$diff[which(df3$vr_cov == "positive")] < lower_bound_pos | df3$diff[which(df3$vr_cov == "negative")] > upper_bound_pos))
# 
# lower_bound_neg <- quantile(df3$diff[which(df3$vr_cov == "negative")], 0.01, na.rm = T)
# upper_bound_neg <- quantile(df3$diff[which(df3$vr_cov == "negative")], 0.99, na.rm = T)
# incl_id_neg <- which(!(df3$diff[which(df3$vr_cov == "negative")] < lower_bound_neg | df3$diff[which(df3$vr_cov == "negative")] > upper_bound_neg))
# 
# 
# diff_pos <- lmer(log(diff) ~ life.expect + iteroparity +
#                       auto_cat + sig.strength + 
#                       auto_cat:life.expect + auto_cat:iteroparity +
#                       (1|lh_id), 
#                     data = df3[incl_id_pos,])
# diff_neg <- lmer(diff ~ life.expect + iteroparity +
#                    auto_cat + sig.strength + 
#                    auto_cat:life.expect + auto_cat:iteroparity +
#                    (1|lh_id), 
#                  data = df3[incl_id_neg,])
# 
# summary(diff_pos)
# summary(diff_neg)
# 
# saveRDS(object = data.frame(Estimate_pos = round(fixef(diff_pos), 4),
#                             Estimate_neg = round(fixef(diff_neg), 4)),
#         file = "results/table_diff.rds")


