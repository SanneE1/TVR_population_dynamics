


output_dir <- "results/05_COMPADRE_studies/"
df <- read.csv(file.path(output_dir, "species_information.csv"))

files <- list.files(file.path(output_dir, "rds"), full.names = T)

x <- files[[1]]

readRDS(x)



##-----------------------------------------------------------------------------------------
t_df <- df %>% select(SpeciesAuthor, MatrixPopulation, R0, La, Gen_time, tmean_auto, contains("Tmean", ignore.case = F)) %>%
  pivot_longer(cols = contains("Tmean", ignore.case = F), values_to = "lambda") %>%
  mutate(CD = stringr::str_split(name, pattern = "_", simplify = T)[,1],
         sim = stringr::str_split(name, pattern = "_", simplify = T)[,2],
         sig = stringr::str_split(name, pattern = "_", simplify = T)[,3]) %>%
  select(-name) %>% 
  pivot_wider(values_from = lambda, names_from = sim) %>%
  rename(control = c, MCD = U) %>%
  mutate(rel_MCD = (MCD-control)/abs(control)) %>%
  filter(!is.na(control))

p_df <- df %>% select(SpeciesAuthor, MatrixPopulation, R0, La, Gen_time, prec_auto, contains("Pmean", ignore.case = F)) %>%
  pivot_longer(cols = contains("Pmean", ignore.case = F), values_to = "lambda") %>%
  mutate(CD = stringr::str_split(name, pattern = "_", simplify = T)[,1],
         sim = stringr::str_split(name, pattern = "_", simplify = T)[,2],
         sig = stringr::str_split(name, pattern = "_", simplify = T)[,3]) %>%
  select(-name) %>% 
  pivot_wider(values_from = lambda, names_from = sim) %>% 
  rename(control = c, MCD = U) %>%
  mutate(rel_MCD = (MCD-control)/abs(control)) %>% 
  filter(!is.na(control) & !is.infinite(rel_MCD))


summary(glm(rel_MCD ~ sig + tmean_auto + R0 + La + Gen_time, data = t_df ))
summary(glm(rel_MCD ~ sig + tmean_auto + La + Gen_time, data = t_df ))
summary(glm(rel_MCD ~ sig + tmean_auto + Gen_time, data = t_df ))
summary(glm(rel_MCD ~ sig + tmean_auto, data = t_df ))


summary(glm(rel_MCD ~ sig + prec_auto + R0 + La + Gen_time, data = p_df ))
summary(glm(rel_MCD ~ sig + R0 + La + Gen_time, data = p_df ))
summary(glm(rel_MCD ~ sig + R0 + La, data = p_df ))
summary(glm(rel_MCD ~ sig + La, data = p_df ))
summary(glm(rel_MCD ~ sig, data = p_df ))



##-----------------------------------------------------------------------------------------
ggplot(t_df) +
  geom_point(aes(x = tmean_auto, y = rel_MCD)) + 
  facet_grid(rows = vars(sig))




##-----------------------------------------------------------------------------------------
ggplot(p_df) +
  geom_point(aes(x = prec_auto, y = rel_MCD, colour = sig)) + 
  scale_colour_manual(values = RColorBrewer::brewer.pal(6, "Blues")[3:6])

ggplot(p_df) +
  geom_point(aes(x = R0, y = rel_MCD, colour = sig)) + 
  scale_colour_manual(values = RColorBrewer::brewer.pal(6, "Blues")[3:6])

ggplot(p_df) +
  geom_point(aes(x = La, y = rel_MCD, colour = sig)) + 
  scale_colour_manual(values = RColorBrewer::brewer.pal(6, "Blues")[3:6])

ggplot(p_df) +
  geom_point(aes(x = Gen_time, y = rel_MCD, colour = sig)) + 
  scale_colour_manual(values = RColorBrewer::brewer.pal(6, "Blues")[3:6])

