library(dplyr)
library(faux)
library(parallel)
library(pbapply)
library(patchwork)

# Transforms all values below/above limits in min/max size


# Survival function
s_z <- function(z, pars, clim) {
  x <- pars$s_int + pars$s_size * log(z) + pars$s_temp * clim$surv
  p <- 1/(1+exp(-(x)))

  return(p)
}

# Growth function
g_z1z <- function(z1, z, pars, clim) {


  mean <- exp(pars$g_int + pars$g_size * log(z) + pars$g_temp * clim$growth)  ## mean size next year
  p_grow <- dnorm(z1, mean = mean, sd = pars$g_sd)   ### probability distribution that you grow to size z1
  ### based on mean size and sd
  return(p_grow)
}


# Flower probability
fp_z <- function(z, pars) {

  x <- pars$fp_int + pars$fp_size * log(z)
  p <- 1/(1+exp(-(x)))

  return(p)
}
# Seed production
seed_z <- function(z, pars, species = "H") {
  
  if(species == "H") {
  n <- exp(pars$seed_int + pars$seed_size * log(z))
  } else if(species == "O") {
    n <- rnorm(1, pars$seed_int, pars$seed_size)
  }

  return(n)
}
# recruit size
fd_z1 <- function(z1, pars) {

  fd <- dnorm(z1, mean = pars$fd_int, sd = pars$fd_sd)    ## probability of size z1 recruits

  return(fd)
}

# Fecundity function
fec <- function(z1, z, pars, species) {
  fec <- fp_z(z, pars) * seed_z(z, pars, species) * pars$germ_int * fd_z1(z1, pars)
  return(fec)
}

# HEQU
# pars <- data.frame(species = "H",
#                    s_int = -0.229,
#                    s_size = 1.077,
#                    s_temp = 1.233,
#                    g_int = 0.424,
#                    g_size = 0.846,
#                    g_temp = -0.066,
#                    g_sd = 1.076,
#                    fp_int = -3.970,
#                    fp_size = 1.719,
#                    fd_int = 1.178749,
#                    fd_sd = 0.76,
#                    seed_int = -0.6652477,
#                    seed_size = 0.8747764,
#                    germ_int = 0.1262968,
#                    germ_sd = 0.2725941,
# 
#                    mat_siz = 100,
#                    L = 1,
#                    U = 115
# 
# )


# MEROW example
pars <- data.frame(species = "O",
                   s_int = -0.364,
                   s_size = 0.632,
                   s_temp = -0.4,
                   g_int = 1.497,
                   g_size = 0.0635,
                   g_temp = -0.746,
                   g_sd = 1.778,
                   fp_int = -5.113,
                   fp_size = 2.017,
                   fd_int = -3.25,
                   fd_sd = 0.77,
                   seed_int = 139.5,
                   seed_size = 10.6,
                   germ_int = -3.43,
                   germ_sd = 0.21,
                   
                   mat_siz = 100,
                   L = 1,
                   U = 115
                   
)


kernel <- function(pars, clim, species){

  n   <- pars$mat_siz
  L   <- pars$L
  U   <- pars$U
  #these are the upper and lower integration limits
  h   <- (U-L)/n                   #Bin size
  b   <- L+c(0:n)*h                #Lower boundaries of bins
  y   <- 0.5*(b[1:n]+b[2:(n+1)])   #Bins' midpoints


  # Fertility matrix
  Fmat <- matrix(0,(n),(n))
  Fmat <- outer(y, y, fec, pars, species) * h


  # Survival matrix
  smat <- matrix(0,(n), (n))
  smat <- s_z(y, pars, clim)

  # Growth matrix
  gmat <- matrix(0, (n), (n))
  gmat <- outer(y,y, g_z1z, pars, clim) * h


  # Growth/survival matrix
  Tmat <- gmat # place holder

  # fix eviction of offspring
  for(i in 1:(pars$mat_siz/2)) {
    gmat[1,i]<- gmat[1,i]+1-sum(gmat[,i])
    Tmat[,i]<-gmat[,i]*smat[i]
  }
  # fix eviction of large adults
  for(i in (pars$mat_siz/2+1):pars$mat_siz) {
    gmat[n,i]<-gmat[n,i]+1-sum(gmat[,i])
    Tmat[,i]<-gmat[,i]*smat[i]
  }

  # Full kernel
  k_yx <- Fmat + Tmat

  return(list(k_yx = k_yx,
              smat = smat,
              gmat = gmat,
              Tmat = Tmat,
              Fmat = Fmat))

}


#----------------------------------------------------------------------------------------------------------------
# Calculate variance of, and correlation between survival and growth matrices
#----------------------------------------------------------------------------------------------------------------

calc_corr_lambda <- function(clim_corr) {
  clim_seq <- rnorm_multi(n = 5000, mu = c(0,0), sd = 2, r = clim_corr, varnames = c("surv", "growth"))
  init_pop_vec <- runif(pars$mat_siz)


  lambdas <- rep(0, length(clim_seq$surv))
  pop_vec <- as.list(c(1:length(clim_seq$surv)))
  smats <- as.list(c(1:length(clim_seq$surv)))
  gmats <- as.list(c(1:length(clim_seq$surv)))

  pb <- txtProgressBar(min = 0, max = length(clim_seq$surv), style = 3)

  for (i in c(1:length(clim_seq$surv))) {

    if(i == 1) {
      pop = init_pop_vec
    } else {pop = pop_vec[[i-1]]}

    k <- kernel(pars, clim_seq[i,], pars$species)
    pop_vec[[i]] <- k$k_yx %*% pop
    smats[[i]] <- k$smat
    gmats[[i]] <- k$gmat

    setTxtProgressBar(pb, i)

  }

  a <- sapply(pop_vec, function(x) sum(x))
  b <- a[-1]

  c <- data.frame(T1 = a,T2 = c(b, NA))
  c$lambda <- log(c$T2 / c$T1)

  d <- c[-c(1:(0.1 * length(c$lambda))),] %>%
    filter(is.finite(lambda))


  correlations <- matrix(0, (pars$mat_siz), (pars$mat_siz))
  variance_s <- matrix(0, (pars$mat_siz), (pars$mat_siz))
  variance_g <- matrix(0, (pars$mat_siz), (pars$mat_siz))
  covariance <- matrix(0, (pars$mat_siz), (pars$mat_siz))


  pb <- txtProgressBar(min = 0, max = length(pars$mat_siz), style = 1)
  for (c in c(1:pars$mat_siz)) {
    s <- sapply(smats, function(x) x[c])
    for (r in c(1:pars$mat_siz)) {
      g <- sapply(gmats, function(x) x[r,c])

      correlations[r,c] <- cor(s,g)
      variance_s[r,c] <- var(s)
      variance_g[r,c] <- var(g)
      covariance[r,c] <- cov(s,g)

    }
    setTxtProgressBar(pb, c)

  }

  return(tibble(clim_corr = clim_corr,
                lambda = d$lambda %>% mean(na.rm = T),
                correlation = list(correlations),
                variance_s = list(variance_s),
                variance_g = list(variance_g),
                covariance = list(covariance)))

}

cores <- detectCores()
cl <- makeCluster(cores-2)
clusterExport(cl = cl, ls())
clusterEvalQ(cl, c(library(faux), library(ipmr), library(dplyr)))

lambdas <- pblapply(as.list(rep(c(-0.9,0,0.9), 15)), calc_corr_lambda, cl = cl) %>% bind_rows

stopCluster(cl)

saveRDS(lambdas, file = "results/simulations/correlation_between_s&g.R")

lambdas <- readRDS("results/simulations/correlation_between_s&g.rds")


mat_cor <- matrix(0, (pars$mat_siz), (pars$mat_siz))
mat_var_s <- matrix(0, (pars$mat_siz), (pars$mat_siz))
mat_var_g <- matrix(0, (pars$mat_siz), (pars$mat_siz))
mat_cov <- matrix(0, (pars$mat_siz), (pars$mat_siz))

correlations <- list()
variance_s <- list()
variance_g <- list()
covariance <- list()

for (i in unique(lambdas$clim_corr)) {
  df <- lambdas[which(lambdas$clim_corr == i),]
  for (c in c(1:pars$mat_siz)) {
    for (r in c(1:pars$mat_siz)) {

      mat_cor[r,c] <- mean(sapply(df$correlation, function(x) x[r,c]), na.rm = T)
      mat_var_s[r,c] <- mean(sapply(df$variance_s, function(x) x[r,c]), na.rm = T)
      mat_var_g[r,c] <- mean(sapply(df$variance_g, function(x) x[r,c]), na.rm = T)
      mat_cov[r,c] <- mean(sapply(df$covariance, function(x) x[r,c]), na.rm = T)

    }

  }

  correlations <- append(correlations, list(mat_cor))
  variance_s <- append(variance_s, list(mat_var_s))
  variance_g <- append(variance_g, list(mat_var_g))
  covariance <- append(covariance, list(mat_cov))

}


names(correlations) <- c("neg", "zero", "pos")
names(variance_s) <- c("neg", "zero", "pos")
names(variance_g) <- c("neg", "zero", "pos")
names(covariance) <- c("neg", "zero", "pos")

saveRDS(correlations, file = "results/simulations/correlation_between_s&g_correlation.rds")
save(variance_s, variance_g, file = "results/simulations/correlation_between_s&g_variance.RData")
saveRDS(covariance, file = "results/simulations/correlation_between_s&g_covariance.rds")


#--------------------------------------------------------
# Plot correlation between survival and growth matrices
#--------------------------------------------------------

df <- data.frame("Climate correlation" = c("-0.9", "0", "0.9"),
                 "mean correlation s/g" = c(mean(correlations[[1]], na.rm = T),
                                            mean(correlations[[2]], na.rm = T),
                                            mean(correlations[[3]], na.rm = T)),
                 "sd" = c(sd(correlations[[1]], na.rm = T),
                          sd(correlations[[2]], na.rm = T),
                          sd(correlations[[3]], na.rm = T)))

a <- lapply(correlations, function(x) x %>% 
              as_tibble %>% 
              tibble::rowid_to_column(var = "X") %>% 
              gather(key = "Y", value = "corr", -1) %>%
              mutate(Y = as.numeric(gsub("V", "", Y))) %>%
              ggplot(aes(X, Y, fill = corr)) + geom_tile() + theme_minimal() + coord_flip() + scale_x_reverse()
)

plot_cor <- ((a$neg + ggtitle("-0.9")) + 
               (a$zero + ggtitle("0")) + 
               (a$pos + ggtitle("0.9"))) / wrap_elements(gridExtra::tableGrob(df)) + 
  plot_annotation(title = "Correlation between survivial and growth kernels",
                  subtitle = "with climate covariance")

ggsave(filename = "results/Correlation_s_g_correlation_OPIM.png", plot_cor, width = 18, height = 7)



#--------------------------------------------------------
# Plot variance of survival and growth matrices
#--------------------------------------------------------

df_s <- data.frame("Climate correlation" = c("-0.9", "0", "0.9"),
                   "mean correlation s/g" = c(mean(variance_s[[1]], na.rm = T),
                                              mean(variance_s[[2]], na.rm = T),
                                              mean(variance_s[[3]], na.rm = T)),
                   "sd" = c(sd(variance_s[[1]], na.rm = T),
                            sd(variance_s[[2]], na.rm = T),
                            sd(variance_s[[3]], na.rm = T)))

df_g <- data.frame("Climate correlation" = c("-0.9", "0", "0.9"),
                   "mean correlation s/g" = c(mean(variance_g[[1]], na.rm = T),
                                              mean(variance_g[[2]], na.rm = T),
                                              mean(variance_g[[3]], na.rm = T)),
                   "sd" = c(sd(variance_g[[1]], na.rm = T),
                            sd(variance_g[[2]], na.rm = T),
                            sd(variance_g[[3]], na.rm = T)))

b <- lapply(variance_s, function(x) x %>% 
              as_tibble %>% 
              tibble::rowid_to_column(var = "X") %>% 
              gather(key = "Y", value = "var", -1) %>%
              mutate(Y = as.numeric(gsub("V", "", Y))) %>%
              ggplot(aes(X, Y, fill = var)) + geom_tile() + theme_minimal() + coord_flip() + scale_x_reverse()
)

c <- lapply(variance_g, function(x) x %>% 
              as_tibble %>% 
              tibble::rowid_to_column(var = "X") %>% 
              gather(key = "Y", value = "var", -1) %>%
              mutate(Y = as.numeric(gsub("V", "", Y))) %>%
              ggplot(aes(X, Y, fill = var)) + geom_tile() + theme_minimal() + coord_flip() + scale_x_reverse()
)


plot_var1 <- ((b$neg + ggtitle("-0.9")) + 
                (b$zero + ggtitle("0")) + 
                (b$pos + ggtitle("0.9"))) +
  plot_annotation(title = "Variance in survivial kernel",
                  subtitle = "with climate covariance")

plot_var2 <- ((c$neg + ggtitle("-0.9")) + 
                (c$zero + ggtitle("0")) + 
                (c$pos + ggtitle("0.9"))) + plot_annotation(title = "Variance in growth kernel",
                                                            subtitle = "with climate covariance")

tables <- (wrap_elements(gridExtra::tableGrob(df_g)) + ggtitle("Survival")) + 
  plot_spacer() +
  (wrap_elements(gridExtra::tableGrob(df_s)) + ggtitle("Growth"))

ggplot2::ggsave(filename = "results/Correlation_s_g_variance_OPIM.png", 
                plot_var1 / plot_var2 / tables + plot_layout(heights = c(3,3,1)), width = 21, height =10)





#--------------------------------------------------------
# Plot covariance of survival and growth matrices
#--------------------------------------------------------

d <- lapply(covariance, function(x) x %>% 
              as_tibble %>% 
              tibble::rowid_to_column(var = "X") %>% 
              gather(key = "Y", value = "cov", -1) %>%
              mutate(Y = as.numeric(gsub("V", "", Y))) %>%
              ggplot(aes(X, Y, fill = cov)) + geom_tile() + theme_minimal() + coord_flip() + scale_x_reverse()
)

plot_cov <- ((d$neg + ggtitle("-0.9")) + 
               (d$zero + ggtitle("0")) + 
               (d$pos + ggtitle("0.9"))) + 
  plot_annotation(title = "Covariance between survivial and growth kernels",
                  subtitle = "with climate covariance")

ggsave(filename = "results/Correlation_s_g_covariance_OPIM.png", plot_cov, width = 18, height = 7)
