
library(ipmr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)


params_list <- list(
  s_int = -0.229,
  s_slope = 1.077,
  # s_temp = 1.233,
  g_int = 0.424,
  g_slope = 0.846,
  # g_temp = -0.066,
  g_sd = 1.076,
  fp_int = -3.970,
  fp_slope = 1.719,
  fpC_int = -3.859385,
  fpC_slope = 1.719768,
  fpC_temp = -0.6169492,
  fn_int = -0.6652477,
  fn_slope = 0.7048809,
  fnC_int = -0.5661762,
  fnC_slope = 0.7048782,
  fnC_temp = -0.3398345,
  germ_mean = 0.1262968,
  germ_sd = 0.2725941,
  fd_mean = 1.178749,
  fd_sd = 0.8747764
)


source("analysis/simulations/ipmr_functions.R")


clim_sd <- rep(seq(from = 0, to = 2, length.out = 10), 30)
clim_corr <- rep(rep(c(-0.9,0,0.9), each = 10), 10)

n_corr_sum <- rep(NA, length(clim_sd))
s_corr_sum <- rep(NA, length(clim_sd))
g_corr_sum <- rep(NA, length(clim_sd))

n_corr_hist <- as.list(rep(NA, length(clim_sd)))
s_corr_hist <- as.list(rep(NA, length(clim_sd)))
g_corr_hist <- as.list(rep(NA, length(clim_sd)))

n_lambda <- rep(NA, length(clim_sd))
s_lambda <- rep(NA, length(clim_sd))
g_lambda <- rep(NA, length(clim_sd))

start <- Sys.time()
start

for (n in c(1:length(clim_sd))) {
  
a <-P_lambdas(n_it = 1000, 
          clim_sd = clim_sd[n], 
          clim_corr = clim_corr[n], 
          params_list = params_list, 
          clim_params = list(s_temp = 1.233, 
                             g_temp = 0.066),
          n_mesh = 100,
          save_K = T)

b <- lapply(a$M_non_lagged[[1]][c((0.8*n_it):n_it)], function(x) as.vector(x)) %>% bind_rows %>% t
c <- corrr::correlate(b)
d <- as.matrix(c[,-1])
n_corr_hist[[n]] <- d
n_corr_sum[n] <- mean(d, na.rm = T)
n_lambda[n] <- a$non_lagged

e <- lapply(a$M_s_lagged[[1]][c((0.8*n_it):n_it)], function(x) as.vector(x)) %>% bind_rows %>% t
f <- corrr::correlate(e)
g <- as.matrix(f[,-1])
s_corr_hist[[n]] <- g
s_corr_sum[n] <- mean(g, na.rm = T)
s_lambda[n] <- a$lagged_s

h <- lapply(a$M_g_lagged[[1]], function(x) as.vector(x)) %>% bind_rows %>% t
i <- corrr::correlate(h)
j <- as.matrix(i[,-1])
g_corr_hist[[n]] <- j
g_corr_sum[n] <- mean(j, na.rm = T)
g_lambda[n] <- a$lagged_g

}

Sys.time()

Sys.time() - start


df <- data.frame(clim_corr, clim_sd,
            n_lambda, s_lambda, g_lambda,
            n_corr_sum, s_corr_sum, g_corr_sum)

a <- df %>% pivot_longer(cols = contains("lambda"), names_to = "type", values_to = "lambda") %>%
  select(clim_corr, clim_sd, type, lambda)
a$type <- gsub("_lambda", "", a$type)

b <- df %>% pivot_longer(cols = contains("corr_sum"), names_to = "type", values_to = "cell_corr_sum") %>%
  select(clim_corr, clim_sd, type, cell_corr_sum)
b$type <- gsub("_corr_sum", "", b$type)

data <- left_join(a,b)



l <- ggplot(data)+ 
  geom_smooth(aes(x = clim_sd,
                  y = lambda, 
                  color = type), 
              size = 1) + 
 facet_grid(rows = vars(clim_corr)) 

c <- ggplot(data) +
  geom_smooth(aes(x = clim_sd,
                  y = cell_corr_sum,
                  color = type),
              size = 1) +
  facet_grid(rows = vars(clim_corr))

l + c

x <- vector()

for (i in c(1:length(n_corr_hist))) {
  if(all(is.na(n_corr_hist[[i]])) == T) {
    x <- append(x, i)
  }
}

n_corr_hist <- n_corr_hist[-x]


n_hists <- as.list(c(1:8)) 
names(n_hists) <- clim_sd[2:9]

for(j in c(1:8)) {
  n_hists[[j]] <- hist(n_corr_hist[[j]])
}     

y <-lapply(n_hists, function(x) ggplot(as.data.frame(x[c("mids", "density")]), aes(x = mids, y = density)) + geom_bar(stat = "identity") + xlab("cell correlations"))
