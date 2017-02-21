library(pracma)
library(MASS)
library(ggplot2)
rm(list=ls()); gc(reset = TRUE)
source("iFormFunctionalMapping_OLS.R")


a <- 180
b <- 3
r <- 0.9
sigma <- 1
rho <- 0.2

m <- 10
t <- 0:9
mu <- a / (1 + b * exp(-r*t))
L <- t(Legendre(t, 4))

  
COVAR <- rho^(outer(seq(0, (m-1)), seq(0, (m-1)), function(a,b) abs(a-b)))
diag(COVAR) <- sigma^2
#COVAR <- kronecker(diag(j), COVAR)
#mvrnorm(n = 1, rep(0, 10), COVAR)

sim_list <- list()
for(i in 1:100){

snps <- matrix(rbinom(10000, 1, 0.5), nrow = 200)
snps <- data.frame(id = 1:200, snps)

time_df <- do.call(rbind, lapply(1:200, function(id) cbind(id, t)))

df <- merge(snps, time_df, by = "id", all = TRUE)

{
# snp1_effect <- rnorm(4, 2, 0.25)
# snp1 <- rep(snp1_effect %*% L, 100)
# snp2_effect <- rnorm(4, 2, 0.25)
# snp2 <- rep(snp2_effect %*% L, 100)
# snp3_effect <- rnorm(4, 2, 0.25)
# snp3 <- rep(snp3_effect %*% L, 100)
# snp4_effect <- rnorm(4, 2, 0.25)
# snp4 <- rep(snp4_effect %*% L, 100)
# snp5_effect <- rnorm(4, 2, 0.25)
# snp5 <- rep(snp5_effect %*% L, 100)
# snp6_effect <- rnorm(4, 2, 0.25)
# snp6 <- rep(snp6_effect %*% L, 100)
# snp7_effect <- rnorm(4, 2, 0.25)
# snp7 <- rep(snp7_effect %*% L, 100)


# snp1_effect <- rnorm(4, 3, 1)
# snp1 <- rep(snp1_effect %*% t(Legendre(t, 4)), 100)
# snp2_effect <- rnorm(2, 3, 1)
# snp2 <- rep(snp2_effect %*% t(Legendre(t, 2)), 100)
# snp3_effect <- rnorm(3, 3, 1)
# snp3 <- rep(snp3_effect %*% t(Legendre(t, 3)), 100)
# snp4_effect <- rnorm(2, 3, 1)
# snp4 <- rep(snp4_effect %*% t(Legendre(t, 2)), 100)
# snp5_effect <- rnorm(4, 3, 1)
# snp5 <- rep(snp5_effect %*% t(Legendre(t, 4)), 100)
# snp6_effect <- rnorm(2, 3, 1)
# snp6 <- rep(snp6_effect %*% t(Legendre(t, 2)), 100)
# snp7_effect <- rnorm(3, 3, 1)
# snp7 <- rep(snp7_effect %*% t(Legendre(t, 3)), 100)
}

snp1_effect <- rnorm(1, 3, 1)
snp1 <- snp1_effect * rep(Legendre(t, 4)[, 4], 100)
snp2_effect <- rnorm(1, 3, 1)
snp2 <- snp2_effect * rep(Legendre(t, 2)[, 2], 100)
snp3_effect <- rnorm(1, 3, 1)
snp3 <- snp3_effect * rep(Legendre(t, 3)[, 3], 100)
snp4_effect <- rnorm(1, 3, 1)
snp4 <- snp4_effect * rep(Legendre(t, 2)[, 2], 100)
snp5_effect <- rnorm(1, 3, 1)
snp5 <- snp5_effect * rep(Legendre(t, 4)[, 4], 100)
snp6_effect <- rnorm(1, 3, 1)
snp6 <- snp6_effect * rep(Legendre(t, 2)[, 2], 100)
snp7_effect <- rnorm(1, 3, 1)
snp7 <- snp7_effect * rep(Legendre(t, 3)[, 3], 100)


y <- as.vector(t(mvrnorm(n = 100, mu, COVAR))) + 
  snp1 * df[, "X1"] + 
  snp2 * df[, "X3"] + 
  snp3 * df[, "X5"] + 
  snp4 * df[, "X8"] + 
  snp5 * df[, "X1"] * df[, "X3"] + 
  snp6 * df[, "X1"] * df[, "X5"] + 
  snp7 * df[, "X3"] * df[, "X5"]

df <- data.frame(y, df)
samp <- sample(200, 50)
test_df <- df[df$id %in% samp, ]
train_df <- df[!{df$id %in% samp}, ]
# ggplot(df, aes(t, y)) + geom_path()
# ggplot(df, aes(t, y)) + geom_point() + geom_smooth()


# practice
{
# a_hat <- max(y)
# form <- y ~ X1 + X3 + X5
# params <- c(a_hat, 1, 1, rep(1, 12))
# 
# system.time({
#   out <- fminsearch(logistic_legendre_fit, params, 
#                     formula = form,
#                     data = df,
#                     time_col = "t")
# })
# 
# out
# 
# sol <- c("X1", "X3")
# cand <- c("X5", "X8")
# candidates <- c("X5", "X8")
# params <- c(a_hat, rep(1, 2 + 4 * (length(sol) + 1)))
# time_col <- "t"
# response <- "y"
# 
# ex_out <-  lapply(cand, function(candidates){
# 
#     var_names <- c(sol, candidates)
# 
#     form <- as.formula(paste(response, "~", paste(var_names, collapse = "+")))
# 
#     out <- fminsearch(logistic_legendre_fit, params,
#                       formula = form,
#                       data = df,
#                       time_col = time_col)
#     names(out$xval)<- c("a", "b", "r", as.vector(t(outer(var_names[1:3], paste0("_P", 0:3), paste0))))
#     out
# 
#   })
# 
# 
# system.time({
#  non_par <- minfval_map_func(C = cand, S = sol, response = "y", data = df, time_col = "t")
# })
#   
# system.time({
#   par_fit <- minfval_map_parallel_func(C = cand, S = sol, response = "y", data = df, time_col = "t", no_cores = 2)
# })

}

system.time({
  
# FuncMap_fit <- iForm_FunctionalMap(formula = y ~ .,
#                     data = df,
#                     id = "id",
#                     time_col = "t",
#                     heredity = "strong",
#                     higher_order = FALSE,
#                     no_cores = 8)
# 
# })

FuncMap_fit <- iForm_FunctionalMap(formula = y ~ .,
                                   data = train_df,
                                   id_col = "id",
                                   time_col = "t",
                                   heredity = "strong",
                                   higher_order = FALSE,
                                   poly_num = 5)

})


# coef_names <- names(FuncMap_fit$fit$coefficients)[-1]
# vars <- do.call(rbind,strsplit(coef_names, split = "_"))[, 1]
# poly_order <- as.numeric(substr(do.call(rbind,strsplit(coef_names, split = "_"))[, 2],2,2))


sim_list[[i]] <- list(fit = FuncMap_fit,
                 train_data = df,
                 snp_effects = list(snp1_effect, snp2_effect, snp3_effect,
                                    snp4_effect, snp5_effect, snp6_effect,
                                    snp7_effect),
                 test_df = test_df)

}

saveRDS(sim_list, "simulation_out_gls.rds")
