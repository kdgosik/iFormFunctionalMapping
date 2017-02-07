library(pracma)
library(MASS)
library(ggplot2)
rm(list=ls()); gc(reset = TRUE)
source("iForm_FunctionalMapping.R")


a <- 180
b <- 3
r <- 0.9
sigma <- 1
rho <- 0.2

m <- 10
t <- 0:9
t_s <- (t - min(t))/(max(t) - min(t))
mu <- a / (1 + b * exp(-r*t))
L <- legendre(3, t_s)


COVAR <- rho^(outer(seq(0, (m-1)), seq(0, (m-1)), function(a,b) abs(a-b)))
diag(COVAR) <- sigma^2
#COVAR <- kronecker(diag(j), COVAR)
mvrnorm(n = 1, rep(0, 10), COVAR)

snps <- matrix(rbinom(5000, 1, 0.5), nrow = 100)
snps <- data.frame(id = 1:100, snps)

time_df <- do.call(rbind, lapply(1:100, function(id) cbind(id, t)))

df <- merge(snps, time_df, by = "id", all = TRUE)

snp1_effect <- rnorm(4, 1, 0.125)
snp1 <- rep(snp1_effect %*% L, 100)
snp2_effect <- rnorm(4, 1, 0.125)
snp2 <- rep(snp2_effect %*% L, 100)
snp3_effect <- rnorm(4, 1, 0.125)
snp3 <- rep(snp3_effect %*% L, 100)
snp4_effect <- rnorm(4, 1, 0.125)
snp4 <- rep(snp4_effect %*% L, 100)
snp5_effect <- rnorm(4, 1, 0.125)
snp5 <- rep(snp5_effect %*% L, 100)
snp6_effect <- rnorm(4, 1, 0.125)
snp6 <- rep(snp6_effect %*% L, 100)
snp7_effect <- rnorm(4, 1, 0.125)
snp7 <- rep(snp7_effect %*% L, 100)

y <- as.vector(t(mvrnorm(n = 100, mu, COVAR))) + 
  snp1 * df[, "X1"] + 
  snp2 * df[, "X3"] + 
  snp3 * df[, "X5"] + 
  snp4 * df[, "X8"] + 
  snp5 * df[, "X1"] * df[, "X3"] + 
  snp6 * df[, "X1"] * df[, "X5"] + 
  snp7 * df[, "X3"] * df[, "X5"]
y <- y + abs(min(y))

df <- data.frame(y, df)
ggplot(df, aes(t, y)) + geom_path()
ggplot(df, aes(t, y)) + geom_point() + geom_smooth()


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
# 
# minfval_map_func(C = cand, S = sol, response = "y", data = df, time_col = "t")
}

system.time({
  
FuncMap_fit <- iForm_FunctionalMap(formula = y ~ .,
                    data = df,
                    id = "id",
                    time_col = "t",
                    heredity = "strong",
                    higher_order = FALSE)

})