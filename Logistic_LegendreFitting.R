library(orthopolynom)
library(pracma)
library(data.table)
library(tidyverse)
library(ggplot2)


(leg4coef <- legendre.polynomials(n=4, normalized=TRUE))

x <- 1:100
leg4 <- as.matrix(
  as.data.frame(
    polynomial.values(polynomials = leg4coef, x = scaleX(x, u = -1, v = 1)
                      )
  )
)

colnames(leg4) <- c("leg0", "leg1", "leg2", "leg3", "leg4")
leg4 <- leg4[, 2:ncol(leg4)]


lm(y ~ leg4)


Legendre <- function(x, n, normalized=TRUE, intercept=FALSE, rescale=TRUE)
{
  ## Create a design matrix for Legendre polynomials
  ## x - numeric
  ## n - see orthopolynom
  ## normalized - logical, see orthopolynom
  ## intercept - logical, add intercept
  tmp <- legendre.polynomials(n=n, normalized=normalized)
  if(!intercept) tmp <- tmp[2:length(tmp)]
  polynomial.values(polynomials=tmp, x=x, matrix=TRUE)
}


polynomial.values <- function(polynomials, x, matrix=FALSE)
{
  ## Changed copy of polynomial.vales from orthopolynom in order
  ## to add matrix argument
  require(polynom)
  n <- length(polynomials)
  if(!matrix) {
    values <- vector(mode="list", length=n)
  } else {
    values <- matrix(ncol=n, nrow=length(x))
  }
  j <- 1
  while(j <= n) {
    if(!matrix) {
      values[[j]] <- predict(polynomials[[j]], x)
    } else {
      values[, j] <- predict(polynomials[[j]], x)
    }
    j <- j + 1
  }
  values
}


lm(y ~ Legendre(x = scaleX(x, u = -1, v = 1), n = 4))




library(LPTime)
data(EyeTrack.sample)
head(LPiTrack(EyeTrack.sample), m = c(3, 5, 15), p = 3)





library(XLConnect)
library(pracma)
library(zoo)
library(tidyr)
library(dplyr)
library(stringr)

# fminsearch()
# legendre()

wb <- loadWorkbook("MeiTreeOriginalData.xlsx")
s1 <- readWorksheet(wb, sheet = 1, header = FALSE, startRow = 1)
colnames(s1) <- s1[2,]
time_period <- s1[1:2,-c(1,2)]
gather(time_period, )

s1 <- s1[-c(1,2), ]
s1$progeny <- na.locf(s1$progeny)
colnames(s1)[2] <- "replicate"

s1_gather <- gather(s1, key, value, -progeny, -replicate)
s1_gather <- na.omit(s1_gather)

ht <- s1_gather %>%
  filter(str_detect(key, 'HT'))


t <- 1:15
legendre(3, t)



x <- seq(-1, 1, len = 50)
L <- legendre(5, x)
plot(x, L[1, ], type = "l", col = 1, ylim = c(-2, 3), ylab = "y",
     main = "Legendre Functions of degree 2")
lines(x, L[2, ], col = 2)
lines(x, L[3, ], col = 3)

grid()


# 
# SAD-1
# 
# sigma_1^2                         
# 
# sigma_1*sigma_2*rho_1               sigma_2^2
# 
# sigma_1*sigm_3*rho_1*rho_2          sigma_2*sigma_3*rho_2         sigma_3^2
# 
# sigma_1*sigma_4*rho_1*rho_2*rho_3   sigma_2*sigm_4*rho_2*rho_3    sigma_3*sigma_4*rho_3   sigma_4^2
# 


rho 
sigma
Sigma <- matrix(c(1, 0.2, 0.2, 1),2,2)
Sigma
var(mvrnorm(n=1000, rep(0, 2), Sigma))

logistic_growth <- function(t, a, b, r){
  sapply(t, function(x) a/(1+b*exp(-r*x)))
}

y <- logistic_growth(t = seq(-15,15,len = 50), a = 1.5, b = 0.5, r = 0.5)

y <- y + rnorm(length(y), 0, 0.3)



df <- gather(s1[1:5,which(sapply(colnames(s1), grepl, "HT")) ]) %>%
separate(col = key, into = c("key", "time")) %>%
  mutate(time = sort(rep(1:11,5)), replicate = rep(1:5, 11))


p <- ggplot(df, aes(x = time, y = value, fill = replicate, group = replicate)) + 
  geom_point() + 
  geom_line()

df_sum <- df %>%
  mutate(value = as.numeric(value))
  group_by(time) %>%
  summarise(mean = mean(value))

  
  DT <- as.data.table(df)
  DT_sum <- DT[, mean(as.numeric(value)), by = time]
  
logi_fit <- function(params, data){
  a <- params[1]
  b <- params[2]
  r <- params[3]
  
  t <- data$time
  y <- data$V1
  
  y_p <- a/(1+b*exp(-r*t))
  sum((y - y_p)^2)
}
  
out <- fminsearch(logi_fit, c(100, 1, 1), data = DT_sum)

logi_predict <- function(out, data){
  a <- out$xval[1]
  b <- out$xval[2]
  r <- out$xval[3]
  
  t <- data$time
  
  y_p <- a/(1+b*exp(-r*t))
  y_p
}
