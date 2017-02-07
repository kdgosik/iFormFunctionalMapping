setwd("/home/dudusdau/fGWAS/ulgwas.final.code/")
#system("R CMD SHLIB dcovxx2.c")

dyn.load("./2HIGWAS/src/dcovxx2.so")
load("2HIGWAS/sim.beta.RData")
library(mvtnorm)
library(MSBVAR)
mu <- c(2.647295,3.841694,4.003850,5.975076,10.157559,14.900797,18.382197,21.927349
        ,23.571484,24.656771,25.337658,26.603671,27.318574,28.376232,29.258628,30.021319
        ,31.369278,33.224948,34.861341,35.639403,36.294607,37.095920,38.286930,39.079150)

source("2HIGWAS/ulgwas.sim.R")
source("2HIGWAS/ulgwas.stage1.R")
source("2HIGWAS/ulgwas.stage2.R")


par <- sim_param(mu,sim.beta=sim.beta,nt=24) 

dat <- simulate_longdt(par)

# dc.sis <- DC.sis(dat)
dc.sis <- list(rr=rr,sisnp=sissnp,iid=ii.id)

load("2HIGWAS/snpdc.RData")

y <- as.matrix(dat$pheno_y)

lambda <- seq(0.1,1.5,0.1)
dd <- Vsugpr(lambda,dc.sis$sisnp,y,method="scad",ti=c(1:24))
length(which(rowMeans(abs(dd$beta.aic))[-1]> 1e-2))
