library(lme4)
library(magrittr)
library(nlme)
library(lattice)
rm(list=ls()); gc(reset = TRUE)

data(sleepstudy); str(sleepstudy)
data(Orange); str(Orange)

fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
fm2 <- lmer(Reaction ~ Days + (Days || Subject), sleepstudy)

data(Orthodont,package="nlme"); str(Orthodont)
Orthodont$nsex <- as.numeric(Orthodont$Sex=="Male")
Orthodont$nsexage <- with(Orthodont, nsex*age)
lmer(distance ~ age + (age|Subject) + (0+nsex|Subject) +
       (0 + nsexage|Subject), data=Orthodont)

lmer(distance ~ age + (age|Subject) + Sex, data = Orthodont)



  # nlmer example with Orange data
startvec <- c(Asym = 200, xmid = 725, scal = 350)
(nm1 <- nlmer(circumference ~ SSlogis(age, Asym, xmid, scal) ~ Asym|Tree,
              Orange, start = startvec))


## selfStart Example ##########

## Using selfStart ##
  # rewriting SSlogis function

SSlogis <- selfStart(~ Asym/(1 + exp((xmid - x)/scal)),
                     function(mCall, data, LHS)
                     {
                       xy <- sortedXyData(mCall[["x"]], LHS, data)
                       if(nrow(xy) < 4) {
                         stop("Too few distinct x values to fit a logistic")
                       }
                       z <- xy[["y"]]
                       if (min(z) <= 0) { z <- z + 0.05 * max(z) } # avoid zeroes
                       z <- z/(1.05 * max(z))              # scale to within unit height
                       xy[["z"]] <- log(z/(1 - z))         # logit transformation
                       aux <- coef(lm(x ~ z, xy))
                       parameters(xy) <- list(xmid = aux[1], scal = aux[2])
                       pars <- as.vector(coef(nls(y ~ 1/(1 + exp((xmid - x)/scal)),
                                                  data = xy, algorithm = "plinear")))
                       setNames(c(pars[3], pars[1], pars[2]),
                                mCall[c("Asym", "xmid", "scal")])
                     }, c("Asym", "xmid", "scal"))

## NLME Examples ################

fm1 <- lme(distance ~ age, data = Orthodont)
fm2 <- lme(distance ~ age + Sex, data = Orthodont, random = ~ 1)


wages <- read.table("http://www.ats.ucla.edu/stat/r/examples/alda/data/wages_pp.txt", header=T, sep=",")
wages$ged.exper <- wages$ged*wages$exper
print(wages[wages$id %in% c(206,2365,4384),c(1:5, 16)])

model.a <- lme(lnw~exper+hgc.9+exper:black+ue.7, wages, random= ~exper | id, method="ML") 
2*model.a$logLik

model.b <- lme(lnw~exper+hgc.9+exper:black+ue.7+ged, wages, random= ~ exper+ged | id, method="ML")
2*model.b$logLik

anova(model.a, model.b)


model.c <- lme(lnw~exper+hgc.9+exper:black+ue.7+ged, wages, random= ~ exper | id, method="ML")
2*model.c$logLik

anova(model.b, model.c)

model.d <- lme(lnw~exper+hgc.9+exper:black+ue.7+postexp, wages, random= ~exper+postexp | id, method="ML")
2*model.d$logLik

anova(model.a, model.d)

model.e <- lme(lnw~exper+hgc.9+exper:black+ue.7+postexp, wages, random= ~exper | id, method="ML")
2*model.e$logLik


## Example 2

#reading in the alcohol data
alcohol <- read.table("http://www.ats.ucla.edu/stat/r/examples/alda/data/alcohol1_pp.txt", header=T, sep=",")

#table 4.1, model e
model.e <- lme(alcuse~coa+peer*age_14 , data=alcohol, random= ~ age_14 | id, method="ML")
#obtaining the fixed effects parameters 
fixef.e <- fixef(model.e)
#obtaining the predicted values and squaring them
fit2.ec0p0 <- (fixef.e[[1]] + .655*fixef.e[[3]] +
                 alcohol$age_14[1:3]*fixef.e[[4]] +
                 .655*alcohol$age_14[1:3]*fixef.e[[5]])^2   
fit2.ec0p1 <- (fixef.e[[1]] + 1.381*fixef.e[[3]] +
                 alcohol$age_14[1:3]*fixef.e[[4]] +
                 1.381*alcohol$age_14[1:3]*fixef.e[[5]] )^2
fit2.ec1p0 <- (fixef.e[[1]] + fixef.e[[2]] + .655*fixef.e[[3]] +
                 alcohol$age_14[1:3]*fixef.e[[4]] +
                 .655*alcohol$age_14[1:3]*fixef.e[[5]] )^2
fit2.ec1p1 <- (fixef.e[[1]] + fixef.e[[2]] + 1.381*fixef.e[[3]] +
                 alcohol$age_14[1:3]*fixef.e[[4]] +
                 1.381*alcohol$age_14[1:3]*fixef.e[[5]])^2

plot(alcohol$age[1:3], fit2.ec0p0, ylim=c(0, 3), type="n", 
     ylab="predicted alcuse squared", xlab="age")
lines(spline(alcohol$age[1:3], fit2.ec0p0), pch=2, type="b")
lines(spline(alcohol$age[1:3], fit2.ec0p1), type="b", pch=0)   
lines(spline(alcohol$age[1:3], fit2.ec1p0), type="b", pch=17)   
lines(spline(alcohol$age[1:3], fit2.ec1p1), type="b", pch=15)   

title("Non-linear Change") 
legend(14, 3, c("COA=0, low peer", "COA=0, high peer", 
                "COA=1, low peer", "COA=1, high peer"))


## Growth of 8 children example

# Reading in the fox and geese data.

fg <- read.table("http://www.ats.ucla.edu/stat/r/examples/alda/data/foxngeese_pp.txt", header=T, sep=",")

# Fig. 6.8, p. 227.
# Empirical growth plots for 8 children in the fox and geese data.

xyplot(nmoves~game | id, data=fg[fg$id %in% c(1, 4, 6, 7, 8, 11, 12, 15), ],
       ylim=c(0, 25), as.table=T)


# Table 6.6, p. 231. Fitting a logistic model to the fox and geese data.
# 
# Notice, this model does not correspond to equation (6.8) in the book. Instead, it corresponds to the following equation:
#   
#   Yij = 1 + 19/(1 + π0*exp(- (π1+u1)*Time - u0)) + εij

# game is time element

model.a <- nlme(nmoves ~ 1 + 19/ (1 + xmid*exp( -scal*game + u)),
                fixed = scal+xmid~1, 
                random= scal+u~1 |id, 
                start = c(scal = 0.2, xmid = 12), data = fg)

summary(model.a)

model.b <- nlme(nmoves ~ 1 + 19/(1+xmid*exp(-scal10*game -scal01*read -scal11*read*game + u)),
                fixed = scal10+scal01+scal11+xmid~1, 
                random = scal10+u~1 |id, 
                start = c(scal10=.12, scal01 = -0.4, scal11 = 0.04, xmid = 12), data = fg)

summary(model.b)

fixef.a <- fixef(model.a)
fit.a <- 1 + 19/(1 + fixef.a[[2]]*exp(-fixef.a[[1]]*fg$game[1:27]))

plot(fg$game[1:27], fit.a, ylim=c(0, 25), type="l", ylab="predicted nmoves", xlab="game")
title("Model A \n Unconditional logistic growth")


fixef.b <- fixef(model.b)
fit.b.high <- 1 + 19/(1+fixef.b[[4]]*exp(-fixef.b[[1]]*fg$game[1:27] - fixef.b[[2]]*1.58 - fixef.b[[3]]*1.58*fg$game[1:27]))

fit.b.low <- 1 + 19/(1+fixef.b[[4]]*exp(-fixef.b[[1]]*fg$game[1:27] - fixef.b[[2]]*(-1.58) - fixef.b[[3]]*(-1.58)*fg$game[1:27]))

plot(fg$game[1:27], fit.b.high, ylim=c(0, 25), type="l", 
     ylab="predicted nmoves", xlab="game")
lines(fg$game[1:27], fit.b.low, lty=3)

title("Model B \n Fitted logistic growth by reading level")
legend(1, 25, c("High Reading level","Low reading level"), lty=c(1, 3))





## nlme ####
  #CO2 example
plot(CO2, outer = ~Treatment * Type, layout = c(4, 1))

x = seq(2.8, 4, 0.01)
f = function(a, b, c, x) {
  return(a * (1 - exp(-exp(b) * (x - c))))
}
y = f(1.5, 2, 3, x)
plot(x, y, type = "l", ylim = c(-2, 3))
grid(5, 5, lwd = 2)

CO2.list <- nlsList(SSasympOff, CO2)
plot(intervals(CO2.list))


dia <- read.csv("L-abr-1_DIA.csv")
dia <- fread("L-abr-1_DIA.csv")
setnames(dia, names(dia), c("Tree", "a", "b", "r"))

dia[, paste(1:12) := lapply(1:12, function(i) a/(1+b*exp(-r*i)))]
dia_melt <- melt(dia[,.SD,.SDcols = c("Tree", paste0(1:12))], "Tree")

xyplot(value~ variable| Tree, data=dia_melt[Tree %in% c(2, 3,  4, 5, 6, 7, 9, 11), ],
       ylim=c(0, 15), as.table=T)



  # need to combine these two formulas and also include selection

  # Using AR(1) model
lme(y ~ time * tx, 
    random = ~ time | subjects, 
    correlation = corAR1(),
    data=data)


# changing the functional form of time
lme(y ~ (time + I(time^2)) * tx,
    random = ~time + I(time^2) | subjects,
    data=data)   


# first attempt
  # a/(1+b*exp(-rt))

lme(y ~ (a/(1 + b*exp(-r * time))) * (Selected SNPs and Interactions go here),
    random = ~(a/(1 + b*exp(-r * time)))|tree,
    correlation = corAR1(),
    data = mei_tree_data)


nlme(model, data, fixed, random, groups, start)



#   Yij = 1 + 19/(1 + π0*exp(- (π1+u1)*Time - u0)) + εij

# game is time element

model.a <- nlme(nmoves ~ 1 + 19/ (1 + xmid*exp( -scal*game + u)),
                fixed = scal+xmid~1, 
                random= scal+u~1 |id, 
                start = c(scal = 0.2, xmid = 12), data = fg)


library(car)
head(Wong)

# Yij = θ1i + θ2ie−θ3iX1ij + εij 
# θ1i = β1 + β2*sqrt(X2i) + δ1i
# θ2i = β3 + β4*sqrt(X2i) + δ2i
# θ3i = β5

form1 <- theta1~1+sqrt(duration)

wong.mod.1 <- nlme(piq ~ theta1 + theta2*exp(-theta3*days), data=Wong,
                   fixed=list(
                    form1,
                    theta2 ~ 1 + sqrt(duration),
                    theta3 ~ 1),
                   random=list(id = list(theta1 ~ 1, theta2 ~ 1)),
                   start=list(fixed=c(100, -2, -10, 0, 0.007)))
summary(wong.mod.1)

mei_model <- nlme(y ~ a / (1 + b*exp( -r * time)),
                fixed = a + b + r ~ SNP, 
                random = s +b + r ~ 1 | Tree, 
                start = c(a = 5, b = 1, r = 1), data = Mei_Trees)


# comparing models
anova(m1, m2)

1 - pchisq(m1$logLik*-2 - m2$logLik*-2, 1)


# Mei Tree practice

dia <- read.csv("L-abr-1_DIA.csv")
dia <- fread("L-abr-1_DIA.csv")
setnames(dia, names(dia), c("Tree", "a", "b", "r"))
setkeyv(dia, "Tree")
dia <- unique(dia)

dia[, paste(1:12) := lapply(1:12, function(i) a/(1+b*exp(-r*i)))]
dia_melt <- melt(dia[,.SD,.SDcols = c("Tree", paste0(1:12))], "Tree")
dia_melt[, Tree := as.character(Tree)]
dia_melt[,variable := as.numeric(variable)]


markers <- fread("markers_impute_two_levels.csv")
setnames(markers, "V1", "Tree")
setkeyv(markers, "Tree")


dat <- markers[dia_melt]
setnames(dat, c("variable", "value"), c("time", "y"))
dat <- dat[!is.na(CATG_lm_ll_1571), ]


mei_model <- nlme(y ~ a / (1 + b * exp(-r*time)),
                  data = dat,
                  fixed = list(a ~ 1 + CATG_lm_ll_1571,
                               b ~ 1 + CATG_lm_ll_1571,
                               r ~ 1 + CATG_lm_ll_1571), 
                  random = list(Tree = list(a ~ 1, 
                                            b ~ 1, 
                                            r ~ 1)), 
                  start = list(fixed = rep(1, 6)), 
                  method = "ML")
summary(mei_model)



## optimx example ####

library(optimx)

## Show multiple outputs of optimx using all.methods
# genrose function code
genrose.f<- function(x, gs=NULL){ # objective function
  ## One generalization of the Rosenbrock banana valley function (n parameters)
  n <- length(x)
  if(is.null(gs)) { gs=100.0 }
  fval<-1.0 + sum (gs*(x[1:(n-1)]^2 - x[2:n])^2 + (x[2:n] - 1)^2)
  return(fval)
}

genrose.g <- function(x, gs=NULL){
  # vectorized gradient for genrose.f
  # Ravi Varadhan 2009-04-03
  n <- length(x)
  if(is.null(gs)) { gs=100.0 }
  gg <- as.vector(rep(0, n))
  tn <- 2:n
  tn1 <- tn - 1
  z1 <- x[tn] - x[tn1]^2
  z2 <- 1 - x[tn]
  gg[tn] <- 2 * (gs * z1 - z2)
  gg[tn1] <- gg[tn1] - 4 * gs * x[tn1] * z1
  return(gg)
}

genrose.h <- function(x, gs=NULL) { ## compute Hessian
  if(is.null(gs)) { gs=100.0 }
  n <- length(x)
  hh<-matrix(rep(0, n*n),n,n)
  
  for (i in 2:n) {
    z1<-x[i]-x[i-1]*x[i-1]
    z2<-1.0-x[i]
    hh[i,i]<-hh[i,i]+2.0*(gs+1.0)
    hh[i-1,i-1]<-hh[i-1,i-1]-4.0*gs*z1-4.0*gs*x[i-1]*(-2.0*x[i-1])
    hh[i,i-1]<-hh[i,i-1]-4.0*gs*x[i-1]
    hh[i-1,i]<-hh[i-1,i]-4.0*gs*x[i-1]
  }
  return(hh)
}

startx <- 4*seq(1:10)/3.
ans8 <- optimx(startx,fn=genrose.f,gr=genrose.g, hess=genrose.h,
             control=list(all.methods=TRUE, save.failures=TRUE, trace=0), gs=10)
ans8
ans8[, "gevals"]
ans8["spg", ]
summary(ans8, par.select = 1:3)
summary(ans8, order = value)[1, ] # show best value
head(summary(ans8, order = value)) # best few
## head(summary(ans8, order = "value")) # best few -- alternative syntax
## order by value. Within those values the same to 3 decimals order by fevals.
## summary(ans8, order = list(round(value, 3), fevals), par.select = FALSE)
summary(ans8, order = "list(round(value, 3), fevals)", par.select = FALSE)
## summary(ans8, order = rownames, par.select = FALSE) # order by method name
summary(ans8, order = "rownames", par.select = FALSE) # same
summary(ans8, order = NULL, par.select = FALSE) # use input order
## summary(ans8, par.select = FALSE) # same



## fminsearch for R #####

library(pracma)

# Rosenbrock function
rosena <- function(x, a) 100*(x[2]-x[1]^2)^2 + (a-x[1])^2 # min: (a, a^2)
fminsearch(rosena, c(-1.2, 1), a = sqrt(2))
# x = (1.414214 2.000010) , fval = 1.239435e-11
fminsearch(rosena, c(-1.2, 1), dfree=FALSE, a = sqrt(2))
# x = (1.414214 2.000000) , fval = 3.844519e-26


