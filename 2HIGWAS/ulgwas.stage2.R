Vsugpr <- function(lambda,x,y,method,ti){

  AIC <- c()
  BIC <- c()
  GCV <- c()
  beta <- array(0,dim=c((dim(x)[2]+1),dim(y)[2],length(lambda)))

  for(i in 1:length(lambda)){
    ll <- lambda[i]
    if(method=="lasso")
      aa <- lqa (x, y, family = gaussian (), penalty = lasso(ll),standardize = F,ti)
    else if(method=="adaptive.lasso")
      aa <- lqa (x, y, family = gaussian (), penalty = adaptive.lasso(ll, al.weights=0.5),standardize = F,ti)
    else
      aa <- lqa (x, y, family = gaussian (), penalty = scad(ll),standardize = F,ti)
    beta[,,i] <- aa$coefficients
    BIC <- c(BIC,aa$bic)
    AIC <- c(AIC,aa$aic)
    GCV <- c(GCV,aa$gcv)
    cat("Lambda=",ll,"\n")
  }
  par(mfrow=c(1,3))
  plot(AIC,type="l",col="red")
  plot(BIC,type="l",col="green")
  plot(GCV,type="l",col="blue")
  
  l.aic <- which(AIC==min(AIC))
  l.bic <- which(BIC==min(BIC))
  l.gcv <- which(GCV==min(GCV))
  beta.aic <- beta[,,l.aic]
  beta.bic <- beta[,,l.bic]
  beta.gcv <- beta[,,l.gcv]
  
  return(list(l.aic=l.aic,l.bic=l.bic,l.gcv=l.gcv,beta.aic=beta.aic,beta.bic=beta.bic,
              beta.gcv=beta.gcv))
  
}




lqa <- function (x, y, family = gaussian (), penalty = NULL, ti,
          method = "lqa.update2", start = NULL, etastart = NULL,mustart = NULL, 
          offset = rep (0, nobs), control = lqa.control (), 
          intercept = TRUE, standardize = TRUE){
  
  if (is.character (family)) 
    family <- get (family, mode = "function", envir = parent.frame ())
  
  if (is.function (family)) 
    family <- family ()
  
  if (is.null (family$family)) 
  {
    print (family)
    stop ("'family' not recognized")
  }
  
  if (is.null (method))
    stop ("method not specified") 
  
  ### Check for column of ones in x if intercept == TRUE:
  ### ---------------------------------------------------
  
  if (intercept & (var (x[,1]) > control$var.eps))
    x <- cbind (1, x)
  
  
  ### Standardization:
  ### ----------------
  
  x <- as.matrix (x)
  xnames <- dimnames (x)[[2L]]
  
  ynames <- if (is.matrix (y)) 
    rownames(y)
  else 
    names(y)
  
  nobs <- nrow (x)   
  nvars <- ncol (x)    # number of coefficients in the predictor (including an intercept, if present)
  ones <- rep (1, nobs)
  mean.x <- drop (ones %*% x) / nobs       # computes the vector of means
  
  if (intercept)    # if an intercept is included in the model its corresponding element of mean.x is set to zero
    mean.x[1] <- 0  # (such that x[,1] is not getting centered (and also not standardized later on ...))
  
  x.std <- scale (x, mean.x, FALSE)   # centers the regressor matrix
  norm.x <- if (standardize)
  {
    norm.x <- sqrt (drop (ones %*% (x.std^2)))   # computes the euclidean norm of the regressors
    nosignal <- apply (x, 2, var) < control$var.eps
    if (any (nosignal))    # modify norm.x for variables with too small a variance (e.g. the intercept)
      norm.x[nosignal] <- 1
    
    norm.x
  }
  else
    rep (1, nvars)
  
  x.std <- scale (x.std, FALSE, norm.x)  # standardizes the centered regressor matrix
  
  fit <- lqa.update2(x = x.std, y = y, family = family, penalty = penalty, intercept = intercept, 
                    ti=ti,control = control) 
  ### Back-Transformation of estimated coefficients:
  ### ----------------------------------------------
  
  coef <- fit$coefficients
  
  if (intercept)
  {
    coef[1,] <- coef[1,] - colSums(mean.x[-1] * coef[-1,] / norm.x[-1])
    coef[-1,] <- coef[-1,] / norm.x[-1]
  }
  else
    coef <- coef / norm.x
  
  
  ### Computation of some important statistics:
  ### -----------------------------------------
  
  # Remark: The predictors are identical no matter whether we use standardized or unstandardized values, hence all statistics
  #         based on the predictor eta are also equal
  
  eta <- drop (x %*% coef)
  mu <- family$linkinv (eta)
  dev <- sum (family$dev.resids (y, mu, 1))
  residuals <- (y - mu)
  
  xnames <- colnames (x)
  ynames <- names (y)
  names (residuals) <- names (mu) <- names (eta)<- ynames
  rownames (coef) <- xnames
  
  if (is.null (fit$tr.H))
    stop ("quadpen.fit: Element 'tr.H' has not been returned from 'method'")
  
  model.aic <- dev + 2 * fit$tr.H
  model.bic <- dev + log (nobs) * fit$tr.H  
  model.gcv <- dev/nobs *((1+2*fit$tr.H/nobs)) 
  fit <- list (coefficients = coef, residuals = residuals, fitted.values = mu, family = family, penalty = penalty, 
               linear.predictors = eta, deviance = dev, aic = model.aic, bic = model.bic, gcv=model.gcv,n.iter = fit$stop.at, 
               best.iter = fit$m.stop,converged = fit$converged, mean.x=mean.x,abr.mat =fit$abr.mat,
               norm.x = norm.x, method = method, rank = fit$tr.H, x = x, y = y, fit.obj = fit)
  
  class (fit) <- c ("lqa")
  fit
}


lqa.update2 <-function (x, y, family = NULL, penalty = NULL, intercept = TRUE, ti,control = lqa.control (), initial.beta, mustart, eta.new, gamma1 = 1, ...)
{
  gamma <- gamma1
  
  if (is.null (family))
    stop ("lqa.update: family not specified")
  
  if (is.null (penalty))
    stop ("lqa.update: penalty not specified")
  
  #if (!is.null (dim (y)))
  #stop ("lqa.update: y must be a vector")
  x <- as.matrix (x)
  converged <- FALSE
  n.iter <- control$max.steps
  eps <- control$conv.eps
  c1 <- control$c1
  ti <<- ti
  t <- dim(y)[2]
  stop.at <- n.iter
  p <- ncol (x)
  nobs <- nrow (x)
  converged <- FALSE
  
  g1 <- function(par){
    if(any(par<0)){
      par[which(par<0)]<- 0.5
    }
    if(par[7]>1)
      par[7] <- 0.5
    a1 <- par[1];b1 <- par[2];r1 <- par[3]
    a2 <- par[4];b2 <- par[5];r2 <- par[6];p <- par[7]
    p*a1/(1+b1*exp(-r1*ti)) + (1-p)*a2/(1+b2*exp(-r2*ti))
  }
  
  
  sumpqg1 <- sumsq <- function(par,yy) { 
    sum( (yy - g1(par) )^2 )
  }
  j <- 4
  abr.mat <- matrix(rep(0,p*18),ncol=18)
  for(i in 1:p){
    if(i==1) 
      abr.mat[i,1:7] <-  c(19.8,5.26,0.26,10,1.5,0.1,0.5) 
    else{    
      abr.mat[i,1:10] <- c(0.1:1)
    }     
  }
  
  
  beta.mat <- array (0, dim=c(nrow = n.iter, ncol = p,t=t))   # to store the coefficient updates
  
  if (missing (initial.beta)){
    initial.beta <- matrix(rep (0,t*dim(abr.mat)[1]),nrow=dim(abr.mat)[1])
    for(i in 1:p){
      if(i==1) 
        initial.beta[i,] <- g1(abr.mat[i,1:7])
      else{
        initial.beta[i,] <- Legendre(ti,np.order=4)%*%abr.mat[i,1:4]
      }     
    }
  }
  #initial.beta <- solve(crossprod (x)) %*% t (x) %*% y
  else
    eta.new <- drop (x %*% initial.beta)
  
  for(j in 4:4){
    
    
    maxit <- 200
    k <- 1
    for(ii in 1:maxit){
      old.mat <- abr.mat
      for(i in 1:p){
        for(i in 1:p){
          if(i==1) 
            initial.beta[i,] <- g1(abr.mat[i,1:7])
          else{
            initial.beta[i,] <- Legendre(ti,np.order=j)%*%abr.mat[i,1:j]
          }     
        }
      }
      if (missing (mustart))
      {
        etastart <- drop (x %*% initial.beta)
        eval (family$initialize)
      }
      
      if (missing (eta.new))
        eta.new <- family$linkfun (mustart)    # predictor
      
      
      for (i in 1 : n.iter)
      {
        beta.mat[i,,] <- initial.beta  
        mu.new <- family$linkinv (eta.new)      # fitted values
        d.new <- matrix(family$mu.eta (eta.new),nrow=nobs)        # derivative of response function
        v.new <- matrix(family$variance (mu.new),nrow= nobs)      # variance function of the response
        x.star <-  x   
        y.tilde.star <- eta.new  + (y - mu.new) / d.new 
        
        A.lambda <- get.Amat (initial.beta = initial.beta, penalty = penalty, intercept = intercept, c1 = c1, x = x) 
        p.imat.new <- crossprod (x.star)/nobs + A.lambda       # penalized information matrix
        
        chol.pimat.new <- chol (p.imat.new)               # applying cholesky decomposition for matrix inversion
        inv.pimat.new <- chol2inv (chol.pimat.new)        # inverted penalized information matrix
        beta.new <- gamma * drop (inv.pimat.new %*% (t(x.star) %*% y.tilde.star)/nobs) + (1 - gamma) * beta.mat[i,,]  # computes the next iterate of the beta vector
        
        #if ((sum (abs (beta.new - initial.beta)) / sum (abs (initial.beta)) <= eps)) 
        if (all((abs (beta.new - initial.beta))<= eps)) # check convergence condition
        {
          converged <- TRUE
          stop.at <- i
          if (i < n.iter)
            break
        } 
        else
        {
          initial.beta <- beta.new    # update beta vector
          eta.new <- drop (x %*% beta.new)      
        }
      }
      #upperD <- c(17,11,0.2,19,15,1.2,17,11,0.2,19,15,1.2,17,11,0.2,19,15,1.2)
      upperd <- c(100,100,1,70,100,1,0.3)
      #upperd <- c(Inf,Inf,Inf,Inf,Inf,Inf,Inf)
      lowerd <- c(0.1,0.1,0.1,0.1,0.1,0.001,0.2)
      for(i in 1:p){        
        if(i==1){
          abr.mat[i,1:7] <- optim(abr.mat[i,1:7],sumpqg1,yy=beta.new[i,],
                                  method="L-BFGS-B",upper=upperd,lower=lowerd)$par
        }
        else{
          XX <- Legendre(ti,np.order=j)
          abr.mat[i,1:j] <- solve(t(XX)%*%XX)%*%t(XX)%*%matrix(beta.new[i,],ncol=1)
        }    
      }
      k <- k + 1
      if(all(abs(abr.mat-old.mat) < 1e-5))
        break;
    }
    
    Hatmat <- x.star %*% inv.pimat.new %*% t (x.star)
    tr.H <- sum (diag (Hatmat))
    dev.m <- sum (family$dev.resids (y, mu.new,1))
  }  
  if (!converged & (stop.at == n.iter))
    cat ("lqa.update with ", penalty$penalty, ": convergence warning! (lambda = ", penalty$lambda, ")\n")
  
  
  fit <- list (coefficients = beta.new, beta.mat = beta.mat[1 : stop.at,,], tr.H = tr.H, fitted.values = mu.new, family = family, Amat = A.lambda, converged = converged, stop.at = stop.at, m.stop = stop.at, linear.predictors = eta.new, 
               k.iter=k,p.imat = p.imat.new, inv.pimat = inv.pimat.new, x.star = x.star, v.new = v.new,Hatmat=Hatmat,abr.mat=abr.mat)
}



lasso <- function (lambda = NULL) 
{
  lambda.check(lambda)
  if (length(lambda) != 1) 
    stop("lambda must be a scalar \n")
  names(lambda) <- "lambda"
  getpenmat <- function(beta = NULL, c1 = lqa.control()$c1 ) {
    if (is.null(beta)) 
      stop("'beta' must be the current coefficient vector \n")
    if (c1 < 0) 
      stop("'c1' must be non-negative \n")
    penmat <- lambda * diag(1/(sqrt(rowSums(beta^2)) + c1)) * as.integer(rowSums(beta) != 
                                                                  0)
    penmat
  }
  first.derivative <- function(beta) {
    p <- length(beta)
    return(rep(lambda, p))
  }
  structure(list(penalty = "lasso", lambda = lambda, getpenmat = getpenmat, 
                 first.derivative = first.derivative), class = "penalty")
}

scad <- function (lambda = NULL) 
{
  lambda.check(lambda)
  if (length(lambda) > 2) 
    stop("The scad penalty must consist on two parameters! \n")
  if (length(lambda) == 1) 
    lambda[2] <- 3.7
  if (lambda[2] <= 2) 
    stop("'lambda[2]' must be '> 2'")
  names(lambda) <- c("lambda", "a")
  first.derivative <- function(beta = NULL) {
    if (is.null(beta)) 
      stop("'beta' must be the current coefficient vector \n")
    lambda1 <- lambda[1]
    a <- lambda[2]
    theta <- sqrt(rowSums(beta^2)) + 1e-06
    p <- length(beta)
    help1 <- sapply(theta, function(theta) {
      max(a * lambda1 - theta, 0)/((a - 1) * lambda1)
    })
    lambda1 * ((theta <= lambda1) + help1 * (theta > lambda1))
  }
  getpenmat <- function(beta = NULL, c1 = lqa.control()$c1) {
    if (is.null(beta)) 
      stop("'beta' must be the current coefficient vector \n")
    if (c1 < 0) 
      stop("'c1' must be non-negative \n")
    A <- sqrt(rowSums(beta^2)) + 1e-06
    penmat <- diag(first.derivative(beta = beta)/(2*A), length(sqrt(rowSums(beta^2)) + A))
    penmat
  }
  structure(list(penalty = "scad", lambda = lambda, first.derivative = first.derivative, 
                 getpenmat = getpenmat), class = "penalty")

}

adaptive.lasso <- function (lambda = NULL, al.weights = NULL) 
{
  lambda.check(lambda)
  if (length(lambda) != 1) 
    stop("lambda must be a scalar \n")
  names(lambda) <- "lambda"
  if (is.null(al.weights)) 
    al.weights <- 1
  getpenmat <- function(beta = NULL, c1 = lqa.control()$c1) {
    if (is.null(beta)) 
      stop("'beta' must be the current coefficient vector \n")
    if (c1 < 0) 
      stop("'c1' must be non-negative \n")
    penmat <- lambda * diag(al.weights/(sqrt(rowSums(beta^2) + c1))) *  as.integer(rowSums(beta) != 0)
    penmat
  }
  first.derivative <- function(beta) {
    p <- length(beta)
    return(rep(lambda * al.weights, p))
  }
  structure(list(penalty = "adaptive.lasso", lambda = lambda, 
                 getpenmat = getpenmat, first.derivative = first.derivative), 
            class = "penalty")
}




lambda.check <-function (lambda) 
{
  if (is.null(lambda)) 
    stop("Tuning parameter lambda missing \n")
  if (any(lambda < 0)) 
    stop("lambda must be non-negative \n")
}


  lqa.control <- function (x = NULL, var.eps = .Machine$double.eps, max.steps = 5000, 
            conv.eps = 0.001, conv.stop = TRUE, c1 = 1e-08, digits = 5) 
  {
    if (!is.numeric(var.eps) || var.eps < 0) 
      stop("value of var.eps must be >= 0")
    max.steps <- as.integer(max.steps)
    digits <- as.integer(digits)
    if (max.steps < 0) 
      stop("max.steps must be positive integer")
    if (!is.numeric(conv.eps) || conv.eps <= 0) 
      stop("value of conv.eps must be > 0")
    if (!is.logical(conv.stop)) 
      stop("conv.stop must be 'TRUE' or 'FALSE'")
    if (!is.numeric(c1) || c1 < 0) 
      stop("value of 'c1' must be >= 0")
    if (!is.numeric(digits) || digits < 0) 
      stop("value of 'digits' must be >= 0")
    list(var.eps = var.eps, max.steps = max.steps, conv.eps = conv.eps, 
         conv.stop = conv.stop, digits = digits, c1 = c1)
  }



get.Amat <- function (initial.beta = NULL, penalty = NULL, intercept = TRUE, 
          c1 = lqa.control()$c1, x = NULL) 
{
  if (intercept) 
    x <- x[, -1]
  if (is.null(initial.beta)) 
    stop("get.Amat: 'initial.beta' is missing.")
  if (is.null(penalty)) 
    stop("get.Amat: 'penalty' is missing.")
  coefm0 <- if (intercept) 
    initial.beta[-1,]
  else initial.beta
  A.lambda <- penalty$getpenmat(beta = coefm0)
  if (intercept) {
    A.lambda <- cbind(0, A.lambda)
    A.lambda <- rbind(0, A.lambda)
  }
  A.lambda
}

Legendre<-function( t, np.order=1,tmin=NULL, tmax=NULL )
{
  u <- -1
  v <- 1
  if (is.null(tmin)) tmin<-min(t)
  if (is.null(tmax)) tmax<-max(t)
  nt <- length(t)
  ti    <- u + ((v-u)*(t-tmin))/(tmax - tmin)
  np.order.mat <- matrix(rep(0,nt*np.order),nrow=nt)
  if(np.order >=1)
    np.order.mat[,1] <- rep(1,nt)
  if (np.order>=2)
    np.order.mat[,2] <- ti
  if (np.order>=3)
    np.order.mat[,3] <- 0.5*(3*ti*ti-1)
  if (np.order>=4)
    np.order.mat[,4] <- 0.5*(5*ti^3-3*ti)
  if (np.order>=5)
    np.order.mat[,5] <- 0.125*(35*ti^4-30*ti^2+3)
  if (np.order>=6)
    np.order.mat[,6] <- 0.125*(63*ti^5-70*ti^3+15*ti)
  if (np.order>=7)
    np.order.mat[,7] <- (1/16)*(231*ti^6-315*ti^4+105*ti^2-5)
  if (np.order>=8)
    np.order.mat[,8] <- (1/16)*(429*ti^7-693*ti^5+315*ti^3-35*ti)
  if (np.order>=9)
    np.order.mat[,9] <- (1/128)*(6435*ti^8-12012*ti^6+6930*ti^4-1260*ti^2+35)
  if (np.order>=10)
    np.order.mat[,10] <- (1/128)*(12155*ti^9-25740*ti^7+18018*ti^5-4620*ti^3+315*ti)
  if (np.order>=11)
    np.order.mat[,11] <- (1/256)*(46189*ti^10-109395*ti^8+90090*ti^6-30030*ti^4+3465*ti^2-63)
  return(np.order.mat)
}