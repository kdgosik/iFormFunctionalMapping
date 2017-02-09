library(foreach)
library(doParallel)
library(iterators)
library(Rcpp)
library(pracma)


# formula <- y ~ .
# data = df[c(1:8,53)]
# id <- "id"
# time_col <- "t"
# heredity = "strong"
# higher_order = FALSE

iForm_FunctionalMap <- function(formula, 
                                data, 
                                id, 
                                time_col, 
                                heredity = "strong", 
                                higher_order = FALSE) {
  
  formula_vars <- all.vars(formula)
  response <- formula_vars[1]
  t <- data[, time_col, drop = FALSE]
  p <- ncol(data) - 3
  n <- nrow(data)
  C <- names(data)[!{names(data) %in% c(id, time_col, formula_vars)}]
  S <- NULL
  M <- NULL
  bic <- NULL
  output_list <- NULL
  step <- 1
  
  a_hat <- max(data[, reponse])
  g_out <- fminsearch(g1, c(a_hat, 1, 1),
                      data = data,
                      response = response,
                      time_col = time_col)
  
  mu_t <- g_out$xval[1]/(1 + (g_out$xval[2]) * exp(-(g_out$xval[3]) * t))
  
  X <- matrix(mu_t, ncol = 1)
  
  
  repeat{
    
    # min_out <- minfval_map_parallel_func(C = C, 
    #                                      S = S, 
    #                                      response = response, 
    #                                      data = data, 
    #                                      time_col = time_col,
    #                                      no_cores = no_cores)
    # output_list[[step]] <- min_out
    
    RSS <-    rss_map_func(C = C,
                           S = S,
                           response = response,
                           data = data,
                           time_col = time_col,
                           design_mat = X)
    
    # sapply(seq_along(min_out), function(i) min_out[[i]]$fval)
    
    ## update names of design matrix?
    form <- as.formula(paste(response, "~0+", C[which.min(unlist(RSS))]))
    d <- drop(model.matrix(form, data))
    X <- cbind(X, Legendre(t, 3) * d)
    
    S <- c(S, C[which.min(unlist(RSS))])
    C <- C[-which.min(unlist(RSS))]
    
    order2 <- switch( heredity,
                      `none` = NULL,
                      `strong` = strong_order2(S = S, data = data),
                      `weak` = weak_order2(S = S, C = C, data = data)
    )
    
    C <- union(C, order2)
    
    if( higher_order ) {
      
      order3 <- switch( heredity,
                        `strong` = strong_order3( S = S, data = data),
                        `weak` = weak_order3(S = S, C = C, data = data)
      )
      
      C <- union(C, order3)
      
    }
    
    bic_val <- log(min(unlist(RSS))/n) + length(S) * (log(n) + 2 * log(p))/n
    bic <- append(bic, bic_val)
    if(length(bic) > 10) break
    step <- step + 1
  }
  
  out <- output_list[[which.min(bic)]]
  rss_out <- sapply(seq_along(out), function(i) out[[i]]$fval)
  out[[which.min(rss_out)]]
  
}


## Asymptotic Growth Equation
g1 <- function(parameters, data, response, time_col){
  y <- data[, response, drop = FALSE]
  ti <- data[, time_col, drop = FALSE]
  
  if(any(parameters < 0)){
    parameters[which(parameters < 0)]<- 0.5
  }
  a1 <- parameters[1];b1 <- parameters[2];r1 <- parameters[3]
  y_hat <- a1/(1+b1*exp(-r1*ti))
  
  sum( (y - y_hat )^2 )
}

## solving the parameters by psuedo least squares
g_out <- fminsearch(g1, c(150, 1, 1),
           data = df,
           response = "y",
           time_col = "t")

# mu_t <- g_out$xval[1]/(1 + (g_out$xval[2]) * exp(-(g_out$xval[3]) * df$t))
# X <- matrix(mu_t, ncol = 1)



## rss mapping function
rss_map_func <- function(C, S, response, data, time_col, design_mat){
  
  ti <- data[, time_col]
  

  sapply(C, function(candidates){

    form <- as.formula(paste(response, "~0+", candidates))
    d <- drop(model.matrix(form, data))
    X <- cbind(design_mat, Legendre(ti, 3) * d)
    
    sum((y - X %*% (solve(t(X) %*% X)) %*% (t(X) %*% y)) ^ 2) 
    
  })
  
}



## Help Functions ############################

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




## Heredity Selection

# strong order 2
strong_order2 <- function(S, data) {
  
  tryCatch({
    
    main_effects <- sort(S[S %in% names(data)])
    combn(main_effects, 2, paste0, collapse = ":")
    
  }, error = function(e) NULL)
  
}


# weak order 2
weak_order2 <- function( S, C, data ) {
  
  tryCatch({
    
    main_effects <- sort(S[S %in% names(data)])
    as.vector(outer(main_effects, C, paste, sep = ":"))
    
  }, error = function(e) NULL)
  
}


# strong order 3
strong_order3 <- function( S, data ) {
  
  tryCatch({
    
    main_effects <- sort(S[S %in% names(data)])
    combn(main_effects, 3, paste0, collapse = ":")
    
  }, error = function(e) NULL)
  
}


# weak order 3
weak_order3 <- function( S, C, data ) {
  
  tryCatch({
    
    interaction_effects <- unlist(
      Map(function(int_term) paste0(int_term, collapse = ":"),
          Filter(function(vec) {length(vec) == 2}, strsplit(S, "[.]|[:]"))
      )
    )
    
    weak_three <- as.vector(
      outer(interaction_effects, C, paste, sep = ":")
    )
    
    as.vector(
      unlist(
        Map(paste0, collapse = ":",
            Map(sort, strsplit(weak_three, ":"))
        )
      )
    )
    
  }, error = function(e) NULL)
}


