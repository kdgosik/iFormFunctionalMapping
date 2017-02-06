library(foreach)
library(doParallel)
library(iterators)
library(Rcpp)

iForm_FunctionalMap <- function(formula, data, heredity = "strong", higher_order = FALSE) {

  dat <- model.frame(formula, data)
  y <- dat[ , 1]
  x <- dat[ , -1]
  p <- ncol(x)
  n <- nrow(x)
  C <- names(x)
  S <- NULL
  M <- NULL
  bic <- NULL

  fit <- iformselect_FunctionalMap(x, y, p, n, C, S, bic, heredity, higher_order)

  y <- fit$y
  S <- fit$S
  bic <- fit$bic

  model_formula <- as.formula(paste("y ~0+", paste(S[1:which.min(bic)], collapse = "+")))
  lm(model_formula, data = x)

}






iformselect_FunctionalMap <- function( x, y, p, n, C, S, bic, heredity, higher_order ) {

  repeat{

    RSS <- rss_map_func(C = C, S = S, y = y, data = x)

    S <- c(S, C[which.min(unlist(RSS))])
    C <- C[-which.min(unlist(RSS))]

      order2 <- switch( heredity,
              `none` = NULL,
              `strong` = strong_order2(S = S, data = x),
              `weak` = weak_order2(S = S, C = C, data = x)
              )

      C <- union(C, order2)

      if( higher_order ) {

        order3 <- switch( heredity,
                `strong` = strong_order3( S = S, data = x ),
                `weak` = weak_order3(S = S, C = C, data = x)
                )

        C <- union(C, order3)

      }

    bic_val <- log(min(unlist(RSS))/n) + length(S) * (log(n) + 2 * log(p))/n
    bic <- append(bic, bic_val)
    if(length(bic) > 20) break

  }

  list(y = y, S = S, bic = bic)
  stopCluster(cl); gc(reset = TRUE)
}





## Help Functions ############################

logistic_legendre_fit <- function(parameters, formula, data, time_col){
  
  a <- parameters[1]
  b <- parameters[2]
  r <- parameters[3]
  
  mf <- model.frame(formula, data)
  y <- mf[, 1, drop = FALSE]
  X <- mf[, -1, drop = FALSE]
  
  p <- ncol(X)
  
  for( i in 1 : p ){
    eval(
      parse(
        text = paste0("beta", i, "<-parameters[", (4*i), ":", (4*i+3), "]")
      )
    )
  }
  
  t <- data[, time_col] # time variable
  t_s <- (t - min(t))/(max(t) - min(t))
  L <- legendre(3, t_s)
  
  eval(
    parse(
      text = paste0(
        c("y_hat <- a/(1+b*exp(-r*t))", paste0("(beta", 1:p, "%*%L)*X[,", 1:p, "]")), collapse = "+")
    )
  )
  
  sum((y - y_hat)^2)
  
}


minfval_map_func <- function(C, S, y, data, time_col){
  params <- rep(1, 3 + 4 * (length(S) + 1))
  lapply(C, function(candidates){
    var_names <- c(S, candidates)
    
    form <- as.formula(paste(y, "~0+", paste(var_names, collapse = "+")))
    
    out <- fminsearch(logistic_legendre_fit, params, 
                      formula = form,
                      data = data,
                      time_col = time_col)
    out
    
  })
}


rss_map_func <- function( C, S, y, data ) {

  sapply(C, function(candidates) {
    var_names <- c(S, candidates)

    X <- model.matrix(as.formula(paste("~0+", paste(var_names, collapse = "+"))), data = data)

    tryCatch({
      sum((y - X %*% (solve(t(X) %*% X)) %*% (t(X) %*% y)) ^ 2)
    }, error = function(e) Inf)


  })

}

## fastRSS in Rcpp

# adapted for just RSS calculation
cppFunction(depends = 'RcppArmadillo',
            'List fastRSS(NumericVector yr, NumericMatrix Xr) {
            int n = Xr.nrow(), k = Xr.ncol();
            arma::mat X(Xr.begin(), n, k, false);
            arma::colvec y(yr.begin(), yr.size(), false);
            arma::colvec coef = arma::solve(X, y);
            arma::colvec resid = y - X*coef;
            double rss = pow(arma::norm(resid), 2);

            return List::create(Named("RSS") = rss);
            }')

rss_map_func_cpp <- function( C, S, y, data ) {

  sapply(C, function(candidates) {
    var_names <- c(S, candidates)

    X <- model.matrix(as.formula(paste("~0+", paste(var_names, collapse = "+"))), data = data)

    tryCatch({
      fastRSS(y, X)
    }, error = function(e) Inf)


  })

}


rss_map_func_parallel <- function( C, S, y, data, no_cores = 2 ) {

  cl <- makeCluster(no_cores)
  registerDoParallel(cl)

  foreach(candidate = C, .combine = "c") %dopar%{
    var_names <- c(S, candidate)

    X <- model.matrix(as.formula(paste("~0+", paste(var_names, collapse = "+"))), data = data)

    tryCatch({
      sum((y - X %*% (solve(t(X) %*% X)) %*% (t(X) %*% y)) ^ 2)
    }, error = function(e) Inf)


  }

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


