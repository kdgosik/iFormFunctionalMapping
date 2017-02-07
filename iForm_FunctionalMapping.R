library(foreach)
library(doParallel)
library(iterators)
library(Rcpp)


# formula <- y ~ .
# data = df[c(1:8,53)]
# id <- "id"
# time_col <- "t"
# heredity = "strong"
# higher_order = FALSE

iForm_FunctionalMap <- function(formula, data, id, time_col, heredity = "strong", higher_order = FALSE) {

  formula_vars <- all.vars(formula)
  response <- formula_vars[1]
  p <- ncol(data) - 3
  n <- nrow(data)
  C <- names(data)[!{names(data) %in% c(id, time_col, formula_vars)}]
  S <- NULL
  M <- NULL
  bic <- NULL
  output_list <- NULL
  step <- 1

  
  repeat{
    
    min_out <- minfval_map_func(C = C, S = S, response = response, data = data, time_col = time_col)
    output_list[[step]] <- min_out
    
    RSS <-sapply(seq_along(min_out), function(i) min_out[[i]]$fval)
    
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
    if(length(bic) > 5) break
    step <- step + 1
  }
  
  out <- output_list[[which.min(bic)]]
  rss_out <- sapply(seq_along(out), function(i) out[[i]]$fval)
  out[[which.min(rss_out)]]

}





# probably delete

# iformselect_FunctionalMap <- function( x, y, p, n, C, S, bic, heredity, higher_order ) {
# 
#   repeat{
# 
#     RSS <- rss_map_func(C = C, S = S, y = y, data = x)
# 
#     S <- c(S, C[which.min(unlist(RSS))])
#     C <- C[-which.min(unlist(RSS))]
# 
#       order2 <- switch( heredity,
#               `none` = NULL,
#               `strong` = strong_order2(S = S, data = x),
#               `weak` = weak_order2(S = S, C = C, data = x)
#               )
# 
#       C <- union(C, order2)
# 
#       if( higher_order ) {
# 
#         order3 <- switch( heredity,
#                 `strong` = strong_order3( S = S, data = x ),
#                 `weak` = weak_order3(S = S, C = C, data = x)
#                 )
# 
#         C <- union(C, order3)
# 
#       }
# 
#     bic_val <- log(min(unlist(RSS))/n) + length(S) * (log(n) + 2 * log(p))/n
#     bic <- append(bic, bic_val)
#     if(length(bic) > 20) break
# 
#   }
# 
#   list(y = y, S = S, bic = bic)
# }





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


minfval_map_func <- function(C, S, response, data, time_col){
  
  a_hat <- max(data[,response])
  params <- c(a_hat, rep(1, 2 + 4 * (length(S) + 1)))
  lapply(C, function(candidates){
    var_names <- c(S, candidates)
    
    form <- as.formula(paste(response, "~0+", paste(var_names, collapse = "+")))
    
    out <- fminsearch(logistic_legendre_fit, params, 
                      formula = form,
                      data = data,
                      time_col = time_col)
    names(out$xval)<- c("a", "b", "r", as.vector(t(outer(var_names, paste0("_P", 0:3), paste0))))
    out
    
  })
  
}



minfval_map_parallel_func <- function(C, S, response, data, time_col, no_cores = 2){
  
  a_hat <- max(data[,response])
  params <- c(a_hat, rep(1, 2 + 4 * (length(S) + 1)))
  
  cl <- makeCluster(no_cores)
  registerDoParallel(cl)
  
  foreach(candidates = C, .packages = "pracma", .export = c("a_hat", "params")) %dopar% {
    var_names <- c(S, candidates)
    
    form <- as.formula(paste(response, "~0+", paste(var_names, collapse = "+")))
    
    out <- fminsearch(logistic_legendre_fit, params, 
                      formula = form,
                      data = data,
                      time_col = time_col)
    names(out$xval)<- c("a", "b", "r", as.vector(t(outer(var_names, paste0("_P", 0:3), paste0))))
    out
    
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


