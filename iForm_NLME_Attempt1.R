require(ff)
require(nlme)
require(lattice)
require(reshape2)

## iForm with ff package ##########
iForm <- function(formula, data, strong = TRUE, higher.order = FALSE){
  
  dat <- model.frame(formula, data)
  y <- dat[ , 1]
  x <- dat[ , -1]
  p <- ncol(x)
  n <- nrow(x)
  candidate <- as.ff(as.matrix(x), colnames = colnames(x), overwrite = TRUE)
  solution <- NULL
  model <- NULL
  step <- 1
  bic <- NULL
  
  fit <- iformselect(x, y, p, n, candidate, solution, model, bic, step, strong, higher.order)
  y <- fit$y
  solution <- fit$solution
  model <- fit$model
  bic <- fit$bic
  
  
  model <- data.frame(solution[ , 1 : which.min(bic), drop = FALSE])
  lm(y ~ ., data = model)
  
}



iformselect <- function(x, y, p, n, candidate, solution, model, bic, step, strong, higher.order){
  
  repeat{
    
    rss <- sapply(colnames(candidate), function(i){
      xx <- as.matrix(cbind(solution, candidate[ , i]))
      tryCatch({
        sum((y - xx %*% (solve(t(xx) %*% xx)) %*% (t(xx) %*% y)) ^ 2)
      }, error = function(e) Inf)
    })
    
    rss <- as.numeric(rss)
    
    solution <- cbind(solution, candidate[ , which.min(rss)])
    colnames(solution)[step] <- colnames(candidate)[which.min(rss)]
    cnames <- setdiff(colnames(candidate), colnames(solution))
    
    candidate[, which.min(rss) : (ncol(candidate) - 1)] <- candidate[, (which.min(rss) + 1) : ncol(candidate)]
    dim(candidate) <- c(n, (ncol(candidate) - 1))
    colnames(candidate) <- cnames
    
    if( ncol(solution) > 1 ) {
      if( strong ){
        interaction.formula <- paste("~0+(", paste(colnames(solution)[colnames(solution) %in% colnames(x)], collapse= "+"),")^2")
      }else{
        interaction.formula <- paste("~0+(", paste(colnames(solution)[colnames(solution) %in% colnames(x)], collapse= "+"),")*(", paste(colnames(x),collapse = "+"),")")
      }
      
      if(higher.order){
        if(sum(!{colnames(solution) %in% colnames(x)}) > 2) {
          
          if(strong){
            tmp <- Filter(function(x) {length(x) == 2}, strsplit(colnames(solution), "[.]|[:]"))
            tmp <- combn(length(tmp), 3, function(x) Reduce(intersect, combn(x, 2, function(y) Reduce(union, tmp[y]), simplify=F)), simplify=F)
            tmp <- Map(function(x) {paste0(x, collapse=":")}, Filter(function(x){length(x)==3}, tmp))
            if(length(tmp) > 0) {interaction.formula <- paste(interaction.formula, paste0(unlist(tmp), collapse = "+"), sep = "+")}
            
          }else{
            tmp <- Filter(function(x) {length(x) == 2}, strsplit(colnames(solution), "[.]|[:]"))
            tmp <- Map(function(x) {paste0(x, collapse=":")}, Filter(function(x){length(x)==2}, tmp))
            if(length(tmp) > 0) {interaction.formula <- paste(interaction.formula, paste("(",paste(unlist(tmp), collapse="+"), ")*(", paste(colnames(x), collapse = "+"), ")"), sep = "+")}
          }
          
        }
      }
      
      interaction.formula <- as.formula(interaction.formula)
      interactions <- model.matrix(interaction.formula, x)
      interactions <- interactions[ , setdiff(colnames(interactions), colnames(solution)), drop = F]
      interactions <- interactions[ , setdiff(colnames(interactions), colnames(candidate)), drop = F]
      
      if( ncol(interactions) > 0 ){
        newdimstart <- ncol(candidate) + 1
        newdimend <- ncol(candidate) + ncol(interactions)
        dim(candidate) <- c(n, newdimend)
        candidate[, newdimstart : newdimend] <- interactions
        colnames(candidate) <- c(cnames, colnames(interactions))
      }
      
    }
    
    bic[step] <- log(as.numeric(min(rss))/n) + (dim(solution)[2]) * (log(n) + 2 * log(p))/n
    if((dim(solution)[2]) > (n / log(n, 2))) break
    step <- step + 1
    
  }
  
  list(y = y, solution = solution, model = model, bic = bic)
}

######



iForm_lme <- function(formula, data, strong = TRUE, higher.order = FALSE){
  
  dat <- model.frame(formula, data)
  y <- dat[ , 1]
  x <- dat[ , -1]
  p <- ncol(x)
  n <- nrow(x)
  candidate <- as.ff(as.matrix(x), colnames = colnames(x), overwrite = TRUE)
  solution <- NULL
  model <- NULL
  step <- 1
  bic <- NULL
  
  fit <- iformselect(x, y, p, n, candidate, solution, model, bic, step, strong, higher.order)
  y <- fit$y
  solution <- fit$solution
  model <- fit$model
  bic <- fit$bic
  
  
  model <- data.frame(solution[ , 1 : which.min(bic), drop = FALSE])
  lm(y ~ ., data = model)
  
}



iformselect_lme <- function(x, y, p, n, candidate, solution, model, bic, step, strong, higher.order){
  
  repeat{
    
    rss <- sapply(colnames(candidate), function(i){
      xx <- as.matrix(cbind(solution, candidate[ , i]))
      tryCatch({
        sum((y - xx %*% (solve(t(xx) %*% xx)) %*% (t(xx) %*% y)) ^ 2)
      }, error = function(e) Inf)
    })
    
    rss <- as.numeric(rss)
    
    solution <- cbind(solution, candidate[ , which.min(rss)])
    colnames(solution)[step] <- colnames(candidate)[which.min(rss)]
    cnames <- setdiff(colnames(candidate), colnames(solution))
    
    candidate[, which.min(rss) : (ncol(candidate) - 1)] <- candidate[, (which.min(rss) + 1) : ncol(candidate)]
    dim(candidate) <- c(n, (ncol(candidate) - 1))
    colnames(candidate) <- cnames
    
    if( ncol(solution) > 1 ) {
      if( strong ){
        interaction.formula <- paste("~0+(", paste(colnames(solution)[colnames(solution) %in% colnames(x)], collapse= "+"),")^2")
      }else{
        interaction.formula <- paste("~0+(", paste(colnames(solution)[colnames(solution) %in% colnames(x)], collapse= "+"),")*(", paste(colnames(x),collapse = "+"),")")
      }
      
      if(higher.order){
        if(sum(!{colnames(solution) %in% colnames(x)}) > 2) {
          
          if(strong){
            tmp <- Filter(function(x) {length(x) == 2}, strsplit(colnames(solution), "[.]|[:]"))
            tmp <- combn(length(tmp), 3, function(x) Reduce(intersect, combn(x, 2, function(y) Reduce(union, tmp[y]), simplify=F)), simplify=F)
            tmp <- Map(function(x) {paste0(x, collapse=":")}, Filter(function(x){length(x)==3}, tmp))
            if(length(tmp) > 0) {interaction.formula <- paste(interaction.formula, paste0(unlist(tmp), collapse = "+"), sep = "+")}
            
          }else{
            tmp <- Filter(function(x) {length(x) == 2}, strsplit(colnames(solution), "[.]|[:]"))
            tmp <- Map(function(x) {paste0(x, collapse=":")}, Filter(function(x){length(x)==2}, tmp))
            if(length(tmp) > 0) {interaction.formula <- paste(interaction.formula, paste("(",paste(unlist(tmp), collapse="+"), ")*(", paste(colnames(x), collapse = "+"), ")"), sep = "+")}
          }
          
        }
      }
      
      interaction.formula <- as.formula(interaction.formula)
      interactions <- model.matrix(interaction.formula, x)
      interactions <- interactions[ , setdiff(colnames(interactions), colnames(solution)), drop = F]
      interactions <- interactions[ , setdiff(colnames(interactions), colnames(candidate)), drop = F]
      
      if( ncol(interactions) > 0 ){
        newdimstart <- ncol(candidate) + 1
        newdimend <- ncol(candidate) + ncol(interactions)
        dim(candidate) <- c(n, newdimend)
        candidate[, newdimstart : newdimend] <- interactions
        colnames(candidate) <- c(cnames, colnames(interactions))
      }
      
    }
    
    bic[step] <- log(as.numeric(min(rss))/n) + (dim(solution)[2]) * (log(n) + 2 * log(p))/n
    if((dim(solution)[2]) > (n / log(n, 2))) break
    step <- step + 1
    
  }
  
  list(y = y, solution = solution, model = model, bic = bic)
}
