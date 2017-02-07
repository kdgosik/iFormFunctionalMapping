DC.sis <- function(dat){
  
  y <- dat$pheno_y#aa$phe0#
  snps <- t(dat$snps)#gene_table#t(dat$snps)
  #snps <- 1 - abs(snps - 1)
  nn <- dim(y)[1]
  #m <- dim(y)[2]
  ################dcov^2(y,y)###################
  #dcovYY_2 <- dcovYY2(n,y)
  snpdc <- c()
  for(j in 1:dim(snps)[2]){
    missing <- which(is.na(snps[,j]))
    if(length(missing)>0){
      yy <- y[-missing,]
      snp <- snps[-missing,j]
    }
    else {
      yy <- y
      snp <- snps[,j]
    }
    n <- dim(yy)[1]
    m <- dim(yy)[2]
    yy <- as.matrix(yy)
    ################dcov^2(y,x_j)###################
    #dcovYX_2 <- dcovYX2(n,y,snp)
    dcovYX_2 <- .C("dcovXY2",as.integer(m),as.integer(n),as.numeric(yy),as.numeric(snp),as.numeric(0))[[5]]
    dcovYY_2 <- .C("dcovYY2",as.integer(m),as.integer(n),as.numeric(yy),as.numeric(0))[[4]]
    ################dcov^2(x_j,x_j)###################
    #dcovXX_2 <- dcovXX2(n,snp)
    dcovXX_2 <- .C("dcovXX2",as.integer(n),as.numeric(snp),as.numeric(0))[[3]]
    
    dcorryx <- sqrt(dcovYX_2)/sqrt(sqrt(dcovYY_2)*sqrt(dcovXX_2))
    snpdc <- c(snpdc,dcorryx)
  }
  save(snpdc,file="snpdc.RData")
  names(snpdc) <- colnames(snps)
  d <- nn/log(nn)
  rr <- sort(snpdc,decreasing=TRUE)[1:round(d)]
  rrnames <- names(rr)
  sissnp <- c()
  ii.id <- c()
  for(i in 1:length(rr)){
    ii <- which(rr[i]==snpdc)
    sissnp <- cbind(sissnp,snps[,ii])
    ii.id <- c(ii.id,ii)
  }
  colnames(sissnp) <- rrnames
  return(list(rr=rr,sisnp=sissnp,iid=ii.id))
}
  

dcovXX2 <- function(n,snp){
  tmps1xx <- c()
  for(i in 1:n){
    tmp1 <- c()
    for(k in 1:n){
      res <- abs(snp[i]-snp[k])^2
      tmp1 <- c(tmp1,res)
    }
    tmps1xx <- c(tmps1xx,sum(tmp1))
  }
  S1xx <- sum(tmps1xx)/n^2
  
  tmps2xx <- c()
  for(i in 1:n){
    tmp2 <- c()
    for(k in 1:n){
      res <- abs(snp[i]-snp[k])
      tmp2 <- c(tmp2,res)
    }
    tmps2xx <- c(tmps2xx,sum(tmp2))
  }  
  S2xx <- (sum(tmps2xx)/n^2)^2
  
  tmps3xx <- c()
  for(i in 1:n){
    tmp4 <- c()
    for(k in 1:n){
      tmp44 <- c()
      for(l in 1:n){
        res <-abs(snp[i]-snp[l])*abs(snp[k]-snp[l])
        tmp44 <- c(tmp44,res)
      }
      tmp4 <- c(tmp4,sum(tmp44))
    }
    tmps3xx <- c(tmps3xx,sum(tmp4))
  }
  S3xx <- sum(tmps3xx)/n^3
  
  dcovXX2 <- S1xx + S2xx -2*S3xx
  return(dcovXX2)
  
}
  
  
dcovYY2 <- function(n,y){
  
  tmps1yy <- c()
  for(i in 1:n){
    tmp1 <- c()
    for(k in 1:n){
      res <- Enorm(y[i,] - y[k,])*Enorm(y[i,] - y[k,])
      tmp1 <- c(tmp1,res)
    }
    tmps1yy <- c(tmps1yy,sum(tmp1))
  }
  S1yy <- sum(tmps1yy)/n^2
  
  tmps2yy <- c()
  for(i in 1:n){
    tmp2 <- c()
    for(k in 1:n){
      res <- Enorm(y[i,] - y[k,])
      tmp2 <- c(tmp2,res)
    }
    tmps2yy <- c(tmps2yy,sum(tmp2))
  }  
  S2yy <- (sum(tmps2yy)/n^2)^2
  
  tmps3yy <- c()
  for(i in 1:n){
    tmp4 <- c()
    for(k in 1:n){
      tmp44 <- c()
      for(l in 1:n){
        res <- Enorm(y[i,] - y[l,])* Enorm(y[k,] - y[l,])
        tmp44 <- c(tmp44,res)
      }
      tmp4 <- c(tmp4,sum(tmp44))
    }
    tmps3yy <- c(tmps3yy,sum(tmp4))
  }
  S3yy <- sum(tmps3yy)/n^3
  
  dcovYY2 <- S1yy + S2yy -2*S3yy
  return(dcovYY2)
}
  
  
dcovYX2 <- function(n,y,snp){
  snp <- snp
  tmps1yx <- c()
  for(i in 1:n){
    tmp1 <- c()
    for(k in 1:n){
      res <- Enorm(y[i,] - y[k,])*abs(snp[i]-snp[k])
      tmp1 <- c(tmp1,res)
    }
    tmps1yx <- c(tmps1yx,sum(tmp1))
  }
  S1yx <- sum(tmps1yx)/n^2
  
  tmps2yx1 <- c()
  for(i in 1:n){
    tmp2 <- c()
    for(k in 1:n){
      res <- Enorm(y[i,] - y[k,])
      tmp2 <- c(tmp2,res)
    }
    tmps2yx1 <- c(tmps2yx1,sum(tmp2))
  }  
  tmps2yx2 <- c()
  for(i in 1:n){
    tmp3 <- c()
    for(k in 1:n){
      res <- abs(snp[i]-snp[k])
      tmp3 <- c(tmp3,res)
    }
    tmps2yx2 <- c(tmps2yx2,sum(tmp3))
  } 
  S2yx <- (sum(tmps2yx1)/n^2)*(sum(tmps2yx2)/n^2) 
  
  tmps3yx <- c()
  for(i in 1:n){
    tmp4 <- c()
    for(k in 1:n){
      tmp44 <- c()
      for(l in 1:n){
        res <- Enorm(y[i,] - y[l,])*abs(snp[k]-snp[l])
        tmp44 <- c(tmp44,res)
      }
      tmp4 <- c(tmp4,sum(tmp44))
    }
    tmps3yx <- c(tmps3yx,sum(tmp4))
  }
  S3yx <- sum(tmps3yx)/n^3
  
  dcovYX2 <- S1yx + S2yx -2*S3yx
  return(dcovYX2)
  
}
  
  
  
  

Enorm <- function(xx){
  
  er <- sqrt(sum(xx^2))
  return(er)
}