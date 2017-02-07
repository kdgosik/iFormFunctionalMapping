require(mvtnorm)
require(MSBVAR)
require(parallel)
require(snow)
sim_param <- function(mu,sim.beta,nt)
{
  par <- list();
  class(par)<- "sim.par";
  
  par$model	 <- "long";
  par$simu_n   <- 200;
  # the number of total predictors
  par$simu_p   <- 500; 
  # pairwise correlation coefficient
  par$simu_rho <- 0; 
  par$simu_chr <- 1;
  
 
  par$simu_times <- c(1:nt); 
  par$simu <- mu;
  par$simu_sigma2 <- 0.01;
  #sim.beta[5,] <- Legendre(ti,np=4)%*%c(2.08, 1.77, -4.11, 1.09)/2
  #sim.beta[6,] <- Legendre(ti,np=4)%*%c(2.34, -0.4, 1.48, -9.43)/2
  #sim.beta[10,] <- Legendre(ti,np=4)%*%c(2.80, -4.5, 2.00,  0.00)/2
  #sim.beta[11,] <- Legendre(ti,np=4)%*%c(2.98, -4.27, 9.64, 2.85)/2
  #sim.beta[14,] <- Legendre(ti,np=4)%*%c(-2.09, -2.64, -3.81,  3.07)/2
  #sim.beta[8,] <- Legendre(ti,np=4)%*%c(2.53, -2.45, 5.42, -3.92)/2;
  #sim.beta[12,] <- Legendre(ti,np=4)%*%c(-3.53, -2.45, 1.42, -4.92)/2;
  # locations of significant SNPs, less than p
  #par$simu_a_posi  <- c(100,300,1000,10000,15000,40000,80000); 
  par$simu_a_posi  <- c(10,30,100,150,300,400,800,820,120,380);
  
  par$simu_a_effect<- sim.beta[c(1,2,3,4,5,7,9,10,13,14),]
  par$simu_d_posi  <- c(200,640,330,420,500);
  par$simu_d_effect<- sim.beta[c(6,8,11,12,15),]
  
  # the number of good predictors
  pos <- c(par$simu_a_posi, par$simu_d_posi);
  par$sig_p    <- length(unique(pos)); 
  par$cluster <- 2
  # times point
  
  return(par);
}







simulate_longdt <- function( par )
{
  dat<-list(
    model        = par$model,
    name         = "longdt",
    chr          = par$simu_chr,
    sample_times = 1,
    sample_n     = par$simu_n,
    sample_p 	 = par$simu_p,
    pheno_file   = "simu.pheno",
    geno_file    = "simu.geno",
    #--the follwoings are vector or matrix data.
    #--SNP raw data(include missing data and rare data)
    raw_snps     = NULL,
    #--SNP data filted by proc_rare_snp;
    snps         = NULL,
    #--individul id list
    subj_ids     = NULL,
    #--phenotype
    pheno_y      = NULL);
  
  class(dat) <- "dat.longdt";
  
  set.seed( 200 );
  
  dat$subj_ids <- paste("sub", c(1:par$simu_n),sep="");
  #--generate snp data
  source('2HIGWAS/ulgwas.sim.R')
  cl <- makeCluster(rep("localhost",2),type="SOCK")
  clusterEvalQ(cl,source("2HIGWAS/ulgwas.sim.R"))
  SNP <- clusterCall(cl, simu_geno,par);
  snpres <- c()
  for(i in 1:par$cluster){
    snpres <- cbind(snpres,SNP[[i]])
  }
  colnames(snpres)<-paste("snp", c(1:(par$simu_p*par$cluster)),sep="");
  rownames(snpres)<-paste("id", c(1:par$simu_n),sep="");
  dat$raw_snps <- t(snpres);
  dat$snps     <- dat$raw_snps;
  
  #-- generate SNP info
  snps.info    <- cbind(chr=1, position=c(1:(par$simu_p*par$cluster))); 
  rownames(snps.info) <- paste("snp", c(1:(par$simu_p*par$cluster)), sep="");
  dat$snps.info <- snps.info
  
  #--generate phenotype data
  r.pheno      <- simu_pheno( par, dat$raw_snps);
  dat$pheno_y <- r.pheno
  dat$sample_times <- par$simu_times
  
  ##output to console
  cat("** Simulation is successful and data object is returned.\n");
  stopCluster(cl)
  return(dat);
}








simu_geno <- function( par ){
  c <- qnorm(1/4);
  
  #--simulating correlated SNPs, needs library(mvtnorm)
  sigma_x <- matrix(seq(par$simu_rho,par$simu_rho,length=par$sig_p*par$sig_p), ncol=par$sig_p)
  sigma_x <- sigma_x + diag(seq(1-par$simu_rho,1-par$simu_rho,length=par$sig_p));
  
  x <- rmvnorm(n=par$simu_n, mean=seq(0, 0, length=par$sig_p), sigma=sigma_x, method="chol")
  #--x has simu_n rows, and sig_p columns, check the mean and var-cov,colMeans(x),var(x)
  
  # simulating the rest of uncorrelated SNPs
  #uncorr_SNPs <- c();
  #for (i in 1:par$simu_n) 
  #	uncorr_SNPs <- rbind(uncorr_SNPs,rnorm(par$simu_p - par$sig_p, mean = 0, sd = 1))
  # bind all SNPs and convert to -1, 0, 1
  #all_SNPs <- cbind(x,uncorr_SNPs);
  
  sig_pos <- unique(c(par$simu_a_posi, par$simu_d_posi));
  all_SNPs <- c();
  for (i in 1:par$simu_p) 
  {
    fi <- which(sig_pos==i);
    if (length(fi)>0 )
      all_SNPs <- cbind(all_SNPs, x[,fi[1]])
    else
      all_SNPs <- cbind(all_SNPs, rnorm(par$simu_n, mean = 0, sd = 1))
  }	
  
  for (i in 1:par$simu_n) 
  {
    for (j in 1:par$simu_p) 
    {
      tmp <- all_SNPs[i,j];
      all_SNPs[i,j] <- -1*(tmp<c)+(tmp>-c);
    }
  }
  
  colnames(all_SNPs)<-paste("snp", c(1:par$simu_p),sep="");
  rownames(all_SNPs)<-paste("id", c(1:par$simu_n),sep="");
  
  ##output to console
  cat("** Genotypical data is simulated successfully.\n");
  
  #--QQ:2, Qq:1, qq:0
  return (all_SNPs+1);
}

simu_pheno <- function( par, snp )
{
  genA <- t(snp)
  genD <- 1 - abs(genA-1);

  
  tmax <- max(par$simu_times)
  tmin <- min(par$simu_times)
  m <- tmax - tmin +1
  sigma0 <- diag(m)
  for (ii in tmin:tmax)
    for (jj in tmin:tmax)
      sigma0[ii,jj] <- par$simu_sigma2^abs( ii - jj );
  
  
  nnt <- length(par$simu_times)
  p <- dim(snp)[1]
  #betat0 <- par$simu_times/m
  #betat1 <- (1 - 2*par$simu_times/m)^2
  #betat3 <- sqrt(betat0)
  #betat7 <- exp(par$simu_times/(2*m))
  
  betaA <- matrix(rep(0,m*(p)),nrow=m)
  betaD <- matrix(rep(0,m*(p)),nrow=m)
  #beta[,1] <- betat0
  #beta[,2] <- betat1
  #beta[,4] <- betat3
  #beta[,8] <- betat7
  # create a vector of addtive effects
  for(i in 1:length(par$simu_a_posi)){
    betaA[,par$simu_a_posi[i]] <- par$simu_a_effect[i,][1:nnt]
  }
  
  for(i in 1:length(par$simu_d_posi)){
    betaD[,par$simu_d_posi[i]] <- par$simu_d_effect[i,][1:nnt]
  }
  
  Beta <- cbind(par$simu[1:nnt],betaA)
  Xa <- cbind(rep(1,m),genA)
  
  Beta1 <- cbind(rep(0,nnt),betaD)
  Xd <- cbind(rep(1,m),genD)

  gene_eff <- Xa%*%t(Beta) + Xd%*%t(Beta1)
  
  simu_Y <- matrix(rep(0,nnt*par$simu_n),nrow=par$simu_n)
  for(i in 1:par$simu_n){
    simu_Y[i,] <-  rmultnorm(1, gene_eff[i,],sigma0) 
  }
 
  ##output to console
  cat("** Phenotypical data is simulated successfully.\n");
  
  return(simu_Y);
}
