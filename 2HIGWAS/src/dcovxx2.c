#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <math.h> 


void dcovXX2(int *n,double *snp,double *tmp){
   double tmps1xx,tmps2xx,tmps3xx;
   double z,z1;
   double S1xx,S2xx,S3xx,dcovXX2;
   tmps1xx=0;
   for(int i=0;i<*n;i++){
     z=0;
     for(int k=0;k<*n;k++){
       z= z + pow(abs(snp[i]-snp[k]),2);
     }
     tmps1xx = tmps1xx + z;
   }
  
  S1xx = tmps1xx/pow(*n,2);
  
  tmps2xx = 0;
  for(int i=0;i<*n;i++){
    z = 0;
    for(int k=0;k<*n;k++){
      z= z + abs(snp[i]-snp[k]);
    }
    tmps2xx = tmps2xx + z;
  }  
  S2xx = (tmps2xx/pow(*n,2))*(tmps2xx/pow(*n,2));
  
  tmps3xx = 0;
  for(int i=0;i<*n;i++){
    z=0;
    for(int k=0;k<*n;k++){
      z1=0;
      for(int l=0;l<*n;l++){
        z1 = z1 + abs(snp[i]-snp[l])*abs(snp[k]-snp[l]);
      }
      z = z + z1;
    }
   tmps3xx = tmps3xx + z;
  }
  S3xx = tmps3xx/pow(*n,3);
  
  dcovXX2 = S1xx + S2xx -2*S3xx;
  
  *tmp=dcovXX2;
}



void dcovYY2(int *m,int *n,double *y,double *tmp){
   double tmps1yy,tmps2yy,tmps3yy;
   double z,z1;
   int mm = *m;
   int nn = *n;
   double S1yy,S2yy,S3yy,dcovYY2;
   double Enorm(int i,int k,int nn,int mm,double y[]);
   tmps1yy=0;
   for(int i=0;i<*n;i++){
     z=0;
     for(int k=0;k<*n;k++){
       z= z + pow(Enorm(i,k,nn,mm,y),2);
     }
     tmps1yy = tmps1yy + z;
   }
  
  S1yy = tmps1yy/pow(*n,2);
  
  tmps2yy =0;
  for(int i=0;i<*n;i++){
    z = 0;
    for(int k=0;k<*n;k++){
      z = z + Enorm(i,k,nn,mm,y);
    }
    tmps2yy = tmps2yy + z;
  }
  
  S2yy = tmps2yy/pow(*n,2)*tmps2yy/pow(*n,2);
  
  tmps3yy =0;
  for(int i=0;i<*n;i++){
    z =0 ;
    for(int k=0;k<*n;k++){
      z1 =0 ;
      for(int l=0;l<*n;l++){
        z1 = z1 + Enorm(i,l,nn,mm,y)* Enorm(k,l,nn,mm,y);
      }
      z = z + z1;
    }
    tmps3yy = tmps3yy + z;
  }
  S3yy = tmps3yy/pow(*n,3);
  dcovYY2 = S1yy + S2yy -2*S3yy;
  *tmp=dcovYY2;
}

void dcovXY2(int *m,int *n,double *y,double *snp, double *tmp){
   double tmps1xy,tmps2xy1,tmps2xy2,tmps3xy;
   double z,z1;
   int mm = *m;
   int nn = *n;
   double S1xy,S2xy,S3xy,dcovXY2;
   double Enorm(int i,int k,int nn,int mm,double y[]);
   tmps1xy=0;
   for(int i=0;i<*n;i++){
     z=0;
     for(int k=0;k<*n;k++){
       z= z + Enorm(i,k,nn,mm,y)*abs(snp[i]-snp[k]);
     }
     tmps1xy = tmps1xy + z;
   }
  
  S1xy = tmps1xy/pow(*n,2);
  
 tmps2xy1 =0;
  for(int i=0;i<*n;i++){
    z =0;
    for(int k=0;k<*n;k++){
     z = z + Enorm(i,k,nn,mm,y);
    }
    tmps2xy1 = tmps2xy1 + z;
  }  
  tmps2xy2 = 0;
  for(int i=0;i<*n;i++){
    z =0;
    for(int k=0;k<*n;k++){
      z = z + abs(snp[i]-snp[k]);
    }
    tmps2xy2 = tmps2xy2 + z;
  } 
  S2xy = (tmps2xy1/pow(*n,2))*(tmps2xy2/pow(*n,2));
  
  tmps3xy =0;
  for(int i=0;i<*n;i++){
    z =0 ;
    for(int k=0;k<*n;k++){
      z1 =0 ;
      for(int l=0;l<*n;l++){
        z1 = z1 + Enorm(i,l,nn,mm,y)*abs(snp[k]-snp[l]);
      }
      z = z + z1;
    }
    tmps3xy = tmps3xy + z;
  }
  S3xy = tmps3xy/pow(*n,3);
  dcovXY2 = S1xy + S2xy -2*S3xy;
  *tmp=dcovXY2;
}


double Enorm(int i,int k,int nn,int mm,double y[]){
  
  double res =0;
  int ii;
  for(ii =0;ii < mm;ii++){
   res = res + pow((y[i+nn*ii]-y[k+nn*ii]),2);
  }
  
  return(sqrt(res));
}


