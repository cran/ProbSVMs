#include "makeKMat.h"
 #include <math.h>

RcppExport
SEXP makeGauKMat(const SEXP traindat_s, const SEXP sigma_s)
{
  const NumericMatrix traindat(traindat_s);    
  const double sigma(as<double>(sigma_s)),sigma2=sigma*sigma;
  const int n=traindat.nrow(), p=traindat.ncol();

  NumericVector KMat(n*(n+1)/2);  
//  arma::vec KMat(n*(n+1)/2);
  for (int col=0,KMind=0;col<n;col++) for (int row=col;row<n;row++,KMind++)
  {
    if (row==col) KMat(KMind) = 1.;
    else {
      double l2norm = 0.;
      for (int j=0;j<p;j++) {
        double dif = traindat(row,j) - traindat(col,j);
        l2norm += dif*dif;
      }
      KMat(KMind) = exp(-l2norm/sigma2);
    }
  }

  return wrap(KMat);
};


RcppExport
SEXP makeGauKMat2(const SEXP traindat_s, const SEXP newdat_s, const SEXP sigma_s)
{
  const NumericMatrix traindat(traindat_s);
  const NumericMatrix newdat(newdat_s);
  const double sigma(as<double>(sigma_s)),sigma2=sigma*sigma;
  const int n=traindat.nrow(), p=traindat.ncol(), nnwdt=newdat.nrow();

  NumericMatrix KMat(nnwdt,n);
//  arma::mat KMat(nnwdt,n);  
  for (int i=0;i<nnwdt;i++) for (int j=0;j<n;j++) {
    double l2norm = 0.;
    for (int k=0;k<p;k++) {
      double dif = traindat(j,k)-newdat(i,k);
      l2norm += dif*dif;
    }
    KMat(i,j) = exp(-l2norm/sigma2);
  }

  return wrap(KMat);
};

RcppExport
SEXP makeLinKMat(const SEXP traindat_s)
{
  const NumericMatrix traindat(traindat_s);
  const int n=traindat.nrow(), p=traindat.ncol();

  NumericVector KMat(n*(n+1)/2);
//  arma::vec KMat(n*(n+1)/2); 
  for (int col=0,KMind=0;col<n;col++) for (int row=col;row<n;row++,KMind++)
  {
    double l2norm = 0.;
    for (int j=0;j<p;j++) l2norm += traindat(row,j) * traindat(col,j);
    KMat(KMind) = l2norm;
  }

  return wrap(KMat);  
};

RcppExport
SEXP makeLinKMat2(const SEXP traindat_s, const SEXP newdat_s)
{
  const NumericMatrix traindat(traindat_s);
  const NumericMatrix newdat(newdat_s);
  const int n=traindat.nrow(), p=traindat.ncol(), nnwdt=newdat.nrow();

  NumericMatrix KMat(nnwdt,n);
//  arma::mat KMat(nnwdt,n);
  for (int i=0;i<nnwdt;i++) for (int j=0;j<n;j++) {
    double l2norm = 0.;
    for (int k=0;k<p;k++)  l2norm += traindat(j,k)*newdat(i,k);
    KMat(i,j) = l2norm;
  }

  return wrap(KMat);
};


RcppExport
SEXP makePolKMat(const SEXP traindat_s, SEXP degree_s, SEXP scale_s, SEXP offset_s)
{
  const NumericMatrix traindat(traindat_s);
  const int n=traindat.nrow(), p=traindat.ncol();
  const int degree(as<int>(degree_s));
  double scale(as<double>(scale_s)),offset(as<double>(offset_s));

  NumericVector KMat(n*(n+1)/2);
//   arma::vec KMat(n*(n+1)/2); 
  for (int col=0,KMind=0;col<n;col++) for (int row=col;row<n;row++,KMind++)
  {
    double l2norm = 0.;
    for (int j=0;j<p;j++) l2norm += traindat(row,j) * traindat(col,j);
    KMat(KMind) = pow(fabs(scale*l2norm+offset),degree);
  }

  return wrap(KMat);
};

RcppExport
SEXP makePolKMat2(const SEXP traindat_s, const SEXP newdat_s, SEXP degree_s, SEXP scale_s, SEXP offset_s)
{
  const NumericMatrix traindat(traindat_s);
  const NumericMatrix newdat(newdat_s);
  const int n=traindat.nrow(), p=traindat.ncol(), nnwdt=newdat.nrow();
  const int degree(as<int>(degree_s));
  double scale(as<double>(scale_s)),offset(as<double>(offset_s));

  NumericMatrix KMat(nnwdt,n);
//  arma::mat KMat(nnwdt,n);
  for (int i=0;i<nnwdt;i++) for (int j=0;j<n;j++) {
    double l2norm = 0.;
    for (int k=0;k<p;k++)  l2norm += traindat(j,k)*newdat(i,k);
    KMat(i,j) = pow(fabs(scale*l2norm+offset),degree);
  }

  return wrap(KMat);
};
