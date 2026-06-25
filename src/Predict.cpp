#include "Predict.h"

template<>
void predict<arma::mat,arma::mat>(const arma::mat& alpha, const arma::mat& K, arma::uvec& classpred)
{
  int n1=K.n_rows,k=alpha.n_cols,bstc;
  double score,maxscor(-infinity);
  arma::mat scores = K * alpha;
  arma::uvec::iterator pos=classpred.begin();
  for (int obs=0;obs<n1;obs++,pos++) {
    maxscor = scores(obs,bstc=0);
    for (int c=1;c<k;c++)  {
      if ( (score=scores(obs,c)) > maxscor )  {
        bstc = c ;
        maxscor = score;
      }
    }
    *pos = bstc;
  }

  return;
}

template<>
void predict<arma::mat,NumericVector>(const arma::mat& alpha, const NumericVector& K, arma::uvec& classpred)
{
  int n=classpred.size(),k=alpha.n_cols,bstc;
  double score,maxscor(-infinity);
  arma::mat scores(n,k,arma::fill::zeros);

  for (int Ki=0,Krow=0,Kcol=0;Ki<K.size();Ki++) {
    for (int c=0;c<k;c++) scores(Krow,c) += K(Ki)*alpha(Kcol,c);
    if (Krow!=Kcol) for (int c=0;c<k;c++) scores(Kcol,c) += K(Ki)*alpha(Krow,c);
    if (Krow==n-1) Krow = ++Kcol;
    else Krow++;
  }
  arma::uvec::iterator pos=classpred.begin();
  for (int obs=0;obs<n;obs++,pos++) {
    maxscor = scores(obs,bstc=0);
    for (int c=1;c<k;c++)  {
      if ( (score=scores(obs,c)) > maxscor )  {
        bstc = c ;
        maxscor = score;
      }
    }
    *pos = bstc;
  }

  return;
}

double CLhitrate(const arma::uvec& orgclass, arma::uvec& predclass)
{
  int n(orgclass.size()),err(0);
  //  assert(predclass,size()==n);
  arma::uvec::iterator j=predclass.begin();
  for (int i=0;i<n;i++,j++) if (orgclass(i)!=*j) err++;
  return 1. - static_cast<double>(err)/n;
}

double CLhitrate(const unsigned orgclass, arma::uvec& predclass)
{
  //  assert(predclass.size()==1);
  if (orgclass!=predclass(0)) return 1.;
  return 0.;
}
