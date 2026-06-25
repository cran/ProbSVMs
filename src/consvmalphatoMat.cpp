#include "consvmalphatoMat.h"

void LLWconsvmalphatoMat(const arma::vec& svmalpha, const arma::uvec& Y, arma::mat& alpha)
{
  unsigned n=alpha.n_rows,k=alpha.n_cols;
  for (unsigned i=0,obs=0;obs<n;obs++) for (unsigned c=0;c<k;c++)
    if (c==Y(obs)) alpha(obs,c) = 0.;
    else alpha(obs,c) = -svmalpha(i++);
  return;
}

void WWconsvmalphatoMat(const arma::vec& svmalpha, const arma::uvec& Y, arma::mat& alpha)
{
  unsigned n=alpha.n_rows,k=alpha.n_cols;
  for (unsigned i=0,obs=0;obs<n;obs++) {
    for (unsigned c=0,c1=0;c<k;c++) {
      if (c==Y(obs)) {
        alpha(obs,c) = svmalpha(i);
        for (unsigned c2=1;c2<k-1;c2++) alpha(obs,c) += svmalpha(i+c2);
        alpha(obs,c) /= 2;
      }
      else alpha(obs,c) = -svmalpha(i+c1++)/2;
    }
    i += k-1;
  }
  return;
}
