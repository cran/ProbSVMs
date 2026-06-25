#ifndef consvmalphatoMat_H
#define consvmalphatoMat_H

#include <RcppArmadillo.h>
using namespace Rcpp;

void LLWconsvmalphatoMat(const arma::vec& svmalpha, const arma::uvec& Y, arma::mat& alpha);
void WWconsvmalphatoMat(const arma::vec& svmalpha, const arma::uvec& Y, arma::mat& alpha);


#endif
