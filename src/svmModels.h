#ifndef SVMMODELS_H
#define SVMMODELS_H

#include "SVMMat.h"

void LLWQ_C(const arma::uvec& grouping, const double lambda, const arma::mat& K, const arma::vec& w, arma::vec& v, ARMASVMMatrix& Q, arma::vec& u);
void LLWQ_C(const arma::uvec& grouping, const double lambda, const NumericVector& K, const arma::vec& w, arma::vec& v, LLWCMatrix& Q, arma::vec& u);
void LLWQ_C(const arma::uvec& grouping, const double lambda, const arma::mat& K, const arma::vec& w, arma::vec& v, LLWCMatrix& Q, arma::vec& u);
void WWQ_C(const arma::uvec& grouping, const double lambda, const arma::mat& K, const arma::vec& w, arma::vec& v, ARMASVMMatrix& Q, arma::vec& u);
void WWQ_C(const arma::uvec& grouping, const double lambda, const NumericVector& K, const arma::vec& w, arma::vec& v, ARMASVMMatrix& Q, arma::vec& u);
void WWQ_C(const arma::uvec& grouping, const double lambda, const NumericVector& K, const arma::vec& w, arma::vec& v, LLWCMatrix& Q, arma::vec& u);
void WWQ_C(const arma::uvec& grouping, const double lambda, const arma::mat& K, const arma::vec& w, arma::vec& v, LLWCMatrix& Q, arma::vec& u);

#endif
