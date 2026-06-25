#ifndef TRAINCANSVMS_H
#define TRAINCANSVMS_H

#include "SVMMat.h"
using namespace Rcpp;

arma::ivec findWrkind(const arma::vec& g, const double epsilon, const double tol, const bool Qnormalized,
                      arma::vec& alpha, const arma::vec& u, SVMCMatrix* Qp,
                      std::list<int>& VatLB, std::list<int>& VatUB, std::list<int>& VatDomInt);

inline double normf2d(const double muB1, const double muB2, const double g1, const double g2, const double Q12)
{
  return (g1-(muB1+2*Q12*muB2)/2)*muB1 + (g2-muB2/2)*muB2;
}
inline double unormf2d(const double muB1, const double muB2, const double g1, const double g2, const double Q12, const double Q11, const double Q22)
{
  return (g1-(Q11*muB1+2*Q12*muB2)/2)*muB1 + (g2-Q22*muB2/2)*muB2;
}
inline double normf1d(const double muB1, const double g1)
{
  return (g1-muB1/2)*muB1;
}
inline double unormf1d(const double muB1, const double g1, const double Q11)
{
  return (g1-Q11*muB1/2)*muB1;
}
inline double normonedimsol(const double muB2, const double alpha1, const double u1, const double g1, const double Q12)
{
  double arg1(g1-Q12*muB2), arg2(arg1>-alpha1 ? arg1:-alpha1), arg3(u1-alpha1);
  return arg2<arg3 ? arg2:arg3;
}
inline double unormonedimsol(const double muB2, const double alpha1, const double u1, const double g1, const double Q12, const double Q11)
{
  double arg1((g1-Q12*muB2)/Q11), arg2(arg1>-alpha1 ? arg1:-alpha1), arg3(u1-alpha1);
  return arg2<arg3 ? arg2:arg3;
}
void normtwocntviol(bool l2scvio, double& muB1, double& muB2,
                    const double alpha1, const double alpha2, const double u1, const double u2,
                    const double g1, const double g2, const double Q12);
void unormtwocntviol(bool l2scvio, double& muB1, double& muB2,
                     const double alpha1, const double alpha2, const double u1, const double u2,
                     const double g1, const double g2, const double Q12, const double Q11, const double Q22);
double solve2dprob(arma::vec& muB, const bool Qnormalized,  const double tol, const double epsilon,
                    const double alpha1, const double alpha2, const double g1, const double g2,
                   const double Q12, const double u1, const double u2, const double Q11=1. ,const  double Q22=1.);
double solve1dprob(arma::vec& muB, const bool Qnormalized,  const double tol, const double epsilon,
                                      const double alpha1, const double g1, const double u1, const double Q11=1.);

RcppExport
SEXP TrainCanSVM(SEXP v_s, SEXP Q_s, SEXP u_s, SEXP alpha_s, SEXP VatLB_s, SEXP VatUB_s, SEXP VatDomInt_s,
                 SEXP epsilon_s, SEXP maxiter_s, SEXP tol_s);

bool TrainCanSVM_C(const arma::vec& v, SVMCMatrix* Qp, const arma::vec& u,
                   std::list<int>& VatLB, std::list<int>& VatUB, std::list<int>& VatDomInt,
                   const bool Qnormalized, const double epsilon, const int maxiter, const double tol,
                   arma::vec& alpha, arma::vec& g, double& fval, int& iter);

#endif
