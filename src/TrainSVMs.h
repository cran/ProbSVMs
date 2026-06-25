#ifndef TRAINSVMS_H
#define TRAINSVMS_H

#include "SVMMat.h"
#include "Losses.h"
#include "StartOpt.h"

using namespace Rcpp;

RcppExport
SEXP TrSVMs(SEXP grouping_s, SEXP lambda_s, SEXP K_s, SEXP W_s, SEXP epsilon_s, SEXP maxiter_s, SEXP tol_s,
            SEXP model_s, SEXP retotpst_s, SEXP start_s , SEXP alpha0_s);

template<class Ktype>
bool TrSVMs_C(const arma::uvec& grouping, const double lambda, const Ktype& K, const arma::mat& W, arma::vec& v, SVMCMatrix* Qp, arma::vec& u,
              const bool normalize, const double epsilon0, const int maxiter, const double tol0,
              arma::vec& alpha, arma::mat& Alpha, double& fval, int& iter, const int model, const bool retotpst, const int start, const arma::vec alpha0)
{
  const int k(W.n_cols),n(grouping.size()),m(n*(k-1)),nsvms(W.n_rows);
//  double vcoef;
  bool succes;
  std::list<int> VatLB,VatUB,VatDomInt;
  bool Qnormalized(true);
  arma::vec dQ(m),g(m);
  arma::vec::iterator dQpos;
  double maxdQ(0.),dQi;

  if (start==WARMSTARTS) {
    double alphai,ui;
    int i(1);
    arma::vec::iterator alphapos=alpha.begin();
    arma::vec::const_iterator upos=u.begin();
    for (arma::vec::const_iterator alpha0pos=alpha0.begin();alpha0pos<alpha0.end();alpha0pos++,upos++,alphapos++,i++) {
      if ( (alphai=*alpha0pos) < 0. || alphai>(ui=*upos)) {
        Rprintf("Error: Initial alpha(%i) value (= %f) outside of its allowable range given by [0. %f]\n",i,alphai,ui);
        return false;
      }
      *alphapos = alphai;
    }
  }

  if (normalize) {
    dQpos=dQ.begin();
    for (int i=0;i<m&&Qnormalized;i++,dQpos++) {
      *dQpos = dQi = sqrt((*Qp)(i,i));
      if (dQi > maxdQ) maxdQ = dQi;
      if (dQi < tol0) Qnormalized = false;
    }
    if (Qnormalized) {
      Qp->cov2cor(dQ);
      dQpos=dQ.begin();
      arma::vec::iterator vpos=v.begin(),alphapos=alpha.begin();
      for (arma::vec::iterator upos=u.begin();upos<u.end();dQpos++,upos++,vpos++,alphapos++) {
          if ( (dQi = *dQpos) !=1. ) {
          if (start==WARMSTARTS) *alphapos *= dQi;
          *vpos /= dQi;
          *upos *= dQi;
        }
      }
    }
  }  else  Qnormalized = false;

//  int totaliter = 0;
  
  for (int i=0;i<nsvms;i++) {

    if (i>0) {
//      if (start==WARMSTARTS) alpha = clone(alpha0);
      if (start==WARMSTARTS) alpha = alpha0;
      arma::vec normw = W.row(i).t() / n;
      if (Qnormalized) dQpos=dQ.begin();
      arma::vec::iterator upos=u.begin();
      for (int j=0;j<n;j++) {
        int y = grouping(j);
        for (int c=0;c<k;c++) if (c!=y) {
          *upos = normw(y);
          if (Qnormalized) {
            if (start==WARMSTARTS) alpha(i) *= *dQpos;
            *upos *= *dQpos++;
          }
          upos++;
        }
      }
    }

    VatLB.clear();
    VatUB.clear();
    VatDomInt.clear();
    double minu,uj;
    if (Qnormalized) minu = maxdQ/n;
    else minu = 1./n;
    for (arma::vec::iterator upos=u.begin();upos<u.end();upos++) if ( (uj=*upos) < minu )  minu = uj;
    double tol=tol0*minu,epsilon=epsilon0*minu;
    fval = 0.;
//  if (start==ALLUB || start==COMPUBWITHLB) alpha = clone(u);
    if (start==ALLUB || start==COMPUBWITHLB) alpha = u;

    if (start!=ALLLB) {
      double vi,gi;
      arma::vec::iterator gpos=g.begin();
      arma::vec::const_iterator vpos=v.begin();
      for (int i=0;i<m;i++,gpos++,vpos++) {
        gi = vi = *vpos;
        arma::vec::iterator alphajpos=alpha.begin();
        for (int j=0;j<m;j++,alphajpos++) gi -= *alphajpos * (*Qp)(j,i);
        fval += (vi+(*gpos=gi)) * alpha(i);
      }
      fval /= 2;
    }

    if (start==ALLUB || (start==COMPUBWITHLB && fval>0.) )  for (int i=0;i<m;i++) VatUB.push_back(i);
    else if ( (start==ALLLB) || (start==COMPUBWITHLB && fval<=0.) ) {
      g = v;
      for (int i=0;i<m;i++) {
         alpha(i) = 0.;
         VatLB.push_back(i);
      }
      fval = 0.;
    } else if (start==WARMSTARTS) {
      double alphai;
      arma::vec::const_iterator upos=u.begin();
      for (int i=0;i<m;i++,upos++) {
        if ( (alphai=alpha(i)) < tol ) VatLB.push_back(i);
        else if ( alphai > *upos-tol ) VatUB.push_back(i);
             else VatDomInt.push_back(i);
      }
    }

    iter = 0;
    succes = TrainCanSVM_C(v, Qp, u, VatLB, VatUB, VatDomInt, Qnormalized, epsilon, maxiter, tol, alpha, g, fval, iter);
    if (!succes) return false;
    if (nsvms==1) {
      if (Qnormalized) alpha = alpha/dQ;
    } else {
      if (Qnormalized) Alpha.row(i) = (alpha / dQ).t();
      else Alpha.row(i) = alpha;
    }
//    totaliter += iter;

  }

  return true;

}

#endif
