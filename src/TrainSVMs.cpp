#include "svmModels.h"
#include "TrainCanSVMs.h"
#include "TrainSVMs.h"

RcppExport
SEXP TrSVMs(SEXP grouping_s, SEXP lambda_s, SEXP K_s, SEXP W_s, SEXP epsilon_s, SEXP maxiter_s, SEXP tol_s,
            SEXP model_s, SEXP retotpst_s, SEXP start_s , SEXP alpha0_s)
{
  const IntegerVector grouping0(grouping_s);
  const arma::uvec grouping(Rcpp::as<arma::uvec>(grouping0));
  const double lambda(as<double>(lambda_s)),epsilon(as<double>(epsilon_s)),tol(as<double>(tol_s));
  const NumericMatrix W0(W_s);
  int n(grouping.size()),k(W0.ncol()),m(n*(k-1)),nsvms(W0.nrow());
  const arma::mat W(W0.begin(),nsvms,k);

  const bool normalize(true),retotpst(as<int>(retotpst_s));
  const int maxiter(as<int>(maxiter_s)),start(as<int>(start_s));
  const int model(as<int>(model_s));
  const NumericVector alpha00(alpha0_s);
  const arma::vec alpha0(Rcpp::as<arma::vec>(alpha00));

  arma::vec alpha(m),v(m),u(m);
  arma::mat Alpha(nsvms,m);
  int iter(0);

  double fval(0.);
  bool succes(false);

  if (model==LLW) {
    const NumericVector K(K_s);
    LLWCMatrix Q(&grouping,n,k);
    LLWQ_C(grouping,lambda,K,W.row(0).t(),v,Q,u);
    succes = TrSVMs_C<NumericVector>(grouping,lambda,K,W,v,&Q,u,normalize,epsilon,maxiter,tol,alpha,Alpha,fval,iter,model,retotpst,start,alpha0);
  } else if (model==WW) {          // Note: this is a temporary code -- the real Q compact implementation of the WW model is yet to be implemented !!!
    const NumericVector K(K_s);
    arma::mat Q0(m,m);
    ARMASVMMatrix Q(Q0,n,k);
    WWQ_C(grouping,lambda,K,W.row(0).t(),v,Q,u);

    succes = TrSVMs_C<NumericVector>(grouping,lambda,K,W,v,&Q,u,normalize,epsilon,maxiter,tol,alpha,Alpha,fval,iter,model,retotpst,start,alpha0);
  }
  if (!succes) return R_NilValue;

  if (nsvms==1)
  {
    if (!retotpst)  return List::create(Named("alpha")=alpha);
    else {
      return List::create(
        Named("alpha")=alpha,
        Named("fval")=fval,
        Named("iter")=iter
      );
    }
  }
  else {
    return wrap(Alpha);
  }
}
