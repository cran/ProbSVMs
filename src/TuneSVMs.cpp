#include "svmModels.h"
#include "TrainCanSVMs.h"
#include "TrainSVMs.h"
#include "TuneSVMs.h"

RcppExport
SEXP TunSVMs(SEXP trgroup_s, SEXP tungroup_s, SEXP trK_s, SEXP tunK_s,  SEXP W_s, SEXP TrTunInd_s,
              SEXP epsilon_s, SEXP maxiter_s, SEXP tol_s, SEXP model_s,
              SEXP Csrchpar_s, SEXP retotpst_s, SEXP start_s , SEXP alpha0_s)
{
  const IntegerVector trgroup0(trgroup_s),tungroup0(tungroup_s);
  const arma::uvec trgroup(Rcpp::as<arma::uvec>(trgroup0)),tungroup(Rcpp::as<arma::uvec>(tungroup0));
  const List Csrchpar(Csrchpar_s),TrTunInd(TrTunInd_s);
  double Cpowerbase(as<double>(Csrchpar["Cpowerbase"]));
  int Cgridinlev(as<int>(Csrchpar["Cgridinlev"])),Cnloops(as<int>(Csrchpar["Cnloops"]));
  const double epsilon(as<double>(epsilon_s)),tol(as<double>(tol_s));
  const bool normalize(true),retotpst(as<int>(retotpst_s));
  const int model(as<int>(model_s)),maxiter(as<int>(maxiter_s)),start(as<int>(start_s));
  const NumericVector alpha00(alpha0_s);
  const arma::vec alpha0(Rcpp::as<arma::vec>(alpha00));

  const NumericMatrix W0(W_s);
  int n(trgroup.size()),n1(tungroup.size()),k(W0.ncol()),m(n*(k-1));
  if (W0.nrow() > 1) {
    Rprintf("Error: TunSVMs only works with one set of weight specifications\n");
    return R_NilValue;
  }
  const arma::mat W(W0.begin(),1,k);
  double bestC;
  arma::vec alpha(m),v(m),u(m);
  arma::mat NullMat;
  int iter(0);
  double fval(0.);

  const NumericMatrix tunK0(tunK_s);
  if (n1!=tunK0.nrow()) {
    Rprintf("Error: size of the tungroup vector does not agree with the number of columns of tunK0\n");
    return R_NilValue;
  }

  const arma::mat tunK(tunK0.begin(),n1,n);
  arma::mat curnewalpha(n,k),bstnewalpha(n,k);
  arma::uvec classpred(n);
  srchres curres={curnewalpha, 0., 0, 0.}, bstres={bstnewalpha, 0., 0, 0.};

  std::vector<arma::uvec>  TrnInd,TunInd;
  if (TrTunInd.size()>0) {
//    List TrIndL(TrTunInd["trIndL"]),TunIndL(TrTunInd["tunIndL"]);
    List TrIndL(as<List>(TrTunInd["trIndL"])),TunIndL(as<List>(TrTunInd["tunIndL"]));
    int nrep=TrIndL.size();
//    IntegerMatrix TrIndIM(as<IntegerMatrix>(TrTunInd["trIndL"])),TunIndIM(as<IntegerMatrix>(TrTunInd["tunIndL"]));
//    int nrep=TrIndIM.nrow();
    TrnInd.resize(nrep);
    TunInd.resize(nrep);
    for (int rep=0;rep<nrep;rep++) {
        IntegerVector TrIndrepasIV = TrIndL(rep);
//        IntegerVector TrIndrepasIV = TrIndIM(rep,_);
        TrnInd[rep] = Rcpp::as<arma::uvec>(TrIndrepasIV);
        IntegerVector TunIndrepasIV = TunIndL(rep);
//        IntegerVector TunIndrepasIV = TunIndIM(rep,_);
        TunInd[rep] = Rcpp::as<arma::uvec>(TunIndrepasIV);
      }
  }

  bool succes(false);
  if (model==LLW) {
    const NumericVector trK(trK_s);
    if (TrTunInd.size()==0) {
      LLWCMatrix Q(&trgroup,n,k);
      EvalKrnSvms<NumericVector,LLWCMatrix>  EvalSvms(trgroup,tungroup,&trK,&tunK,W,normalize,start,&alpha0,k,n,maxiter,model,epsilon,tol,
                                                      &v,&Q,u,alpha,iter,&classpred);
      bestC = CSearch< EvalKrnSvms<NumericVector,LLWCMatrix>, srchres >(Cgridinlev, Cnloops, Cpowerbase, EvalSvms, &curres, &bstres);
      double lambda = 1/(2*n*bestC);
      iter = 0; fval = 0.;
      succes = TrSVMs_C<NumericVector>(trgroup,lambda,trK,W,v,&Q,u,normalize,epsilon,maxiter,tol,alpha,NullMat,fval,iter,model,retotpst,start,alpha0);
    } else {
      LLWCMatrix Q(&trgroup,n,k),Q1(&trgroup,n,k);
      arma::mat trK1(n,n);
      cnvVct2armamat(trK,trK1);
      CrossVEvalKrnSvms<LLWCMatrix>  EvalSvms(trgroup,&trK1,W,&TrnInd,&TunInd,
                                              normalize,start,k,n,maxiter,model,epsilon,tol,
                                              &v,&Q,&Q1,u,alpha,iter,&classpred);
      bestC = CSearch< CrossVEvalKrnSvms<LLWCMatrix>, srchres >(Cgridinlev, Cnloops, Cpowerbase, EvalSvms, &curres, &bstres);
      double lambda = 1/(2*n*bestC);
      LLWQ_C(trgroup,lambda,trK,W.row(0).t(),v,Q,u);
      iter = 0; fval = 0.;
      succes = TrSVMs_C<NumericVector>(trgroup,lambda,trK,W,v,&Q,u,normalize,epsilon,maxiter,tol,alpha,NullMat,fval,iter,model,retotpst,start,alpha0);
    }
  } else if (model==WW) {   // Note: this is a temporary code -- the real Q compact implementation of the WW model is yet to be implemented !!!
    const NumericVector trK(trK_s);
    arma::mat Q0(m,m);
    ARMASVMMatrix Q(Q0,n,k);
    if (TrTunInd.size()==0) {
      EvalKrnSvms<NumericVector,ARMASVMMatrix>  EvalSvms(trgroup,tungroup,&trK,&tunK,W,normalize,start,&alpha0,k,n,maxiter,model,epsilon,tol,
                                                         &v,&Q,u,alpha,iter,&classpred);
      bestC = CSearch< EvalKrnSvms<NumericVector,ARMASVMMatrix>, srchres >(Cgridinlev, Cnloops, Cpowerbase, EvalSvms, &curres, &bstres);
// Rprintf("bestC = %f\n",bestC);      
      double lambda = 1/(2*n*bestC);
      iter = 0; fval = 0.;
      succes = TrSVMs_C<NumericVector>(trgroup,lambda,trK,W,v,&Q,u,normalize,epsilon,maxiter,tol,alpha,NullMat,fval,iter,model,retotpst,start,alpha0);
    } else {
      arma::mat Q10(m,m);
      ARMASVMMatrix Q1(Q10,n,k);
      arma::mat trK1(n,n);
      cnvVct2armamat(trK,trK1);
      CrossVEvalKrnSvms<ARMASVMMatrix>  EvalSvms(trgroup,&trK1,W,&TrnInd,&TunInd,
                                                    normalize,start,k,n,maxiter,model,epsilon,tol,
                                                    &v,&Q,&Q1,u,alpha,iter,&classpred);
      bestC = CSearch< CrossVEvalKrnSvms<ARMASVMMatrix>, srchres >(Cgridinlev, Cnloops, Cpowerbase, EvalSvms, &curres, &bstres);
      double lambda = 1/(2*n*bestC);
      WWQ_C(trgroup,lambda,trK,W.row(0).t(),v,Q,u);
      iter = 0; fval = 0.;
      succes = TrSVMs_C<NumericVector>(trgroup,lambda,trK,W,v,&Q,u,normalize,epsilon,maxiter,tol,alpha,NullMat,fval,iter,model,retotpst,start,alpha0);
      if (model==LLW)  LLWconsvmalphatoMat(alpha,trgroup,bstnewalpha);
      else if (model==WW)  WWconsvmalphatoMat(alpha,trgroup,bstnewalpha);
    }

  }

  if (!succes) return R_NilValue;
  
  if (TrTunInd.size()==0) {
//    if (!retotpst)  return List::create(Named("alpha")=bstres.alpha);
    if (!retotpst)  return List::create(
        Named("bestC") = bestC,
        Named("alpha")=bstres.alpha
       );
    else {
      return List::create(
        Named("bestC") = bestC,
        Named("alpha") = bstres.alpha,
        Named("fval") = bstres.fval,
        Named("iter") = bstres.iter,
        Named("hitrate") = bstres.hitrate
      );
    }
  } else {
//    if (!retotpst)  return List::create(Named("alpha")=alpha);
    if (!retotpst)  return List::create(
        Named("bestC") = bestC,
//        Named("alpha")=alpha
        Named("alpha")=bstnewalpha
       );
    else {
      return List::create(
        Named("bestC") = bestC,
//        Named("alpha")=alpha,
        Named("alpha")=bstnewalpha,
        Named("fval")=fval,
        Named("iter")=iter,
        Named("hitrate") = bstres.hitrate
      );
    }
  }
}

void cnvVct2armamat(const arma::vec& v, arma::mat& M)
  {
  int n(M.n_rows);
  double val;
  for (int vind=0,col=0;col<n;col++) {
    for (int row=col;row<n;row++,vind++)  {
      M(row,col) = val = v(vind);
      if (row>col) M(col,row) = val;
    }
  }
  return;
}
