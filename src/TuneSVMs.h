#ifndef TUNESVM_H
#define TUNESVM_H

#include <assert.h>
#include "GridSearch.h"
#include "Predict.h"
#include "TrainSVMs.h"
#include "TuneAux.h"

RcppExport
SEXP TunSVMs(SEXP trgroup_s, SEXP tungroup_s, SEXP trK_s, SEXP tunK_s,  SEXP W_s, SEXP TrTunInd_s,
             SEXP epsilon_s, SEXP maxiter_s, SEXP tol_s, SEXP model_s,
             SEXP Csrchpar_s, SEXP retotpst_s, SEXP start_s, SEXP alpha0_s);

template<class KType, class SVMMType>
class EvalKrnSvms
{
public:

      EvalKrnSvms(const arma::uvec& trgroup, const arma::uvec& tungroup,
                  const KType* trKp, const arma::mat* tunKp, const arma::mat& W,
                  const bool normalize, const int  start, const arma::vec* alpha0p,
                  const int k,  const int n, const int maxiter, const int model, const double epsilon, const double tol,
                  arma::vec* vp, SVMMType* Qp, arma::vec& u, arma::vec& alpha, int iter, arma::uvec* classpredp) :
                    trgroup_(trgroup), tungroup_(tungroup), trKp_(trKp), tunKp_(tunKp), W_(W),
                    normalize_(normalize), start_(start), alpha0p_(alpha0p),
                    k_(k), n_(n), maxiter_(maxiter), model_(model), epsilon_(epsilon), tol_(tol),
                    vp_(vp), Qp_(Qp), u_(u), alpha_(alpha), iter_(iter), classpredp_(classpredp)
  {  }

        ~EvalKrnSvms() { }

        double operator () (double C, srchres* resp)  {

          double fval(0.),lambda(1./(2*n_*C)),hitrate(0.);
          if (model_==LLW) LLWQ_C(trgroup_,lambda,*trKp_,W_.row(0).t(),*vp_,*Qp_,u_);
          else if (model_==WW) WWQ_C(trgroup_,lambda,*trKp_,W_.row(0).t(),*vp_,*Qp_,u_);
          arma::mat NullMat;

          bool succes = TrSVMs_C<KType>(trgroup_,lambda,*trKp_,W_,*vp_,Qp_,u_,normalize_,epsilon_,maxiter_,tol_,alpha_,NullMat,fval,iter_,
                                        model_,false,start_,*alpha0p_);
          if (!succes) return -infinity;
          (*resp).fval = fval;
          (*resp).iter = iter_;
          if (tungroup_.size() > 0) {
            predict<arma::mat>(alpha_, trgroup_, *tunKp_, (*resp).alpha,  *classpredp_, model_);
            hitrate = (*resp).hitrate = CLhitrate(tungroup_, *classpredp_);
          } else {
            predict<KType>(alpha_, trgroup_, *trKp_, (*resp).alpha,  *classpredp_, model_);
            hitrate = (*resp).hitrate = CLhitrate(trgroup_, *classpredp_);
          }
          return hitrate;

        }

private:

      const arma::uvec trgroup_, tungroup_;
      const KType *trKp_;
      const arma::mat* tunKp_;
      const arma::mat W_;
      const bool normalize_;
      const int start_ ;
      const arma::vec* alpha0p_;
      const int k_, n_, maxiter_, model_;
      const double epsilon_, tol_;
      arma::vec* vp_;
      SVMMType* Qp_;
      arma::vec u_;
      arma::vec alpha_;
      int iter_;
      arma::uvec* classpredp_;

};

template<class SVMMType>
class CrossVEvalKrnSvms
{
public:

    CrossVEvalKrnSvms(const arma::uvec& trgroup, const arma::mat* trKp, const arma::mat& W,
                      const std::vector<arma::uvec>* trIndLp, const std::vector<arma::uvec>* tunIndLp,
                      const bool normalize, const int  start,
                      const int k,  const int n, const int maxiter, const int model, const double epsilon, const double tol,
                      arma::vec* vp, SVMMType* Qp, SVMMType* Q1p, arma::vec& u, arma::vec& alpha, int iter, arma::uvec* classpredp) :
                      
                        trgroup_(trgroup), trKp_(trKp), W_(W), trIndLp_(trIndLp), tunIndLp_(tunIndLp),
                        normalize_(normalize), start_(start),
                        k_(k), n_(n), maxiter_(maxiter), model_(model), epsilon_(epsilon), tol_(tol),
                        vp_(vp), Qp_(Qp), Q1p_(Q1p), u_(u), alpha_(alpha), iter_(iter), classpredp_(classpredp)
                    { if (start_ == WARMSTARTS) start_ = ALLUB; }

        ~CrossVEvalKrnSvms() { }

        double operator () (double C, srchres* resp)  {

          double fval(0.),lambda(1./(2*n_*C)),hitratessum(0.);
          int nrep(trIndLp_->size());
          for (int rep=0;rep<nrep;rep++) {
            arma::uvec curtrInd = (*trIndLp_)[rep]-1, curtunInd = (*tunIndLp_)[rep]-1;
            arma::uvec curtrgroup = trgroup_.elem(curtrInd), curtungroup = trgroup_.elem(curtunInd);
            const unsigned cutTrsize(curtrInd.size()),m1((k_-1)*cutTrsize);
            if (v1_.size() != m1) {
              v1_.resize(m1);
              u1_.resize(m1);
              alpha1_.resize(m1);
              alpha1asmat_.resize(cutTrsize,k_);
              classpred1_.resize(curtunInd.size());
              Q1p_->resize(&curtrgroup,m1,k_);
            }
            arma::mat curtrK = trKp_->submat(curtrInd,curtrInd), curtunK = trKp_->submat(curtunInd,curtrInd);
            if (model_==LLW) LLWQ_C(curtrgroup,lambda,curtrK,W_.row(0).t(),v1_,*Q1p_,u1_);
            else if (model_==WW) WWQ_C(curtrgroup,lambda,curtrK,W_.row(0).t(),v1_,*Q1p_,u1_);
            arma::mat NullMat;
            arma::vec NullVec;
            bool succes = TrSVMs_C<arma::mat>(curtrgroup,lambda,curtrK,W_,v1_,Q1p_,u1_,normalize_,epsilon_,maxiter_,tol_,
                                               alpha1_,NullMat,fval,iter_,model_,false,start_,NullVec);
             if (!succes) return -infinity;
             predict<arma::mat>(alpha1_, curtrgroup, curtunK, alpha1asmat_,  classpred1_, model_);
             double hitrate = CLhitrate(curtungroup, classpred1_);
            hitratessum += hitrate;
          }

          return (*resp).hitrate=hitratessum/nrep;
        }

private:

      const arma::uvec trgroup_;
      const arma::mat *trKp_;
      const arma::mat W_;
      const std::vector<arma::uvec> *trIndLp_,*tunIndLp_;
      const bool normalize_;
      int start_ ;
      const int k_, n_, maxiter_, model_;
      const double epsilon_, tol_;
      arma::vec* vp_;
      SVMMType *Qp_, *Q1p_;
      arma::vec u_,alpha_;
      int iter_;
      arma::uvec* classpredp_;
      
      arma::vec v1_,u1_,alpha1_;
      arma::mat alpha1asmat_;
      arma::uvec classpred1_;
};

#endif
