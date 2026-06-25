#include "TrainCanSVMs.h"
//#include <assert.h>

const double NOINC = 999e99;
const double NONPOSQMAT = 998e99;
const double SEMIPOSQMAT = 997e99;

RcppExport
SEXP TrainCanSVM(SEXP v_s, SEXP Q_s, SEXP u_s, SEXP alpha_s, SEXP VatLB_s, SEXP VatUB_s, SEXP VatDomInt_s,
                 SEXP epsilon_s, SEXP maxiter_s, SEXP tol_s)
{

    const NumericVector v_r(v_s), u_r(u_s);
    const IntegerVector VatLB0(VatLB_s), VatUB0(VatUB_s), VatDomInt0(VatDomInt_s);
    const bool Qnormalized(true);
    const int maxiter(as<int>(maxiter_s));

    const double epsilon(as<double>(epsilon_s)), tol(as<double>(tol_s));

    NumericVector alpha_r(alpha_s);

    const int m(alpha_r.size());
    arma::vec alpha(alpha_r.begin(),m),g(m);

    const arma::vec v(v_r.begin(),m);
    const arma::vec u(u_r.begin(),m);
    NumericMatrix Q0(Q_s);
    arma::mat Q1(Q0.begin(),m,m);
    ARMASVMMatrix Q(Q1,m);

    int iter(0);
    double fval(0.);

    std::list<int> VatLB,VatUB,VatDomInt;
    if (VatLB0(0)!=-1) for (int i=0;i<VatLB0.size();i++) VatLB.push_back(VatLB0(i));
    if (VatUB0(0)!=-1) for (int i=0;i<VatUB0.size();i++) VatUB.push_back(VatUB0(i));
    if (VatDomInt0(0)!=-1) for (int i=0;i<VatDomInt0.size();i++) VatDomInt.push_back(VatDomInt0(i));

    double vi,gi;
    arma::vec::iterator gpos=g.begin();
    arma::vec::const_iterator vpos=v.begin();
    for (int i=0;i<m;i++,gpos++,vpos++) {
      gi = vi = *vpos;
      arma::vec::iterator alphajpos=alpha.begin();
      for (int j=0;j<m;j++,alphajpos++) gi -= *alphajpos * Q(j,i);
      fval += (vi+(*gpos=gi)) * alpha[i];
    }
    fval /= 2;

    bool succes = TrainCanSVM_C(v, &Q, u, VatLB, VatUB, VatDomInt, Qnormalized, epsilon, maxiter, tol, alpha, g, fval, iter);
    if (!succes) return R_NilValue;

    return List::create(
      Named("alpha")=alpha,
      Named("fval")=fval,
      Named("iter")=iter
    );
}


bool TrainCanSVM_C(const arma::vec& v, SVMCMatrix* Qp, const arma::vec& u,
                   std::list<int>& VatLB, std::list<int>& VatUB, std::list<int>& VatDomInt,
                   const bool Qnormalized, const double epsilon, const int maxiter, const double tol,
                   arma::vec& alpha, arma::vec& g, double& fval, int& iter)
{
//    const int m=alpha.size();
    arma::vec mu(2);
    double inc;

    bool stop;
    if (maxiter>0) stop = false;
    else stop = true;
    while (!stop) {
      arma::ivec B = findWrkind(g,epsilon,tol,Qnormalized,alpha,u,Qp,VatLB,VatUB,VatDomInt);
      if (B[0]==-1) stop = true;
      else {
        iter++;
        if (B[1]!=-1) {
          if (Qnormalized) inc = solve2dprob(mu,true,tol,epsilon,alpha(B[0]),alpha(B[1]),g(B[0]),g(B[1]),(*Qp)(B[0],B[1]),u[B[0]],u[B[1]]);
          else inc = solve2dprob(mu,false,tol,epsilon,alpha(B[0]),alpha(B[1]),g(B[0]),g(B[1]),(*Qp)(B[0],B[1]),u[B[0]],u[B[1]],(*Qp)(B[0],B[0]),(*Qp)(B[1],B[1]));
//          if (inc==NONPOSQMAT || inc==SEMIPOSQMAT) return false;
//          if (inc==NOINC) stop = true;
          if (inc==NOINC || inc==NONPOSQMAT || inc==SEMIPOSQMAT ) stop = true;
          else {
            fval += inc;
            alpha(B[0]) += mu[0];
            alpha(B[1])  += mu[1];
            if (iter==maxiter) stop = true;
            else {
              int b(B[1]);
              double newalpha2(alpha(b)),u2(u[b]);
              if (newalpha2 <= tol) VatLB.push_back(b);
              else if (newalpha2 >= u2-tol) VatUB.push_back(b);
              else VatDomInt.push_back(b);
              const arma::vec B0mu0(Qp->colscalarprod(B[0],mu[0])),B1mu1(Qp->colscalarprod(B[1],mu[1]));
              g -=  (B0mu0 + B1mu1);
            }
          }
        } else {
          if (Qnormalized) inc = solve1dprob(mu,true,tol,epsilon,alpha(B[0]),g(B[0]),u[B[0]]);
          else inc = solve1dprob(mu,false,tol,epsilon,alpha(B[0]),g(B[0]),u[B[0]],(*Qp)(B[0],B[0]));
//          if (inc==SEMIPOSQMAT) return false;
//          if (inc==NOINC) stop = true;
          if (inc==NOINC || inc==NONPOSQMAT || inc==SEMIPOSQMAT ) stop = true;
          else {
            fval += inc;
            alpha(B[0]) += mu[0];
            if (iter==maxiter) stop = true;
            else {
              const arma::vec B0mu0 = Qp->colscalarprod(B[0],mu[0]);
              g -=  B0mu0;
            }
          }
        }
        if (!stop) {
          int a(B[0]);
          double newalpha1(alpha(a)),u1(u[a]);
          if (newalpha1 <= tol) VatLB.push_back(a);
          else if (newalpha1 >= u1-tol) VatUB.push_back(a);
          else VatDomInt.push_back(a);
        }
      }

    }

    return true;
}

arma::ivec findWrkind(const arma::vec& g, const double epsilon, const double tol, const bool Qnormalized,
                             arma::vec& alpha, const arma::vec& u, SVMCMatrix* Qp,
                             std::list<int>& VatLB, std::list<int>& VatUB, std::list<int>& VatDomInt)
{
  int m(g.size());
  arma::ivec B = {-1,-1};
  double maxg(0.);
  int a,VatLele;
  bool aattUB(false),aatBnd(false);

  if (VatUB.size() > 0) {
    std::list<int>::iterator it=VatUB.begin();
    for (int i=0; i<VatUB.size(); it++, i++) {
      double gi = -g(VatLele=*it);
      if (gi > maxg) {
        a = VatLele;
        maxg = gi;
      }
    }
  }

  if (maxg > epsilon)  aattUB = aatBnd = true;

  if (VatLB.size() > 0) {
    std::list<int>::iterator it=VatLB.begin();
    for (int i=0; i<VatLB.size(); it++, i++) {
      double gi = g(VatLele=*it);
      if (gi > maxg) {
        a = VatLele;
        maxg = gi;
        if (aattUB) aattUB = false;
      }
    }
  }
  if (!aatBnd && maxg > epsilon) aatBnd = true;

  if (VatDomInt.size() > 0) {
    std::list<int>::iterator it=VatDomInt.begin();
    for (int i=0; i<VatDomInt.size(); it++, i++) {
      double gi = fabs(g(VatLele=*it));
      if (gi > maxg) {
        a = VatLele;
        maxg = gi;
        if (aattUB) aattUB = false;
        if (aatBnd) aatBnd = false;
      }
    }
  }

  if (maxg <= epsilon) return B;
  // assert(a >= 0 && a < m);

  if (aattUB) VatUB.remove(a);
  else if (aatBnd) VatLB.remove(a);
       else VatDomInt.remove(a);
  B[0] = a;

  double maxtwicegain(0.), ga(g(a)), alphaa(alpha(a)), ua(u[a]);
//  double Qaa, Qbb, Qab, detQB, alphab, gb, ub;
  double Qaa, Qbb, Qab, alphab, gb, ub;
  if (!Qnormalized) Qaa = (*Qp)(a,a);
  for (int b=0;b<m;b++) if (b!=a)
  {
    alphab = alpha(b);  gb = g(b); ub = u[b];
    Qab = (*Qp)(a,b);
/*
    if (Qnormalized)  detQB = 1. - Qab*Qab;
    else {
      Qbb = (*Qp)(b,b);
      detQB = Qaa*Qbb - Qab*Qab;
    }
*/
    if (!Qnormalized)  Qbb = (*Qp)(b,b);
    double mua;
    if (alphaa <= tol) mua = ua;
    else if (alphaa >= ua-tol) mua = -ua;
    else {
      if (Qnormalized) {
        if (ga > Qab*gb) mua = ua-alphaa;
        else mua = -alphaa;
      } else {
        if (Qbb*ga > Qab*gb) mua = ua-alphaa;
        else mua = -alphaa;
      }
    }

    double posdir2gain(0.),negdir2gain(0.);
    if (Qnormalized) {
      if (alphab <= ub-tol) posdir2gain = (2*(gb+ub)-Qab*mua)*ub;
      if (alphab >= tol) negdir2gain = (2*(ub-gb)+Qab*mua)*ub;
    } else {
      if (alphab <= ub-tol) posdir2gain = (2*(gb+Qbb*ub)-Qab*mua)*ub;
      if (alphab >= tol) negdir2gain = (2*(Qbb*ub-gb)+Qab*mua)*ub;
    }
    double twicegain = posdir2gain > negdir2gain ? posdir2gain : negdir2gain;
    if (twicegain > maxtwicegain) {
      maxtwicegain = twicegain;
      B[1] = b;
    }
  }

  if (maxtwicegain < epsilon) return B;
  int bestb = B[1];
  double alphabb = alpha(bestb);
  if (alphabb <= tol) VatLB.remove(bestb);
  else if (alphabb >= ub - tol) VatUB.remove(bestb);
       else VatDomInt.remove(bestb);

  return B;
}

void normtwocntviol(bool l2scvio, double& muB1, double& muB2,
                    const double alpha1, const double alpha2, const double u1, const double u2,
                    const double g1, const double g2, const double Q12)
{
  double f1muB1(muB1),f1muB2,f1,f2;
  f1muB2 = muB2 = normonedimsol(muB1,alpha2,u2,g2,Q12);
  f1 = normf2d(muB1,muB2,g1,g2,Q12);
  if (l2scvio) muB2 = -alpha2;
  else muB2 = u2-alpha2;
  muB1 = normonedimsol(muB2,alpha1,u1,g1,Q12);
  f2 = normf2d(muB1,muB2,g1,g2,Q12);
  if (f1>f2) {
    muB1 = f1muB1;
    muB2 = f1muB2;
  }
}

void unormtwocntviol(bool l2scvio, double& muB1, double& muB2,
                     const double alpha1, const double alpha2, const double u1, const double u2,
                     const double g1, const double g2, const double Q12, const double Q11, const double Q22)
{
  double f1muB1(muB1),f1muB2,f1,f2;
  f1muB2 = muB2 = unormonedimsol(muB1,alpha2,u2,g2,Q12,Q22);
  f1 = unormf2d(muB1,muB2,g1,g2,Q12,Q11,Q22);
  if (l2scvio) muB2 = -alpha2;
  else muB2 = u2-alpha2;
  muB1 = unormonedimsol(muB2,alpha1,u1,g1,Q12,Q11);
  f2 = unormf2d(muB1,muB2,g1,g2,Q12,Q11,Q22);
  if (f1>f2) {
    muB1 = f1muB1;
    muB2 = f1muB2;
  }
}

double solve2dprob(arma::vec& muB, const bool Qnormalized, const double tol, const double epsilon,
                   const double alpha1, const double alpha2, const double g1, const double g2,
                   const double Q12, const double u1, const double u2, const double Q11, const double Q22)
{
  double detQB,muB1,muB2,inc;
  if (Qnormalized) detQB = 1. - Q12*Q12;
  else detQB = Q11*Q22 - Q12*Q12;
  if (detQB < -tol) return NONPOSQMAT;
  if (detQB <= tol) return SEMIPOSQMAT;
  if (Qnormalized) {
    muB1 = (g1-Q12*g2) / detQB;
    muB2 = (g2-Q12*g1) / detQB;
  } else {
    muB1 = (Q22*g1-Q12*g2) / detQB;
    muB2 = (Q11*g2-Q12*g1) / detQB;
  }
  double newalpha1 = alpha1 + muB1;
  double newalpha2 = alpha2 + muB2;
  if (newalpha1 > tol) {
    if (newalpha1 < u1-tol) {
      if (newalpha2 > tol) {
        if (newalpha2 >= u2-tol) {
          muB2 = u2-alpha2;
          if (Qnormalized) muB1 = normonedimsol(muB2,alpha1,u1,g1,Q12);
          else muB1 = unormonedimsol(muB2,alpha1,u1,g1,Q12,Q11);
        }
      } else {
        muB2 = -alpha2;
        if (Qnormalized) muB1 = normonedimsol(muB2,alpha1,u1,g1,Q12);
        else muB1 = unormonedimsol(muB2,alpha1,u1,g1,Q12,Q11);
      }
    } else {
      muB1 = u1-alpha1;
      if (muB2 > tol) {
        if (newalpha2 < u2-tol)  {
          if (Qnormalized) muB2 = normonedimsol(muB1,alpha2,u2,g2,Q12);
          else muB2 = unormonedimsol(muB1,alpha2,u2,g2,Q12,Q22);
        } else {
          if (Qnormalized) normtwocntviol(false,muB1,muB2,alpha1,alpha2,u1,u2,g1,g2,Q12);
          else  unormtwocntviol(false,muB1,muB2,alpha1,alpha2,u1,u2,g1,g2,Q12,Q11,Q22);
        }
      } else {
        if (Qnormalized) normtwocntviol(true,muB1,muB2,alpha1,alpha2,u1,u2,g1,g2,Q12);
        else  unormtwocntviol(true,muB1,muB2,alpha1,alpha2,u1,u2,g1,g2,Q12,Q11,Q22);
      }
    }
  } else {
    muB1 = -alpha1;
    if (newalpha2 > tol) {
      if (newalpha2 < u2-tol)  {
        if (Qnormalized) muB2 = normonedimsol(muB1,alpha2,u2,g2,Q12);
        else muB2 = unormonedimsol(muB1,alpha2,u2,g2,Q12,Q22);
      } else {
        if (Qnormalized) normtwocntviol(false,muB1,muB2,alpha1,alpha2,u1,u2,g1,g2,Q12);
        else  unormtwocntviol(false,muB1,muB2,alpha1,alpha2,u1,u2,g1,g2,Q12,Q11,Q22);
      }
    } else {
      if (Qnormalized) normtwocntviol(true,muB1,muB2,alpha1,alpha2,u1,u2,g1,g2,Q12);
      else  unormtwocntviol(true,muB1,muB2,alpha1,alpha2,u1,u2,g1,g2,Q12,Q11,Q22);
    }
  }

  if (Qnormalized) inc = normf2d(muB1,muB2,g1,g2,Q12);
  else inc = unormf2d(muB1,muB2,g1,g2,Q12,Q11,Q22);

  if (inc < epsilon) return NOINC;
  muB[0] = muB1; muB[1] = muB2;

  return inc;
}

double solve1dprob(arma::vec& muB, const bool Qnormalized,  const double tol, const double epsilon,
                    const double alpha1, const double  g1, const double u1, const double Q11)
{
  double muB1,inc;

  if (Qnormalized) muB1 = g1;
  else muB1 = g1/Q11;
  double newalpha1 = alpha1 + muB1;
  if (newalpha1 < tol)  muB1 = -alpha1;
  else if (newalpha1 > u1-tol) muB1 = u1-alpha1;
  if (Qnormalized) inc = normf1d(muB1,g1);
  else inc = unormf1d(muB1,g1,Q11);
  if (inc < epsilon) return NOINC;
  muB[0] = muB1;
  return inc;
}

