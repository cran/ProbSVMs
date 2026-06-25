#include "svmModels.h"


void LLWQ_C(const arma::uvec& grouping, const double lambda, const arma::mat& K, const arma::vec& w, arma::vec& v, ARMASVMMatrix& Q, arma::vec& u)
{
  int n=grouping.size(),k=w.size(),m=n*(k-1);
  double lambdainvover2=1./(2*lambda);

  double vcoef = 1./(k-1);
  for (arma::vec::iterator pos=v.begin();pos<v.end();pos++) *pos = vcoef;
  Q.fill(-lambdainvover2/k);
  for (int i=0;i<n;i++) {
    int yi = grouping(i);
    int rowind0 = i*(k-1);
    for (int j=0;j<=i;j++) {
      int yj = grouping(j);
      int colind0 = j*(k-1);
        for (int c=0;c<k;c++) if(c!=yi && c!=yj) {
          int hi,ri,hj,cj;
          if (c<yi) hi=c;  else hi = c-1;
          ri = rowind0 + hi;
          if (c<yj) hj=c;  else hj = c-1;
          cj = colind0 + hj;
          Q(ri,cj) += lambdainvover2;
        }
        for (int c1=0;c1<k-1;c1++) {
          if (j<i) for (int c2=0;c2<k-1;c2++) Q(rowind0+c1,colind0+c2) *= K(i,j);
          else for (int c2=0;c2<=c1;c2++) Q(rowind0+c1,colind0+c2) *= K(i,j);
        }
    }
  }
  for (int row=0;row<m;row++) for (int col=0;col<row;col++) Q(col,row) = Q(row,col);

  arma::vec normw(k);
  for (int c=0;c<k;c++) normw(c) = w(c) / n;
  arma::vec::iterator pos=u.begin();
  for (int i=0;i<n;i++) {
    int y = grouping(i);
    for (int c=0;c<k;c++) if (c!=y) *pos++ = normw(y);
  }

  return;
}

void LLWQ_C(const arma::uvec& grouping, const double lambda, const NumericVector& K, const arma::vec& w, arma::vec& v, LLWCMatrix& Q, arma::vec& u)
{
  int n=grouping.size(),k=w.size();
  double lambdainvover2=1./(2*lambda),lambdainvover2k=lambdainvover2/k;

  double vcoef = 1./(k-1);
  for (arma::vec::iterator pos=v.begin();pos<v.end();pos++) *pos = vcoef;

  if (Q.grouping()==0) Q.setgrp(&grouping);
  for (int datai=0;datai<K.size();datai++) {
    double Kij = K(datai);
    double const1 = -Kij*lambdainvover2k;
    Q.setdatael(0,datai,const1);
    Q.setdatael(1,datai,const1+Kij*lambdainvover2);
  }
  arma::vec normw = w / n;
  arma::vec::iterator pos=u.begin();
  for (int i=0;i<n;i++) {
    int y = grouping(i);
    for (int c=0;c<k;c++) if (c!=y) *pos++ = normw(y);
  }

  return;
}

void LLWQ_C(const arma::uvec& grouping, const double lambda, const arma::mat& K, const arma::vec& w, arma::vec& v, LLWCMatrix& Q, arma::vec& u)
{
  int n=grouping.size(),k=w.size();
  double lambdainvover2=1./(2*lambda),lambdainvover2k=lambdainvover2/k;

  double vcoef = 1./(k-1);
  for (arma::vec::iterator pos=v.begin();pos<v.end();pos++) *pos = vcoef;

  if (Q.grouping()==0) Q.setgrp(&grouping);
  for (int datai=0,col=0;col<n;col++) for (int row=col;row<n;row++,datai++)
  {
    double Kij = K(row,col);
    double const1 = -Kij*lambdainvover2k;
    Q.setdatael(0,datai,const1);
    Q.setdatael(1,datai,const1+Kij*lambdainvover2);
  }
  arma::vec normw = w / n;
  arma::vec::iterator pos=u.begin();
  for (int i=0;i<n;i++) {
    int y = grouping(i);
    for (int c=0;c<k;c++) if (c!=y) *pos++ = normw(y);
  }

  return;
}

void WWQ_C(const arma::uvec& grouping, const double lambda, const arma::mat& K, const arma::vec& w, arma::vec& v, ARMASVMMatrix& Q, arma::vec& u)
{
  int n=grouping.size(),k=w.size(),m=n*(k-1);
  double lambdainvover8=1./(8*lambda);

  for (int i=0;i<m;i++) v(i) = 1.;
  Q.fill(0.);
  for (int i=0;i<n;i++) {
    int yi = grouping(i);
    int rowind0 = i*(k-1);
    for (int j=0;j<=i;j++) {
      int yj = grouping(j);
      int colind0 = j*(k-1);
      for (int c=0;c<k;c++) if(c!=yi && c!=yj) {
        int hi,ri,hj,cj;
        if (c<yi) hi=c;  else hi = c-1;
        ri = rowind0 + hi;
        if (c<yj) hj=c;  else hj = c-1;
        cj = colind0 + hj;
        Q(ri,cj) += lambdainvover8;
      }
      if (yi==yj) {
        for (int ri=rowind0;ri<rowind0+k-1;ri++) for (int cj=colind0;cj<colind0+k-1;cj++)
          Q(ri,cj) += lambdainvover8;
      } else {
        int hyij,cyij;
        if (yi<yj) hyij=yi;  else hyij = yi-1;
        cyij= colind0 + hyij;
        for (int ri=rowind0;ri<rowind0+k-1;ri++) Q(ri,cyij) -= lambdainvover8;
        int hyji,ryji;
        if (yj<yi) hyji=yj;  else hyji = yj-1;
        ryji= rowind0 + hyji;
        for (int cj=colind0;cj<colind0+k-1;cj++) Q(ryji,cj) -= lambdainvover8;
      }
      for (int c1=0;c1<k-1;c1++) {
        if (j<i) for (int c2=0;c2<k-1;c2++) Q(rowind0+c1,colind0+c2) *= K(i,j);
        else for (int c2=0;c2<=c1;c2++) Q(rowind0+c1,colind0+c2) *= K(i,j);
      }

    }
  }
  for (int row=0;row<m;row++) for (int col=0;col<row;col++) Q(col,row) = Q(row,col);

  arma::vec normw(k);
  for (int c=0;c<k;c++) normw(c) = w(c) / n;
  arma::vec::iterator pos=u.begin();
  for (int i=0;i<n;i++) {
    int y = grouping(i);
    for (int c=0;c<k;c++) if (c!=y) *pos++ = normw(y);
  }

  return;
}

void WWQ_C(const arma::uvec& grouping, const double lambda, const NumericVector& K, const arma::vec& w, arma::vec& v, ARMASVMMatrix& Q, arma::vec& u)
{
//  int n=grouping.size(),k=w.size(),m=n*(k-1);
  int n=grouping.size(),k=w.size();

  arma::mat K1(n,n);
  for (int Kpos=0,row=0;row<n;row++) for (int col=row;col<n;col++,Kpos++)
  {
    if (row==col) K1(row,col) = K(Kpos);
    else K1(row,col) = K1(col,row) = K(Kpos);
  }
  WWQ_C(grouping,lambda,K1,w,v,Q,u);
}

void WWQ_C(const arma::uvec& grouping, const double lambda, const NumericVector& K, const arma::vec& w, arma::vec& v, LLWCMatrix& Q, arma::vec& u)
{
  //  assert(0==1);
  Rprintf("Error:  This is a shallow function only used for compiling compatibility that should never be called\n");
  return;
}


void WWQ_C(const arma::uvec& grouping, const double lambda, const arma::mat& K, const arma::vec& w, arma::vec& v, LLWCMatrix& Q, arma::vec& u)
{
  //  assert(0==1);
  Rprintf("Error:  This is a shallow function only used for compiling compatibility that should never be called\n");
  return;
}

