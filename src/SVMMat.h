#ifndef SVMMAT_H
#define SVMMAT_H

// #include <vector>
#include <RcppArmadillo.h>
using namespace Rcpp;

class SVMCMatrix {
  public:
    virtual unsigned dim(void) = 0;
    virtual unsigned n(void) = 0;
    virtual unsigned k(void) = 0;
    virtual void fill(double scl) = 0;
    virtual void setgrp(const arma::uvec* igrps) = 0;
    virtual void resize(const arma::uvec* igrps, const unsigned m, const unsigned k) = 0;
    virtual double& operator()(unsigned row, unsigned col) = 0;
//    virtual const double operator()(unsigned row, unsigned col) = 0;
    virtual arma::vec colscalarprod(unsigned col, double scl) = 0;
    virtual void cov2cor(arma::vec& dQ) = 0;
    void cov2cor(void) {
      unsigned m = this->dim();
      arma::vec dQ(m);
      for (unsigned i=0;i<m;i++) dQ(i) = sqrt((*this)(i,i));
      this->cov2cor(dQ);
    }
    virtual void cor2cov(arma::vec& dQ) = 0;
};

class ARMASVMMatrix : public SVMCMatrix
{
  private:
     arma::mat* data_;
     unsigned n_;
     unsigned k_;
     unsigned m_;

  public:
    ARMASVMMatrix(arma::mat& idata, unsigned m) : data_(&idata), n_(0), k_(0), m_(m) { }
    ARMASVMMatrix(arma::mat& idata, unsigned n, unsigned k) : data_(&idata), n_(n), k_(k), m_(n*(k-1))  
      { if ( idata.n_rows!=m_ || idata.n_cols!=m_ ) Rf_error("Wrong dimensions in ARMASVMMatrix construction\n");  }
    virtual unsigned dim(void) { return m_; }
    virtual unsigned n(void) { return n_; }
    virtual unsigned k(void) { return k_; }
    arma::mat* data(void) { return data_; }
    virtual void resize(const arma::uvec* igrps, const unsigned m, const unsigned k) { m_ = m; k_ = k;  n_ = m/(k-1); data_->resize(m_,m_); }
    void resize(const unsigned m) { m_ = m; k_ = 0;  n_ = 0; ; data_->resize(m_,m_); }
    virtual void setgrp(const arma::uvec* igrps) {  }
    virtual double& operator()(unsigned row, unsigned col) { return (*data_)(row,col); }
//    virtual const double operator()(unsigned row, unsigned col) { return (*data_)(row,col); }
    virtual void fill(double scl) { data_->fill(scl); }
    virtual arma::vec colscalarprod(unsigned col, double scl) { return scl * data_->col(col); }
    virtual void cov2cor(arma::vec& dQ) {
      double dQidQj;
      for (unsigned i=0;i<m_;i++) {
        if (dQ(i)!=1.) (*data_)(i,i) = 1.;
        for (unsigned j=i+1;j<m_;j++) {
          if (dQ(i)!=1.) dQidQj = dQ(i)*dQ(j);
          else dQidQj = dQ(j);
          if (dQidQj!=1) (*data_)(i,j) = (*data_)(j,i) = (*data_)(i,j)/dQidQj;
        }
      }
    }
    virtual void cor2cov(arma::vec& dQ) {
      double dQi,dQj,dQidQj;
      for (unsigned i=0;i<m_;i++) {
        dQi = dQ(i);
        if (dQi!=1.) (*data_)(i,i) = dQi*dQi;
        for (unsigned j=i+1;j<m_;j++) {
          dQj = dQ(j);
          if (dQi!=1.) dQidQj = dQi*dQj;
          else dQidQj = dQj;
          if (dQidQj!=1) (*data_)(i,j) = (*data_)(j,i) = (*data_)(i,j)*dQidQj;
        }
      }
    }

};

class LLWCMatrix : public SVMCMatrix
{
  private:
     std::vector<arma::vec> data_;
     const arma::uvec* grps_;
     unsigned n_;
     unsigned k_;
     unsigned m_;
  public:
    LLWCMatrix(std::vector<arma::vec> idata, arma::uvec* igrps, unsigned n, unsigned k) : data_(idata), grps_(igrps), n_(n), k_(k), m_(n*(k-1)) { }
    LLWCMatrix(const arma::uvec* igrps, unsigned n, unsigned k) :  grps_(igrps), n_(n), k_(k), m_(n*(k-1)) { data_.resize(2); unsigned dtdim(n*(n+1)/2); data_[0].resize(dtdim); data_[1].resize(dtdim); }
    virtual void resize(const arma::uvec* igrps, const unsigned m, const unsigned k) { grps_ = igrps,  m_ = m; k_ = k; n_ = m/(k-1); unsigned dtdim(n_*(n_+1)/2); data_[0].resize(dtdim); data_[1].resize(dtdim); }
    virtual unsigned dim(void) { return m_; }
    virtual unsigned n(void) { return n_; }
    virtual unsigned k(void) { return k_; }
    const arma::uvec* grouping(void) { return grps_; }
    const std::vector<arma::vec> data(void) { return data_; }
    virtual void setgrp(const arma::uvec* igrps) { grps_ = igrps; }
    virtual double& operator()(unsigned row, unsigned col) { unsigned i,j; getind(row,col,i,j); return data_[i](j); }
    virtual void fill(double scl) { std::fill(data_[0].begin(),data_[0].end(),scl); std::fill(data_[1].begin(),data_[1].end(),scl);  }
    virtual arma::vec colscalarprod(unsigned col, double scl) {
      arma::vec res(m_);
      unsigned dataind, obsc(col/(k_-1)), grpc(col-obsc*(k_-1)), yj((*grps_)(obsc));
      if (grpc>=yj) grpc++;
      for (unsigned obsr=0;obsr<n_;obsr++) {
        unsigned j0,j1;
        if (obsr<=obsc) { j0 = obsr; j1 = obsc; }
        else { j0 = obsc ; j1 = obsr; }
        if (j0==0) dataind = j1;
        else if (j0==1) dataind = n_ + j1-j0 ;
        else dataind = j0*n_ - j0*(j0-1)/2 + j1-j0;
        double scres0(scl*data_[0](dataind)), scres1(scl*data_[1](dataind));
        unsigned baseind(obsr*(k_-1)), yi((*grps_)(obsr));
        for (unsigned gindr=0;gindr<k_-1;gindr++) {
          unsigned grpr(gindr);
          if (gindr>=yi) grpr++;
          if (grpr!=grpc) res(baseind+gindr) = scres0;
          else res(baseind+gindr) = scres1;
        }
      }
      return res;
    }
    virtual void getind(unsigned row, unsigned col, unsigned& i, unsigned& j) {
      unsigned obsr(row/(k_-1)), grpr(row-obsr*(k_-1)), yi((*grps_)(obsr));
      if (grpr>=yi) grpr++;
      unsigned obsc(col/(k_-1)), grpc(col-obsc*(k_-1)), yj((*grps_)(obsc));
      if (grpc>=yj) grpc++;
      if (grpr!=grpc) i=0;  else  i=1;
      unsigned j0,j1;
      if (obsr<=obsc) { j0 = obsr; j1 = obsc; }
      else { j0 = obsc ; j1 = obsr; }
      if (j0==0) j = j1;
      else if (j0==1) j = n_ + j1-j0 ;
      else j = j0*n_ - j0*(j0-1)/2 + j1-j0;
      return;
    }
    void setdatael(unsigned i, unsigned j, double val) { data_[i](j) = val; }
    virtual void cov2cor(arma::vec& dQ) {
      unsigned datai,dataij;
      double dQi,dQj,dQidQj;
      for (unsigned i=0;i<n_;i++) {
        dQi = dQ(i*(k_-1));
        if (i==0) datai = 0;
        else if (i==1) datai = n_;
             else datai = i*n_ - i*(i-1)/2;
        for (unsigned j=i;j<n_;j++) {
          dQj = dQ(j*(k_-1));
          dQidQj = dQi*dQj;
          dataij = datai + j-i;
          data_[0](dataij) /= dQidQj;
          data_[1](dataij) /= dQidQj;
        }
      }
    }
    virtual void cor2cov(arma::vec& dQ) {
      unsigned datai,dataij;
      double dQi,dQj,dQidQj;
      for (unsigned i=0;i<n_;i++) {
        dQi = dQ(i*(k_-1));
        if (i==0) datai = 0;
        else if (i==1) datai = n_;
             else datai = i*n_ - i*(i-1)/2;
        for (unsigned j=i;j<n_;j++) {
          dQj = dQ(j*(k_-1));
          dQidQj = dQi*dQj;
          dataij = datai + j-i;
          data_[0](dataij) *= dQidQj;
          data_[1](dataij) *= dQidQj;
        }
      }

    }


};

#endif
