#ifndef Predict_H
#define Predict_H

#include <limits>
#include "consvmalphatoMat.h"
#include "Losses.h"

const double infinity = std::numeric_limits<double>::max();

template<class ClassifierType, class Ktype>
void predict(const ClassifierType& svm, const Ktype& K, arma::uvec& classpred);

template<class Ktype>
void predict(const arma::vec& svmalpha0, const arma::uvec& Y, const Ktype& K, arma::mat& svmalpha, arma::uvec& classpred,const int model)
{
  if (model==LLW)  LLWconsvmalphatoMat(svmalpha0,Y,svmalpha);
  else if (model==WW)  WWconsvmalphatoMat(svmalpha0,Y,svmalpha);
  predict<arma::mat,Ktype>(svmalpha,K,classpred);
}

double CLhitrate(const arma::uvec& orgclass, arma::uvec& predclass);
double CLhitrate(const unsigned orgclass, arma::uvec& predclass);

#endif
