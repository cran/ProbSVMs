#ifndef makeKMat_H
#define makeKMat_H

#include <RcppArmadillo.h>
using namespace Rcpp;

RcppExport
SEXP makeGauKMat(const SEXP traindat_s, const SEXP sigma_s);
RcppExport
SEXP makeGauKMat2(const SEXP traindat_s, const SEXP newdat_s, const SEXP sigma_s);
RcppExport
SEXP makeLinKMat(const SEXP traindat_s);
RcppExport
SEXP makeLinKMat2(const SEXP traindat_s, const SEXP newdat_s);
RcppExport
SEXP makePolKMat(const SEXP traindat_s, SEXP degree_s, SEXP scale_s, SEXP offset_s);
RcppExport
SEXP makePolKMat2(const SEXP traindat_s, const SEXP newdat_s, SEXP degree_s, SEXP scale_s, SEXP offset_s);

#endif
