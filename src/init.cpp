#include "makeKMat.h"
#include "svmModels.h"
#include "TrainCanSVMs.h"
#include "TrainSVMs.h"
#include "TuneSVMs.h"

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

static const R_CallMethodDef CallEntries[] = {
  {"makeGauKMat",   (DL_FUNC) &makeGauKMat,   2},
  {"makeGauKMat2", (DL_FUNC) &makeGauKMat2,    3},
  {"TrainCanSVM", (DL_FUNC) &TrainCanSVM,    10},
  {"TrSVMs", (DL_FUNC) &TrSVMs,    11},
  {"TunSVMs", (DL_FUNC) &TunSVMs,    14},
  {"makeLinKMat",   (DL_FUNC) &makeLinKMat,   1},
  {"makeLinKMat2",   (DL_FUNC) &makeLinKMat2,  2},
  {"makePolKMat",   (DL_FUNC) &makePolKMat,   4},
  {"makePolKMat2",   (DL_FUNC) &makePolKMat2,  5},
  {NULL, NULL, 0}
};

void R_init_ProbSVMs(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
