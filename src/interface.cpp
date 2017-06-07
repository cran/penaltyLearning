/* -*- compile-command: "R CMD INSTALL .." -*- */

#include "modelSelection.h"
#include "largestContinuousMinimum.h"
#include <R.h>
#include <R_ext/Rdynload.h>

void largestContinuousMinimum_interface
(int *n_data, double *cost_vec, double *size_vec, int *index_vec){
  int status = largestContinuousMinimum
    (*n_data, cost_vec, size_vec, index_vec);
  if(status==ERROR_SIZES_MUST_BE_POSITIVE){
    error("sizes must be positive");
  }
  if(status != 0){
    error("error code %d", status);
  }
}

void modelSelection_interface
(double *loss, double *complexity, int *n_models,
 int *before, double *lambda
 ){
  int status = modelSelection(loss, complexity, *n_models, before, lambda);
  if(status == ERROR_LOSS_NOT_DECREASING){
    error("loss not decreasing");
  }
  if(status == ERROR_COMPLEXITY_NOT_INCREASING){
    error("complexity not increasing");
  }
  if(status != 0){
    error("error code %d", status);
  }
}
  
R_CMethodDef cMethods[] = {
  {"modelSelection_interface",
   (DL_FUNC) &modelSelection_interface, 5
   //,{REALSXP, REALSXP, INTSXP, INTSXP, REALSXP}
  },
  {"largestContinuousMinimum_interface",
   (DL_FUNC) &largestContinuousMinimum_interface, 4
   //,{INTSXP, REALSXP, REALSXP, INTSXP}
  },
  {NULL, NULL, 0}
};

extern "C" {
  void R_init_penaltyLearning(DllInfo *info) {
    R_registerRoutines(info, cMethods, NULL, NULL, NULL);
    //R_useDynamicSymbols call says the DLL is not to be searched for
    //entry points specified by character strings so .C etc calls will
    //only find registered symbols.
    R_useDynamicSymbols(info, FALSE);
  }
}
