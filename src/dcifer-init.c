#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP llik(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP llikEqr(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"llik",    (DL_FUNC) &llik,    11},
    {"llikEqr", (DL_FUNC) &llikEqr, 11},
    {NULL, NULL, 0}
};

void R_init_unfazeD(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
