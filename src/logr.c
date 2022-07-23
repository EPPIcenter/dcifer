#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>

/******************************************************************************/
/*                                                                            */
/*                         Calculate logs for reval                           */
/*        logr, log1r, number of 1's, number of non-01, sum(log1r)            */
/*                                                                            */ 
/******************************************************************************/

SEXP logReval(SEXP Rreval, SEXP Rneval, SEXP Rnm)
{
  int neval, nm, nn, i, j, mtemp, ikeep, shift;
  double rtemp, stemp;
  int *m1, *nmid;
  double *reval, *logr, *log1r, *sum1r;
  SEXP Rres, Rlogr, Rlog1r, Rm1, Rnmid, Rsum1r;

  reval = REAL(Rreval);
  neval = INTEGER(Rneval)[0];
  nm    = INTEGER(Rnm)[0];
  nn    = neval*nm;
  
  Rres   = PROTECT(allocVector(VECSXP,  5));
  Rlogr  = PROTECT(allocVector(REALSXP, nn));
  Rlog1r = PROTECT(allocVector(REALSXP, nn));  
  Rm1    = PROTECT(allocVector(INTSXP,  neval));
  Rnmid  = PROTECT(allocVector(INTSXP,  neval));
  Rsum1r = PROTECT(allocVector(REALSXP, neval));
  logr  = REAL(Rlogr);
  log1r = REAL(Rlog1r);
  m1    = INTEGER(Rm1);
  nmid  = INTEGER(Rnmid);
  sum1r = REAL(Rsum1r);

  shift = 0;
  for (j = 0; j < neval; j++) {
    mtemp = 0;
    stemp = 0;
    ikeep = 0;
    for (i = 0; i < nm; i++) {
      rtemp = reval[shift + i];
      if (rtemp == 1) {
	mtemp++;
      } else if (rtemp != 0) { 
	logr[ shift + ikeep] = log(rtemp);  
	log1r[shift + ikeep] = log(1 - rtemp); 
	stemp       += log1r[shift + ikeep];  
	ikeep++;
      }
    }
    m1[   j] = mtemp;
    nmid[ j] = ikeep;
    sum1r[j] = stemp;
    shift += nm;
  }

  SET_VECTOR_ELT(Rres, 0, Rlogr);
  SET_VECTOR_ELT(Rres, 1, Rlog1r);
  SET_VECTOR_ELT(Rres, 2, Rm1);
  SET_VECTOR_ELT(Rres, 3, Rnmid);
  SET_VECTOR_ELT(Rres, 4, Rsum1r);  
  UNPROTECT(6);
  return Rres;
}
    

  
  
