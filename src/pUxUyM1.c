#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>
#include "llik.h"

/******************************************************************************/
/*                                                                            */ 
/*                                   M = 1                                   */
/*                      likelihood for a pair of samples                      */
/*   fix probs for X; subsets of Sxy: subtract from Y update multinom coefs   */
/*                                                                            */ 
/******************************************************************************/

/* ixy - indices of which ux are in uy (or uxy), sorted; likewise for iyx */
/* logj(starts with 1), factj (starts with 0)                             */
SEXP llikM1(SEXP Rux, SEXP Ruy, SEXP Rixy, SEXP Riyx, SEXP Rnx, SEXP Rny,
	    SEXP Rlogp, SEXP Rlogj, SEXP Rfactj, SEXP Rr, SEXP Rneval)
{
  int nux, nuy, nuxy, nx, ny, ax, ay, bx, by, i, ix, iy, mmax, neval;
  int *ux1, *uy1, *ixy1, *iyx1;  // to not rewrite R obj  
  double combx, comby, sump0, sump1, dif10;
  double *logp, *logj, *factj, *r, *lik;
  SEXP Rlik; 

  ux1   = INTEGER(Rux);   // one-based indices (passed from R)
  uy1   = INTEGER(Ruy);   // one-based indices (passed from R)
  ixy1  = INTEGER(Rixy);  // one-based indices (passed from R)
  iyx1  = INTEGER(Riyx);  // one-based indices (passed from R)
  nx    = INTEGER(Rnx)[0];
  ny    = INTEGER(Rny)[0];
  neval = INTEGER(Rneval)[0];
  logp  = REAL(Rlogp);
  logj  = REAL(Rlogj);
  factj = REAL(Rfactj);
  r     = REAL(Rr);  
  nux   = length(Rux);         
  nuy   = length(Ruy);
  nuxy  = length(Rixy);   // length(Rixy) = length(Riyx)

  Rlik  = PROTECT(allocVector(REALSXP, neval));
  lik   = REAL(Rlik);
  
  for (i = 0; i < neval; i++) {
    lik[i] = 0;
  }

  ax = nux - 1;  /* can     be 0 */
  ay = nuy - 1;
  bx = nx - ax;  /* can not be 0 */
  by = ny - ay;

  int ux[nux], uy[nuy], ixy[nuxy], iyx[nuxy], vx[nux], vy[nuy],
    vmaxx[nux], vmaxy[nuy];  
  double logpx[nux], logpy[nuy], logpxy[nuxy], ppx[ax], ppy[ay], pplastx[bx],
    pplasty[by], sprob[2];

  sump0 = 0;
  sump1 = 0;
  sprob[1] = 0;  /* in case nuxy = 0 */

  /* ux, uy, ixy, and iyx: zero-based indices for convenience (optional) */
  for (i = 0; i < nux;  i++) ux[i] = ux1[i] - 1;
  for (i = 0; i < nuy;  i++) uy[i] = uy1[i] - 1;
  for (i = 0; i < nuxy; i++) {
    ixy[i] = ixy1[i] - 1;
    iyx[i] = iyx1[i] - 1;
  }
  
  /* logpx, logpy, logpxy */
  for (i = 0; i < nux;  i++) logpx[ i] = logp[ux[i]];
  for (i = 0; i < nuy;  i++) logpy[ i] = logp[uy[i]];
  for (i = 0; i < nuxy; i++) logpxy[i] = logp[ux[ixy[i]]];

  /* vmaxx, vmaxy */
  vmaxx[0] = bx;
  vmaxy[0] = by;
  for (i = 1; i < nux; i++) vmaxx[i] = 1;
  for (i = 1; i < nuy; i++) vmaxy[i] = 1;  

  /* pplastx, pplasty - partial probs for last category */
  pplastx[0] = logpx[ax];
  for (i = 1; i < bx; i++) {
    pplastx[i] = pplastx[i - 1] + pplastx[0] - logj[i];
  }
  pplasty[0] = logpy[ay];
  for (i = 1; i < by; i++) {
    pplasty[i] = pplasty[i - 1] + pplasty[0] - logj[i];
  }

  /* initialize vx, nxlast, ppx with nonexistent "pre-first" combination */
  for (i = 0; i < ax; i++) {
    vx[i]  = 1;
    ppx[i] = logpx[i];
  }
  if (ax > 0) {
    vx[ ax - 1] = 0;
    ppx[ax - 1] = 0;
  }
  vx[ax] = bx + 1;

  /* subsequent combinations */ 
  ix = ax - 1;
  while (!equalArr(vx, vmaxx, nux)) {       
    if (ax == 0) {  
      vx[ax]--;     
    } else if (vx[ax] == 1 || vx[ix] == bx) {
      vx[ax] += vx[ix] - 1;
      vx[ix] = 1;
      ppx[ix] = logpx[ix];
      ix--;
      continue;
    } else {
      vx[ix]++;
      vx[ax]--;
      ppx[ix] += logpx[ix] - logj[vx[ix] - 1];
      ix = ax - 1;
    }
    combx = pplastx[vx[ax] - 1]; 
    for (i = 0; i < ax; i++) {
      combx += ppx[i];
    }
    combx += factj[nx];

    /* initialize vy, nylast, ppy */
    for (i = 0; i < ay; i++) {
      vy[i]  = 1;
      ppy[i] = logpy[i];
    }
    if (ay > 0) {
      vy[ ay - 1] = 0;
      ppy[ay - 1] = 0;
    }
    vy[ay] = by + 1;

    /* subsequent combinations for y */
    iy = ay - 1;
    while (!equalArr(vy, vmaxy, nuy)) {         
      if (ay == 0) {
	vy[ay]--;
      } else if (vy[ay] == 1 || vy[iy] == by) {
	vy[ay] += vy[iy] - 1; 
	vy[iy] = 1;
        ppy[iy] = logpy[iy];
	iy--;
	continue;
      } else {
	vy[iy]++;
	vy[ay]--;
        ppy[iy] += logpy[iy] - logj[vy[iy] - 1];
	iy = ay - 1;
      } 
      comby = pplasty[vy[ay] - 1];
      for (i = 0; i < ay; i++) {
        comby += ppy[i];
      }
      comby += factj[ny]; 

      /* calculate conditional probabilities P(Sx, Sy | m) */
      probSxSyCond(vx, vy, logpxy, logj, factj, nx, ny, nux, nuy, nuxy, ixy, iyx,
		   combx, comby, sprob, &mmax, 1);
      sump0 += sprob[0];
      sump1 += sprob[1];
    }
  }

  for (i = 0; i < neval; i++) {
    dif10 = sump1 - sump0;
    lik[i] = log(r[i]*dif10 + sump0); 
  }

  UNPROTECT(1);
  return Rlik;
}

/******************************************************************************/
/*                                                                            */ 
/*                                   M = 1                                    */
/*           sum(P(Sx, Sy | IBD = 0)) and sum(P(Sx, Sy | IBD = 1))            */
/*                                                                            */ 
/******************************************************************************/

/* return P1, P0 (can be used for llik and derivatives)                       */
SEXP p0p1(SEXP Rux, SEXP Ruy, SEXP Rixy, SEXP Riyx, SEXP Rnx, SEXP Rny,
	    SEXP Rlogp, SEXP Rlogj, SEXP Rfactj)
{
  int nux, nuy, nuxy, nx, ny, ax, ay, bx, by, i, ix, iy, mmax;
  int *ux1, *uy1, *ixy1, *iyx1;  // to not rewrite R obj  
  double combx, comby, sump0, sump1;
  double *logp, *logj, *factj, *res;
  SEXP Rres;

  ux1   = INTEGER(Rux);   // one-based indices (passed from R)
  uy1   = INTEGER(Ruy);   // one-based indices (passed from R)
  ixy1  = INTEGER(Rixy);  // one-based indices (passed from R)
  iyx1  = INTEGER(Riyx);  // one-based indices (passed from R)
  nx    = INTEGER(Rnx)[0];
  ny    = INTEGER(Rny)[0];
  logp  = REAL(Rlogp);
  logj  = REAL(Rlogj);
  factj = REAL(Rfactj);
  nux   = length(Rux);         
  nuy   = length(Ruy);
  nuxy  = length(Rixy);   // length(Rixy) = length(Riyx)

  Rres  = PROTECT(allocVector(REALSXP, 2));
  res   = REAL(Rres);
  
  ax = nux - 1;  /* can     be 0 */
  ay = nuy - 1;
  bx = nx - ax;  /* can not be 0 */
  by = ny - ay;

  int ux[nux], uy[nuy], ixy[nuxy], iyx[nuxy], vx[nux], vy[nuy],
    vmaxx[nux], vmaxy[nuy];  
  double logpx[nux], logpy[nuy], logpxy[nuxy], ppx[ax], ppy[ay], pplastx[bx],
    pplasty[by], sprob[2];
  
  sump0 = 0;
  sump1 = 0;  
  sprob[1] = 0;  /* in case nuxy = 0 */

  /* ux, uy, ixy, and iyx: zero-based indices for convenience (optional) */
  for (i = 0; i < nux;  i++) ux[i] = ux1[i] - 1;
  for (i = 0; i < nuy;  i++) uy[i] = uy1[i] - 1;
  for (i = 0; i < nuxy; i++) {
    ixy[i] = ixy1[i] - 1;
    iyx[i] = iyx1[i] - 1;
  }
  
  /* logpx, logpy, logpxy */
  for (i = 0; i < nux;  i++) logpx[ i] = logp[ux[i]];
  for (i = 0; i < nuy;  i++) logpy[ i] = logp[uy[i]];
  for (i = 0; i < nuxy; i++) logpxy[i] = logp[ux[ixy[i]]];

  /* vmaxx, vmaxy */
  vmaxx[0] = bx;
  vmaxy[0] = by;
  for (i = 1; i < nux; i++) vmaxx[i] = 1;
  for (i = 1; i < nuy; i++) vmaxy[i] = 1;  

  /* pplastx, pplasty - partial probs for last category */
  pplastx[0] = logpx[ax];
  for (i = 1; i < bx; i++) {
    pplastx[i] = pplastx[i - 1] + pplastx[0] - logj[i];
  }
  pplasty[0] = logpy[ay];
  for (i = 1; i < by; i++) {
    pplasty[i] = pplasty[i - 1] + pplasty[0] - logj[i];
  }

  /* initialize vx, nxlast, ppx with nonexistent "pre-first" combination */
  for (i = 0; i < ax; i++) {
    vx[i]  = 1;
    ppx[i] = logpx[i];
  }
  if (ax > 0) {
    vx[ ax - 1] = 0;
    ppx[ax - 1] = 0;
  }
  vx[ax] = bx + 1;

  /* subsequent combinations */ 
  ix = ax - 1;
  while (!equalArr(vx, vmaxx, nux)) {       
    if (ax == 0) {  
      vx[ax]--;     
    } else if (vx[ax] == 1 || vx[ix] == bx) {
      vx[ax] += vx[ix] - 1;
      vx[ix] = 1;
      ppx[ix] = logpx[ix];
      ix--;
      continue;
    } else {
      vx[ix]++;
      vx[ax]--;
      ppx[ix] += logpx[ix] - logj[vx[ix] - 1];
      ix = ax - 1;
    }
    combx = pplastx[vx[ax] - 1]; 
    for (i = 0; i < ax; i++) {
      combx += ppx[i];
    }
    combx += factj[nx];

    /* initialize vy, nylast, ppy */
    for (i = 0; i < ay; i++) {
      vy[i]  = 1;
      ppy[i] = logpy[i];
    }
    if (ay > 0) {
      vy[ ay - 1] = 0;
      ppy[ay - 1] = 0;
    }
    vy[ay] = by + 1;

    /* subsequent combinations for y */
    iy = ay - 1;
    while (!equalArr(vy, vmaxy, nuy)) {         
      if (ay == 0) {
	vy[ay]--;
      } else if (vy[ay] == 1 || vy[iy] == by) {
	vy[ay] += vy[iy] - 1; 
	vy[iy] = 1;
        ppy[iy] = logpy[iy];
	iy--;
	continue;
      } else {
	vy[iy]++;
	vy[ay]--;
        ppy[iy] += logpy[iy] - logj[vy[iy] - 1];
	iy = ay - 1;
      } 
      comby = pplasty[vy[ay] - 1];
      for (i = 0; i < ay; i++) {
        comby += ppy[i];
      }
      comby += factj[ny]; 

      /* calculate conditional probabilities P(Sx, Sy | m) */
      probSxSyCond(vx, vy, logpxy, logj, factj, nx, ny, nux, nuy, nuxy, ixy, iyx,
		   combx, comby, sprob, &mmax, 1);
      sump0 += sprob[0]; 
      sump1 += sprob[1];
    }
  }

  res[0] = sump0;
  res[1] = sump1;
  UNPROTECT(1);
  return Rres;
}
