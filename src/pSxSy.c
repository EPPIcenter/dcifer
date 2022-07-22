#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>
#include "llik.h"

/******************************************************************************/
/*                                                                            */
/*                                 P(Sx, Sy)                                  */
/*                           MIRSA implementation                             */
/*                                                                            */ 
/******************************************************************************/

/* ixy: 0-based sorted indices - which ux are in uy (or uxy); likewise for iyx */
/* mmax = min(|Sxy|, nm) [mmax <- min(sum(vxy), nm)], r[i] = 0 for i > nm     */
/* sprob is on log scale, so only non-zero ones, can't go beyond mmax */
// prob is sum1r

double probSxSy(double *logr, double *log1r, double prob, double *sprob,
		int nm, int mmax, int nmid, int m1) 
{
  if (mmax < m1) {  // number of r = 1 is greater than number of shared alleles
    return 0;
  }
  if (nmid == 0) {
    return exp(prob + sprob[m1]);
  }

  int i;
  double sum;
  int ibd[nmid], ibdmax[nmid];

  for (i = 0; i < nmid; i++) {
    ibd[i] = 0;
  }
  /* MIRSA */
  for (i = 0; i < nmid; i++) {
    ibdmax[i] = 1;
  }
  for (i = mmax - m1; i < nmid; i++) {  /* if mmax < nmid */ 
    ibdmax[i] = 0;
  }

  sum = exp(prob + sprob[m1]);  
  i = nmid - 1; 
  while (!equalArr(ibd, ibdmax, nmid)) {  // equalArr() returns 1 if nmid = 0
    if (ibd[i] == 1 || m1 == mmax) {
      if (ibd[i] == 0) {                /* m1 == mmax */
	i--;
	continue; 
      }
      m1 -= ibd[i];  /* ibd[i] can be 0 */     
      ibd[i] = 0;
      prob += log1r[i] - logr[i]; 
      i--;
    } else {
      ibd[i] = 1;
      m1++;
      prob += logr[i] - log1r[i];
      sum += exp(prob + sprob[m1]);  // that's where we add 0's and 1's
      i = nmid - 1;
    }
  } 
  return sum;
}

double probSxSyEqr(double logr, double log1r, double *sprob, double *sproblog,
		   int nm, int mmax, double *cbinom, int nmid, int m1) 
{
  int m;
  double sum;

  /* edge cases */
  if (nmid == 0) {
    if (m1 == 0) {              /* r = 0 */
      return sprob[0];
    }
    if (mmax < nm) {            /* r = 1 and not enough shared alleles */
      return 0;
    }
    return sprob[mmax];         /* r = 1 */
  }

  /* 0 < r < 1 */
  sum = 0;
  for (m = 0; m <= mmax; m++) {
    sum += exp(cbinom[m] + m*logr + (nm - m)*log1r + sproblog[m]);
  }
  return sum;
}

/******************************************************************************/
/*                                                                            */
/*                              P(Sx, Sy | m)                                 */
/*                                                                            */ 
/******************************************************************************/

void probSxSyCond(int *vx, int *vy, double *logpxy, double *logj, double *factj,
		  int numx, int numy, int nux, int nuy, int nuxy, int *ixy,
		  int *iyx, double combx, double comby, double *sprob,
		  int *mmax, int nm)
{
  int vxy[nuxy], s[nuxy], vxc[nuxy], vyc[nuxy], vmax[nuxy];
  int i, nums, nsi, mm;  // numerators: numx = nx; numy = ny
  double prob; 

  prob = combx + comby;
  nums = 0;

  /* find vxy (for Sxy), initiate s with 0's, calculate mmax, vmax */
  *mmax = 0;
  for (i = 0; i < nuxy; i++) {
    s[i]   = 0;                       /* first subset - empty */
    vxc[i] = vx[ixy[i]];
    vyc[i] = vy[iyx[i]];
    vxy[i] = MIN(vxc[i], vyc[i]);
    *mmax += vxy[i];
  }
  *mmax = MIN(*mmax, nm);  
  vMax(vxy, nuxy, *mmax, vmax);  /* vmax equal vxy if *mmax = nm */
  mm = *mmax;  

  /* initialize sprob with 0's, fill in the independence case */
  sprob[0] = exp(prob);    /* independence, Sxy is an empty set */  
  for (i = 1; i <= mm; i++) sprob[i] = 0; 

  /* go through the subsets */
  i = nuxy - 1;

  // !!! note factj starts with 0 (but not logj) !!! 
  while (!equalArr(s, vmax, nuxy)) {
    if (s[i] == vxy[i] || nums == mm) {
      nsi = s[i]; 
      if (nsi == 0) {  // small gain in efficiency
	i--;
	continue;
      }
      prob -= factj[nums] + factj[numx]   + factj[numy]  // subtract old
            - factj[nsi]  - factj[vxc[i]] - factj[vyc[i]];
      s[i] = 0;
      vxc[i] += nsi;
      vyc[i] += nsi;
      nums -= nsi;
      numx += nsi;
      numy += nsi;
      prob += factj[nums] + factj[numx]   + factj[numy]  // add new
                          - factj[vxc[i]] - factj[vyc[i]]
            + logpxy[i]*nsi;
      i--;
    } else {
      s[i]++;
      vxc[i]--;
      vyc[i]--;
      nums++;
      numx--;
      numy--; 
      prob += logj[nums - 1] - logj[numx]   - logj[numy]
	    - logj[s[i] - 1] + logj[vxc[i]] + logj[vyc[i]]
	    - logpxy[i]; 
      sprob[nums] += exp(prob);    
      i = nuxy - 1;
    }
  }
}

void vMax(int *base, int nbase, int sumlim, int *vmax)
{
  int i, rem;
  for (i = 0; i < nbase; i++) {
    vmax[i] = 0;
  }
  rem = sumlim;
  i = 0;
  while (rem > 0) {
    vmax[i] = MIN(base[i], rem);
    rem -= base[i];
    i++;
  }
}

int equalArr(int *arr1, int *arr2, int narr)
{
  int i;
  for (i = 0; i < narr; i++) {
    if (arr1[i] != arr2[i]) {
      return 0;
    }
  }
  return 1;
}


  
