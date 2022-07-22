#define MAX(a, b) ((a) > (b)  ? (a) : (b))
#define MIN(a, b) ((a) < (b)  ? (a) : (b))
#define EPS pow(10, -9)

void probSxSyCond(int *vx, int *vy, double *logpxy, double *logj, double *factj,
		  int numx, int numy, int nux, int nuy, int nuxy, int *ixy,
		  int *iyx, double combx, double comby, double *sprob,
		  int *mmax, int nm);
double probSxSy(double *logr, double *log1r, double prob, double *sprob, int nm,
		int mmax, int nmid, int m1);
double probSxSyEqr(double logr, double log1r, double *sprob, double *sproblog,
		   int nm, int mmax, double *cbinom, int nmid, int m1);
int equalArr(int *arr1, int *arr2, int narr);
void vMax(int *base, int nbase, int sumlim, int *vmax);

