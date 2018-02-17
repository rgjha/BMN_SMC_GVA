#include "utilities.h"
#include "action.h"

void obs(const Gauge_Field &U, const Site_Field phi[NSCALAR], const Site_Field
F[NFERMION],
double &act_s, double &act_F,double eigenvals[T][NSCALAR][NCOLOR]);

extern "C" void zgeev_( char*, char*, int*, double at[], int *, double b[], double dummy[],
                       int *, double dummy2[], int*, double work[], int *, double work2[], int *);
