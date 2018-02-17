#include "utilities.h"

void eigenvalues(Complex M[FERMIONSIZE][FERMIONSIZE],double &MIN, double&MAX);

extern "C" void zgeev_( char*, char*, int*, double at[], int *, double b[], double dummy[],
                       int *, double dummy2[], int*, double work[], int *, double work2[], int *);


