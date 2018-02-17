#include "utilities.h"

#ifndef CG_RESIDUAL_H
#define CG_RESIDUAL_H
const double CG_RESIDUAL = 0.00000000000001;
#endif


void MCG_solver(const Gauge_Field &, const Site_Field phi[], const Site_Field rhs[NFERMION], double shift[], 
Site_Field sol[NFERMION][DEGREE], 
Site_Field psol[NFERMION][DEGREE]);
