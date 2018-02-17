#include "utilities.h"
#include "MCG_solver.h"
#include "fermion_forces.h"

void force(const Gauge_Field &U, Gauge_Field &f_U, const Site_Field phi[],
           Site_Field f_phi[],
           const Site_Field F[], Site_Field f_F[]);
