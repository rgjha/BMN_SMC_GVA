#include "kinetic_energy.h"

double kinetic_energy(const Gauge_Field &p_U, const Site_Field p_phi[],
                      const Site_Field p_F[]) {

  int sites, mu, i;
  Complex dum = Complex();
  Lattice_Vector x;

  // note minus sign to take into account AH Lambdas
  sites = 0;
  while (loop_over_lattice(x,sites)){
    for (mu = 0; mu < D; mu++)
      dum = dum + 0.5 * Tr(Adj(p_U.get(x, mu)) * p_U.get(x, mu));
  }

  for (i = 0; i < NSCALAR; i++)
    dum = dum + 0.5 * Tr(Adj(p_phi[i]) * p_phi[i]);

  for (i = 0; i < NFERMION; i++)
    dum = dum + Tr(Adj(p_F[i]) * p_F[i]);

  return dum.real();

}
