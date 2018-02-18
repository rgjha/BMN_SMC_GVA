#include "evolve_fields.h"

void evolve_fields(Gauge_Field &U, Gauge_Field &p_U, Gauge_Field &f_U,
                   Site_Field phi[], Site_Field p_phi[], Site_Field f_phi[],
                   Site_Field F[], Site_Field p_F[], Site_Field f_F[]) {

  Lattice_Vector x, y;
  Gauge_Field nf_U;
  Site_Field nf_phi[NSCALAR], nf_F[NFERMION];
  int mu, i, site;

  // Leap frog algorithm
  // Update fields
  site = 0;
  while(loop_over_lattice(x, site)) {
    for (mu = 0; mu < D; mu++)
      U.set(x, mu, exp(DT * p_U.get(x, mu) + 0.5 * DT * DT * f_U.get(x, mu)) * U.get(x, mu));
  }

  for (i = 0; i < NSCALAR; i++)
    phi[i] = phi[i] + DT * p_phi[i] + 0.5 * DT * DT * f_phi[i];

  for (i = 0; i < NFERMION; i++)
    F[i] = F[i] + DT * p_F[i] + 0.5 * DT * DT * f_F[i];

  // Update forces
  force(U, nf_U, phi, nf_phi, F, nf_F);

  // Update momenta
  site = 0;
  while(loop_over_lattice(x, site)) {
    for (mu = 0; mu < D; mu++)
      p_U.set(x,mu,p_U.get(x,mu)+0.5*DT*(nf_U.get(x,mu)+f_U.get(x,mu)));
  }

  for (i = 0; i < NSCALAR; i++)
    p_phi[i]=p_phi[i]+0.5*DT*(nf_phi[i]+f_phi[i]);

  for (i = 0; i < NFERMION; i++)
    p_F[i]=p_F[i]+0.5*DT*(nf_F[i]+f_F[i]);

  // Store final forces for next iteration
  f_U = nf_U;
  for (i = 0; i < NSCALAR; i++)
    f_phi[i]=nf_phi[i];
  for (i = 0; i < NFERMION; i++)
    f_F[i]=nf_F[i];

  return;
}
