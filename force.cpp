#include "force.h"

void force(const Gauge_Field &U, Gauge_Field &f_U, const Site_Field phi[],
           Site_Field f_phi[], const Site_Field F[], Site_Field f_F[]) {

  Lattice_Vector x, e_mu, xm1;
  int sites, mu, nu, i, j, k, n;
  double td;
  Complex trace;
  Site_Field sol[NFERMION][DEGREE], psol[NFERMION][DEGREE];
  Gauge_Field Udag, utmp;
  Site_Field ptmp[NFERMION],stmp[NFERMION];
  Site_Field phitmp[NSCALAR];
  Umatrix tmp, tU;

  // Gauge force -- scalar contribution
  Udag = Adj(U);
  sites = 0;
  while (loop_over_lattice(x, sites)) {
    for (mu = 0; mu < D; mu++) {
      e_mu = Lattice_Vector(mu);
      tmp = Umatrix();

      for (i = 0; i < NSCALAR; i++) {
        tU = U.get(x, mu) * phi[i].get(x + e_mu) * Udag.get(x, mu);
        tmp = tmp + phi[i].get(x) * tU;
      }

      // Make traceless anti-hermitian
      tmp = tmp - Adj(tmp);
      trace = Tr(tmp);
      if (trace.norm() > TRACETOL)
        tmp = tmp - (1.0 / NCOLOR) * trace * Umatrix(1);

      f_U.set(x, mu, -2.0 * BETA * tmp);
    }
  }

  // Scalar force contributions from scalar--scalar self-interactions
  for (i = 0; i < NSCALAR; i++) {
    sites = 0;
    while (loop_over_lattice(x, sites)) {
      tmp = Umatrix();
      for  (j = 0; j < NSCALAR; j++) {
        if (i == j)
          continue;
        tmp = tmp + comm(phi[j].get(x), comm(phi[i].get(x), phi[j].get(x)));
      }

      trace = Tr(tmp);
      if (trace.norm() > TRACETOL)
        tmp = tmp - (1.0 / NCOLOR) * trace * Umatrix(1);

      f_phi[i].set(x, -2.0 * BETA * tmp);
    }
  }

  // Scalar force contributions from mass terms
  td = 2.0 * BETA * MU * MU;
  for (i = 0; i < 3; i++) {
    sites = 0;
    while (loop_over_lattice(x, sites)) {
      f_phi[i].set(x, f_phi[i].get(x) - td * phi[i].get(x));
    }
  }

  td = 0.5 * BETA * MU * MU;
  for (i = 3; i < NSCALAR; i++) {
    sites = 0;
    while (loop_over_lattice(x, sites))
      f_phi[i].set(x, f_phi[i].get(x) - td * phi[i].get(x));
  }

  td = BETA * 6.0 * sqrt(2.0) * MU;
  for (i = 0; i < 3; i++) {
    sites = 0;
    while (loop_over_lattice(x, sites)) {
      tmp = Umatrix();
      for (j = 0; j < 3; j++) {
        if (i == j)
          continue;
        for (k = 0; k < 3; k++) {
          if ((i==k) || (j==k))
            continue;
          tmp = tmp + phi[j].get(x) * phi[k].get(x) * epsilon[i][j][k];
        }
      }
      f_phi[i].set(x, f_phi[i].get(x) - td * tmp);
    }
  }

  // Scalar force contributions from kinetic term
  td = 2.0 * BETA;
  for (i = 0; i < NSCALAR; i++) {
    sites = 0;
    while (loop_over_lattice(x, sites)) {
      tmp = Umatrix();
      for (mu = 0; mu < D; mu++) {
        e_mu = Lattice_Vector(mu);
        xm1 = x - e_mu;
        tU = Udag.get(xm1, mu) * phi[i].get(xm1) * U.get(xm1, mu);
        tmp = tmp + U.get(x, mu) * phi[i].get(x + e_mu) * Udag.get(x, mu)
                  + tU - 2.0 * phi[i].get(x);
      }

      trace = Tr(tmp);
      if (trace.norm() > TRACETOL)
        tmp = tmp - (1.0 / NCOLOR) * trace * Umatrix(1);

      f_phi[i].set(x, f_phi[i].get(x) + td * tmp);
    }
  }

  // Have checked that scalar forces are already traceless
//  for (i = 0; i < NSCALAR; i++) {
//    sites = 0;
//    while (loop_over_lattice(x, sites)) {
//      trace = Tr(f_phi[i].get(x));
//      if (trace.norm() > TRACETOL)
//        f_phi[i].set(x, f_phi[i].get(x) - (1.0 / NCOLOR) * trace * Umatrix(1));
//    }
//  }

  // Add contributions from pseudofermions
  // Use partial fraction approx to inverse square root of operator
  for (i = 0; i < NFERMION; i++)
    f_F[i] = Site_Field();

  if (FERMIONS) {
    for (i = 0; i < NFERMION; i++) {
      f_F[i] = ampdeg * Adj(F[i]);}

    MCG_solver(U, phi, F, shift, sol, psol);
    for (n = 0; n < DEGREE; n++) {
      for (i = 0; i < NFERMION; i++) {
        stmp[i] = sol[i][n];
        ptmp[i] = psol[i][n];
      }

      fermion_forces(U, utmp, phi, phitmp, stmp, ptmp);
      for (i = 0; i < NFERMION; i++)
        f_F[i] = f_F[i] + amp[n] * Adj(sol[i][n]);

      // Add in kick from fermion effective action
      sites = 0;
      while (loop_over_lattice(x, sites)) {
        for (mu = 0; mu < D; mu++)
          f_U.set(x, mu, f_U.get(x, mu) + amp[n] * utmp.get(x, mu));
      }

      for (i = 0; i < NSCALAR; i++)
        f_phi[i] = f_phi[i] + amp[n] * phitmp[i];
    }
  }

  for (i = 0; i < NFERMION; i++)
    f_F[i] = Adj(-1.0 * f_F[i]);

  return;
}
