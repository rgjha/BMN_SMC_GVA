 #include "action.h"

double action(const Gauge_Field &U, const Site_Field phi[NSCALAR],
              const Site_Field F[NFERMION]) {

  int mu, nu, site, i, j, k, n;
  double act_s = 0.0, act_F = 0.0, td, td2;
  double Myers = 0.0, so3 = 0.0, so6 = 0.0;
  Lattice_Vector x, e_mu;
  Umatrix tU;
  Gauge_Field Udag;
  Site_Field sol[NFERMION][DEGREE],psol[NFERMION][DEGREE];    // NFERMION = 16

  Udag = Adj(U);

  // Scalar kinetic term
  //   Tr[phi_i(x) {  U_mu(x) phi_i(x + mu) * Udag_mu(x)
  //                + Udag_mu(x - mu) * phi_i(x - e_mu) * U_mu(x - mu)
  //                - 2phi_i(x)}]
  // Equivalently (by cyclicity and sum over x)
  //   2Tr[phi_i(x) {U_mu(x) phi_i(x + mu) Udag_mu(x) - phi_i(x)}]
  td = 2.0 * BETA;
  for (i = 0; i < NSCALAR; i++) {           // NSCALAR = 9
    site = 0;
    while (loop_over_lattice(x, site)) {
      for (mu = 0; mu < D; mu++) {         // D = 1
        e_mu = Lattice_Vector(mu);

        tU = U.get(x, mu) * phi[i].get(x + e_mu) * Udag.get(x,mu);
        tU = tU - phi[i].get(x);
        act_s = act_s + td * Tr(phi[i].get(x) * tU).real();
      }
    }
  }

  // Scalar--scalar potential -sum_{i<j} Tr[phi_i, phi_j]^2
  for (i = 0; i < NSCALAR; i++) {
    for (j = i + 1; j < NSCALAR; j++) {
      site = 0;
      while (loop_over_lattice(x, site)) {
        tU = comm(phi[i].get(x), phi[j].get(x));
        act_s = act_s - BETA * Tr(tU * tU).real();
      }
    }
  }

  // BMN mass terms -mu^2/9  Tr[phi_i phi_i] for i = 0, 1, 2
  //                -mu^2/36 Tr[phi_i phi_i] for i = 3, 4, 5, 6, 7, 8
  td = BETA * MU * MU / 9.0;
  for (i = 0; i < 3; i++) {
    site = 0;
    while (loop_over_lattice(x, site)) {
      td2 = td * Tr(phi[i].get(x) * phi[i].get(x)).real();
      act_s = act_s - td2;
      so3 = so3 - td2;
    }
  }

  td = BETA * MU * MU / 36.0;
  for (i = 3; i < NSCALAR; i++) {
    site = 0;
    while (loop_over_lattice(x, site)) {
      td2 = td * Tr(phi[i].get(x) * phi[i].get(x)).real();
      act_s = act_s - td2;
      so6 = so6 - td2;
    }
  }

  // Cubic term (-sqrt(8) mu / 3) epsilon_ijk phi_i phi_j phi_k
  td = BETA * sqrt(2.0) * MU / 1.5;
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
        if (i == j)
          continue;
      for (k = 0; k < 3; k++) {
        if ((i == k) || (j == k))
          continue;
        site = 0;
        while (loop_over_lattice(x, site)) {
          tU = phi[i].get(x) * phi[j].get(x) * phi[k].get(x);
          td2 = td * epsilon[i][j][k] * Tr(tU).real();
          act_s = act_s - td2;
          Myers = Myers - td2;
        }
      }
    }
  }
  cout << "so3 " << so3 << " so6 " << so6 << " Myers " << Myers
       << " boson " << act_s << endl;

  // Pseudofermion contribution
  if (FERMIONS) {
    for (i = 0; i < NFERMION; i++)
      act_F = act_F + ampdeg * Tr(Adj(F[i]) * F[i]).real();

    MCG_solver(U,phi,F, shift, sol,psol);

    for (n = 0; n < DEGREE; n++) {
      for (i = 0; i < NFERMION; i++)
        act_F = act_F + amp[n] * Tr(Adj(F[i]) * sol[i][n]).real();
    }
  }

  return (act_s + act_F);
}
