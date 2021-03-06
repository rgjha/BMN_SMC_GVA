#include "sym.h"

int SWEEPS, GAP, THERM, SEED, READIN;
double BETA, DT, MASS, TIME, MU;

double amp[DEGREE], shift[DEGREE], ampdeg;
Umatrix Lambda[RANK];
Gamma_Matrix Gamma[NSCALAR - 2], Gam123;
double f[RANK][RANK][RANK];
double LARGECUT, SMALLCUT;
int TRAJECTORY_LENGTH;
int epsilon[3][3][3];

int main(void) {
  int i, sweep, sites;
  double WIDTH = 0.1;
  Gauge_Field U;
  Site_Field phi[NSCALAR];
  Site_Field F[NFERMION];

  SEED = -1;
  read_param();

  if (READIN) {
    read_in(U, phi, F);
    check(U);
  }
  else {
    U = Gauge_Field(1);
    for (i = 0; i < NSCALAR; i++)
      phi[i] = WIDTH * Site_Field(1);
    if (FERMIONS) {
      for (i = 0; i < NFERMION; i++)
        F[i] = Site_Field(2);
    }
  }

  cout << "Warming up" << "\n" << flush;

  DT = 0.5 * DT;
  cout << "using thermalization DT " << DT << "\n" << flush;
  for (sweep = 1; sweep <= THERM; sweep++) {
    update(U, phi, F);
    write_out(U, phi, F);
  }

  cout << "Commencing measurement sweeps" << "\n" << flush;
  DT = DT * 2;
  cout << "Updated DT is " << DT << "\n" << flush;
  for (sweep = 1; sweep <= SWEEPS; sweep++) {
    update(U, phi, F);
    write_out(U, phi, F);

    // debug: set all matrices in vacuum
#if 0
    Umatrix vac=Umatrix();
    for (int i = 0; i < NCOLOR - 1; i++)
      vac.set(i, i, Complex(0.0,1.0));
    vac.set(NCOLOR - 1, NCOLOR - 1, Complex(0.0, 1.0 - NCOLOR));

    int sites = 0,mu;
    Lattice_Vector x;

    sites = 0;
    while (loop_over_lattice(x,sites)) {
      for (i = 0; i < NSCALAR; i++)
        phi[i].set(x,myrandom() * vac);

      for (mu = 0; mu < D; mu++)
        U.set(x, mu, exp(0.1 * vac));
    }
#endif
    check(U);

    // Measure config
    cout << "sweep no. " << sweep << "\n" << flush;
    if (sweep % GAP == 0) {
      measure(U, phi, F);
      write_out(U, phi, F);
    }
  }

  write_out(U, phi, F);
  return 0;
}
