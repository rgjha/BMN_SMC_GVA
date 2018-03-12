#include "read_in.h"

void read_in(Gauge_Field &u, Site_Field phi[], Site_Field F[]) {
  int j, site, LDUM, TDUM;
  double BETADUM, NCOLORDUM, DTDUM;
  Umatrix dummy;
  Lattice_Vector x;
  Site_Field s;
  ifstream f_read;

  f_read.open("config");
  if (f_read.bad())
    cout << "error opening config file to read\n" << flush;

  f_read >> LDUM >> TDUM >> BETADUM >> DTDUM >> NCOLORDUM; 
  if ((LDUM != L) || (TDUM != T))
    cout << "wrong size lattice read in - abort\n";

  cout << "config coupling is " << BETADUM << "\n";
  cout << "time step used for input config " << DTDUM << "\n" << flush;
  cout << "number of colors " << NCOLORDUM << "\n" << flush;

  site = 0;
  while (loop_over_lattice(x, site)) {
    for (j = 0; j < D; j++) {
      f_read >> dummy;
      u.set(x, j, dummy);
    }
  }

  for (j = 0; j < NSCALAR; j++) {
    site = 0;
    while (loop_over_lattice(x, site)) {
      f_read >> dummy;
      s.set(x, dummy);
    }
    phi[j] = s;
  }

  if (FERMIONS) {
    for (j = 0; j < NFERMION; j++) {
      site = 0;
      while (loop_over_lattice(x, site)) {
        f_read >> dummy;
        s.set(x,dummy);
      }
      F[j] = s;
    }
  }

  f_read.close();
  return;
}
