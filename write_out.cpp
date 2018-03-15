#include "write_out.h"

void write_out(const Gauge_Field &u, const Site_Field phi[],
               const Site_Field F[]) {

  static ofstream f_write;
  int j, mu, site;
  Lattice_Vector x;

  f_write.open("/tmp/dump");   // !!! Needs absolute path
  if (f_write.bad())
    cout << "error opening config file\n" << flush;

  f_write << L << "\t" << T << "\t" << BETA << "\t" << DT << "\t" << NCOLOR << "\n";

  site = 0;
  while (loop_over_lattice(x, site)) {
    for (mu = 0; mu < D; mu++)
      f_write << u.get(x, mu) << "\n";
  }
  f_write << "\n" << flush;

  for (j = 0; j < NSCALAR;j++) {
    site = 0;
    while (loop_over_lattice(x, site))
      f_write << phi[j].get(x) << "\n";
  }
  f_write << "\n" << flush;

  if (FERMIONS) {
    for (j = 0; j < NFERMION; j++) {
      site = 0;
      while (loop_over_lattice(x, site))
        f_write << F[j].get(x) << "\n";
    }
  }

  f_write.close();
  return;
}
