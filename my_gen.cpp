#include "my_gen.h"

const int POS1 = (NCOLOR*(NCOLOR-1)/2);
const int POS2 = (NCOLOR*(NCOLOR-1));

const double ROOT2 = 1.41421356237309504880168872421;

int posmat[NCOLOR][NCOLOR];
Complex genmat[RANK][NCOLOR][NCOLOR];
Complex tr[RANK][RANK], strconst[RANK][RANK][RANK];

/* SU(N) routines */
void initgen() {
  int i, j, a = 0;

  for (i = 0; i < NCOLOR; i++) {
    for (j = i + 1; j < NCOLOR; j++) {
      posmat[i][j] = a;
      a += 1;
    }
  }
}

/* Adjoint to fundamental */
void vectomat(Complex vec[RANK], Complex mat[NCOLOR][NCOLOR]) {
  int i, j;
  Complex mult;

  for (i = 0; i < NCOLOR; i++) {
    for (j = i + 1; j < NCOLOR; j++) {
      mat[i][j] = 0.5 * (vec[posmat[i][j]] - Complex(0.0, 1.0) * vec[POS1 + posmat[i][j]]);
      mat[j][i] = 0.5 * (vec[posmat[i][j]] + Complex(0.0, 1.0) * vec[POS1 + posmat[i][j]]);
    }
  }

  for (i = 0; i < NCOLOR; i++)
    mat[i][i] = Complex(0.0, 0.0);

  for (i = 0; i < NCOLOR - 1; i++) {
    mult = vec[POS2 + i] * (1.0 / sqrt(2.0 + 2.0 / (1.0 + i)) / (1.0 + i));
    for (j = 0; j < i + 1; j++)
      mat[j][j] = mat[j][j] + mult;

    mat[i + 1][i + 1] = mat[i + 1][i + 1] - 1.0 * Complex(1.0 + i, 0.0) * mult;
  }
}

/* Compute generator matrices */
void computegen() {
  int a, i, j;
  Complex adj[RANK];
  Complex fund[NCOLOR][NCOLOR];

  for (a = 0; a < RANK; a++) {
    for (i = 0; i < RANK; i++)
      adj[i] = Complex(0.0, 0.0);
    adj[a] = Complex(1.0, 0.0);

    vectomat(adj, fund);

    for (i = 0; i < NCOLOR; i++) {
      for (j = 0; j < NCOLOR; j++)
        genmat[a][i][j] = fund[i][j];
    }
  }
}

int my_gen() {
  int a, i, j;

  initgen();
  computegen();

  /* map to data structures in code */
  for (a = 0; a < RANK; a++) {
    for (i = 0; i < NCOLOR; i++) {
      for (j = 0; j < NCOLOR; j++)
        Lambda[a].set(i, j, (Complex(0.0, 1.0) * sqrt(2.0)) * genmat[a][i][j]);
    }
//    Lambda[a].print();
  }
  return 1;
}
