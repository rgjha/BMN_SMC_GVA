#include "myrandom.h"
#include "utilities.h"

#include <stdio.h>
#include <sys/time.h>

unsigned int random_seed() {
  unsigned int seed;
  struct timeval tv;
  FILE *devrandom;

  // This should let us specify the seed by inputting a positive number
  if (SEED > 0) {
    printf("Got seed %d from user input\n", SEED);
    return SEED;
  }

  if ((devrandom = fopen("/dev/random","r")) == NULL) {
    gettimeofday(&tv,0);
    seed = tv.tv_sec + tv.tv_usec;
    printf("Got seed %u from gettimeofday()\n",seed);
  }
  else {
    fread(&seed,sizeof(seed),1,devrandom);
    printf("Got seed %u from /dev/random\n",seed);
    fclose(devrandom);
  }

  return(seed);
}


#include <stdlib.h>

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

double ran3(long *idum) {
  static int inext,inextp;
  static long ma[56];
  static int iff=0;
  long mj,mk;
  int i,ii,k;

  if (*idum < 0 || iff == 0) {
    iff=1;
    mj=labs(MSEED-labs(*idum));
    mj %= MBIG;
    ma[55]=mj;
    mk=1;
    for (i=1;i<=54;i++) {
      ii=(21*i) % 55;
      ma[ii]=mk;
      mk=mj-mk;
      if (mk < MZ) mk += MBIG;
      mj=ma[ii];}
    for (k=1;k<=4;k++)
      for (i=1;i<=55;i++) {
        ma[i] -= ma[1+(i+30) % 55];
        if (ma[i] < MZ) ma[i] += MBIG;}
    inext=0;
    inextp=31;
    *idum=1;}
  if (++inext == 56) inext=1;
  if (++inextp == 56) inextp=1;
  mj=ma[inext]-ma[inextp];
  if (mj < MZ) mj += MBIG;
  ma[inext]=mj;
  return mj*FAC;
}

#undef MBIG
#undef MSEED
#undef MZ
#undef FAC

double myrandom(void) {
  static long idum;
  static int first_time=1;

  if(first_time) {
    idum=(-random_seed());
    first_time=0;
  }

  return(ran3(&idum));
}

double gasdev(void) {
  static int iset=0;
  static double gset;
  double fac,rsq,v1,v2;
  if(iset==0) {
    do{
      v1=2.0*myrandom()-1.0;
      v2=2.0*myrandom()-1.0;
      rsq=v1*v1+v2*v2;
    }
    while(rsq>=1.0 || rsq == 0.0);
    fac=sqrt(-2.0*log(rsq)/rsq);
    gset=v1*fac;
    iset=1;
    return(v2*fac);}
  else {
    iset=0;
    return(gset);}
}
