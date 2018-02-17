#include "eigenvalues.h"

void eigenvalues(Complex M[FERMIONSIZE][FERMIONSIZE],double &MIN, double &MAX){
int i,j,ok,c1,c2,c3;
char c4;
double b[2*FERMIONSIZE],dummy[2],work[4*FERMIONSIZE];
double at[2*FERMIONSIZE*FERMIONSIZE];


for(i=0;i<FERMIONSIZE;i++){
for(j=0;j<FERMIONSIZE;j++){
at[2*(j+FERMIONSIZE*i)]=M[j][i].real();
at[2*(j+FERMIONSIZE*i)+1]=M[i][j].imag();
}
}

c1=FERMIONSIZE;
c2=2*FERMIONSIZE;
c3=1;
c4='N';

zgeev_(&c4,&c4,&c1,at,&c1,b,dummy,&c3,dummy,&c3,work,&c2,work,&ok);

MIN=1000.0;
MAX=0.0;
for(j=0;j<2*FERMIONSIZE;j=j+2){
//cout << "( " << b[j] << "," << b[j+1] << ")" << "\n";
if ((b[j]*b[j]+b[j+1]*b[j+1])<MIN) MIN=b[j]*b[j]+b[j+1]*b[j+1];
if ((b[j]*b[j]+b[j+1]*b[j+1])>MAX) MAX=b[j]*b[j]+b[j+1]*b[j+1];
}

MIN=NORM*NORM*MIN;
MAX=NORM*NORM*MAX;
return;
}
