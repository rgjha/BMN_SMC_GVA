#include "Gamma_Matrices.h"
#define B 2

#define NR_END 1
#define FREE_ARG char*

double **matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
long i,nrow=nrh-nrl+1,ncol=nch-ncl+1;
double **m;

/* allocate pointers to rows */

m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double *)));
if(!m) printf("allocation failure 1 in matrix()");
m+= NR_END;
m-= nrl;

/* allocate rows and set pointers to them */
m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
if(!m[nrl]) printf("allocation failure 2 in matrix()");
m[nrl]+= NR_END;
m[nrl]-= ncl;

for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

/* return pointer to array of pointers to rows */
return m;
}

void free_matrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by matrix () */
{
free((FREE_ARG) (m[nrl]+ncl-NR_END));
free((FREE_ARG) (m+nrl-NR_END));
}

void direct_prod(double **A, double **D, int na, int nb, double **C)
{
int i,j,k,l,I,J;
//printf("in direct_prod with dimension %d\n",na);
for(i=0;i<na;i++){
for(j=0;j<na;j++){
for(k=0;k<nb;k++){
for(l=0;l<nb;l++){
I=nb*i+k;
J=nb*j+l;
C[I][J]=A[i][j]*D[k][l];
}}}}

return;
}



void Gamma_Matrices(void){
//int main(){

double **sigma1,**sigma2,**sigma3,**unit;
double **b1,**b2;
int i,j,k;

//printf("beginning\n");
cout << "Computing Gamma matrices\n";

sigma1=matrix(0,B-1,0,B-1);
sigma2=matrix(0,B-1,0,B-1);
sigma3=matrix(0,B-1,0,B-1);
unit=matrix(0,B-1,0,B-1);
b1=matrix(0,2*B-1,0,2*B-1);
b2=matrix(0,4*B-1,0,4*B-1);

//printf("assigning sigma etc\n");
sigma1[0][0]=0.0;
sigma1[0][1]=1.0;
sigma1[1][0]=1.0;
sigma1[1][1]=0.0;
    
sigma2[0][0]=0.0;
sigma2[0][1]=-1.0;
sigma2[1][0]=1.0;
sigma2[1][1]=0.0;
    
sigma3[0][0]=1.0;
sigma3[0][1]=0.0;
sigma3[1][0]=0.0;
sigma3[1][1]=-1.0;
    
unit[0][0]=1.0;
unit[0][1]=0.0;
unit[1][0]=0.0;
unit[1][1]=1.0;

// compute 16 dim Majorana rep in terms of (real) 8 dim chiral rep.

direct_prod(sigma2,sigma2,B,B,b1);
direct_prod(b1,sigma2,(2*B),B,b2);


for(i=0;i<4*B;i++){
for(j=0;j<4*B;j++){
Gamma[0].set(i,j,b2[i][j]);
}}
    
    

direct_prod(sigma1,sigma2,B,B,b1);
direct_prod(b1,unit,2*B,B,b2);


for(i=0;i<4*B;i++){
for(j=0;j<4*B;j++){
Gamma[1].set(i,j,b2[i][j]);
}}

direct_prod(sigma3,sigma2,B,B,b1);
direct_prod(b1,unit,2*B,B,b2);

for(i=0;i<4*B;i++){
for(j=0;j<4*B;j++){
Gamma[2].set(i,j,b2[i][j]);
}}

direct_prod(sigma2,unit,B,B,b1);
direct_prod(b1,sigma1,2*B,B,b2);


for(i=0;i<4*B;i++){
for(j=0;j<4*B;j++){
Gamma[3].set(i,j,b2[i][j]);
}}

direct_prod(sigma2,unit,B,B,b1);
direct_prod(b1,sigma3,2*B,B,b2);


for(i=0;i<4*B;i++){
for(j=0;j<4*B;j++){
Gamma[4].set(i,j,b2[i][j]);
}}

direct_prod(unit,sigma1,B,B,b1);
direct_prod(b1,sigma2,2*B,B,b2);


for(i=0;i<4*B;i++){
for(j=0;j<4*B;j++){
Gamma[5].set(i,j,b2[i][j]);
}}

direct_prod(unit,sigma3,B,B,b1);
direct_prod(b1,sigma2,2*B,B,b2);


for(i=0;i<4*B;i++){
for(j=0;j<4*B;j++){
Gamma[6].set(i,j,b2[i][j]);
}}
    
    
    

Gam123=Gamma_Matrix();
for(i=0;i<3;i++){
for(j=0;j<3;j++){
for(k=0;k<3;k++){
if((i==j) || (i==k) || (j==k)) continue;
Gam123=Gam123+(epsilon[i][j][k]*1.0/6.0)*(Gamma[i]*Gamma[j]*Gamma[k]);
}}}



free_matrix(sigma1,0,B-1,0,B-1);
free_matrix(sigma2,0,B-1,0,B-1);
free_matrix(sigma3,0,B-1,0,B-1);
free_matrix(unit,0,B-1,0,B-1);
free_matrix(b1,0,2*B-1,0,2*B-1);
free_matrix(b2,0,4*B-1,0,4*B-1);

//return(0);
return;
}
