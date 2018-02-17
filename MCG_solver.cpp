#include "MCG_solver.h"
// multimass CG solver 
// can return solution to (M^daggerM +beta_i) x=b for all shifts beta_i
// for price of smallest shift 

void MCG_solver(const Gauge_Field &U, const Site_Field phi[],
const Site_Field b[NFERMION], double shift[], Site_Field sol[NFERMION][DEGREE], 
Site_Field psol[NFERMION][DEGREE]){

Site_Field Scalar[KDFERMION][KDFERMION];

Site_Field r[NFERMION],p[NFERMION];
double alpha1,alpha2,beta0,beta1,rrtmp,rrtmp2,resid,psdot;

Site_Field plist[NFERMION][DEGREE];

Site_Field s[NFERMION],t[NFERMION];
double alphalist[DEGREE],betalist[DEGREE],xi2[DEGREE];
double xi0[DEGREE],xi1[DEGREE];


int count2;
static int av_count2=0,no_calls=0;
int i,n;
double dummy;

no_calls++;
//cout << "in MCG-solver\n" << flush;

count2=0;
Build_Scalar(phi,Scalar);

// initialize solver

beta0=1;
alpha1=0.0;

for(n=0;n<DEGREE;n++){

alphalist[n]=alpha1;
betalist[n]=beta0;
xi0[n]=beta0;
xi1[n]=beta0;
xi2[n]=beta0;

for(i=0;i<NFERMION;i++){
sol[i][n]=Site_Field();
plist[i][n]=b[i];}

}

for(i=0;i<NFERMION;i++){
r[i]=b[i];
p[i]=b[i];}


do{

Fermion_operator(U,Scalar,phi,p,t,1);

Fermion_operator(U,Scalar,phi,t,s,-1);

rrtmp=psdot=0.0;
for(i=0;i<NFERMION;i++){
rrtmp=rrtmp+Tr(Adj(r[i])*r[i]).real();
psdot=psdot+Tr(Adj(p[i])*s[i]).real();}

beta1=-rrtmp/psdot;

for(n=0;n<DEGREE;n++){
xi2[n]=(xi1[n]*xi0[n]*beta0)/
 (
beta1*alpha1*(xi0[n]-xi1[n])+
xi0[n]*beta0*(1-shift[n]*beta1)
 );
 
if(fabs(xi2[n])<1.0e-50){xi2[n]=1.0e-50;}

betalist[n]=beta1*xi2[n]/xi1[n];
}

for(i=0;i<NFERMION;i++){
r[i]=r[i]+beta1*s[i];}

for(n=0;n<DEGREE;n++){

for(i=0;i<NFERMION;i++){
sol[i][n]=sol[i][n]-betalist[n]*plist[i][n];}

}

rrtmp2=0.0;
for(i=0;i<NFERMION;i++){
rrtmp2=rrtmp2+Tr(Adj(r[i])*r[i]).real();}

alpha2=rrtmp2/rrtmp;

for(n=0;n<DEGREE;n++){
alphalist[n]=alpha2*(xi2[n]/xi1[n])*(betalist[n]/beta1);
}

for(i=0;i<NFERMION;i++){
p[i]=r[i]+alpha2*p[i];
}

for(n=0;n<DEGREE;n++){
for(i=0;i<NFERMION;i++){
plist[i][n]=xi2[n]*r[i]+alphalist[n]*plist[i][n];
}}


resid=sqrt(rrtmp2/(2*FERMIONSIZE));


for(n=0;n<DEGREE;n++){
xi0[n]=xi1[n];
xi1[n]=xi2[n];
}

alpha1=alpha2;
beta0=beta1;


count2++;
//cout << "residual is " << resid << "\n" << flush;
}
while((resid>CG_RESIDUAL)&&(count2<(10*FERMIONSIZE)));

av_count2+=count2;


if(0){
cout<<"exiting residual is " << resid << " at " 
<< count2 << " iterations\n" <<
flush;}



if(no_calls%100==0){
cout << "average number of CG iterations " <<
(double)av_count2/(no_calls) << "\n" << flush;

no_calls=0;

av_count2=0;
}



// check solutions

for(int n=0;n<DEGREE;n++){

for(i=0;i<NFERMION;i++){
p[i]=sol[i][n];}

Fermion_operator(U,Scalar,phi,p,t,1);
for(i=0;i<NFERMION;i++){
psol[i][n]=t[i];}

Fermion_operator(U,Scalar,phi,t,s,-1);
for(i=0;i<NFERMION;i++){
s[i]=s[i]+shift[n]*p[i];
t[i]=b[i]-s[i];}

dummy=0.0;
for(i=0;i<NFERMION;i++){
dummy=dummy+Tr(Adj(t[i])*t[i]).real();}

if(sqrt(dummy/(2*FERMIONSIZE))>10000*CG_RESIDUAL){
cout << "poor " << n << "soln: " << sqrt(dummy)/(2*FERMIONSIZE) << "\n" <<
flush;
}
}


return;
}
