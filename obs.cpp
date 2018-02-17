#include "obs.h"

void obs(const Gauge_Field &U, const Site_Field phi[], const Site_Field
F[NFERMION],double &act_s, double &act_F,double eigenvals[T][NSCALAR][NCOLOR]){
Lattice_Vector x,e_mu;
int mu,nu,site,i,j,k,s;
Gauge_Field Udag;
Umatrix av;
Site_Field sol[NFERMION][DEGREE],psol[NFERMION][DEGREE]; 

double b[2*NCOLOR],dummy[2],work[4*NCOLOR];
double at[2*NCOLOR*NCOLOR];

act_F=act_s=0.0;

//scalar action
  
Udag=Adj(U);
  
//scalar action
  
for(i=0;i<NSCALAR;i++){
  
site=0;
while(loop_over_lattice(x,site)){				
for(mu=0;mu<D;mu++){	     
e_mu=Lattice_Vector(mu);

act_s=act_s+BETA*Tr(phi[i].get(x)*(
U.get(x,mu)*phi[i].get(x+e_mu)*Udag.get(x,mu)+                                  
Udag.get(x-e_mu,mu)*phi[i].get(x-e_mu)*U.get(x-e_mu,mu)-
2.0*phi[i].get(x))).real();
 
}
}  
}

for(i=0;i<NSCALAR;i++){
for(j=i+1;j<NSCALAR;j++){

site=0;
while(loop_over_lattice(x,site)){
act_s=act_s-BETA*Tr(comm(phi[i].get(x),phi[j].get(x))*
                    comm(phi[i].get(x),phi[j].get(x))).real();
}
}
}
// BMN  mass terms

for(i=0;i<3;i++){
site=0;
while(loop_over_lattice(x,site)){
act_s=act_s-BETA*MU*MU*Tr(phi[i].get(x)*phi[i].get(x)).real();
}
}

for(i=3;i<NSCALAR;i++){
site=0;
while(loop_over_lattice(x,site)){
act_s=act_s-BETA*0.25*MU*MU*Tr(phi[i].get(x)*phi[i].get(x)).real();
}}

// cubic term

for(i=0;i<3;i++){
for(j=0;j<3;j++){
for(k=0;k<3;k++){
if((i==j) || (i==k) || (j==k)) continue;
site=0;
while(loop_over_lattice(x,site)){
act_s=act_s-BETA*2*sqrt(2.0)*MU*epsilon[i][j][k]*Tr(phi[i].get(x)*phi[j].get(x)*phi[k].get(x)).real();
}}}}


// pseudofermion contribution
if(FERMIONS){

for(i=0;i<NFERMION;i++){
act_F=act_F+ampdeg*Tr(Adj(F[i])*F[i]).real();}


MCG_solver(U,phi,F,shift,sol,psol);

for(int n=0;n<DEGREE;n++){
for(i=0;i<NFERMION;i++){
act_F=act_F+amp[n]*Tr(Adj(F[i])*sol[i][n]).real();}
}


}


site=0;
while(loop_over_lattice(x,site)){
for(s=0;s<NSCALAR;s++){

av=phi[s].get(x);


int i,j,ok,c1,c2,c3;
char c4;


for(i=0;i<NCOLOR;i++){
for(j=0;j<NCOLOR;j++){
at[2*(j+NCOLOR*i)]=av.get(j,i).real();
at[2*(j+NCOLOR*i)+1]=av.get(j,i).imag();
}
}

c1=NCOLOR;
c2=2*NCOLOR;
c3=1;
c4='N';

zgeev_(&c4,&c4,&c1,at,&c1,b,dummy,&c3,dummy,&c3,work,&c2,work,&ok);

for(j=0;j<2*NCOLOR;j=j+2){
eigenvals[x.get(D-1)][s][j/2]=b[j+1];}

}
}

return;
}
