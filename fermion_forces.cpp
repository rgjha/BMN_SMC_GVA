#include "fermion_forces.h"

void fermion_forces(const Gauge_Field &U, Gauge_Field &f_U, 
const Site_Field phi[], Site_Field f_phi[],
const Site_Field s[], const Site_Field p[]){
Lattice_Vector x,e_mu;
int i,a,b,sites,mu;
Umatrix tmp;
Gauge_Field Udag;
Site_Field pdag[NFERMION];

Udag=Adj(U);
for(i=0;i<NFERMION;i++){
pdag[i]=Adj(p[i]);}

f_U=Gauge_Field();
for(i=0;i<NSCALAR;i++){
f_phi[i]=Site_Field();}

// compute fermion kick to gauge link force
// f_U=Cjg(A).D_U M(U) sol + h.c
// equals D_U (Cjg(A).sol) if A and sol taken U indep.
 


// Dplus term
 
sites=0;
while(loop_over_lattice(x,sites)){
for(mu=0;mu<D;mu++){
e_mu=Lattice_Vector(mu);
tmp=Umatrix();


for(i=0;i<KDFERMION;i++){
tmp=tmp-
U.get(x,mu)*s[i+KDFERMION].get(x+e_mu)*Udag.get(x,mu)*pdag[i].get(x)*BC(x,e_mu);
tmp=tmp+
pdag[i].get(x)*U.get(x,mu)*s[i+KDFERMION].get(x+e_mu)*Udag.get(x,mu)*BC(x,e_mu);}
	 
tmp=tmp-Adj(tmp);

f_U.set(x,mu,f_U.get(x,mu)+tmp);
}}

// Dminus term

sites=0;
while(loop_over_lattice(x,sites)){
for(mu=0;mu<D;mu++){
tmp=Umatrix();
e_mu=Lattice_Vector(mu);

for(i=0;i<KDFERMION;i++){
tmp=tmp+
U.get(x,mu)*pdag[i+KDFERMION].get(x+e_mu)*Udag.get(x,mu)*s[i].get(x)*BC(x,e_mu);
tmp=tmp-
s[i].get(x)*U.get(x,mu)*pdag[i+KDFERMION].get(x+e_mu)*Udag.get(x,mu)*BC(x,e_mu);
}

tmp=tmp-Adj(tmp);

f_U.set(x,mu,f_U.get(x,mu)+tmp);
}
}

// SU(N) case

if(RANK==(NCOLOR*NCOLOR-1)){
sites=0;
while(loop_over_lattice(x,sites)){
for(mu=0;mu<D;mu++){
f_U.set(x,mu,f_U.get(x,mu)-(1.0/NCOLOR)*Tr(f_U.get(x,mu))*Umatrix(1));
}
}
}

//Yukawas ...

for(i=0;i<NSCALAR-2;i++){
for(a=0;a<KDFERMION;a++){
for(b=0;b<KDFERMION;b++){
f_phi[i]=f_phi[i]+Gamma[i].get(a,b)*comm(pdag[a],s[b+KDFERMION]);
f_phi[i]=f_phi[i]+Gamma[i].get(b,a)*comm(pdag[a+KDFERMION],s[b]);
}}}

for(a=0;a<KDFERMION;a++){
f_phi[NSCALAR-2]=f_phi[NSCALAR-2]+Complex(-1.0,0.0)*comm(pdag[a],s[a]);
f_phi[NSCALAR-2]=f_phi[NSCALAR-2]+Complex(1.0,0.0)*
comm(pdag[a+KDFERMION],s[a+KDFERMION]);
}

for(a=0;a<NFERMION;a++){
f_phi[NSCALAR-1]=f_phi[NSCALAR-1]+Complex(0.0,1.0)*comm(pdag[a],s[a]);
}

for(i=0;i<NSCALAR;i++){
sites=0;
while(loop_over_lattice(x,sites)){
f_phi[i].set(x,(f_phi[i].get(x)-Adj(f_phi[i].get(x))));
}
}


if(RANK==(NCOLOR*NCOLOR-1)){

for(i=0;i<NSCALAR;i++){
sites=0;
while(loop_over_lattice(x,sites)){
f_phi[i].set(x,f_phi[i].get(x)-(1.0/NCOLOR)*Tr(f_phi[i].get(x))*Umatrix(1));
}
}


}

sites=0;
while(loop_over_lattice(x,sites)){
for(mu=0;mu<D;mu++){
f_U.set(x,mu,NORM*f_U.get(x,mu));
}
}

for(i=0;i<NSCALAR;i++){
f_phi[i]=NORM*f_phi[i];}


return;
}
