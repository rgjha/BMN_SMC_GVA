 #include "action.h"

double action(const Gauge_Field &U, const Site_Field phi[NSCALAR], const 
Site_Field F[NFERMION]){
Lattice_Vector x,e_mu;
Gauge_Field Udag;

double act_s=0.0, act_F=0.0;
int mu,nu,site,i,j,k; 
Site_Field sol[NFERMION][DEGREE],psol[NFERMION][DEGREE];
  
Udag=Adj(U);
  
//scalar kinetic term
  
for(i=0;i<NSCALAR;i++){     // NSCALAR = 9, NFERMION = 16
  
site=0;
while(loop_over_lattice(x,site)){				
for(mu=0;mu<D;mu++){	       // D is = 1
e_mu=Lattice_Vector(mu);

act_s=act_s-BETA*Tr(phi[i].get(x)*(
U.get(x,mu)*phi[i].get(x+e_mu)*Udag.get(x,mu)+                                  
Udag.get(x-e_mu,mu)*phi[i].get(x-e_mu)*U.get(x-e_mu,mu)-
2.0*phi[i].get(x))).real();
    
// Tr [φ_i(x) { U_μ(x) * φ_i(x + e_μ) * Udag_μ(x) + Udag_μ(x - e_μ) * φ_i(x - e_μ) * Udag_μ(x - e_μ) - 2φ_i(x) } ]
 
}
}
   
}
// scalar-scalar potential

for(i=0;i<NSCALAR;i++){
for(j=i+1;j<NSCALAR;j++){

site=0;
while(loop_over_lattice(x,site)){
act_s=act_s-BETA*Tr(comm(phi[i].get(x),phi[j].get(x))*
                    comm(phi[i].get(x),phi[j].get(x))).real();
// - Tr [φ_i,φ_j]^2
}
}
}

// BMN  mass terms

for(i=0;i<3;i++){
site=0;
while(loop_over_lattice(x,site)){
act_s=act_s-BETA*MU*MU*Tr(phi[i].get(x)*phi[i].get(x)).real();
// - μ^2 * Tr [φ_i * φ_i]   ; i = 0,1,2
}
}

for(i=3;i<NSCALAR;i++){
site=0;
while(loop_over_lattice(x,site)){
act_s=act_s-BETA*0.25*MU*MU*Tr(phi[i].get(x)*phi[i].get(x)).real();
// - 0.25 * μ^2 * Tr [φ_i * φ_i]   ; i = 3,....8
}}

// cubic term

for(i=0;i<3;i++){
for(j=0;j<3;j++){
for(k=0;k<3;k++){
if((i==j) || (i==k) || (j==k)) continue;
site=0;
while(loop_over_lattice(x,site)){
act_s=act_s-BETA*2*sqrt(2.0)*MU*epsilon[i][j][k]*Tr(phi[i].get(x)*phi[j].get(x)*phi[k].get(x)).real();
// - 2 √ 2 μ ε_ijk * φ_i * φ_j * φ_k
}}}}
    

//pseudofermion contribution
if(FERMIONS){

for(i=0;i<NFERMION;i++){
act_F=act_F+ampdeg*Tr(Adj(F[i])*F[i]).real();}


MCG_solver(U,phi,F,shift,sol,psol);

for(int n=0;n<DEGREE;n++){
for(i=0;i<NFERMION;i++){
act_F=act_F+amp[n]*Tr(Adj(F[i])*sol[i][n]).real();}

}

}
return(act_s+act_F);
}
