
#include "update.h"

void update(Gauge_Field &U, Site_Field phi[NSCALAR],
Site_Field F[NFERMION]){

Gauge_Field p_U,old_U,old_f_U;
Site_Field p_phi[NSCALAR],old_phi[NSCALAR],old_f_phi[NSCALAR];
Site_Field p_F[NFERMION],old_F[NFERMION],old_f_F[NFERMION];

double H_old,K_old,K_new,S_new,H_new,hmc_test;

static Gauge_Field f_U;
static Site_Field f_phi[NSCALAR];
static Site_Field f_F[NFERMION];

static double S_old=0.0;
static int first_time=1,accept=0,no_calls=0;
static ofstream f_hmc;

int i;

no_calls++;

// refresh momenta
p_U=Gauge_Field(2);

for(i=0;i<NSCALAR;i++){
p_phi[i]=Site_Field(1);}

for(i=0;i<NFERMION;i++){
p_F[i]=Site_Field(2);}


K_old=kinetic_energy(p_U,p_phi,p_F);

//cout << "kold  is " << K_old << "\n";


if(first_time){
f_hmc.open("hmc_test");
if(f_hmc.bad()){
cout << "failed to open hmc_test file\n" << flush;}

S_old=action(U,phi,F);
//cout << "computed action\n" << flush;

force(U,f_U,phi,f_phi,F,f_F);
first_time=0;
}	

//cout << "sold is " << S_old << "\n";

if((no_calls%100==0)&&(!first_time)){
cout << "acceptance rate " << (double)accept/(double)no_calls << "\n" <<
flush;
no_calls=0;
accept=0;
}


 

//save copies of fields
old_U=U;
for(i=0;i<NSCALAR;i++){
old_phi[i]=phi[i];}
for(i=0;i<NFERMION;i++){
old_F[i]=F[i];}

H_old=S_old+K_old;

old_f_U=f_U;
for(i=0;i<NSCALAR;i++){
old_f_phi[i]=f_phi[i];}
for(i=0;i<NFERMION;i++){
old_f_F[i]=f_F[i];}

cout << "-------------*****--------------" << endl ;
cout << "H at beginning traj " << H_old << "\n" << flush;

//classical evolution

for(int i=0;i<TRAJECTORY_LENGTH;i++){
evolve_fields(U,p_U,f_U,phi,p_phi,f_phi,F,p_F,f_F);
}

S_new=action(U,phi,F);
K_new=kinetic_energy(p_U,p_phi,p_F);

H_new=S_new+K_new;

cout << "H at end of traj " << H_new << "\n" << flush;

hmc_test=exp(-H_new+H_old);
f_hmc << hmc_test << "\n" << flush;
      	    
//metropolis test
	if(myrandom()<exp(H_old-H_new)){
		S_old=S_new;
		accept++;
        cout << "ACCEPTED with deltaH of " << H_new - H_old << endl ;
		return;
	}
	else{
	cout << "hmc_test " << hmc_test << " failed\n" << flush;
    cout << "REJECTED with deltaH of " << H_new - H_old << endl ;
	// if fails copy back fields
	U=old_U;
        for(i=0;i<NSCALAR;i++){
	phi[i]=old_phi[i];
        f_phi[i]=old_f_phi[i];}
        for(i=0;i<NFERMION;i++){
	F[i]=old_F[i];
        f_F[i]=old_f_F[i];}
	f_U=old_f_U;
	return;
	}	   



return;	    
}	    
	    
	    
	    
	    
	    
