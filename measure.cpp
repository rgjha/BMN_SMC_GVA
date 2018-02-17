#include "measure.h"  

void measure(const Gauge_Field &U, const Site_Field phi[],
const Site_Field F[NFERMION]){

	static int first_time=1;
	static ofstream f_data,f_av,f_eigen,f_pfaff;
        double act_s,act_F,MAX,MIN;
        double eigenvals[T][NSCALAR][NCOLOR];
	int sites,i,j,k;
	Lattice_Vector x;
	Site_Field G,phitmp[NSCALAR];
	Gauge_Field Utmp;
        Complex dum;
	
	if(first_time){
	f_data.open("data",ios::app);
	if(f_data.bad()){ 
	cout << "failed to open action file\n" << flush ;}
	
        f_av.open("scalar",ios::app);
        if(f_av.bad()){
        cout << "failed to open scalar eigenvalues file\n" << flush;}
	
	f_eigen.open("eigen",ios::app);
	if(f_eigen.bad()){ 
	cout << "failed to open eigenvalues file\n" << flush ;}

	f_pfaff.open("pfaffian",ios::app);
	if(f_pfaff.bad()){ 
	cout << "failed to open pfaffian file\n" << flush ;}
	first_time=0;
	}
	  
        line(U);	   
	obs(U,phi,F,act_s,act_F,eigenvals);
        f_data << act_s/(4.50*NCOLOR*NCOLOR*SITES) <<  "\t" << act_F/(16.0*NCOLOR*NCOLOR*SITES) << "\n" << flush;

       if(NCOLOR<5){
        Complex M[FERMIONSIZE][FERMIONSIZE];
	full_fermion_op(U,phi,M);
	eigenvalues(M,MIN,MAX);
	dum=Pfaffian(M);
	f_pfaff << log(dum.norm()) << "\t" << dum.real()/dum.norm() << "\t" << dum.imag()/dum.norm() << "\n"
	<< flush;
	
	f_eigen << MIN << "\t" << MAX << "\n" << flush;

	if(MAX>LARGECUT){cout << "outside max eigenvalue " << MAX << "\n";}
        if(MIN<SMALLCUT){cout << "outside min eigenvalue " << MIN << "\n";}
        }
	
	for(i=0;i<T;i++){
	for(j=0;j<NSCALAR;j++){
	for(k=0;k<NCOLOR;k++){
        f_av << eigenvals[i][j][k] << "\t" ;}
	f_av << "\n";}}
	f_av << flush;

      
	return;
}
