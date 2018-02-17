#include "write_out.h"

void write_out(const Gauge_Field &u, const Site_Field phi[],const Site_Field F[]){
static ofstream f_write;
static int first_time=1;
Lattice_Vector x;
int site,j,mu;

f_write.open("dump");
if(f_write.bad()){
cout << "error opening config file\n" << flush;}
f_write << L << "\t" << T << "\t" << BETA << "\t" << DT << "\t" << NCOLOR << "\n";

site=0;
while(loop_over_lattice(x,site)){
for(int mu=0;mu<D;mu++){
f_write << u.get(x,mu) << "\n";
}}
f_write << "\n" << flush;

for(j=0;j<NSCALAR;j++){
site=0;
while(loop_over_lattice(x,site)){
f_write << phi[j].get(x) << "\n";
}
}
f_write << "\n" << flush;

for(j=0;j<NFERMION;j++){
site=0;
while(loop_over_lattice(x,site)){
f_write << F[j].get(x)<< "\n";
}
}


f_write.close();
return;
}

