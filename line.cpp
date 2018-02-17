#include "line.h"

// computes Polyakov line
void line(const Gauge_Field  &U){
int site=0,t;
Lattice_Vector x,y,e_mu;
Umatrix prod;
static int first_time=1,count=0;
static ofstream f_line;
Complex poly=Complex();

count++;

if(first_time){
f_line.open("lines",ios::app);
if(f_line.bad()){
cerr << "failed to open lines file" << "\n";exit(1);}

first_time=0;
}


site=0;
while(loop_over_lattice(x,site)){

e_mu=Lattice_Vector(D-1);

prod=Umatrix(1);

y=x;
for(t=1;t<=T;t++){
prod=prod*U.get(y,D-1);
y=y+e_mu;
}

poly=poly+(1.0/NCOLOR)*Tr(prod);
}

 
f_line  << poly*(1.0/SITES) << "\n" << flush;


return;
}
