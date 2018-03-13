#include "read_param.h"

void read_param(void){
double LAMBDA;
ifstream f_in("parameters");
if(!f_in.good()) {
  cout << "\ncan't open file parameters to read data!\n";
  exit(1);
}
f_in>>SWEEPS>>THERM>>GAP>>LAMBDA>>MU>>DT>>READIN>>SEED;

BETA=(NCOLOR*0.25)/LAMBDA; 

TRAJECTORY_LENGTH=(int)(1/DT);

if(FERMIONS){cout << "16 supercharge theory in 1D\n" ;}
else
{cout << "0 supercharge theory in 1D\n";}

cout << "--------------------------------" << endl ;
cout << "Number of colors " << NCOLOR <<  "\n";
cout << "Temporal extent " << T << "\n";
cout << "Inverse temperature " << LAMBDA << "\n";
cout << "Mass parameter " << MU << "\n";
cout << "Lattice Coupling " << BETA << "\n";
cout << "Thermalization sweeps " << THERM << "\n";
cout << "Number of sweeps " << SWEEPS << "\n";
cout << "Gap between measurements " << GAP << "\n";
cout << "Time step in leapfrog eqs " << DT << "\n";
cout << "Trajectory length " << TRAJECTORY_LENGTH << "\n";
cout << "Minimax approx degree " << DEGREE << "\n";
cout << "Reading initial config: (1 for yes, 0 for no) " << READIN << "\n";

if (PBC==1.0) {cout << "periodic temporal bc for fermions" << "\n";}
else{cout << "antiperiodic temporal bc for fermions" << "\n";}

cout << "--------------------------------" << endl ;

setup();

return;
}
