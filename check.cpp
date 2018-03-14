#include "check.h"

void check(Gauge_Field &U){
Lattice_Vector x;
int mu,site;
double dummy;

site=0;
while(loop_over_lattice(x,site)){
for(mu=0;mu<D;mu++){

dummy=fabs((1.0/NCOLOR)*Tr(U.get(x,mu)*Adj(U.get(x,mu))).real()-1.0);

if(dummy>GAUGETOL)
{cout << "U field not unitary - 1/NTr(U*Adj(U))-1 is " << dummy << "\n";
dummy=sqrt(1.0/NCOLOR*(Tr(U.get(x,mu)*Adj(U.get(x,mu))).real()));
U.set(x,mu,U.get(x,mu)*(1.0/dummy)); }

}
}

return;
}
