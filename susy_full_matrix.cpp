#include "susy_full_matrix.h"

int lat_pack(const Lattice_Vector &x){
int i,site=0;
static int Lattice_Map[D];
static int first_time=1;
if(first_time){
for(i=0;i<D;i++)
Lattice_Map[i]=(int)pow((double)L,(double)i);
first_time=0;
}
for(i=0;i<D;i++){
site=site+Lattice_Map[i]*x.get(i);
}

return(site);
}

void lat_unpack(const int n, Lattice_Vector &x){
static int Lattice_Map[D];
static int first_time=1;
int current,i;

if(first_time){
for(i=0;i<D;i++)
Lattice_Map[i]=(int)pow((double)L,(double)i);
first_time=0;
}

current=n;
for(i=D-1;i>=0;i--){
x.set(i,current/Lattice_Map[i]);
current=current-Lattice_Map[i]*x.get(i);
}

return;
}

int pack(const int flavor, const Lattice_Vector &x, const int color){
int i;
i=(SITES*RANK)*flavor+RANK*lat_pack(x)+color;

return(i);
}

void unpack(const int i, int &flavor, Lattice_Vector &x, int &color){
int tmp;
flavor=i/(SITES*RANK);
tmp=i-flavor*(SITES*RANK);
lat_unpack(tmp/RANK,x);
color=tmp-(tmp/RANK)*RANK;
return;
}

Adjoint_Matrix::Adjoint_Matrix(void){
for(int i=0;i<RANK;i++)
for(int j=0;j<RANK;j++){
amat[i][j]=Complex();}
return;
}

Complex Adjoint_Matrix::get(const int i, const int j) const{
return(amat[i][j]);}

void Adjoint_Matrix::set(const int i, const int j, const Complex &c){
amat[i][j]=c;
return;
}

Adjoint_Matrix::Adjoint_Matrix(Complex m[RANK][RANK]){
for(int i=0;i<RANK;i++){
for(int j=0;j<RANK;j++){
amat[i][j]= m[i][j];
}}
return;
}

ostream& operator<<(ostream& out, Adjoint_Matrix s){
for(int i=0;i<RANK;i++){
for(int j=0;j<RANK;j++){
	out<<s.get(i,j)<<'\t';}}
	return out;}
			
istream& operator>>(istream& in, Adjoint_Matrix & s){
	Complex v[RANK][RANK];
	for(int j=0;j<RANK;j++){
	for(int i=0;i<RANK;i++){
	in>>v[j][i];}}
	s=Adjoint_Matrix(v);
	return in;
	}		


Adjoint_Links::Adjoint_Links(void){
for(int i=0;i<SITES;i++)
for(int j=0;j<D;j++){
alinks[i][j]=Adjoint_Matrix();
}
return;
} 

Adjoint_Matrix Adjoint_Links::get(const Lattice_Vector &x, const int mu) const{
int site=0,i;
static int first_time=1;
static int Lattice_Map[D];

if(first_time){
for(i=0;i<D;i++)
Lattice_Map[i]=(int)pow((double)L,(double)i);
first_time=0;
}

for(i=0;i<D;i++)
site=site+x.get(i)*Lattice_Map[i];

return(alinks[site][mu]);
}


void Adjoint_Links::set(const Lattice_Vector &x, const int mu,
const Adjoint_Matrix &u){
int site=0,i;
static int first_time=1;
static int Lattice_Map[D];

if(first_time){
for(i=0;i<D;i++)
Lattice_Map[i]=(int)pow((double)L,(double)i);
first_time=0;
}

for(i=0;i<D;i++)
site=site+x.get(i)*Lattice_Map[i];

alinks[site][mu]=u;
return;
}

void Adjoint_Links::print(void){
cout << "adjoint links values\n" << flush;
for(int i=0;i<SITES;i++){
cout << "site= " << i << "\n" << flush;
for(int j=0;j<D;j++){
cout << "j= " << j << ":" << alinks[i][j]<< "\n" << flush;}}
return;
}


Adjoint_Sites::Adjoint_Sites(void){
for(int i=0;i<SITES;i++){
asites[i]=Adjoint_Matrix();
}
return;
} 

Adjoint_Matrix Adjoint_Sites::get(const Lattice_Vector &x) const{
int site=0,i;
static int first_time=1;
static int Lattice_Map[D];

if(first_time){
for(i=0;i<D;i++)
Lattice_Map[i]=(int)pow((double)L,(double)i);
first_time=0;
}

for(i=0;i<D;i++)
site=site+x.get(i)*Lattice_Map[i];

return(asites[site]);
}


void Adjoint_Sites::set(const Lattice_Vector &x,
const Adjoint_Matrix &u){
int site=0,i;
static int first_time=1;
static int Lattice_Map[D];

if(first_time){
for(i=0;i<D;i++)
Lattice_Map[i]=(int)pow((double)L,(double)i);
first_time=0;
}

for(i=0;i<D;i++)
site=site+x.get(i)*Lattice_Map[i];

asites[site]=u;
return;
}

void Adjoint_Sites::print(void){
cout << "adjoint links values\n" << flush;
for(int i=0;i<SITES;i++){
cout << "site= " << i << "\n" << flush;}
return;
}

Adjoint_Sites operator*(const double d, const Adjoint_Sites &A){
Adjoint_Sites dum;
Adjoint_Matrix aa;
int sites=0;
Lattice_Vector x;
Complex temp;

while(loop_over_lattice(x,sites)){
for(int a=0;a<RANK;a++){
for(int b=0;b<RANK;b++){
temp=d*A.get(x).get(a,b);
aa.set(a,b,temp);
}
}
dum.set(x,aa);
}

return(dum);
}

void compute_Adjoints(const Gauge_Field &U, Adjoint_Links &V, const Site_Field
phi[], Adjoint_Sites aphi[]){

int mu,a,b,sites;
Lattice_Vector x;
Adjoint_Matrix tmp;
Umatrix dum;
//cout << "in compute_Adjoint_Links\n" << flush;

sites=0;
while(loop_over_lattice(x,sites)){
for(mu=0;mu<D;mu++){

for(a=0;a<RANK;a++){
for(b=0;b<RANK;b++){
tmp.set(a,b,Tr(Lambda[a]*U.get(x,mu)*Lambda[b]*Adj(U.get(x,mu))));
}
}
V.set(x,mu,tmp);
}
}

for(int i=0;i<NSCALAR;i++){
sites=0;
while(loop_over_lattice(x,sites)){
for(a=0;a<RANK;a++){
for(b=0;b<RANK;b++){
tmp.set(a,b,Tr(Lambda[a]*phi[i].get(x)*Lambda[b])-Tr(Lambda[b]*phi[i].get(x)*Lambda[a]));
}}
aphi[i].set(x,tmp);
}
}

return;
}

Complex delta(const int a, const int b){
if(a==b) return(Complex(1.0,0.0));
else
return(Complex(0.0,0.0));
}

void full_fermion_op(const Gauge_Field &U, const Site_Field phi[], Complex
M[FERMIONSIZE][FERMIONSIZE]){
int sites,a,b,c,j1,j2,i,mu,nu,f;
Lattice_Vector x,e_mu,e_nu;
Adjoint_Links V;
Adjoint_Sites aphi[NSCALAR];

for(j1=0;j1<FERMIONSIZE;j1++){
for(j2=0;j2<FERMIONSIZE;j2++){
M[j1][j2]=Complex();
}}

// compute adjoint links ....


compute_Adjoints(U, V,phi,aphi);

// kinetic terms 

for(f=0;f<KDFERMION;f++){
sites=0;
while(loop_over_lattice(x,sites)){
for(a=0;a<RANK;a++){
i=pack(f,x,a);

for(mu=0;mu<D;mu++){
e_mu=Lattice_Vector(mu);
for(b=0;b<RANK;b++){

j1=pack(f+KDFERMION,x,b);
j2=pack(f+KDFERMION,x+e_mu,b);
M[i][j1]=M[i][j1]+delta(a,b);
M[i][j2]=M[i][j2]+V.get(x,mu).get(a,b)*BC(x,e_mu);
}
}

}
}
}


for(f=0;f<KDFERMION;f++){
sites=0;
while(loop_over_lattice(x,sites)){
for(mu=0;mu<D;mu++){
e_mu=Lattice_Vector(mu);

for(a=0;a<RANK;a++){
i=pack(f+KDFERMION,x,a);

for(b=0;b<RANK;b++){

j1=pack(f,x,b);
j2=pack(f,x-e_mu,b);
M[i][j2]=M[i][j2]-V.get(x-e_mu,mu).get(b,a)*BC(x,-e_mu);
M[i][j1]=M[i][j1]-delta(a,b);
}}

}
}
}
// BMN mass

for(a=0;a<RANK;a++){
sites=0;
while(loop_over_lattice(x,sites)){

for(f=0;f<KDFERMION;f++){
i=pack(f,x,a);
for(b=0;b<KDFERMION;b++){

j1=pack(b+KDFERMION,x,a);
M[i][j1]=M[i][j1]+Gam123.get(f,b)*MU*3.0/4.0*Complex(1.0,0.0);

}}}}

for(a=0;a<RANK;a++){
sites=0;
while(loop_over_lattice(x,sites)){

for(f=0;f<KDFERMION;f++){
i=pack(f+KDFERMION,x,a);
for(b=0;b<KDFERMION;b++){

j1=pack(b,x,a);
M[i][j1]=M[i][j1]-Gam123.get(b,f)*MU*3.0/4.0*Complex(1.0,0.0);

}}}}

// Yukawas

//cout << "computing Yukawas\n" << flush;

sites=0;
while(loop_over_lattice(x,sites)){

for(f=0;f<KDFERMION;f++){
for(a=0;a<RANK;a++){
for(c=0;c<KDFERMION;c++){
for(b=0;b<RANK;b++){

j1=pack(f,x,a);
j2=pack(c+KDFERMION,x,b);

for(i=0;i<NSCALAR-2;i++){
M[j1][j2]=M[j1][j2]+Gamma[i].get(f,c)*aphi[i].get(x).get(a,b);
}

j1=pack(f+KDFERMION,x,a);
j2=pack(c,x,b);

for(i=0;i<NSCALAR-2;i++){
M[j1][j2]=M[j1][j2]+Gamma[i].get(c,f)*aphi[i].get(x).get(a,b);
}

j1=pack(f,x,a);
j2=pack(c,x,b);

M[j1][j2]=M[j1][j2]+delta(c,f)*Complex(-1.0,0.0)*aphi[NSCALAR-2].get(x).get(a,b);
M[j1][j2]=M[j1][j2]+delta(c,f)*Complex(0.0,1.0)*aphi[NSCALAR-1].get(x).get(a,b);

j1=pack(f+KDFERMION,x,a);
j2=pack(c+KDFERMION,x,b);

M[j1][j2]=M[j1][j2]+delta(c,f)*Complex(1.0,0.0)*aphi[NSCALAR-2].get(x).get(a,b);
M[j1][j2]=M[j1][j2]+delta(c,f)*Complex(0.0,1.0)*aphi[NSCALAR-1].get(x).get(a,b);

}}}}
}

//test AH nature
for(j1=0;j1<FERMIONSIZE;j1++){
for(j2=j1+1;j2<FERMIONSIZE;j2++){
//cout << "M[j1][j2] is " << M[j1][j2] << "\n" <<flush;
Complex dum=M[j1][j2]+M[j2][j1];
if(dum.norm()>0.00000001) {cout << "problem " << j1 << j2 << "\n" << flush;}
}
}


return;
}


double cnorm(const Complex &c){
double dum;
dum=c.real()*c.real()+c.imag()*c.imag();
return(sqrt(dum));
}

Complex Pfaffian(Complex M[FERMIONSIZE][FERMIONSIZE])
{
	int i,j,k,jpiv,interchange=1;
        static int firsttime=1;
	double pivot,totalangle,angle,cosine,sine,mag;
	Complex dum[FERMIONSIZE],scale,f;
//        cout << "in Pfaffian\n" << flush;
        static ofstream f_pfaff;
       
        if(firsttime){ 
        f_pfaff.open("pfaffian",ios::app);
        if(f_pfaff.bad()){
        cout << "failed to open pfaffian file\n" << flush ;}
        firsttime=0;
}
// loop over all rows in steps of 2
	for(i=0;i<FERMIONSIZE-2;i+=2){

// first row i:	
// find col whose ith component is biggest to use as pivot
	pivot=cnorm(M[i][i+1]);
	jpiv=i+1;
	for(j=i+2;j<FERMIONSIZE;j++){
	if(cnorm(M[i][j])>pivot){pivot=cnorm(M[i][j]);jpiv=j;}}

// interchange col(i+1) with col(jpiv)
        for(j=0;j<FERMIONSIZE;j++){
	dum[j]=M[j][i+1];}
	for(j=0;j<FERMIONSIZE;j++){
	M[j][i+1]=M[j][jpiv];
	M[j][jpiv]=dum[j];}
// interchange row(i+1) with row(jpiv)
	for(j=0;j<FERMIONSIZE;j++){
	dum[j]=M[i+1][j];}
	for(j=0;j<FERMIONSIZE;j++){
	M[i+1][j]=M[jpiv][j];
	M[jpiv][j]=dum[j];}
	
	if(jpiv!=i+1) interchange*=(-1);
	
// using this zero progressively elements of row M[i][j], j=i+2...FERMIONSIZE-1
        for(j=i+2;j<FERMIONSIZE;j++){
	scale=M[i][j]/M[i][i+1];
	
	for(k=0;k<FERMIONSIZE;k++){
	M[k][j]=M[k][j]-scale*M[k][i+1];}
// zero out elements along corresponding column M[j][i] too
	for(k=0;k<FERMIONSIZE;k++){
	M[j][k]=M[j][k]-scale*M[i+1][k];}
	}

	
// next row i+1;
        
	// using this zero progressively elements M[i][j], j=i+2...FERMIONSIZE-1
        for(j=i+2;j<FERMIONSIZE;j++){
	scale=M[i+1][j]/M[i+1][i];
	
	for(k=0;k<FERMIONSIZE;k++){
	M[k][j]=M[k][j]-scale*M[k][i];}
// zero out elements along corresponding column too
	for(k=0;k<FERMIONSIZE;k++){
	M[j][k]=M[j][k]-scale*M[i][k];}
	}
	
	}
	
	
	f=Complex(1.0,0.0);
        mag=0.0;
	totalangle=0.0;
	for(i=0;i<FERMIONSIZE;i+=2){
        cosine=M[i][i+1].real()/M[i][i+1].norm();
        sine=M[i][i+1].imag()/M[i][i+1].norm();
        if ((cosine>0.0) && (sine>0.0)) angle=acos(cosine);
        if ((cosine<0.0) && (sine>0.0)) angle=acos(cosine);
        if ((cosine>0.0) && (sine<0.0)) angle=4*acos(0.0)-1.0*acos(cosine);
        if ((cosine<0.0) && (sine<0.0)) angle=4*acos(0.0)-1.0*acos(cosine);
        mag+=log(M[i][i+1].norm());
        totalangle+=angle;
	f=f*M[i][i+1];
	}

        f_pfaff << mag << "\t" << cos(totalangle)*interchange << "\t" 
                << sin(totalangle)*interchange << "\n";
return (f*interchange);
}

