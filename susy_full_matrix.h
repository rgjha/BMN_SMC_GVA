#include "utilities.h"

int lat_pack(const Lattice_Vector &x);
void lat_unpack(const int n, Lattice_Vector &x);
int pack(const int flavor,const int type, const Lattice_Vector &x, const int color);
void unpack(const int i, int &flavor,int &type, Lattice_Vector &x, int &color);

class Adjoint_Matrix{
private:
       Complex amat[RANK][RANK];
public:
       Adjoint_Matrix(void);
       Adjoint_Matrix(Complex m[RANK][RANK]);
       Complex get(const int, const int) const;
       void set(const int, const int, const Complex &);
       friend ostream& operator<<(ostream &, Adjoint_Matrix);
       friend istream& operator>>(istream &, Adjoint_Matrix &);
       };
       

class Adjoint_Links{
private:
       Adjoint_Matrix alinks[SITES][D];
       
public:
       Adjoint_Links(void);
       Adjoint_Matrix get(const Lattice_Vector &, const int) const;
       void set(const Lattice_Vector &, const int, const Adjoint_Matrix &);
       void print(void);
       };

class Adjoint_Sites{
private: 
	Adjoint_Matrix asites[SITES];
	
public:
	Adjoint_Sites(void);
	Adjoint_Matrix get(const Lattice_Vector &) const;
	void set(const Lattice_Vector &, const Adjoint_Matrix &);
	void print(void);
	};

Adjoint_Sites operator*(const double, const Adjoint_Sites &);
	
void compute_Adjoints(const Gauge_Field &, Adjoint_Links &, const Site_Field
phi[], Adjoint_Sites aphi[]);

void full_fermion_op(const Gauge_Field &U, const Site_Field phi[],Complex
M[FERMIONSIZE][FERMIONSIZE]);

double cnorm(const Complex &c);

Complex Pfaffian(Complex M[FERMIONSIZE][FERMIONSIZE]);
