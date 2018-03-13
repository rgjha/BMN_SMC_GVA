#ifndef UTILITIES_H
#define UTILITIES_H

#include <iostream>
#include <fstream>
using namespace std;

#include <math.h>
#include <stdlib.h>

// 16-supercharge SYM in 1 dim
// Compile with g++ -O3 *.cpp -llapack -lblas -lg2c
// Reads parameters from file parameters

const int FERMIONS = 0;
const int L = 1;
const int T = 4;
const int D = 1;
const int NSCALAR = 9;
const int NFERMION = 16;
const int KDFERMION = 8;
const int SITES = L * T;
const int NCOLOR = 2;
const int RANK = NCOLOR * NCOLOR - 1;
const int FERMIONSIZE = RANK * SITES * NFERMION;
const double GAUGETOL = 0.00000000000001;
const double TRACETOL = 1e-8;       // Will be used to test norm()
const int DEGREE = 15;
const double NORM = 0.5 / sqrt(D);
const double PBC = 1.0;

extern double ampdeg,amp[DEGREE],shift[DEGREE];
extern double BETA, DT, MASS, TIME, MU;
extern int SWEEPS, GAP, START, THERM, READIN, SEED;
extern double f[RANK][RANK][RANK];
extern double SMALLCUT, LARGECUT;
extern int TRAJECTORY_LENGTH;
extern int epsilon[3][3][3];

class Complex{
  private:
    double re, im;
  public:
    Complex();
    Complex(double, double);
    double real(void) const;
    double imag(void) const;
    double norm(void);
    void print(void) const;
    friend ostream& operator<<(ostream&,Complex);
    friend istream& operator>>(istream&,Complex &);};

inline Complex conjug(const Complex &o1){return(Complex(o1.real(),-o1.imag()));}
inline Complex operator +(const Complex &o1, const Complex &o2){
  return(Complex(o1.real()+o2.real(),o1.imag()+o2.imag()));}
inline Complex operator -(const Complex &o1, const Complex &o2){
  return(Complex(o1.real()-o2.real(),o1.imag()-o2.imag()));}
inline Complex operator *(const Complex &o1, const Complex &o2){
  return(Complex(o1.real()*o2.real()-o1.imag()*o2.imag(),
    o1.real()*o2.imag()+o1.imag()*o2.real()));}
inline Complex operator *(const Complex &o1, const double o2){
  return(Complex(o1.real()*o2,o1.imag()*o2));}
inline Complex operator *(const double o1, const Complex &o2){
  return(Complex(o2.real()*o1,o2.imag()*o1));}

Complex operator /(const Complex &, const Complex &);
Complex pow(const Complex &, const int);

class Umatrix{
  private:
    Complex mat[NCOLOR][NCOLOR];
  public:
    Umatrix();
    Umatrix(int);
    Umatrix(Complex [NCOLOR][NCOLOR]);
    Complex get(int,int) const;
    void set(int,int,const Complex);
    void print(void);
    friend ostream& operator<<(ostream &, Umatrix);
    friend istream& operator>>(istream &, Umatrix &);};

Umatrix operator +(const Umatrix &o1, const Umatrix &o2);
Umatrix operator -(const Umatrix &o1, const Umatrix &o2);
Umatrix operator *(const Umatrix &, const Umatrix &);
Umatrix operator *(const Umatrix &, const Complex &);
Umatrix operator *(const Complex &, const Umatrix &);
Umatrix operator *(const Umatrix &, const double);
Umatrix operator *(const double, const Umatrix &);
Umatrix comm(const Umatrix &, const Umatrix &);
Umatrix exp(const Umatrix &u);
Umatrix Adj(const Umatrix &u);
Complex Tr(const Umatrix &);
Umatrix real_gaussian_Umatrix(void);

class Gamma_Matrix{
  private:
    double gam[KDFERMION][KDFERMION];
  public:
    Gamma_Matrix();
    Gamma_Matrix(int);
    double get(int,int) const;
    void set(int,int,const double);
    void print(void);};

Gamma_Matrix operator *(const Gamma_Matrix &, const Gamma_Matrix &);
Gamma_Matrix operator *(const double k, const Gamma_Matrix &o);
Gamma_Matrix operator +(const Gamma_Matrix &, const Gamma_Matrix &);

class Lattice_Vector{
private:
  int coords[D];
public:
  Lattice_Vector(void);
  Lattice_Vector(int);
  void set(int, int);
  int get(int) const;
  void print(void) const;
};

Lattice_Vector operator +(const Lattice_Vector &x, const Lattice_Vector &y);
Lattice_Vector operator -(const Lattice_Vector &x, const Lattice_Vector &y);
Lattice_Vector operator -(const Lattice_Vector &x);
double BC(const Lattice_Vector &x, const Lattice_Vector &y);
int loop_over_lattice(Lattice_Vector &, int &);


class Gauge_Field{
private:
  Umatrix link[SITES][D];

public:
        Gauge_Field(void);
  Gauge_Field(int);
  Umatrix get(const Lattice_Vector &, const int) const;
  void set(const Lattice_Vector &, const int, const Umatrix &);
  };

Gauge_Field Adj(const Gauge_Field &);

class Site_Field{
private:
  Umatrix points[SITES];
public:
  Site_Field(void);
  Site_Field(int);
  Umatrix get(const Lattice_Vector &) const;
  void set(const Lattice_Vector &, const Umatrix &);
  void print(void);
  };

Site_Field Adj(const Site_Field &);

Site_Field operator +(const Site_Field &, const Site_Field &);
Site_Field operator -(const Site_Field &, const Site_Field &);
Site_Field operator *(const double, const Site_Field &);
Site_Field operator *(const Complex &, const Site_Field &);
Umatrix operator *(const Site_Field &, const Site_Field &);
Site_Field comm(const Site_Field &, const Site_Field &);

Site_Field Dplus(const Gauge_Field &, const Site_Field &);
Site_Field Dminus(const Gauge_Field &, const Site_Field &);

void Build_Scalar(const Site_Field phi[NSCALAR], Site_Field Scalar[KDFERMION][KDFERMION]);
void Fermion_operator(const Gauge_Field &,
const Site_Field Scalar[KDFERMION][KDFERMION],
const Site_Field phi[NSCALAR],
const Site_Field K[NFERMION], Site_Field K2[NFERMION],const int);

extern Umatrix Lambda[RANK];
extern Gamma_Matrix Gamma[NSCALAR-2],Gam123;

double gasdev(void);
#endif
