#include "utilities.h"

// Complex type
Complex::Complex(void) {re = 0.0; im = 0.0;}
Complex::Complex(double x, double y) {re = x; im = y;}
double Complex::real(void) const{return(re);}
double Complex::imag(void) const{return(im);}
double Complex::norm(void) {return(sqrt(re*re+im*im));}
void Complex::print(void) const {cout << "("<< re << ", " << im << ")";}

ostream& operator<<(ostream& out,Complex c) {
  out<<c.real()<<"\t"<<c.imag();
  return out;
}
istream& operator>>(istream& in,Complex & c) {
  double x,y;
  in>>x>>y;
  c=Complex(x,y);
  return in;
}

Complex operator /(const Complex &o1, const Complex &o2) {
  Complex dum;
  double norm;
  norm=o2.real()*o2.real()+o2.imag()*o2.imag();
  dum=Complex((o1.real()*o2.real()+o1.imag()*o2.imag())/norm,
      (o1.imag()*o2.real()-o1.real()*o2.imag())/norm);
  return(dum);
}

Complex pow(const Complex &o1, const int o2) {
  Complex c(1,0);
  for (int i=0;i<o2;i++) c=c*o1;
  return c;
}

// unitary matrix type Umatrix
Umatrix::Umatrix(void) {
  for (int i=0;i<NCOLOR;i++) {
    for (int j=0;j<NCOLOR;j++)
      mat[i][j]=Complex();
  }
}

Umatrix::Umatrix(int k) {
  if (k==1) {
    for (int i=0;i<NCOLOR;i++) {
      for (int j=0;j<NCOLOR;j++) {mat[i][j]=Complex();}
    }
    for (int n=0;n<NCOLOR;n++) {
      mat[n][n]=Complex(1.0,0.0);}
  }
  else {cout << "wrong Umatrix constructor\n" << flush; }

}

Umatrix::Umatrix(Complex m[NCOLOR][NCOLOR]) {
  for (int i=0;i<NCOLOR;i++)
    for (int j=0;j<NCOLOR;j++) mat[i][j]=m[i][j];}


    Complex Umatrix::get(int i, int j) const {return(mat[i][j]);}
    void Umatrix::set(int i, int j, const Complex o) {mat[i][j]=o;}
    void Umatrix::print(void) {
      for (int i=0;i<NCOLOR;i++) {
        for (int j=0;j<NCOLOR;j++) {mat[i][j].print();cout << "\t";}
        cout << "\n";}
    }
Umatrix Adj(const Umatrix &u) {
  Umatrix res;
  for (int i=0;i<NCOLOR;i++)
    for (int j=0;j<NCOLOR;j++) {
      res.set(i,j,conjug(u.get(j,i)));}
  return(res);}             

  ostream& operator<<(ostream& out,Umatrix s) {
    for (int i=0;i<NCOLOR;i++)
      for (int j=0;j<NCOLOR;j++) {
        out<<s.get(i,j)<<'\t';}
    return out;}             
    istream& operator>>(istream& in, Umatrix & s) {
      Complex v[NCOLOR][NCOLOR];
      for (int i=0;i<NCOLOR;i++)
        for (int j=0;j<NCOLOR;j++) {
          in>>v[i][j];}
      s=Umatrix(v);
      return in;
    }     


Umatrix operator *(const Umatrix &o1, const Umatrix &o2) {
  Umatrix r;
  Complex dum;
  for (int i=0;i<NCOLOR;i++)
    for (int j=0;j<NCOLOR;j++) {
      dum=Complex();
      for (int k=0;k<NCOLOR;k++) {
        dum=dum+o1.get(i,k)*o2.get(k,j);}
      r.set(i,j,dum);}
  return(r);}      

  Umatrix operator *(const Umatrix &o1, const Complex &o2) {
    Umatrix dum;
    for (int i=0;i<NCOLOR;i++)
      for (int j=0;j<NCOLOR;j++) {
        dum.set(i,j,o1.get(i,j)*o2);}
    return(dum);}
    Umatrix operator *(const Complex &o2, const Umatrix &o1) {
      Umatrix dum;
      for (int i=0;i<NCOLOR;i++)
        for (int j=0;j<NCOLOR;j++) {
          dum.set(i,j,o1.get(i,j)*o2);}
      return(dum);}
      Umatrix operator *(const Umatrix &o1, const double o2) {
        Umatrix dum;
        for (int i=0;i<NCOLOR;i++)
          for (int j=0;j<NCOLOR;j++) {
            dum.set(i,j,o1.get(i,j)*o2);}
        return(dum);}
        Umatrix operator *(const double o2, const Umatrix &o1) {
          Umatrix dum;
          for (int i=0;i<NCOLOR;i++)
            for (int j=0;j<NCOLOR;j++) {
              dum.set(i,j,o1.get(i,j)*o2);}
          return(dum);}
          Umatrix operator +(const Umatrix &x, const Umatrix &y) {
            Umatrix dum;
            for (int i=0;i<NCOLOR;i++)
              for (int j=0;j<NCOLOR;j++)
                dum.set(i,j,x.get(i,j)+y.get(i,j));
            return(dum);
          }
Umatrix operator -(const Umatrix &x, const Umatrix &y) {
  Umatrix dum;
  for (int i=0;i<NCOLOR;i++)
    for (int j=0;j<NCOLOR;j++)
      dum.set(i,j,x.get(i,j)-y.get(i,j));
  return(dum);
}


Umatrix comm(const Umatrix &o1, const Umatrix &o2) {
  return(o1*o2-o2*o1);
}

Umatrix exp(const Umatrix &u) {
  Umatrix c,del,prod;
  double fac=1.0;
  int i=1;
  prod=Umatrix(1);
  c=Umatrix(1);
  static int sum=0,counter=0;

  do{
    fac=fac*(double)i;
    prod=prod*u;
    del=prod*(1.0/fac);
    c=c+del;
    i++;}
  while(sqrt(Tr(del*Adj(del)).real())>GAUGETOL);

  sum+=i;
  counter++;
  if (counter==1000) {
    cout << "mean no. of terms in exp() "
      << (double)sum/counter << "\n" << flush;
    counter=0;sum=0;}

  return(c);}



  Complex Tr(const Umatrix &o) {
    Complex dum=Complex();
    for (int i=0;i<NCOLOR;i++)
      dum=dum+o.get(i,i);
    return(dum);}      


    Umatrix real_gaussian_Umatrix(void) {
      Umatrix dum=Umatrix();
      for (int a=0;a<RANK;a++) {
        dum=dum+gasdev()*Lambda[a];
      }
      return(dum);
    }


Gamma_Matrix::Gamma_Matrix(void) {
  for (int i=0;i<KDFERMION;i++) {
    for (int j=0;j<KDFERMION;j++) {gam[i][j]=0.0;}}
}

Gamma_Matrix::Gamma_Matrix(int k) {
  if (k==1) {
    for (int i=0;i<KDFERMION;i++) {
      for (int j=0;j<KDFERMION;j++) {gam[i][j]=0.0;}
    }
    for (int n=0;n<KDFERMION;n++) {
      gam[n][n]=1.0;}
  }
  else {cout << "wrong Gamma_Matrix constructor\n" << flush; }

}



double Gamma_Matrix::get(int i, int j) const {return(gam[i][j]);}
void Gamma_Matrix::set(int i, int j, const double o) {gam[i][j]=o;}
void Gamma_Matrix::print(void) {
  for (int i=0;i<KDFERMION;i++) {
    for (int j=0;j<KDFERMION;j++) {cout << gam[i][j];cout << "\t";}
    cout << "\n";}
}


Gamma_Matrix operator *(const Gamma_Matrix &o1, const Gamma_Matrix &o2) {
  Gamma_Matrix r;
  double dum;
  for (int i=0;i<KDFERMION;i++)
    for (int j=0;j<KDFERMION;j++) {
      dum=0.0;
      for (int k=0;k<KDFERMION;k++) {
        dum=dum+o1.get(i,k)*o2.get(k,j);}
      r.set(i,j,dum);}
  return(r);}

  Gamma_Matrix operator *(const double k, const Gamma_Matrix &o) {
    Gamma_Matrix r;
    for (int i=0;i<KDFERMION;i++) {
      for (int j=0;j<KDFERMION;j++) {
        r.set(i,j,k*o.get(i,j));
      }}
    return(r);
  }



Gamma_Matrix operator +(const Gamma_Matrix &o1, const Gamma_Matrix &o2) {
  Gamma_Matrix r;
  for (int i=0;i<KDFERMION;i++) {
    for (int j=0;j<KDFERMION;j++) {
      r.set(i,j,o1.get(i,j)+o2.get(i,j));
    }}
  return(r);
}   

Lattice_Vector::Lattice_Vector(void) {for (int i=0;i<D;i++)coords[i]=0;}
Lattice_Vector::Lattice_Vector(int mu) {for (int i=0;i<D;i++){
  coords[i]=0;}coords[mu]=1;}
  void Lattice_Vector::set(int i, int a) {
    coords[i]=a;
    return;}
    int Lattice_Vector::get(int i) const{
      return(coords[i]);}
      void Lattice_Vector::print(void) const {
        for (int i=0;i<D;i++) {cout << coords[i] << "\t";} cout << "\n";}

        Lattice_Vector operator +(const Lattice_Vector &x, const Lattice_Vector &y) {
          Lattice_Vector dum;
          for (int i=0;i<(D-1);i++)
            dum.set(i,(x.get(i)+y.get(i))%L);
          dum.set(D-1,(x.get(D-1)+y.get(D-1))%T);
          return(dum);}

          Lattice_Vector operator -(const Lattice_Vector &x, const Lattice_Vector &y) {
            Lattice_Vector dum;
            for (int i=0;i<(D-1);i++)
              dum.set(i,(x.get(i)-y.get(i)+L)%L);
            dum.set(D-1,(x.get(D-1)-y.get(D-1)+T)%T);
            return(dum);}

            Lattice_Vector operator -(const Lattice_Vector &x) {
              Lattice_Vector dum;
              for (int i=0;i<D;i++)
                dum.set(i,-1*x.get(i));
              return(dum);
            }

double BC(const Lattice_Vector &x, const Lattice_Vector &y) {

  if (x.get(D-1)+y.get(D-1)<0) return(PBC);
  if (x.get(D-1)+y.get(D-1)>(T-1))return(PBC);

  return(1.0);
}

int loop_over_lattice(Lattice_Vector &x, int &site) {
  int test,i,current;
  static int Lattice_Map[D];
  static int first_time=1;

  if (first_time) {
    for (i=0;i<D;i++)
      Lattice_Map[i]=(int)pow((double)L,(double)i);
    first_time=0;
  }


  current=site;
  for (i=D-1;i>=0;i--) {
    x.set(i,current/Lattice_Map[i]);
    current=current-Lattice_Map[i]*x.get(i);
  }

  if (current!=0) {cout << "error in loop_over_lattice" << "\n";}

  if (site==SITES)
    test=1;
  else
    test=0;

  site++;
  return(!test);
}



Gauge_Field::Gauge_Field(void) {
  for (int i=0;i<SITES;i++)
    for (int j=0;j<D;j++) {
      link[i][j]=Umatrix();}
}

Gauge_Field::Gauge_Field(int hot) {
  if (hot==2) {
    for (int i=0;i<SITES;i++) {
      for (int j=0;j<D;j++) {
        link[i][j]=real_gaussian_Umatrix();}}
    return;
  }
  if (hot==1) {
    for (int i=0;i<SITES;i++) {
      for (int j=0;j<D;j++) {
        link[i][j]=exp(real_gaussian_Umatrix());}}
    return;
  }
  if (hot==0) {
    for (int i=0;i<SITES;i++) {
      for (int j=0;j<D;j++) {
        link[i][j]=Umatrix(1);}}
    return;
  }
  cout << "error in gauge field constructor " << "\n" << flush;

  return;
}


Umatrix Gauge_Field::get(const Lattice_Vector &x, const int mu) const{
  int site=0,i;
  static int first_time=1;
  static int Lattice_Map[D];

  if (first_time) {
    for (i=0;i<D;i++)
      Lattice_Map[i]=(int)pow((double)L,(double)i);
    first_time=0;
  }

  for (i=0;i<D;i++)
    site=site+x.get(i)*Lattice_Map[i];

  return(link[site][mu]);
}

void Gauge_Field::set(const Lattice_Vector &x, const int mu, const Umatrix &u) {
  int site=0,i;
  static int first_time=1;
  static int Lattice_Map[D];

  if (first_time) {
    for (i=0;i<D;i++)
      Lattice_Map[i]=(int)pow((double)L,(double)i);
    first_time=0;
  }

  for (i=0;i<D;i++)
    site=site+x.get(i)*Lattice_Map[i];

  link[site][mu]=u;
  return;
}

Gauge_Field Adj(const Gauge_Field & U) {
  int site;
  Gauge_Field W;
  Lattice_Vector x;
  for (int mu=0;mu<D;mu++) {
    site=0;
    while(loop_over_lattice(x,site))
      W.set(x,mu,Adj(U.get(x,mu)));
  }
  return W;
}

Site_Field::Site_Field(void) {
  for (int i=0;i<SITES;i++) {
    points[i]=Umatrix();}
  return;
}

Site_Field::Site_Field(int c) {
  if (c==2) {
    for (int i=0;i<SITES;i++) {
      points[i]=real_gaussian_Umatrix()+Complex(0.0,1.0)*real_gaussian_Umatrix();
      points[i]=(1.0/sqrt(2.0))*points[i];}
  }
  if (c==1) {
    for (int i=0;i<SITES;i++) {
      points[i]=real_gaussian_Umatrix();}
    return;
  }
}

Umatrix Site_Field::get(const Lattice_Vector &x) const{
  int site=0,i;
  static int first_time=1;
  static int Lattice_Map[D];

  if (first_time) {
    for (i=0;i<D;i++)
      Lattice_Map[i]=(int)pow((double)L,(double)i);
    first_time=0;
  }

  for (i=0;i<D;i++)
    site=site+x.get(i)*Lattice_Map[i];

  return(points[site]);
}

void Site_Field::set(const Lattice_Vector &x, const Umatrix &u) {
  int site=0,i;
  static int first_time=1;
  static int Lattice_Map[D];

  if (first_time) {
    for (i=0;i<D;i++)
      Lattice_Map[i]=(int)pow((double)L,(double)i);
    first_time=0;
  }

  for (i=0;i<D;i++)
    site=site+x.get(i)*Lattice_Map[i];

  points[site]=u;
  return;
}

void Site_Field::print(void) {
  cout << "site field values\n" << flush;
  for (int i=0;i<SITES;i++) {
    cout << "site= " << i << "\n" << flush;
    cout << points[i] << "\n" << flush;}
  return;
}


Site_Field operator +(const Site_Field &s1, const Site_Field &s2) {
  int sites=0;
  Lattice_Vector x;
  Site_Field dum=Site_Field();

  while(loop_over_lattice(x,sites)) {
    dum.set(x,s1.get(x)+s2.get(x));
  }

  return(dum);
}


Site_Field operator -(const Site_Field &s1, const Site_Field &s2) {
  int sites=0;
  Lattice_Vector x;
  Site_Field dum=Site_Field();

  while(loop_over_lattice(x,sites)) {
    dum.set(x,s1.get(x)-s2.get(x));
  }

  return(dum);
}

Site_Field operator *(const double o, const Site_Field &s) {
  int sites=0;
  Lattice_Vector x;
  Site_Field dum=Site_Field();

  while(loop_over_lattice(x,sites)) {
    dum.set(x,o*s.get(x));
  }

  return(dum);
}


Site_Field operator *(const Complex &o, const Site_Field &s) {
  int sites=0;
  Lattice_Vector x;
  Site_Field dum=Site_Field();

  while(loop_over_lattice(x,sites)) {
    dum.set(x,o*s.get(x));
  }

  return(dum);
}

Umatrix operator *(const Site_Field &s1, const Site_Field &s2) {
  int sites=0;
  Lattice_Vector x;
  Umatrix dum=Umatrix();

  while(loop_over_lattice(x,sites)) {
    dum=dum+s1.get(x)*s2.get(x);
  }
  return(dum);
}

Site_Field Adj(const Site_Field &l) {
  int sites;
  Lattice_Vector x;
  Site_Field dum;

  sites=0;
  while(loop_over_lattice(x,sites)) {
    dum.set(x,Adj(l.get(x)));}

  return(dum);
}


Site_Field comm(const Site_Field &phi, const Site_Field &S) {
  Lattice_Vector x;
  int sites;
  Site_Field dum=Site_Field();
  sites=0;
  while(loop_over_lattice(x,sites)) {
    dum.set(x,phi.get(x)*S.get(x)-S.get(x)*phi.get(x));
  }
  return(dum);
}





Site_Field Dplus(const Gauge_Field &U, const Site_Field &s) {
  Lattice_Vector x,e_mu;
  int mu,sites;
  Site_Field dum=Site_Field();
  Umatrix tmp;
  Gauge_Field Udag;
  Udag=Adj(U);

  // computes covariant finite difference on adjoint site field

  sites=0;
  while(loop_over_lattice(x,sites)) {
    for (mu=0;mu<D;mu++) {
      e_mu=Lattice_Vector(mu);
      tmp=U.get(x,mu)*s.get(x+e_mu)*Udag.get(x,mu)*BC(x,e_mu)-s.get(x);
      dum.set(x,tmp);
    }}


  return(dum);
}

Site_Field Dminus(const Gauge_Field &U, const Site_Field &L) {
  Lattice_Vector x,e_mu;
  int mu,sites;
  Site_Field dum=Site_Field();
  Umatrix tmp;
  Gauge_Field Udag;

  Udag=Adj(U);

  sites=0;
  while(loop_over_lattice(x,sites)) {
    tmp=Umatrix();

    for (mu=0;mu<D;mu++) {
      e_mu=Lattice_Vector(mu);

      tmp=tmp+L.get(x)-
        Udag.get(x-e_mu,mu)*L.get(x-e_mu)*U.get(x-e_mu,mu)*BC(x,-e_mu);
    }

    dum.set(x,tmp);
  }
  return(dum);
}


void Build_Scalar(const Site_Field phi[NSCALAR], Site_Field Scalar[KDFERMION][KDFERMION]) {
  int a,b,i;

  for (a = 0; a < KDFERMION; a++) {
    for (b=0;b<KDFERMION;b++) {
      Scalar[a][b]=Site_Field();

      for (i=0;i<NSCALAR-2;i++)
        Scalar[a][b]=Scalar[a][b]+Gamma[i].get(a,b)*phi[i];
    }
  }

  return;
}

void Fermion_operator(const Gauge_Field &U,
    const Site_Field Scalar[KDFERMION][KDFERMION],
    const Site_Field phi[NSCALAR],
    const Site_Field K[NFERMION], Site_Field MK[NFERMION], const int sign) {

  int i,a,b;
  Site_Field Y[NFERMION],B[NFERMION];

  // if sign==1 operator returned. If sign==-1 its adjoint

  // kinetic terms
  for (i = 0; i < KDFERMION; i++) {
    MK[i] = Dplus(U, K[i + KDFERMION]);
    MK[i + KDFERMION] = Dminus(U, K[i]);
  }

  if (sign==(-1)) {
    for (i = 0; i < NFERMION; i++)
      MK[i] = -1.0 * MK[i];
  }

  // BMN mass
  for (a = 0; a < KDFERMION; a++) {
    B[a] = Site_Field();
    B[a + KDFERMION] = Site_Field();
    for (b = 0; b < KDFERMION; b++) {
      B[a] = B[a] + 3 * MU * 0.25 * Gam123.get(a,b) * K[b + KDFERMION];
      B[a + KDFERMION] = B[a + KDFERMION] - 3 * MU * 0.25 * Gam123.get(b, a) * K[b];
    }
  }

  if (sign == 1) {
    for (a = 0; a < NFERMION; a++)
      MK[a] = MK[a] + B[a];
  }
  else {
    for (a = 0; a < NFERMION; a++)
      MK[a] = MK[a] - B[a];
  }

  // Yukawas
  for (a = 0; a < KDFERMION; a++) {
    Y[a] = Site_Field();
    Y[a + KDFERMION] = Site_Field();
    for (b = 0; b < KDFERMION; b++) {
      Y[a] = Y[a] + comm(Scalar[a][b], K[b + KDFERMION]);
      Y[a+KDFERMION] = Y[a + KDFERMION] + comm(Scalar[b][a], K[b]);
    }
  }

  // Last 2 gamma's diagonal in this basis. Gamma[7]=(-I0//0I) Gamma[8]=(i0//0i)
  for (a = 0; a < KDFERMION; a++) {
    Y[a] = Y[a] + Complex(-1.0, 0.0) * comm(phi[NSCALAR - 2], K[a]);
    Y[a+KDFERMION] = Y[a+KDFERMION]
                   + Complex(1.0, 0.0) * comm(phi[NSCALAR - 2], K[a + KDFERMION]);
  }

  if (sign==(-1)) {
    for (a = 0; a < NFERMION; a++) {
      Y[a]=-1.0*Y[a];}}

  for (a = 0; a < NFERMION; a++) {
    Y[a] = Y[a] + Complex(0.0, 1.0) * comm(phi[NSCALAR - 1], K[a]);
  }

  for (i = 0; i < NFERMION; i++)
    MK[i]=NORM*(MK[i]+Y[i]);

  return;
}
