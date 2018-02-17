#include "my_gen.h"

const int POS1 = (NCOLOR*(NCOLOR-1)/2);
const int POS2 = (NCOLOR*(NCOLOR-1));

const double ROOT2 = 1.41421356237309504880168872421;

int posmat[NCOLOR][NCOLOR] ;
Complex genmat[RANK][NCOLOR][NCOLOR] ;

Complex tr[RANK][RANK], strconst[RANK][RANK][RANK] ;


void initgen(void) ;
void vectomat(Complex [RANK], Complex [NCOLOR][NCOLOR]) ;


void computegen(void) ;


int my_gen(void)
{
  int aa, ii, jj ;
  

  initgen() ;

  computegen() ;

  /* map to data structures in code */
  
  for(aa=0;aa<RANK;aa++) {

    for(ii=0;ii<NCOLOR;ii++) {
      for(jj=0;jj<NCOLOR;jj++) {
	
Lambda[aa].set(ii,jj,(Complex(0.0,1.0)*sqrt(2.0))*genmat[aa][ii][jj]);
      } 
    }

  }


  return(1) ;

}



/* Compute generator matrices */

void computegen(void) 
{
  int aa, ii, jj ;
Complex adj[RANK] ;
Complex fund[NCOLOR][NCOLOR] ; 

  for(aa=0;aa<RANK;aa++) {

    for(ii=0;ii<RANK;ii++) adj[ii] = Complex(0.0,0.0) ;
    adj[aa] = Complex(1.0,0.0) ;

    vectomat(adj,fund) ;

    for(ii=0;ii<NCOLOR;ii++) {
      for(jj=0;jj<NCOLOR;jj++) {
	genmat[aa][ii][jj] = fund[ii][jj] ;
      } 
    }

  }

}



/* SU(N) routines */

void initgen(void)
{
  int ii,jj,aa ;

  aa=0 ;

  for(ii=0;ii<NCOLOR;ii++) {
    for(jj=ii+1;jj<NCOLOR;jj++) {

      posmat[ii][jj] = aa ;

      aa+=1 ;

    }
  }

}


/* Adjoint to fundamental */

void vectomat(Complex vec[RANK], Complex mat[NCOLOR][NCOLOR])
{
int ii,jj ;
Complex mult ;

  for(ii=0;ii<NCOLOR;ii++) {
    for(jj=ii+1;jj<NCOLOR;jj++) {
      mat[ii][jj] = 0.5*(vec[posmat[ii][jj]]-Complex(0.0,1.0)*vec[POS1+posmat[ii][jj]]) ;
      mat[jj][ii] = 0.5*(vec[posmat[ii][jj]]+Complex(0.0,1.0)*vec[POS1+posmat[ii][jj]]) ;
    }
  }
    
  for(ii=0;ii<NCOLOR;ii++) 
    mat[ii][ii] = Complex(0.0,0.0) ; 
    
  for(ii=0;ii<NCOLOR-1;ii++) {
    mult = vec[POS2+ii]*(1.0/sqrt(2.+2./(1.+ii))/(1.+ii)) ;
    for(jj=0;jj<ii+1;jj++) {    
      mat[jj][jj] = mat[jj][jj]+mult ;
    }
    mat[ii+1][ii+1] = mat[ii+1][ii+1]- 1.0*Complex(1.0+ii,0.0)*mult ; 
  }
    
}


