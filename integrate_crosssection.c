// 
//  compile with "gcc -c integrate_crosssection.c"
//  link    with "gcc -lm integrate_crosssection.o nvegas.o genps.o -o Integrate"
//  execute ./Integrate
// 
// to link MG use: "icc integrate_crosssection.o nvegas.o genps.o emep_taptam.o MGWrapper.a -lifcore -limf -o Integrate"
// 
// to link PDFs use: "gcc integrate_crosssection.c -lm nvegas.o genps.o ./PDFs/mstwpdf.o ./PDFs/alphaS.o -lgfortran -o Integrate"


#include <math.h>
#include <stdio.h>
#include "vegas.h"


#define DIMENSION 2
#define FUNCTIONS 1




//----------------------------------------------------------------------------------
double Dot(double *pa,double *pb) // Minkowski dot product
{
 
 return pa[0]*pb[0]-pa[1]*pb[1]-pa[2]*pb[2]-pa[3]*pb[3];
 
}
//----------------------------------------------------------------------------------







//----------------------------------------------------------------------------------

double MEsq(double *p1,double *p2,double *p3,double *p4) // dummy matrix element squared
{
 
 double PreFactor = 4;
 double MEsq;
 
    MEsq = PreFactor * Dot(p1,p3)*Dot(p2,p4);
 
 return MEsq;
 
}
//----------------------------------------------------------------------------------





//----------------------------------------------------------------------------------
void TestCrossSection(double x[DIMENSION], double f[FUNCTIONS])
{ 
  
  int NPart=2;
  double CMSEnergy = 300.0;
  double Mass[2];
  Mass[0] = 1.777;
  Mass[1] = 1.777;
  double pOut[2][4];
  double Jacobian;
  genps_(&NPart, &CMSEnergy, x, Mass, pOut, &Jacobian );    // this function input:  number of particles, Center-of-mass energy, integration variables x=0..1, masses
                                                            // this function output: momenta, jacobian factor
  
  double pIn[2][4];
  pIn[0][0] =+CMSEnergy/2.0;
  pIn[0][1] = 0.0;
  pIn[0][2] = 0.0;
  pIn[0][3] =+CMSEnergy/2.0;

  pIn[1][0] =+CMSEnergy/2.0;
  pIn[1][1] = 0.0;
  pIn[1][2] = 0.0;
  pIn[1][3] =-CMSEnergy/2.0;
  
  
  
//------- calling pdf -----------------------------  
  
  char path[] = "mstw2008lo.00.dat";
  int iParton;
  int iSet;
  double x1,q,pdf_u;
  
// initializing  
  q = 500.0;
  iSet=0;

// momentum fraction x1
  x1 = 0.000123;
  
// iParton =   -6,  -5,  -4,  -3,  -2,  -1,0,1,2,3,4,5,6
//         = tbar,bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b,t.
  iParton = 3;  // selecting a strange quark
  
// calling the pdf function
  getonepdf_(&iSet,&x1,&q,&iParton,&pdf_u);

// printing the result
  printf("%f %f %f \n",x1,q,pdf_u);
  
// -------------------------------------------------  
  
  
  
  // check the generated momenta
   
//   printf("x : %10.4f %10.4f \n",x[0],x[1]);;
//   printf("p1: %10.4f %10.4f %10.4f %10.4f \n",pIn[0][0],pIn[0][1],pIn[0][2],pIn[0][3]);
//   printf("p2: %10.4f %10.4f %10.4f %10.4f \n",pIn[0][0],pIn[0][1],pIn[0][2],pIn[0][3]);
//   printf("p3: %10.4f %10.4f %10.4f %10.4f \n",pOut[0][0],pOut[0][1],pOut[0][2],pOut[0][3]);
//   printf("p4: %10.4f %10.4f %10.4f %10.4f \n",pOut[1][0],pOut[1][1],pOut[1][2],pOut[1][3]);
//   printf("check on-shell ness sqrt(p1^2): %10.4f \n", sqrt(Dot(pOut[0],pOut[0])));
//   printf("check on-shell ness sqrt(p2^2): %10.4f \n", sqrt(Dot(pOut[1],pOut[1])));
//   printf("check energy-momentum conservation p1+p2: %10.4f \n",pOut[0][0]+pOut[1][0]);
//   printf("check energy-momentum conservation p1+p2: %10.4f \n",pOut[0][1]+pOut[1][1]);
//   printf("check energy-momentum conservation p1+p2: %10.4f \n",pOut[0][2]+pOut[1][2]);
//   printf("check energy-momentum conservation p1+p2: %10.4f \n",pOut[0][3]+pOut[1][3]);


 f[0] = Jacobian * MEsq(pIn[0],pIn[1],pOut[0],pOut[1]);


//  double MGRes;
//  double MomExtMG[4][4];
//  
//  MomExtMG[0][0]=pIn[0][0];
//  MomExtMG[0][1]=pIn[0][1];
//  MomExtMG[0][2]=pIn[0][2];
//  MomExtMG[0][3]=pIn[0][3];
//  
//  MomExtMG[1][0]=pIn[1][0];
//  MomExtMG[1][1]=pIn[1][1];
//  MomExtMG[1][2]=pIn[1][2];
//  MomExtMG[1][3]=pIn[1][3];
//  
//  MomExtMG[2][0]=pOut[0][0];
//  MomExtMG[2][1]=pOut[0][1];
//  MomExtMG[2][2]=pOut[0][2];
//  MomExtMG[2][3]=pOut[0][3];
//  
//  MomExtMG[3][0]=pOut[1][0];
//  MomExtMG[3][1]=pOut[1][1];
//  MomExtMG[3][2]=pOut[1][2];
//  MomExtMG[3][3]=pOut[1][3];
//  
//  int set=4;
//  mgwrapper(&set,MomExtMG,&MGRes);
//  printf("MG Result %f \n",MGRes);
 


  return;
}
//----------------------------------------------------------------------------------









// MAIN PROGRAM
int main(int argc, char **argv)
{ /* declare variables */
  int i;
  double estim[FUNCTIONS];   /* estimators for integrals                     */
  double std_dev[FUNCTIONS]; /* standard deviations                          */
  double chi2a[FUNCTIONS];   /* chi^2/n                                      */
  double reg[2*DIMENSION];   /* integration domain                           */



  //  initializing the integration range, always from 0.0 to 1.0
  for (i=0; i<DIMENSION; i++) {
    reg[i] = 0.0;
    reg[i+DIMENSION] = 1.0;
  }

  // setting parameters for vegas integrator
  int init=0;
  unsigned long ncall=500000;
  int itmx=10;
  
  // calling vegas integrator
  vegas(reg, DIMENSION, TestCrossSection, init, ncall, itmx, NPRN_INPUT | NPRN_RESULT, 1, 0, 1, estim, std_dev, chi2a);
  

  return(0);
}


