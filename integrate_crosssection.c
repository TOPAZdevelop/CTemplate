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


#define DIMENSION 4

#define FUNCTIONS 2



//----------------------------------------------------------------------------------
double Dot(double *pa,double *pb) // Minkowski dot product
{
 
 return pa[0]*pb[0]-pa[1]*pb[1]-pa[2]*pb[2]-pa[3]*pb[3];
 
}
//----------------------------------------------------------------------------------







//----------------------------------------------------------------------------------

double ME1sq(double *p1,double *p2,double *p3,double *p4) // spin-1 matrix element squared
{
 double Pi=3.14159265;
 double alpha=1.0/128.0;
 double alpha_s=0.11;
 double C_DMspin1=0.1;
 double PreFactor = 8.0*(4.0*Pi*alpha_s)*C_DMspin1*C_DMspin1;
 double MEsq,s,t,u,M2;
 
    M2= 2*Dot(p3,p3);
    s = 2*Dot(p1,p2);
    t = 2*Dot(p1,p3)+M2;
    u = 2*Dot(p2,p3)+M2;

    MEsq = PreFactor * (2*s*s + t*t + u*u + 2*s*(t+u))/(t*u);
 
 return MEsq;
 
}



double ME0sq(double *p1,double *p2,double *p3,double *p4) // spin-0 matrix element squared
{
 double Pi=3.14159265;
 double alpha=1.0/128.0;
 double alpha_s=0.11;
 double C_DMspin0=0.1;
 double PreFactor = 4.0*(4.0*Pi*alpha_s)*C_DMspin0*C_DMspin0;
 double MEsq,s,t,u,M2;
 
    M2= 2*Dot(p3,p3);
    s = 2*Dot(p1,p2);
    t = 2*Dot(p1,p3)+M2;
    u = 2*Dot(p2,p3)+M2;

    MEsq = PreFactor * (s*s + M2*M2)/(t*u);
 
 return MEsq;
 
}

//----------------------------------------------------------------------------------





//----------------------------------------------------------------------------------
void CrossSection(double x[DIMENSION], double f[FUNCTIONS])
{ 
// declare variables   
  int NPart,iSet,iParton;
  double x1,x2,q,Flux,pTg,pTg_cut;
  double Mass[2],pOut[2][4],pIn[2][4];
  double sigma_spin0_qq,sigmahat_spin0_qq;
  double sigma_spin1_qq,sigmahat_spin1_qq;
  double Jacobian,fbGeV2,shat,CMSEnergy,ColliderEnergy,ColAvg_qq;
  double pdf_u,pdf_d,pdf_c,pdf_s,pdf_b,pdf_g;
  double pdf_ub,pdf_db,pdf_cb,pdf_sb,pdf_bb;
  double vev=1.0;
  double mov_u=(2e-3)/vev,mov_d=(5e-3)/vev,mov_c=(1.3)/vev,mov_s=(95e-3)/vev,mov_b=(4.2)/vev;


// initialize variables 
  f[0] = 0.0;
  f[1] = 0.0;
  ColliderEnergy = 13000;
  NPart=2;
  Mass[0] = 50.0; // spin-0,1
  Mass[1] = 0.0; // gluon 
  pTg_cut = 20.0;
  q = 100.0; // Mass[0]+Mass[1];
  iSet=0;
  x1 = x[2];
  x2 = x[3];

  shat = x1*x2*ColliderEnergy*ColliderEnergy;
  CMSEnergy = sqrt(shat);
  Flux = 1.0/(2.0*shat);
  ColAvg_qq = 1.0/3.0/3.0;
  fbGeV2 = 0.389379*1e12;

  if ( CMSEnergy < q ){
    return;
  };

  genps_(&NPart, &CMSEnergy, x, Mass, pOut, &Jacobian );    // this function input:  number of particles, Center-of-mass energy, integration variables x=0..1, masses
                                                            // this function output: momenta, jacobian factor
  pIn[0][0] =+CMSEnergy/2.0;
  pIn[0][1] = 0.0;
  pIn[0][2] = 0.0;
  pIn[0][3] =+CMSEnergy/2.0;

  pIn[1][0] =+CMSEnergy/2.0;
  pIn[1][1] = 0.0;
  pIn[1][2] = 0.0;
  pIn[1][3] =-CMSEnergy/2.0;
  
  pTg = sqrt(pOut[1][1]*pOut[1][1]+pOut[1][2]*pOut[1][2]);
  if ( pTg < pTg_cut ){
    return;
  };  
  
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



  
  

// iParton =   -6,  -5,  -4,  -3,  -2,  -1,0,1,2,3,4,5,6
//         = tbar,bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b,t.

  // P_left( q(x1) ) * P_right( qbar(x2 )

  iParton = 0;
  getonepdf_(&iSet,&x1,&q,&iParton,&pdf_g);
  iParton = 1;
  getonepdf_(&iSet,&x1,&q,&iParton,&pdf_d);
  iParton = 2;
  getonepdf_(&iSet,&x1,&q,&iParton,&pdf_u);
  iParton = 3;
  getonepdf_(&iSet,&x1,&q,&iParton,&pdf_s);
  iParton = 4;
  getonepdf_(&iSet,&x1,&q,&iParton,&pdf_c);
  iParton = 5;
  getonepdf_(&iSet,&x1,&q,&iParton,&pdf_b);

  iParton = -1;
  getonepdf_(&iSet,&x2,&q,&iParton,&pdf_db);
  iParton = -2;
  getonepdf_(&iSet,&x2,&q,&iParton,&pdf_ub);
  iParton = -3;
  getonepdf_(&iSet,&x2,&q,&iParton,&pdf_sb);
  iParton = -4;
  getonepdf_(&iSet,&x2,&q,&iParton,&pdf_cb);
  iParton = -5;
  getonepdf_(&iSet,&x2,&q,&iParton,&pdf_bb);

  sigmahat_spin0_qq = ColAvg_qq * Flux * Jacobian * ME0sq(pIn[0],pIn[1],pOut[0],pOut[1]);
  sigmahat_spin1_qq = ColAvg_qq * Flux * Jacobian * ME1sq(pIn[0],pIn[1],pOut[0],pOut[1]);
  sigma_spin0_qq    = sigmahat_spin0_qq * (pdf_u*pdf_ub*mov_u*mov_u + pdf_d*pdf_db*mov_d*mov_d + pdf_c*pdf_cb*mov_c*mov_c + pdf_s*pdf_sb*mov_s*mov_s + pdf_b*pdf_bb*mov_b*mov_b);
  sigma_spin1_qq    = sigmahat_spin1_qq * (pdf_u*pdf_ub + pdf_d*pdf_db + pdf_c*pdf_cb + pdf_s*pdf_sb + pdf_b*pdf_bb);



  // P_left( qbar(x1) ) * P_right( q(x2 )

  iParton = 0;
  getonepdf_(&iSet,&x2,&q,&iParton,&pdf_g);
  iParton = 1;
  getonepdf_(&iSet,&x2,&q,&iParton,&pdf_d);
  iParton = 2;
  getonepdf_(&iSet,&x2,&q,&iParton,&pdf_u);
  iParton = 3;
  getonepdf_(&iSet,&x2,&q,&iParton,&pdf_s);
  iParton = 4;
  getonepdf_(&iSet,&x2,&q,&iParton,&pdf_c);
  iParton = 5;
  getonepdf_(&iSet,&x2,&q,&iParton,&pdf_b);

  iParton = -1;
  getonepdf_(&iSet,&x1,&q,&iParton,&pdf_db);
  iParton = -2;
  getonepdf_(&iSet,&x1,&q,&iParton,&pdf_ub);
  iParton = -3;
  getonepdf_(&iSet,&x1,&q,&iParton,&pdf_sb);
  iParton = -4;
  getonepdf_(&iSet,&x1,&q,&iParton,&pdf_cb);
  iParton = -5;
  getonepdf_(&iSet,&x1,&q,&iParton,&pdf_bb);

  sigmahat_spin0_qq = ColAvg_qq * Flux * Jacobian * ME0sq(pIn[1],pIn[0],pOut[0],pOut[1]);
  sigmahat_spin1_qq = ColAvg_qq * Flux * Jacobian * ME1sq(pIn[1],pIn[0],pOut[0],pOut[1]);
  sigma_spin0_qq   += sigmahat_spin0_qq * (pdf_u*pdf_ub*mov_u*mov_u + pdf_d*pdf_db*mov_d*mov_d + pdf_c*pdf_cb*mov_c*mov_c + pdf_s*pdf_sb*mov_s*mov_s + pdf_b*pdf_bb*mov_b*mov_b);
  sigma_spin1_qq   += sigmahat_spin1_qq * (pdf_ub*pdf_u + pdf_db*pdf_d + pdf_cb*pdf_c + pdf_sb*pdf_s + pdf_bb*pdf_b);


 f[0] = fbGeV2 * sigma_spin0_qq;
 f[1] = fbGeV2 * sigma_spin1_qq;


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
  vegas(reg, DIMENSION, CrossSection, init, ncall, itmx,0x0001 | 0x0002 | 0x0004, FUNCTIONS, 0, 1, estim, std_dev, chi2a);
  	

  return(0);
}


