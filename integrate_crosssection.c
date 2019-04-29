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
#include <stdlib.h>
#include "vegas.h"


#define DIMENSION 4

#define FUNCTIONS 2

#ifndef max
	#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

const double Pi=3.14159265359;
      double alpha_s;
const double alpha=0.00729927007299; // 0.0075624689028; //
const double GeV=1.0;

const double C_DMspin0=1.0;
const double C_DMspin1=0.1;


#define NUMHIST 2
#define NUMBINS 20

double Histogram_spin0[NUMHIST][NUMBINS];
double Histogram_spin1[NUMHIST][NUMBINS];


//----------------------------------------------------------------------------------
double Dot(double *pa,double *pb) // Minkowski dot product
{
 
 return pa[0]*pb[0]-pa[1]*pb[1]-pa[2]*pb[2]-pa[3]*pb[3];
 
}
//----------------------------------------------------------------------------------







//----------------------------------------------------------------------------------

double ME1sq_qqb(double *p1,double *p2,double *p3,double *p4) // spin-1 matrix element squared
{
 double PreFactor = 8.0*(4.0*Pi*alpha_s)*C_DMspin1*C_DMspin1;
 double MEsq,s,t,u,M2;
 
    M2=   Dot(p3,p3);
    s = 2*Dot(p1,p2);
    t =-2*Dot(p1,p3)+M2;
    u =-2*Dot(p2,p3)+M2;

    MEsq = PreFactor * (2*s*s + t*t + u*u + 2*s*(t+u))/(t*u);

 return MEsq;
 
}



double ME0sq_qqb(double *p1,double *p2,double *p3,double *p4) // spin-0 matrix element squared
{
 double PreFactor = 4.0*(4.0*Pi*alpha_s)*C_DMspin0*C_DMspin0;
 double MEsq,s,t,u,M2;
 
    M2=   Dot(p3,p3);
    s = 2*Dot(p1,p2);
    t =-2*Dot(p1,p3)+M2;
    u =-2*Dot(p2,p3)+M2;

    MEsq = PreFactor * (s*s + M2*M2)/(t*u);

 return MEsq;
 
}

double ME1sq_gq(double *p1,double *p2,double *p3,double *p4) // spin-1 matrix element squared
{
 double PreFactor = 8.0*(4.0*Pi*alpha_s)*C_DMspin1*C_DMspin1;
 double MEsq,s,t,u,M2;
 
    M2=   Dot(p3,p3);
    t = 2*Dot(p1,p2);
    s =-2*Dot(p1,p3)+M2;
    u =-2*Dot(p2,p3)+M2;

    MEsq = -PreFactor * (2*s*s + t*t + u*u + 2*s*(t+u))/(t*u);

 return MEsq;
 
}



double ME0sq_gq(double *p1,double *p2,double *p3,double *p4) // spin-0 matrix element squared
{
 double PreFactor = 4.0*(4.0*Pi*alpha_s)*C_DMspin0*C_DMspin0;
 double MEsq,s,t,u,M2;
 
    M2=   Dot(p3,p3);
    t = 2*Dot(p1,p2);
    s =-2*Dot(p1,p3)+M2;
    u =-2*Dot(p2,p3)+M2;

    MEsq = -PreFactor * (s*s + M2*M2)/(t*u);

 return MEsq;
 
}



double ME1sq_qbg(double *p1,double *p2,double *p3,double *p4) // spin-1 matrix element squared
{
 double PreFactor = 8.0*(4.0*Pi*alpha_s)*C_DMspin1*C_DMspin1;
 double MEsq,s,t,u,M2;
 
    M2=   Dot(p3,p3);
    u = 2*Dot(p1,p2);
    t =-2*Dot(p1,p3)+M2;
    s =-2*Dot(p2,p3)+M2;

    MEsq = -PreFactor * (2*s*s + t*t + u*u + 2*s*(t+u))/(t*u);

 return MEsq;
 
}



double ME0sq_qbg(double *p1,double *p2,double *p3,double *p4) // spin-0 matrix element squared
{
 double PreFactor = 4.0*(4.0*Pi*alpha_s)*C_DMspin0*C_DMspin0;
 double MEsq,s,t,u,M2;
 
    M2=   Dot(p3,p3);
    u = 2*Dot(p1,p2);
    t =-2*Dot(p1,p3)+M2;
    s =-2*Dot(p2,p3)+M2;

    MEsq = -PreFactor * (s*s + M2*M2)/(t*u);

 return MEsq;
 
}



int WhichBin(int iBin, double ObsVal)
{ 
  double LowVal,BinSize;

  switch(iBin) {
    case 1: LowVal=20.0*GeV; BinSize=10*GeV; break;  // pT_glu
    case 2: LowVal=-5.0; BinSize=0.5;        break;  // y_glu
  };

  int TheBin = (ObsVal-LowVal)/BinSize;
  if( TheBin >= NUMBINS ) TheBin=NUMBINS-1;
  if( TheBin < 0 ) TheBin=0;
  return TheBin;
}

//----------------------------------------------------------------------------------


// qqbar+qbarq --> photon+glu
// Value of final        lo integral is   13.6682     +/-  0.25195E-02 nb
// threads:   1 time:    10.53    x-section: 0.136681846494E+08 fb

// check:      1.additional integral=   1.364921e+07+/-  2.9e+04  chi^2/IT n =      0.79


//----------------------------------------------------------------------------------
void CrossSection(double x[DIMENSION], double *weight, double f[FUNCTIONS])
{ 
// declare variables   
  int NPart,iSet,iParton;
  double x1,x2,q,Flux,pTg,yg,pTg_cut;
  double Mass[2],pOut[2][4],pIn[2][4];
  double sigma_spin0_qqb,sigmahat_spin0_qqb;
  double sigma_spin1_qqb,sigmahat_spin1_qqb;
  double sigma_spin0_qbq,sigmahat_spin0_qbq;
  double sigma_spin1_qbq,sigmahat_spin1_qbq;
  double sigma_spin0_qg,sigmahat_spin0_qg;
  double sigma_spin1_qg,sigmahat_spin1_qg;
  double sigma_spin0_qbg,sigmahat_spin0_qbg;
  double sigma_spin1_qbg,sigmahat_spin1_qbg;
  double sigma_spin0_gq,sigmahat_spin0_gq;
  double sigma_spin1_gq,sigmahat_spin1_gq;
  double sigma_spin0_gqb,sigmahat_spin0_gqb;
  double sigma_spin1_gqb,sigmahat_spin1_gqb;
  double Jacobian,shat,CMSEnergy,ColliderEnergy;
  double pdf1_u,pdf1_d,pdf1_c,pdf1_s,pdf1_b,pdf1_g;
  double pdf1_ub,pdf1_db,pdf1_cb,pdf1_sb,pdf1_bb;
  double pdf2_u,pdf2_d,pdf2_c,pdf2_s,pdf2_b,pdf2_g;
  double pdf2_ub,pdf2_db,pdf2_cb,pdf2_sb,pdf2_bb;
  double vev=1*GeV;
  double mov2_u=((2e-3)*GeV/vev)*((2e-3)*GeV/vev);
  double mov2_d=((5e-3)*GeV/vev)*((5e-3)*GeV/vev);
  double mov2_c=((1.3)*GeV/vev)*((1.3)*GeV/vev);
  double mov2_s=((95e-3)*GeV/vev)*((95e-3)*GeV/vev);
  double mov2_b=((4.2)*GeV/vev)*((4.2)*GeV/vev);
  double Qu2=(2.0/3.0)*(2.0/3.0);
  double Qd2=(1.0/3.0)*(1.0/3.0);
  const double PiWgt2 = 1.0/Pi;
  const double ColAvg_qq = 1.0/3.0/3.0;
  const double ColAvg_qg = 1.0/3.0/8.0;
  const double fbGeV2 = 0.389379*1e12;

//x[0]=0.15;
//x[1]=0.76;
//x[2]=0.47;
//x[3]=0.22;

// initialize variables 
  f[0] = 0.0;
  f[1] = 0.0;
  ColliderEnergy = 13000*GeV;
  NPart=2;
  Mass[0] = 10.0*GeV; // spin-0,1 DM
  Mass[1] = 0.0*GeV; // gluon 
  pTg_cut = 20.0*GeV;
  iSet=0;
  x1 = x[2];
  x2 = x[3];

  shat = x1*x2*ColliderEnergy*ColliderEnergy;
  CMSEnergy = sqrt(shat);
  q = CMSEnergy; 
  Flux = 1.0/(2.0*shat);

  if ( CMSEnergy < Mass[0]+Mass[1]+pTg_cut ){
    return;
  };

  genps_(&NPart, &CMSEnergy, x, Mass, pOut, &Jacobian );    // this function input:  number of particles, Center-of-mass energy, integration variables x=0..1, masses
  Jacobian *= PiWgt2;                                       // this function output: momenta, jacobian factor

  pIn[0][0] =+CMSEnergy/2.0;
  pIn[0][1] = 0.0;
  pIn[0][2] = 0.0;
  pIn[0][3] =+CMSEnergy/2.0;

  pIn[1][0] =+CMSEnergy/2.0;
  pIn[1][1] = 0.0;
  pIn[1][2] = 0.0;
  pIn[1][3] =-CMSEnergy/2.0;
  
  pTg = sqrt(pOut[1][1]*pOut[1][1]+pOut[1][2]*pOut[1][2]);
  yg = 0.5*log((pOut[1][0]+pOut[1][3])/(pOut[1][0]-pOut[1][3]));
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
  getonepdf_(&iSet,&x1,&q,&iParton,&pdf1_g);
  iParton = 1;
  getonepdf_(&iSet,&x1,&q,&iParton,&pdf1_d);
  iParton = 2;
  getonepdf_(&iSet,&x1,&q,&iParton,&pdf1_u);
  iParton = 3;
  getonepdf_(&iSet,&x1,&q,&iParton,&pdf1_s);
  iParton = 4;
  getonepdf_(&iSet,&x1,&q,&iParton,&pdf1_c);
  iParton = 5;
  getonepdf_(&iSet,&x1,&q,&iParton,&pdf1_b);

  iParton = -1;
  getonepdf_(&iSet,&x1,&q,&iParton,&pdf1_db);
  iParton = -2;
  getonepdf_(&iSet,&x1,&q,&iParton,&pdf1_ub);
  iParton = -3;
  getonepdf_(&iSet,&x1,&q,&iParton,&pdf1_sb);
  iParton = -4;
  getonepdf_(&iSet,&x1,&q,&iParton,&pdf1_cb);
  iParton = -5;
  getonepdf_(&iSet,&x1,&q,&iParton,&pdf1_bb);    
  
  iParton = -1;
  getonepdf_(&iSet,&x2,&q,&iParton,&pdf2_db);
  iParton = -2;
  getonepdf_(&iSet,&x2,&q,&iParton,&pdf2_ub);
  iParton = -3;
  getonepdf_(&iSet,&x2,&q,&iParton,&pdf2_sb);
  iParton = -4;
  getonepdf_(&iSet,&x2,&q,&iParton,&pdf2_cb);
  iParton = -5;
  getonepdf_(&iSet,&x2,&q,&iParton,&pdf2_bb);

  iParton = 0;
  getonepdf_(&iSet,&x2,&q,&iParton,&pdf2_g);
  iParton = 1;
  getonepdf_(&iSet,&x2,&q,&iParton,&pdf2_d);
  iParton = 2;
  getonepdf_(&iSet,&x2,&q,&iParton,&pdf2_u);
  iParton = 3;
  getonepdf_(&iSet,&x2,&q,&iParton,&pdf2_s);
  iParton = 4;
  getonepdf_(&iSet,&x2,&q,&iParton,&pdf2_c);
  iParton = 5;
  getonepdf_(&iSet,&x2,&q,&iParton,&pdf2_b);

  Jacobian *=1.0/x1/x2;

  
  alphas_(&q,&alpha_s);

  sigmahat_spin0_qqb= ColAvg_qq * Flux * Jacobian * ME0sq_qqb(pIn[0],pIn[1],pOut[0],pOut[1]);
  sigmahat_spin1_qqb= ColAvg_qq * Flux * Jacobian * ME1sq_qqb(pIn[0],pIn[1],pOut[0],pOut[1]);
  sigma_spin0_qqb   = sigmahat_spin0_qqb* (pdf1_u*pdf2_ub*mov2_u + pdf1_d*pdf2_db*mov2_d+ pdf1_c*pdf2_cb*mov2_c + pdf1_s*pdf2_sb*mov2_s + pdf1_b*pdf2_bb*mov2_b);
  sigma_spin1_qqb   = sigmahat_spin1_qqb* (pdf1_u*pdf2_ub + pdf1_d*pdf2_db + pdf1_c*pdf2_cb + pdf1_s*pdf2_sb + pdf1_b*pdf2_bb);
//  sigma_spin1_qq    = sigmahat_spin1_qq* (pdf_u*pdf_ub*Qu2 + pdf_d*pdf_db*Qd2 + pdf_c*pdf_cb*Qu2 + pdf_s*pdf_sb*Qd2 + pdf_b*pdf_bb*Qd2);


  sigmahat_spin0_qbq = ColAvg_qq * Flux * Jacobian * ME0sq_qqb(pIn[1],pIn[0],pOut[0],pOut[1]);
  sigmahat_spin1_qbq = ColAvg_qq * Flux * Jacobian * ME1sq_qqb(pIn[1],pIn[0],pOut[0],pOut[1]);
  sigma_spin0_qbq    = sigmahat_spin0_qbq* (pdf1_ub*pdf2_u*mov2_u + pdf1_db*pdf2_d*mov2_d+ pdf1_cb*pdf2_c*mov2_c + pdf1_sb*pdf2_s*mov2_s + pdf1_bb*pdf2_b*mov2_b);
  sigma_spin1_qbq    = sigmahat_spin1_qbq* (pdf1_ub*pdf2_u + pdf1_db*pdf2_d + pdf1_cb*pdf2_c + pdf1_sb*pdf2_s + pdf1_bb*pdf2_b);
//  sigma_spin1_qq   += sigmahat_spin1_qq* (pdf_u*pdf_ub*Qu2 + pdf_d*pdf_db*Qd2 + pdf_c*pdf_cb*Qu2 + pdf_s*pdf_sb*Qd2 + pdf_b*pdf_bb*Qd2);
  
  
  
  sigmahat_spin0_gq = ColAvg_qg * Flux * Jacobian * ME0sq_gq(pIn[0],pIn[1],pOut[0],pOut[1]);
  sigmahat_spin1_gq = ColAvg_qg * Flux * Jacobian * ME1sq_gq(pIn[0],pIn[1],pOut[0],pOut[1]);
  sigma_spin0_gq    = sigmahat_spin0_gq* (pdf2_u*pdf1_g*mov2_u + pdf2_d*pdf1_g*mov2_d+ pdf2_c*pdf1_g*mov2_c + pdf2_s*pdf1_g*mov2_s + pdf2_b*pdf1_g*mov2_b);
  sigma_spin1_gq    = sigmahat_spin1_gq* (pdf2_u*pdf1_g + pdf2_d*pdf1_g + pdf2_c*pdf1_g + pdf2_s*pdf1_g + pdf2_b*pdf1_g);
  
  sigmahat_spin0_qg = ColAvg_qg * Flux * Jacobian * ME0sq_gq(pIn[1],pIn[0],pOut[0],pOut[1]);
  sigmahat_spin1_qg = ColAvg_qg * Flux * Jacobian * ME1sq_gq(pIn[1],pIn[0],pOut[0],pOut[1]);
  sigma_spin0_qg    = sigmahat_spin0_qg* (pdf1_u*pdf2_g*mov2_u + pdf1_d*pdf2_g*mov2_d+ pdf1_c*pdf2_g*mov2_c + pdf1_s*pdf2_g*mov2_s + pdf1_b*pdf2_g*mov2_b);
  sigma_spin1_qg    = sigmahat_spin1_qg* (pdf1_u*pdf2_g + pdf1_d*pdf2_g + pdf1_c*pdf2_g + pdf1_s*pdf2_g + pdf1_b*pdf2_g);
  
  
  
  sigmahat_spin0_qbg= ColAvg_qg * Flux * Jacobian * ME0sq_qbg(pIn[0],pIn[1],pOut[0],pOut[1]);
  sigmahat_spin1_qbg= ColAvg_qg * Flux * Jacobian * ME1sq_qbg(pIn[0],pIn[1],pOut[0],pOut[1]);
  sigma_spin0_qbg   = sigmahat_spin0_qbg* (pdf2_g*pdf1_ub*mov2_u + pdf2_g*pdf1_db*mov2_d+ pdf2_g*pdf1_cb*mov2_c + pdf2_g*pdf1_sb*mov2_s + pdf2_g*pdf1_bb*mov2_b);
  sigma_spin1_qbg   = sigmahat_spin1_qbg* (pdf2_g*pdf1_ub + pdf2_g*pdf1_db + pdf2_g*pdf1_cb + pdf2_g*pdf1_sb + pdf2_g*pdf1_bb);
 
  sigmahat_spin0_gqb= ColAvg_qg * Flux * Jacobian * ME0sq_qbg(pIn[1],pIn[0],pOut[0],pOut[1]);
  sigmahat_spin1_gqb= ColAvg_qg * Flux * Jacobian * ME1sq_qbg(pIn[1],pIn[0],pOut[0],pOut[1]);
  sigma_spin0_gqb   = sigmahat_spin0_gqb* (pdf1_g*pdf2_ub*mov2_u + pdf1_g*pdf2_db*mov2_d+ pdf1_g*pdf2_cb*mov2_c + pdf1_g*pdf2_sb*mov2_s + pdf1_g*pdf2_bb*mov2_b);
  sigma_spin1_gqb   = sigmahat_spin1_gqb* (pdf1_g*pdf2_ub + pdf1_g*pdf2_db + pdf1_g*pdf2_cb + pdf1_g*pdf2_sb + pdf1_g*pdf2_bb);
  
  
//printf(" ehat %e %e \n",CMSEnergy,q);
//printf(" flux %e \n",Flux);
//printf(" jaco %e \n",Jacobian);
//printf(" ptgl %e \n",pTg);
//printf(" m0sq %e \n",ME0sq(pIn[0],pIn[1],pOut[0],pOut[1])*mov2_u*ColAvg_qq);
//printf(" m1sq %e \n",ME1sq(pIn[0],pIn[1],pOut[0],pOut[1])*ColAvg_qq);
//printf(" pdfu d %e %e\n",pdf_u/x1,pdf_d/x1);
//printf(" pdfb   %e %e\n",pdf_ub/x2,pdf_db/x2);
//printf(" p0gl %e \n",pOut[1][0]);
//printf(" pzgl %e \n",pOut[1][3]);
//exit(0);
  // P_left( qbar(x1) ) * P_right( q(x2 )



 f[0] = fbGeV2 * (sigma_spin0_qqb+sigma_spin0_qbq);
 f[1] = fbGeV2 * (sigma_spin1_qqb+sigma_spin1_qbq);

 f[0] += fbGeV2 * (sigma_spin0_qg+sigma_spin0_gq);
 f[1] += fbGeV2 * (sigma_spin1_qg+sigma_spin1_gq);

 f[0] += fbGeV2 * (sigma_spin0_gqb+sigma_spin0_qbg);
 f[1] += fbGeV2 * (sigma_spin1_gqb+sigma_spin1_qbg);


 int iBin1= WhichBin(1,pTg);
 int iBin2= WhichBin(2,yg);
 
 Histogram_spin0[0][iBin1] += f[0] * (*weight);
 Histogram_spin1[0][iBin1] += f[1] * (*weight);
 Histogram_spin0[1][iBin2] += f[0] * (*weight);
 Histogram_spin1[1][iBin2] += f[1] * (*weight);


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
  int i,j;
  double estim[FUNCTIONS];   /* estimators for integrals                     */
  double std_dev[FUNCTIONS]; /* standard deviations                          */
  double chi2a[FUNCTIONS];   /* chi^2/n                                      */
  double reg[2*DIMENSION];   /* integration domain                           */


   int iord=2;
   double one=1.0;
   double MZ=91.19*GeV, as_mz=0.13;
   double m_c=(1.3)*GeV, m_b=(4.2)*GeV, m_t=172.0*GeV;
   initalphasr0_(&iord, &one, &MZ, &as_mz, &m_c, &m_b, &m_t);
 

  //  initializing the integration range, always from 0.0 to 1.0
  for (i=0; i<DIMENSION; i++) {
    reg[i] = 0.0;
    reg[i+DIMENSION] = 1.0;
  }

  // setting parameters for vegas integrator
  unsigned long ncall=100000;


  int init=0;
  int itmx=5;
  // calling vegas integrator
  vegas(reg, DIMENSION, CrossSection, init, ncall/2, itmx,0x0001 | 0x0002 | 0x0004, FUNCTIONS, 0, 1, estim, std_dev, chi2a);
   
  for (j=0; j<NUMHIST; j++) {
   for (i=0; i<NUMBINS; i++) { 
      Histogram_spin0[j][i]=0; 
      Histogram_spin1[j][i]=0;
   };
  };

  init=1;
  itmx=10;
  // calling vegas integrator
  vegas(reg, DIMENSION, CrossSection, init, ncall, itmx,0x0001 | 0x0002 | 0x0004, FUNCTIONS, 0, 1, estim, std_dev, chi2a);
  


  double LowVal,BinSize;

  for (j=0; j<NUMHIST; j++) {
   for (i=0; i<NUMBINS; i++) { 
	  switch(j) {
	    case 0: LowVal=20.0*GeV; BinSize=10*GeV; break;  // pT_glu
	    case 1: LowVal=-5.0; BinSize=0.5;        break;  // y_glu
	  };
      printf("0  %2i  %6.2f   %e \n",j+1,(i)*BinSize+LowVal,Histogram_spin0[j][i]/itmx);
   };
  };
  printf("\n\n");
  for (j=0; j<NUMHIST; j++) {
   for (i=0; i<NUMBINS; i++) { 
	  switch(j) {
	    case 0: LowVal=20.0*GeV; BinSize=10*GeV; break;  // pT_glu
	    case 1: LowVal=-5.0; BinSize=0.5;        break;  // y_glu
	  };
      printf("1  %2i  %6.2f   %e \n",j+1,(i)*BinSize+LowVal,Histogram_spin1[j][i]/itmx);
   };
  };
	

  return(0);
}


