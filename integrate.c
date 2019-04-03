
#include <math.h>
#include <stdio.h>
#include "vegas.h"


#define DIMENSION 2
#define FUNCTIONS 1



// test integrand: 2-dimensional integration over da1 and da2 with boundaries a1=[-10,15] and a2=[-1,1]
//                 integrand is a1^2+3*a2 
//                 
//  checking the result in mathematica: 
//  
// Integrate[ a1^2+3*a2, {a1,-10,15},{a2,-1,1} ] // N                                                                                                                            
// Out[1]= 2916.67
// 
void TestIntegrand(double x[DIMENSION], double f[FUNCTIONS])
{  double a1,a2,jac_a1,jac_a2;


// x[0] and x[1] run from 0..1
// mapping for the correct a1,a2 boundaries

   a1 = 25.0*x[0] - 10.0;
   a2 = 2*x[1] - 1.0;

// jacobians for this transformation
   jac_a1 = 25.0;
   jac_a2 = 2.0;


//  defining the integrand function
    f[0] = a1*a1 + 3*a2;


//  multiplying the jacobians
    f[0] = f[0] * jac_a1 * jac_a2;


  return;
}




// main program
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
  int itmx=5;
  
  // calling vegas integrator
  vegas(reg, DIMENSION, TestIntegrand, init, ncall, itmx, NPRN_INPUT | NPRN_RESULT, 1, 0, 1, estim, std_dev, chi2a);
  

  return(0);
}


