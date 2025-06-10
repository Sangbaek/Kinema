//                                                                    -*-c++-*-
// 
//  Lab2CM.cc
//

#include <iostream>
#include <stdio.h>
#include <math.h>
#include "const.h"
#include "Lab2CM.h"
#include "relativistic.h"

//
// Function name : VPLab2CM(double  beta, double  q,    double  omega 
//                                        double &q_cm, double &omega_cm)
//
// Description   : calculate q and omega in CM
// Input         : beta, q, omega
// Return        : q_cm, omega_cm
//
int 
VPLab2CM(double beta, double q, double omega, double &q_cm, double &omega_cm) {

  double gamma = calcgamma(beta);

  q_cm         = gamma*(q - beta*omega);
  omega_cm     = gamma*(omega - beta*q);
  
  return 0;
}

void
printVPLab2CM(double q_cm, double omega_cm){

  printf("\n");
  printf(" Center of Mass:\n");
  printf("\t Momentum Transfer    q: %4.1f [MeV/c]\t %2.3f [fm^-1]\n:", 
	 q_cm, q_cm/hc); 
  printf("\t Energy Transfer  omega: %4.1f [MeV/c]                \n:",   
	 omega_cm); 
  printf("\n");

}



//
// Function name : HadronEnergyCM(double W, M1, M2)
//
// Description   : calculate hadron (Mass1) energy in CM
// Input         : W, M1, M2
// Return        : E1_cm
//
double 
HadronEnergyCM(double W, double M1, double M2) {

  return (W*W + M1*M1 - M2*M2)/2/W;

    } //end of HadronEnergyCM




//
// Function name : calcPhotonEquivalentEnergyCM(double W, double M)
//
// Description   : calculate photon equivalent energy in CM
// Input         : W, M
// Return        : kgamma_cm
//
double 
calcPhotonEquivalentEnergyCM(double W, double Mt){	

  return (W*W - Mt*Mt)/2/W ;

}

//
// Function name : calcbetaCM(double q, double omega, double M)
//
// Description   : calculate beta in CM
// Input         : q, W
// Return        : betaCM
//
double 
calcbetaCM(double q, double W){

  return q/W;

}


//
// Function name : calcbetaCM(double q, double omega, double M)
//
// Description   : calculate beta in CM
// Input         : q, omega, M
// Return        : betaCM
//
double 
calcbetaCM(double q, double omega, double M){

  return q/(omega + M);

}


//
// Function name : calcgammaCM(double omega, double W, double M);
//
// Description   : calculate gamma in CM
// Input         : omega, W, M
// Return        : gammaCM
//
double 
calcgammaCM(double omega, double W, double M){

   return (omega + M)/W;

}



//
// Function name : calcJacobianLab
//
// Description   : calculate Jacobian from given theta_pq_Lab angle
// Input         : gammaCM, pLab, pCM, betaCM, ELab, tpqLab
// Return        : Jacobian(Omega_pCM/Omega_pLab)
//
double 
calcJacobianLab(double gammaCM, double pLab, double pCM, 
		double betaCM,  double ELab, double tpqLab){

  double ctpqLab = cos(tpqLab);

  return 1/gammaCM * pLab/pCM * 1/(1 - betaCM*ELab/pLab*ctpqLab) ;

}

//
// Function name : calcJacobianCM
//
// Description   : calculate Jacobian from given theta_pq_CM angle
// Input         : gammaCM, pLab, pCM, betaCM, ECM, tpqCM
// Return        : Jacobian(Omega_pCM/Omega_pLab)
//
double 
calcJacobianCM(double gammaCM, double pLab, double pCM, 
	       double betaCM,  double ECM,  double tpqCM){

  double ctpqCM  = cos(tpqCM);
  double pLabCM3 = pLab/pCM * pLab/pCM * pLab/pCM ;
  
  printf("theta_pq in CM: %g, P Ener. in CM: %g\n",tpqCM, ECM);
  printf("P mom. in Lab: %g, P mom. in CM: %g\n",pLab,pCM);
  printf("gammaCM: %g, betaCM: %g\n",gammaCM,betaCM);

  return 1/gammaCM * pLabCM3 * 1/(1 + betaCM*ECM/pCM*ctpqCM) ;

}



//
// Function name : calcJacobianLab2CM
//
// Description   : calculate Jacobian from Lab to CM frame
// Reference     : E.Byckling and K.Kajantie, "Particle Kinematics", John Willy & Sons
//               : Section III
// Input         : double tCM, double betaCM, double beta1CM
// Return        : Jacobian = Omega_CM/Omega_Lab
//
double 
calcJacobianLab2CM(double tCM, double betaCM, double beta1CM){


  double stCM  = sin(tCM);
  double ctCM  = cos(tCM);
  
  double gammaCM = calcgamma(betaCM);
  double gCM   = betaCM/beta1CM;

  // PPCM2 = P_Lab*P_Lab/(P_CM*P_CM) 
  // PPCM  = P_Lab/P_CM = sin(tCM)/sin(tLab)
  double PPCM2 = gammaCM*gammaCM*(ctCM+gCM)*(ctCM+gCM)+stCM*stCM;

  // dtLab_dtCM = dtLab/dtCM
  double dtLab_dtCM = gammaCM*(1+gCM*ctCM)/PPCM2;

  // J = sin(tCM)/sin(tLab) * dtLab/dtCM 
  double J = sqrt(PPCM2) * dtLab_dtCM ;
   
  return J ;

}





















