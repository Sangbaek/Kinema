//                                                                    -*-c++-*-
// 
//  xs.cc
//

#include <iostream>
#include <stdio.h>
#include <math.h>
#include "const.h"
#include "xs.h"
#include "calc.h"
#include "Lab2CM.h"
#include "relativistic.h"

using namespace std;

//
// calculate virtural photon variables
//
int  
VP::calcVirtualPhoton(double Ee, double Eep, double q, double Q2, 
		      double  te, double W, double Mt, double M1, double M2){ 

  me = MASS_ELECTRON;
  k  = Energy2Mom(Ee,  me);
  kp = Energy2Mom(Eep, me);

  if (te*rad2deg < 5) {
    kgamma   = calcPhotonEquivalentEnergy(k, kp, te);
    epsln    = calcVPhPolarization(k, kp, te, Q2, kgamma);
    
    cerr << " " << endl;
    cerr << " The current program uses different formalism for te < 5 deg.\a "
	 << endl;
    cerr << " " << endl;

  }else{
    kgamma   = calcPhotonEquivalentEnergy(W, Mt);
    epsln    = calcVPhPolarization(q, Q2, te);
  }

  /*  get phase space factor Ppi/k_gamma  */
  Ppikgamma = calcPpikgamma(W, Mt, M1, M2);

  /*  get the fraction Q/omega_cm  */
  Qomega_cm = calcQomegacm(Q2, q, W, Mt); 
  Qomega_cm2= Qomega_cm* Qomega_cm;

  /*  get the fraction q_cm/omega_cm  */
  qCMomegaCM = calcqcmomegacm(Q2, q, W, Mt); 

  /*  get the fraction Q/q_cm  */
  QqCM       = calcQqcm(Q2, q, W, Mt);

  epslnTL  = calcVPhKinVarTL(epsln);
  epslnTT  = epsln;
  epslnTLp = calcVPhKinVarTLp(epsln);
  epslnTTp = calcVPhKinVarTTp(epsln);

  /* rhoTL = sqrt( epslnL*(1+epsln) ) */
  epslnL   = Qomega_cm2 * epsln;
  rhoTL    = Qomega_cm  * epslnTL;
  rhoTT    = epslnTT;
  rhoTLp   = Qomega_cm  * epslnTLp;
  rhoTTp   = epslnTTp;

  Gamma    = calcVPhFlux(Ee, Eep, kgamma, Q2, epsln);

  Jacobian = calcJacobian(Mt, Ee, Eep, W);

  //  printVPh();

  return 0;
}

//
// Class name  : 
// Method name : calcPpikgamma(double W, double Mt, double M1, double M2)
//
// Description : calculate phase space factor: fraction of outgoing pion 
//             : momentum and the virtual photon equivalent energy in CM
// Input       : double W, double Mt, double M1, double M2
// Return      : p1_cm/kgamma_cm
//

double
VP::calcPpikgamma(double W, double Mt, double M1, double M2){

  /* calculate CM variables */
  double E1_cm    = HadronEnergyCM(W, M1, M2);
  double p1_cm    = Energy2Mom(E1_cm, M1);

  kgamma_cm   = calcPhotonEquivalentEnergyCM(W, Mt);
  
  return p1_cm/kgamma_cm;

}


//
// Class name  : 
// Method name : calcQomegacm(double Q2, double q, double omega, double Mt)
//
// Description : calculate the fraction of Q2 and omega_cm
//
// Input       : double Q2, double q, double omega, double Mt
// Return      : Q2/omega_cm
//

double
VP::calcQomegacm(double Q2, double q, double W, double Mt){

  double q_cm, omega_cm;
  double omega = sqrt(W*W + q*q) - Mt;
  double beta     = calcBeta(q, omega, Mt);
  VPLab2CM(beta, q, omega, q_cm, omega_cm);

  //  printVPLab2CM( q_cm, omega_cm);
  return sqrt(Q2)/omega_cm;

}

//
// Class name  : 
// Method name : calcqcmomegacm(double Q2, double q, double omega, double Mt)
//
// Description : calculate the fraction of q_cm and omega_cm
//
// Input       : double Q2, double q, double omega, double Mt
// Return      : q_cm/omega_cm
//

double
VP::calcqcmomegacm(double Q2, double q, double W, double Mt){

  double q_cm, omega_cm;
  double omega = sqrt(W*W + q*q) - Mt;
  double beta     = calcBeta(q, omega, Mt);
  VPLab2CM(beta, q, omega, q_cm, omega_cm);

  //  printVPLab2CM( q_cm, omega_cm);
  return q_cm/omega_cm;

}

//
// Class name  : 
// Method name : calcqcmomegacm(double Q2, double q, double omega, double Mt)
//
// Description : calculate the fraction of q_cm and omega_cm
//
// Input       : double Q2, double q, double omega, double Mt
// Return      : Q/q_cm
//

double
VP::calcQqcm(double Q2, double q, double W, double Mt){

  double q_cm, omega_cm;
  double omega = sqrt(W*W + q*q) - Mt;
  double beta     = calcBeta(q, omega, Mt);
  VPLab2CM(beta, q, omega, q_cm, omega_cm);

  return sqrt(Q2)/q_cm;

}



// 
// Vertual Photon Polarization
//
double 
VP::calcVPhPolarization(double q, double Q2, double te) {

  if ( Q2 ) 
    return  1/( 1 + 2*q*q/Q2*tan(te/2)*tan(te/2) );
  else 
    cerr << "Error in calcVPPolarization() : zero value in Q2" << endl;
    return 0;
}


//
// Class name  : VP
// Method name : calcVPhPolarization(double, double, double, double, double) 
//
// Description : calculate virtual photon polarization without me->0 apprx.
//             : this is for forward angle scattering
// Input       : k, kp 
//             : te, kgamma, Q2
// Return      : kgamma
//
double 
VP::calcVPhPolarization(double k, double kp, double te, 
			double Q2, double kgamma) {

  double x,y;
  if ( te*k*kp ) {
    y = Q2*kgamma*kgamma ;
    x = 2*k*k*kp*kp*sin(te)*sin(te);
    return  1/( 1 + y/x);
  }  else {
    cerr << "Error in calcVPPolarization() : zero value either in k, k', te" 
	 << endl;
    return 0;
  }

}


// 
// kinemtic variables for TL
//
double 
VP::calcVPhKinVarTL(double epsln){

  return sqrt( 2*epsln*(1+epsln) );

}

// 
// kinemtic variables for TL'
//
double 
VP::calcVPhKinVarTLp(double epsln){

  return sqrt( 2*epsln*(1-epsln) );

}


// 
// kinemtic variables for TT'
//
double 
VP::calcVPhKinVarTTp(double epsln){

  return sqrt( 1 - epsln*epsln );

}


//
// Class name  : VP
// Method name : calcPhotonEquivalentEnergy(double W, double Mt)
//
// Description : calculate photon equivalent energy, using me->0 apprx.
//             : this cannot be used for the forward angle scattering
// Input       : double W, double Mt
// Return      : kgamma
//

double 
VP::calcPhotonEquivalentEnergy(double W, double Mt){	

  return (W*W - Mt*Mt)/2/Mt ;
}

//
// Class name  : 
// Method name : calcPhotonEquivalentEnergy(double Ee, double Eep, double te)
//
// Description : calculate photon equivalent energy, not using me->0 apprx.
//             : this is for forward angle scattering
// Input       : double k, double kp, double te
// Return      : kgamma
//

double 
VP::calcPhotonEquivalentEnergy(double k, double kp, double te){	

  return sqrt(k*k + kp*kp - 2*k*kp*cos(te)) ;
}

//
// calculate virtual photon flux in unit of [MeV]^-1
// 
double 
VP::calcVPhFlux(double Ee, double Eep, double kgamma, double Q2, double epsln){

  
  if ( Ee*Q2*(1-epsln) ) 
    return ALPHA/2/M_PI/M_PI * Eep/Ee * kgamma/Q2 * 1/(1-epsln);
  else {
    cerr << "Error in calcVPhFlux() : devited by zero" << endl;
    return 0;
  }
}

//
// calculate Transformation Jacobian from (te,Eep) -> (W,Q2)
//
double 
VP::calcJacobian(double Mt, double Ee, double Eep, double W){

  if ( W ) 
    return 2*Mt*Ee*Eep/W ;
  else 
    cerr << "Error in calcJacobian() : zero value in W" << endl;
  return 0;

}

//
// returns Gamma
//
double
VP::getGamma(){

  return Gamma;
}

//
// print Virtual Phoson calcuated results
//
void   
VP::printVPh(){

  printf("\n");
  printf(" Virtual Photon:\n");
  printf("\t Ppi/kgamma       : %2.3f \n"           ,        Ppikgamma);
  printf("\t q_cm/omega_cm    : %2.3f "             ,       qCMomegaCM);
  printf("\t Q/q_cm           : %2.3f \n"           ,             QqCM);
  printf("\t Q/omega_cm       : %2.3f "             ,        Qomega_cm);
  printf("\t (Q/omega_cm)^2   : %2.3f \n"           ,       Qomega_cm2);
  printf("\t kgamma      [MeV]: %4.1f "             ,           kgamma);
  printf("\t kgamma_cm   [MeV]: %4.1f \n"           ,        kgamma_cm);
  printf("\t Gamma  [MeV/c]^-1: %2.3e \n"           ,            Gamma);
  printf("\t epsilon          : %2.3f "             ,            epsln); 
  printf("\t epsilonL (=rhoL) : %2.3f \n"           ,           epslnL); 
  printf("\t v_TL             : %2.3f "             ,          epslnTL);
  printf("\t rhoTL            : %2.3f \n"           ,            rhoTL); 
  printf("\t v_TT             : %2.3f "             ,          epslnTT);
  printf("\t rhoTT            : %2.3f \n"           ,            rhoTT);
  printf("\t v_TL'            : %2.3f "             ,         epslnTLp);
  printf("\t rhoTL'           : %2.3f \n"           ,           rhoTLp);
  printf("\t v_TT'            : %2.3f "             ,         epslnTTp);
  printf("\t rhoTT'           : %2.3f \n"           ,           rhoTTp);
  printf("\t Jacobian[MeV/c]^2: %8.1f \n"           ,         Jacobian);
  printf("\n");

} // end-of-printVPh(){

//
// returns Virtual Photon variables
//
void
VP::getVP(double &Gamma_, double &epsln_, double &epslnTT_, 
	  double &epslnTL_, double &epslnTLp_, double &Jacobian_){

  Gamma_    = Gamma;
  epsln_    = epsln;
  epslnTT_  = epslnTT;
  epslnTL_  = epslnTL;
  epslnTLp_ = epslnTLp;
  Jacobian_ = Jacobian;

  return ;

} // end-of-getVP()


//
// print Virtual Phoson calcuated results for flux 
//
void   
VP::printVPhFlux(double te){

  printf("\n");
  printf("\t Scattering Angle         : %2.3f [deg.]\n", te*rad2deg);
  printf("\t Gamma                    : %2.3e [MeV/c]^-1\n", Gamma);
  printf("\t Gamma*pi^2/90*sin(the_e) : %2.3e [GeV/c*deg]^-1\n", 
	 Gamma*M_PI*M_PI/90*sin(te)*1e3);
  printf("\n");

}

