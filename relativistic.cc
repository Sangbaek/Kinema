#include <math.h>
#include <iostream>
#include "const.h"
#include "relativistic.h"

//
// Function name : Energy2Mom(double E, double M)
//
// Description   : calculate momentum from energy
// Input         : E, M
// Return        : p
//
double 
Energy2Mom(double E, double M){

  return sqrt(E*E - M*M);

}

//
// Function name : Mom2Energy(double p, double M)
//
// Description   : calculate total energy from momentum
// Input         : p, M
// Return        : E
//
double 
Mom2Energy(double p, double M){

  return sqrt(M*M + p*p);

}


//
// Function name : beta(double p, double E)
//
// Description   : calculate beta
// Input         : Momentum, Total Energy
// Return        : beta
//
double
calcBeta(double p, double E) {

 return p/E ; 

}


//
// Function name : beta(double p, double T, double M)
//
// Description   : calculate beta
// Input         : Momentum, Kinetic Energy, Mass
// Return        : beta
//
double
calcBeta(double p, double T, double M) {

 return p/(T+M) ; 

}


//
// Function name : calcgamma(double beta)
//
// Description   : calculate gamma
// Input         : beta
// Return        : gamma
//
double
calcgamma(double beta) {

 return 1/sqrt(1-beta*beta) ; 

}


//
// Function name : Mandelstam_t
//
// Description   : calculate Mandelstam t
// Input         : Ma, M1, Ea, E1, theta
// Return        : t
//
double
Mandelstam_t(double Ma, double M1, double Ea, double E1, double theta) {

  double pa = Energy2Mom(Ea, Ma);
  double p1 = Energy2Mom(E1, M1);

  double t = Ma*Ma + M1*M1 - 2*Ea*E1 + 2*pa*p1*cos(theta);

  return t ;

}

//
// Function name : Mandelstam_s
//
// Description   : calculate Mandelstam s
// Input         : Ma, Mb, Ea, Eb, theta
// Return        : s
//
double
Mandelstam_s(double Ma, double Mb, double Ea, double Eb, double theta) {

  double pa = Energy2Mom(Ea, Ma);
  double pb = Energy2Mom(Eb, Mb);

  double s = Ma*Ma + Mb*Mb + 2*Ea*Eb - 2*pa*pb*cos(theta);

  return s;

}
