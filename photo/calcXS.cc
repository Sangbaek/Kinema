#include <iostream.h>
#include <fstream.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "const.h"
#include "relativistic.h"
#include "calcXS.h"


//
// Classx name  : 
// Method name : NuclearCoherentFF(double A, int Z, double q, double tq){
//
// Description : calculate nuclear coherent cross section using Cornell's parameters
//             : 
// Input       : double A, double t1, double E0, double Fnc
// Return      : xs [micro barn]
//
double 
NuclearCoherentXS(double A, double t1, double E0, double Fnc){

  // constants obtained from PRL33, 1400 (1974) //
  double bn  = CORNELL_CONST_bn;
  double L2  = 100 * E0*M2G * E0*M2G ; // [microbarn]
  double C   = bn * L2 ;
  double st1 = sin(t1);


  double xs = C * A*A * Fnc * Fnc * st1*st1 ; //[microbarn]

  return xs;

}





//
// Classx name  : 
// Method name : PrimakoffXS
//
// Description : calculate nuclear coherent cross section using Cornell's parameters
//             : 
// Input       : double Fc, double p1, double M1, int Z, double E0, double q, double Fc, double t1
// Return      : xs [mcb/sr]
//
double 
PrimakoffXS(double p1, double M1, int Z, double E0, double q, double Fc, double t1){


  double Gamma = CORNELL_CONST_bc ; // [MeV]
  double st1 = sin(t1);
  double E1 = Mom2Energy(p1,M1);
  double beta = calcBeta(p1, E1);

  double xs  = Gamma * 8 * ALPHA * Z*Z / (M1*M1*M1) ;
  xs *= beta*beta*beta * E0*E0*E0*E0 ;
  xs /= (q*q*q*q) ;
  xs *= Fc*Fc * st1*st1 ; // 1/[MeV]^2/[sr]

  xs *= hc*hc * fm2mcb  ; // [mcb/sr]

  return xs;

}


//
// Classx name  : 
// Method name : InterferenceXS
//
// Description : calculate interference between nuclear coherent and Primakoff
//             : 
// Input       : (double Xc, double Xnc){
// Return      : xs [mcb/sr]
//
double 
InterferenceXS(double Xc, double Xnc){

  // nuclear coherent and Coulomb phase shift [rad]
  double phi = CORNELL_CONST_phi ;

  // unit : [mcb/sr]
  double xs = 2 * sqrt(Xc*Xnc) * cos(phi) ;

  return xs;

}





//
// Classx name  : 
// Method name : NuclearIncoherentXS(double Finc);
//
// Description : calculate nuclear incoherent cross section
//             : 
// Input       : double Finc
// Return      : xs [mcb/sr]
//
double 
NuclearIncoherentXS(double Finc){

  // nuclear incoherent constant from Cornell Paper
  double bb = CORNELL_CONST_bb ;
 
  // unit : [mcb/sr]
  double xs = bb * Finc * Finc ;

  return xs;

}


//
// Classx name  : 
// Method name : TotalXS
//
// Description : calculate total cross section
//             : 
// Input       : double Xnc, double Xc, double Xint, double Xinc
// Return      : xs [mcb/sr]
//
double 
TotalXS(double Xnc, double Xc, double Xint, double Xinc){

  double xs = Xnc + Xc + Xint + Xinc ;

  return xs;

}







