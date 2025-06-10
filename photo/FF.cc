#include <iostream.h>
#include <fstream.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "const.h"
#include "FF.h"
#include "Math.h"

//
// Classx name  : Cornell
// Method name : NuclearCoherentFF(double A, int Z, double q, double tq){
//
// Description : calculate Cornell's nuclear coherent form factor
//             : See documentation
//             : " Parametrize Sergey Gevorkyan's Nuclear Coherent Form Factor 
//             :   using Cornell Parametrization" , I. Nakagawa, June 29, 2004
// Input       : double A, int Z, double q, double tq
// Return      : Fnc
//

double 
Cornell::NuclearCoherentFF(double A, int Z, double q, double tq){


  double a, b;
  a=b=0;
  switch (Z){
  case 92 : // Uranium
      a = 0.42 ;
      b = -320 ;
      break;
  case 82 :// Lead substituted by Uranium
      a = 0.42 ;
      b = -320 ;
      break;
  case 50 :// Sn
      a = 0.49 ;
      b = -212 ;
  case 29 ://Cupper
      a = 0.56 ;
      b = -150 ;
      break ;
  case 4 ://Beryllium
      a = 0.77 ;
      b = -60 ;
      break;
  case 6 ://Carbon substituted by Beryllium
      a = 0.77 ;
      b = -60 ;
      break;
  }
  if (a*b==0) {
    cerr << "NuculearCoherent: Error, A="<<A<<" and Z=" << Z << " is not in database" << endl;
    exit(-1);
  }

  double qt   = q*sin(tq);
  double Fnc2 = a * exp ( b*qt*M2G*qt*M2G );
  double Fnc  = sqrt(Fnc2);

  return Fnc;

}



//
// Classx name : Cornell
// Method name : CoulombFF(double A, int Z, double q, double tq){
//
// Description : calculate Cornell's Coulomb form factor
//             : 
// Input       : double A, int Z, double q, double tq
// Return      : Fc
//

double 
Cornell::CoulombFF(double A, int Z, double q, double tq){

  double b=0;
  switch (Z){
  case 92 : // Uranium
      b = -520 ;
      break;
  case 82 :// Lead substited by Uranium
      b = -520 ;
      break;
  case 50 :// Sn
      b = -316 ;
      break;
  case 29 ://Cupper
      b = -210 ;
      break ;
  case 4 ://Beryllium
      b = -72 ;
      break;
  case 6 ://Carbon substituted by Beryllium
      b = -72 ;
      break;
  }
  if (b==0) {
    cerr << "CoulombFF: Error, A="<<A<<" and Z=" << Z << " is not in database" << endl;
    exit(-1);
  }


  double qt  = q*sin(tq);
  double Fc2 = exp ( b*qt*M2G*qt*M2G );
  double Fc  = sqrt(Fc2);

  return Fc;

}



//
// Classx name : Cornell
// Method name : Cornell::NuclearIncoherentFF(double A, int Z, double theta)
//
// Description : calculate Cornell's Nuclear Incoherent form factor
//             : 
// Input       : double A, int Z, double theta
// Return      : Fi
//

double 
Cornell::NuclearIncoherentFF(double A, int Z, double theta){

  double Fi2  = 1.0 * pow(A,0.75); 
  double Fi   = sqrt(Fi2);

  return Fi;

}


//
// Classx name : Sergey
// Method name : NuclearCoherentFF(double A, int Z, double q, double tq){
//
// Description : calculate Sergey's nuclear coherent form factor
//             : Parametrized by Cornell's parametrization
//             : " Parametrize Sergey Gevorkyan's Nuclear Coherent Form Factor 
//             :   using Cornell Parametrization" , I. Nakagawa, June 29, 2004
// Input       : double A, int Z, double q, double tq
// Return      : Fnc
//

double 
Sergey::NuclearCoherentFF(double A, int Z, double q, double tq){

  double a, b;
  a=b=0;
  switch (Z){
  case 82 :// Lead 
      a = 0.497 ;
      b = -302.7 ;
      break;
  case 50 :// Sn
      a = 0.751 ;
      b = -201.1 ;
      break;
  case 6 :// Carbon
      a = 0.768 ;
      b = -55.20 ;
      break;
  }
  if (a*b==0) {
    cerr << "NuculearCoherent: Error, A="<<A<<" and Z=" << Z << " is not in Sergey's db" << endl;
    exit(-1);
  }

  double qt   = q*sin(tq);
  double Fnc2 = a * exp ( b*qt*M2G*qt*M2G );
  double Fnc  = sqrt(Fnc2);

  return Fnc;

}


//
// Classx name : Cascade
// Method name : Cascade::NuclearIncoherentFF(double A)
//
// Description : calculate Cascade Model Nuclear Incoherent form factor
//             : 
// Input       : double A, int Z, double tq
// Return      : Fi
//

double 
Cascade::NuclearIncoherentFF(double A, int Z, double theta){

  // File Operation 
  char * dir = "/home/itaru/PrimEx/Model/Cascade/dat/";

  char File[255];
  switch (Z){
  case 6 :// 12C
    sprintf(File,"Incoherent_12C_5.0GeV_FF.dat");
    break;
  case 82 :// Pb
    sprintf(File,"Incoherent_Pb_5.6GeV_FF.dat");
    break;
  }

  char DataFile[100] = "";
  strcat(DataFile,     dir);
  strcat(DataFile,    File);
  // cerr << "Incoherent Model = Cascade, database=" << DataFile << endl;

  ifstream in_file;
  in_file.open(DataFile);
  if (in_file.fail())
    {
      cerr <<"Cascade::NuclearIncoherentFF(): File opening Error. \n "<< DataFile
           <<" is not found."<<endl;
    exit(-1);
   }


  // Scan Model database until finding interporation point of theta_pi angle
  int ch;
  int j=0;
  int OutOfRange=1;
  double x[MAXDATA],y[MAXDATA];
  while ( ( ch = in_file.peek()) != EOF )
    {
      in_file >> x[j] >> y[j] ;
      if (theta*r2d<x[j]) {
	OutOfRange = 0;
	break;
      }
      ++j;
    }

  // Check the model range
  if (OutOfRange) {
    cerr << "Cascade::NuclearIncoherentFF: Error. theta is out of database coverage " << endl;
    exit(-1);
  }


  // Perform Linear Interporation 
  double Fi = LinearInterporate(x[j], y[j], x[j-1], y[j-1], theta*r2d);

  return Fi;

}
