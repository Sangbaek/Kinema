//
//  kinema.cc
//

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include "kinema.h"
#include "const.h"
#include "calc.h"
#include "xs.h"
#include "relativistic.h"
#include "Lab2CM.h"

using namespace std;


// Initialize Print Level 
//PrintLeve=0 : Prints Output of Each Functions
//
int PrintLevel = 0;

void
Usage(){
//int InteractiveInput(); 

      cout << "\n Usage: kinema [options]" << endl;
      cout << " " << endl;
      cout << "\t -h \t Show this help    "<< endl;
      cout << "\t -P \t detect Proton" << endl;
      cout << "\t -p \t detect piplus" << endl;
      cout << "\t -e \t Incident electron [e]nergy in [MeV/c]" << endl;
      cout << "\t -f \t final electron energy    in [MeV/c]" << endl;
      cout << "\t -m \t elastic recoil momentum in [MeV/c]" << endl;
      cout << "\t -L \t Print Level=0 : Normal Output(def)" << endl;
      cout << "\t    \t       Level=1 : data only         " << endl;
      cout << "\t    \t       Level=2 : Virtual Photon Gamma data" << endl;
      cout << "\t -t \t Lab theta_e Angle         in [deg.]" << endl;
      cout << "\t -q \t Lab theta_q Angle         in [deg.]" << endl;
      cout << "\t -T \t Lab theta_pq Angle        in [deg.]" << endl;
      cout << "\t -C \t CM  theta_pq Angle        in [deg.]" << endl;
      cout << "\t -Q \t |4-Momentum Transfer|    in [GeV/c]^2" << endl;
      cout << "\t -W \t CM Energy Transfer         in [MeV/c]" << endl;
      cout << "\t -X \t Show Expamples                       " << endl;
      cout << "\t -v \t calculate virtual photon polarization" << endl;
      cout << "\t -G \t calculate VCS kinematics" << endl;
      cout << endl;

} //end-of-Usage

void
Example(){

      cout << endl;
      cout << " Examples:" << endl;
      cout << " \t elastic  :\t kinema -e 950 -t 28.74" << endl;
      cout << " \t elastic  :\t kinema -e 950 -q 62.72" << endl;
      cout << " \t elastic  :\t kinema -e 950 -m 500.00" << endl;
      cout << " \t (in)elastic:\t kinema -e 950 -t 28.74 -f 600 -v -L 2" 
	   << endl;
      cout << " \t inelastic:\t kinema -e 950 -t 28.74 -W 1232" << endl;
      cout << " \t inelastic:\t kinema -e 950 -Q 0.127 -W 1232" << endl;
      cout << " \t exclusive:\t kinema -e 950 -Q 0.127 -W 1232 -T 14.8 -P"<<endl; 
      cout << endl;

} // end-of-Example



int 
main(int argc, char *argv[]) {

  int opt;
  int vp=0;
  int Mode=0;
  int tpqLab=0;

  double q2, Rec, ERec, PRec;  // returns from calceKinema(elastic)
  double Mt, M1, M2, Ee, te, tpq, tpqCM, Q2, W;
  double omega, Eep, q, tq;   // returns from calcKinema(inelastic)
  //double p1[2], E1[2], p2[2], E2[2];   // returns from calcKinema(exclusive)

  Rec = ERec = PRec = 0;
  Mt = M1 = M2 = Ee = te = tpq = tpqCM = Q2 = W = 0;
  omega = Eep = q = q2 = tq = 0;
  //p1[0] = E1[0] = p2[0] = E2[0] = p1[1] = E1[1] = p2[1] = E2[1] = 0;

  // Default Particle 
//  Mt = M1 = 12.01*931.49432 ;
//  Mt = M1 = 1.00794*931.49432 ;
  Mt = M1 = MASS_PROTON ;
  M2 =MASS_PI0    ;

  // parse command line options
  while (EOF != (opt = getopt(argc, argv, "PGph?Xve:f:m:t:q:T:C:Q:W:L:"))) {
    switch (opt) {
    case 'e':
      Ee = atof(optarg);
      break;
    case 'f':
      Eep= atof(optarg);
      break;
    case 'm':
      PRec = atof(optarg);
      break;
    case 'L':
      PrintLevel = atoi(optarg);
      break;
    case 't':
      te = atof(optarg) * deg2rad;
      break;
    case 'q':
      tq = atof(optarg) * deg2rad;
      break;
    case 'T':
      tpq = atof(optarg) * deg2rad; 
      tpqLab = 1;
      break;
    case 'C':
      tpqCM = atof(optarg) * deg2rad;
      break;
    case 'p':
      Mt = MASS_PROTON;
      M1 = MASS_PIPLUS;
      M2 = MASS_NEUTRON;
      break;
    case 'G':
      M2=MASS_PHOTON;
      break;
    case 'Q':
      Q2 = atof(optarg) * 1e6;  // [MeV/c]^2
      break;
    case 'W':
      W = atof(optarg) ;
      break;
    case 'v':
      vp = 1;
      break;
    case 'X':
      Example();
      exit(0);
    case 'h':
    case '?':
    case '*':
      Usage();
      exit(0);
    }
  }

  VP kinema;

  if (!W) {

    if (!Eep) {

        cout <<"elastic"<<endl;
    // caluculate elastic kinematics
      if (te) 
	CalceKinema(Mt, Ee, te, Rec, Eep, q2, Q2, tq, ERec, PRec);
      else if (tq) 
	CalcpKinema(Mt, Ee, tq, Rec, Eep, q2, Q2, te, ERec, PRec);
      else if (PRec) 
	CalcPpKinema(Mt, Ee, PRec, Rec, Eep, q2, Q2, te, ERec, tq);
      else 
	cout << "Error: Input arguments are not proper" << endl;
      if (!PrintLevel) 
	PrinteKinema(Mt, Ee, te, Rec, Eep, q2, Q2, tq, ERec, PRec);
    }else{
    // caluculate kinematics: (Mt,Ee,te,Eep) -> (omega,W,q,Q2,tq)
    CalcWQKinema(Mt, Ee, te, Eep, omega, W, q, Q2, tq);
    if (!PrintLevel) 
	PrintKinema(Mt, Ee, Q2, W, omega, Eep, q, te, tq); 
    } // if (!Eep) 

  } else {

    // calculate inelastic kinematics
    if (te) {
cout<<"11"<<endl;
      // (Mt, Ee, te, W) -> (omega, Eep, q, Q2, tq)
      CalcteKinema(Mt, Ee, te, W, omega, Eep, q, Q2, tq);
    } else {
cout<<"12"<<endl;
      // ( Mt, Ee, Q2, W)  -> (omega, Eep, q, te, tq)
      CalcKinema(Mt, Ee, Q2, W, omega, Eep, q, te, tq); 
    }
      if (!PrintLevel) {
	PrintKinema(Mt, Ee, Q2, W, omega, Eep, q, te, tq); 
	PlainPrint(Ee, te, Eep, omega, W, tq, q, Q2, Mt, M1, M2) ;
      } else {
	PlainPrint(Ee, te, Eep, omega, W, tq, q, Q2, Mt, M1, M2) ;
      } //end PrintLevel

    // calculate hadron arm kinematics
    Hadron PiProd;
      if (tpqLab) {
          cout <<"tpq =="<<tpq<<endl;
        PiProd.CalcHKinema(q, omega, tpq, Mt, M1, M2);
	PiProd.CalcHKinemaCM(q, omega, tpq, W, Mt, M1, M2);
      } else {
cout<<"01 "<<tpqCM<<endl;
	PiProd.CalcHKinemaLab(q, omega, tpqCM, W, Mt, M1, M2);
	//PiProd.CalcHKinema(W, q, omega, tpqCM, Mt, M1, M2);
      }
      if (!PrintLevel)
	PiProd.PrintHKinema(Mode);
  } // end if (te)
      
  if (vp) {
    // calculate virtual photon variables
    kinema.calcVirtualPhoton(Ee, Eep, q, Q2, te, W, Mt, M1, M2); 
    if  (!PrintLevel) {
      kinema.printVPh();
    } else if (PrintLevel == 2) {
      kinema.printVPhFlux(te);
    }
  } // if(vp)

  cout <<"tp0 : "<<(tq-tpq)*rad2deg<<" tp180: "<<(tq+tpq)*rad2deg <<endl;
  return 0;
}

//
// Function name : PlainPrint(double*11)
//
// Description   : printout plain output for kinematic variables needed for xs.
// Input         : Ee, te, Eep, omega, W, tq, q, Q2, Mt, M1, M2
// Return        : 
//
int 
PlainPrint(double Ee, double te, double Eep, double omega, double W, double tq,
	   double  q, double Q2, double  Mt, double    M1, double M2) {

  VP kinema;

  double epsln    = kinema.calcVPhPolarization(q, Q2, te);

  /* calculate CM variables */
  double E1_cm    = HadronEnergyCM(W, M1, M2);
  double p1_cm    = Energy2Mom(E1_cm, M1);
  double kgamma_cm   = calcPhotonEquivalentEnergyCM(W, Mt);
  double pikgamma = p1_cm/kgamma_cm;

  double q_cm, omega_cm;
  double beta     = calcBeta(q, omega, Mt);
  VPLab2CM(beta, q, omega, q_cm, omega_cm);
  //printVPLab2CM( q_cm, omega_cm);
  double Qomega_cm  = sqrt(Q2)/omega_cm;

cout<<" epsilon: "<<epsln<<endl;
/*  printf("  %4.1f  %2.2f  %4.1f  %4.0f  %2.3f  %4.1f  %2.3f  %2.2f  %2.3f\
  %2.3f  %2.3f\n", 
	 Ee, te*rad2deg, Eep, W, Q2*1e-6, omega, q/hc, tq*rad2deg, 
	 epsln, pikgamma, Qomega_cm); 
hamza commented out 
*/

  return 0;

} // end of PlainPrint






