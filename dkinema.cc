// 
// dkinema.cc
//

#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include "calc.h"
#include "const.h"
#include "dkinema.h"

using namespace std;

// Initialize Print Level 
//PrintLeve=0 : Prints Output of Each Functions
//
int PrintLevel = 0;


void
Usage(){
      
      cout << "\n Usage: dkinema [options]" << endl;
      cout << " " << endl;
      cout << "\t -e \t Incident Energy   in MeV/c"       << endl;
      cout << "\t -t \t Lab theta_e Angle in deg "        << endl;
      cout << "\t -Q \t |Momentum Transfer| in (GeV/c)^2" << endl;
      cout << "\t -W \t CM Energy Transfer in MeV/c"      << endl;
      cout << "\t -p \t Momentum Transfer Error [%]"      << endl;
      cout << "\t -s \t Scatterning Angle Error [rad]"   << endl;
      cout << "\t -D \t Readin errors from DKINEMA.DAT [def]: off" << endl;
      cout << "\t -L \t print [L]evel=1: plain data "     << endl;
      cout << "\t -u \t Show units of output data (used with option L)" 
	   << endl;
      cout << "\t -h \t Show this help    "               << endl;
      cout << " "                                         << endl;
      cout << " Examples:"                                << endl;
      cout << " \t dkinema -e 950 -Q 0.127 -W 1232"       << endl;
      cout << " \t dkinema -e 950 -Q 0.127 -W 1232 -L >> outfile" << endl;
      cout << " \t dkinema -e 950 -Q 0.127 -W 1232 -L -u" << endl;
      cout << " \t dkinema -e 950 -Q 0.127 -W 1232 -p 0.1 -s 3e-3" << endl;
      cout << " "                                         << endl;
      exit(0);
}



// Function name : getParameters(double &p1, &p2)
//
// Description   : get parameters for acceptance delta efficiency correction
// Input         : 
// Return        : double p1, double p2
//
void
getParameters(double &p1, double &p2){

  ifstream in_file;
  in_file.open("DKINEMA.DAT");

  if (in_file.fail()) {
    cout << " Error: Opening EffCorr.dat file " << endl;
    exit(-1) ;
      }

  in_file >> p1 >> p2 ;
  in_file.close();

  return ;
}



int 
main(int argc, char *argv[]) {

  int opt;
  int ShowUnits = 0;
  int fromFile = 0;
  double Mt, Ee, te, Q2, W;
  double omega, Eep, q, tq;   // returns from calcKinema(inelastic)
  double dE, dang, dEep, dte;
  Mt = Ee = te = Q2 = W = 0;
  omega = Eep = q = tq = 0;
  dE = dang = dEep = dte = 0;

  // Default Target
  Mt = MASS_PROTON ;

  // parse command line options
  while (EOF != (opt = getopt(argc, argv, "DLhue:t:Q:W:s:p:?"))) {
    switch (opt) {
    case 'e':
      Ee = atof(optarg);
      break;
    case 'L':
      PrintLevel = 1;
      break;
    case 't':
      te = atof(optarg) * deg2rad;
      break;
    case 'u':
      ShowUnits = 1;
      break;
    case 'Q':
      Q2 = atof(optarg) * 1e6;  // [MeV/c]^2
      break;
    case 'W':
      W = atof(optarg) ;
      break;
    case 'p':
      dE = atof(optarg) ;
      break;
    case 's':
      dte = atof(optarg) ;
      break;
    case 'D':
      fromFile = 1;
      break;
    case 'h':
    case '?':
    case '*':
      Usage();
    } // switch 
  }  // end-of-while loop

  if (W) {
    if (Q2) {
      CalcKinema(Mt, Ee, Q2, W, omega, Eep, q, te, tq); 
    } else if (te) {
      CalcteKinema(Mt, Ee, te, W, omega, Eep, q, Q2, tq);
    } else {
    cerr << "dkinema: (-e, -Q, -W) or (-e -f -W) inputs required" << endl;
    exit (-1);
    }
  } else {
    cerr << "dkinema: (-e, -Q, -W) or (-e -f -W) inputs required" << endl;
    exit (-1);
  }

  if (!PrintLevel)
    PrintKinema(Mt, Ee, Q2, W, omega, Eep, q, te, tq); 

  if (fromFile) {
  // getting errors from "DKINEMA.DAT" in the unit of [%]
    getParameters(dE, dang);
    dEep = Eep * dE / 100;
    dte = te * dang /100;
  } else {
    if ( dte*dE == 0 ) { 
      cerr << "Error: Either -p or -s or both options are missing" << endl;
      Usage();
    } else {
      dEep = Eep * dE / 100;
    }
    // Scattered Electron Momentum Resolution
    //dEep = Eep * 1e-2 ; /* BLAST */
    //dEep = Eep * 1.2e-3 ; /* OHIPS */
    //dEep = Eep * 1e-4 ; /* MAINZ */
    //Angular Resolution of Scattered Eelectrons = 5 mrad
    //dte  = 5e-3    ;  /* BLAST */
    //dte  = 3.5e-3  ;  /* OHIPS */
    //dte  = 3.0e-3  ;  /* MAINZ */
  }

  // calculate dtheta_q in quadratic sum of dtheta_e and dEe'
  double dtq = CalcDtq(Ee, Eep, te, q, tq, dte, dEep);

  // Printing Routine
  if (!PrintLevel)
    printDkinema(te, dte, Eep, dEep, tq, dtq);
  else {  // output plain data (no texts)
    if (ShowUnits){
      /*
      printf("\n");
      printf("   W     Q^2     Ee'   te     q  ");
      printf("   tq   dEe'  dte   dtq  \n");
      printf("[MeV] [GeV/c]^2 [MeV] [deg] [fm-1]");
      printf("[deg] [MeV] [mrad] [mrad]\n");
      printf("------------------------------------------------------------\n");
      */
      printf("\n");
      printf("   W     Q^2     Ee'   te     q  ");
      printf("   tq   dEe'   dte     dtq  \n");
      printf("[MeV] [GeV/c]^2 [MeV] [deg] [fm-1]");
      printf("[deg]  [%]  [mrad] [mrad] [%]\n");
      printf("----------------------------------------------------------------\n");
    }
    /*
    printf("  %4.0f  %1.3f  %4.1f  %4.1f  %2.3f  %4.1f  %2.2f  %2.2f  %2.2f\n", 
              W, Q2*1e-6, Eep, te*rad2deg, q/hc, tq*rad2deg, 
	      dEep, dte*1e3, dtq*1e3);
    */
    printf("  %4.0f  %1.3f  %4.1f  %4.1f  %2.3f  %4.1f  %2.2f  %2.2f  %2.2f %2.2f\n", 
              W, Q2*1e-6, Eep, te*rad2deg, q/hc, tq*rad2deg, 
	      dEep/Eep*100, dte*1e3, dtq*1e3, dtq/tq*100);
  }

} // end of main()


//
// Function name : CalcDtq( double*7 )
//
// Description   : calculate angular error of q-vector from dtheta_e and dEe'
//               : See BLAST logbook No.1 P.21 
// Input         : Ee, Eep, te, q, tq, dte, dEep (= dEe')
// Return        : dtq
//
double 
CalcDtq(double Ee, double  Eep, double te, double   q, double   tq,
                                           double dte, double dEep) {

  double cte = cos(te) ;
  double ste = sin(te) ;
  double ctq = cos(tq) ;
  double q3  = q*q*q   ;

  // ceDte: coefficiency Dtheta_e 
  double ceDte   = Eep/q*cte - Ee*Eep*Eep/q3*ste*ste ;
  // ceDEep: coefficiency DEe' 
  double ceDEep  = ste/q - Eep*ste/q3*(Eep - Ee*cte) ;

  double ctqdtq2 =  ceDte*ceDte  *  dte*dte
                 + ceDEep*ceDEep * dEep*dEep;

  double ctqdtq  = sqrt(ctqdtq2) ;

 return ctqdtq/ctq ;

} // end of CalcDtq


//
// Function name : printDkinema(double dte, double dEep, double dtq) 
//
// Description   : printout errors of theta_e, Ee', theta_q
// Input         : te, dte, Eep, dEep, tq, dtq
// Return        : 
//
void 
printDkinema(double te, double dte, double Eep, double dEep, 
	     double tq, double dtq) {

  printf("\n");
  printf(" Errors:\n");
  printf("\t dEe'      : %5.2f [MeV]  (%5.2f [%])\n", dEep, dEep/Eep*100);
  printf("\t dtheta_e  : %5.2f [mrad]  %2.2f [deg]", dte*1e3, dte*rad2deg);
  printf("(%5.2f [%])\n",dte/te*100);
  printf("\t dtheta_q  : %5.2f [mrad]  %2.2f [deg]", dtq*1e3, dtq*rad2deg);
  printf("(%5.2f [%])\n",dtq/tq*100);
  printf("\n");

}

