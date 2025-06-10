//
//  smite.cc
//

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include <tgmath.h>
#include "kinema.h"
#include "const.h"
#include "calc.h"
#include "xs.h"
#include "relativistic.h"
#include "Lab2CM.h"


using namespace std;

double LinearIntegration(double, double, double, double);

// Initialize Print Level 
//PrintLeve=0 : Prints Output of Each Functions
//
int PrintLevel = 0;

//int InteractiveInput(); 

int 
main(int argc, char *argv[]) {

  int opt;
  int vp=0;

  double q2, Rec, ERec, PRec;  // returns from calceKinema(elastic)
  double Mt, M1, M2, Ee, te, tpq, Q2, W;
  double te1,te2;
  double omega, Eep, q, tq;   // returns from calcKinema(inelastic)
  double p1, E1, p2, E2;      // returns from calcKinema(exclusive)

  Rec = ERec = PRec =0;
  Mt = M1 = M2 = Ee = te = tpq = Q2 = W = 0;
  omega = Eep = q = q2 = tq = 0;
  p1 = E1 = p2 = E2 = 0;

  // Default Particle 
  Mt = M1 = MASS_PROTON ;
  M2 = MASS_PI0    ;

  // parse command line options
  while (EOF != (opt = getopt(argc, argv, "ih?ve:f:t:L:T:"))) {
    switch (opt) {
    case 'e':
      Ee = atof(optarg);
      break;
    case 'f':
      Eep= atof(optarg);
      break;
    case 'L':
      PrintLevel = atoi(optarg);
      break;
    case 't':
      te1 = atof(optarg) * deg2rad;
      break;
    case 'T':
      te2 = atof(optarg) * deg2rad;
      break;
    case 'v':
      // calculate vertual photon variables ;
      vp = 1;
      break;
    case 'h':
    case '?':
    case '*':
      cout << "\n Usage: smite [options]" << endl;
      cout << " " << endl;
      cout << "\t -h \t Show this [h]elp    "<< endl;
      cout << "\t -e \t Incident electron [e]nergy in [MeV/c]" << endl;
      cout << "\t -f \t [f]inal electron energy    in [MeV/c]" << endl;
      /*      cout << "\t -L \t Print [L]evel=0 : Normal Output(def)" << endl;
      cout << "\t    \t       [L]evel=1 : data only         " << endl;
      cout << "\t    \t       [L]evel=2 : Virtual Photon Gamma data" << endl;*/
      cout << "\t -t \t Lab [t]heta_e Angle         in [deg.]" << endl;
      //      cout << "\t -v \t calculate [v]irtual photon polarization" << endl;
      cout << " " << endl;
      cout << " Examples:" << endl;
      cout << " \t (in)elastic:\t smite -e 950 -t 28.74 -f 600" << endl;
      cout << " " << endl;
      exit(0);
    }
  }

  VP kinema;

  int    ndiv = 1000;
  double x1, x2, Gamma1, Gamma2, Flux;
  double dte  = (te2-te1)/ndiv;
  Flux=0;

  te=te1;
  int i;
  for (i=0; i<=ndiv; ++i){

    CalcWQKinema(Mt, Ee, te, Eep, omega, W, q, Q2, tq);
    kinema.calcVirtualPhoton(Ee, Eep, q, Q2, te, W, Mt, M1, M2); 
    Gamma1=kinema.getGamma();
    x1=te;

    te=te+dte;
    CalcWQKinema(Mt, Ee, te, Eep, omega, W, q, Q2, tq);
    kinema.calcVirtualPhoton(Ee, Eep, q, Q2, te, W, Mt, M1, M2); 
    Gamma2=kinema.getGamma();
    x2=te;

    Flux += LinearIntegration(x1, Gamma1, x2, Gamma2);

  }

  cout << "Total Flux = " << Flux << endl;

  return 0;
}



//
// Class name  : 
// Method name : LinearIntegration(float dx, float y1, float y2)
//
// Description : Linearly Integrete the Area trapezoid surrounded by
//             : (x1,y1), (x1,0), (x2,0), and (x2,y2)
// Input       : float x1, float y1, float x2, float y2
// Return      : float (y2 - y1)*(x2 - x1)/2
//

double 
LinearIntegration(double x1, double y1, double x2, double y2){

  double Trapezoid = (y2 + y1)*(x2 - x1)/2;
    
  return abs(Trapezoid);

}//end of LinearIntegration()




