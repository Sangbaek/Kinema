#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include "calc.h"
#include "kinema.h"
#include "const.h"
#include "target.h"

using namespace std;

int
Usage(char *argv[]){

  cout << "\n Usage:" << argv[0] << "[-shX] [-e <engy>] [-t <angle>] [-T <target>]" << endl;
  cout << "\n Description: calculate real photon kinematics" << endl;
  cout << "\n Options:" << endl;
  cout << "\t -e <engy>   incident photon energy [GeV]" << endl;
  cout << "\t -t <angle>  pion production angle [deg.]" << endl;
  cout << "\t -T <target> target (He,Be,12C,Sn,Pb) [def:Pb]" << endl; 
  cout << "\t -p \t     parse output to PiProd" << endl;
  cout << "\t -s \t     show unit         " << endl;
  cout << "\t -h \t     show this help    " << endl;
  cout << "\t -X \t     show example      " << endl;
  cout << endl;
  exit(0);

}



int
Example(char *argv[]){

  cout << "\n Exapmle: " << endl;
  cout << "\t" << argv[0] << " -e 5 -t 0.6" << endl;
  cout << endl;
  exit(0);

}



int
main(int argc, char *argv[]) {

  int opt;
  double p1, t2, p2, Q, omega, tq;
  //char *target = "Pb";
  char *target = "H";

  //------------------------ Defults ---------------------------------//
  int PrintMode = 8;   //  bin(1000)
  int CMS    = 0 ;    // output in CMS
  double E0  = 4.4*G2M ; // 5[GeV]
  double t1  = 1*r2d ; // 1[deg]
  double A   = 1.00794; // 208Pb
  int    Z   = 1;    
  //int    Z   = 82;    
  //double A   = 207.19; // 208Pb
  double M1  = MASS_PI0;
  double Mt  = MASS_PROTON;//MassFormula(A,Z);
  double M2  = M1; 

  while (EOF != (opt = getopt(argc, argv, "h?XCspt:e:T:"))) {
    switch (opt) {
    case 'e':
      E0 = atof(optarg)*G2M;
      break;
    case 't':
      t1 = atof(optarg)*d2r;
      break;
    case 'T':
      target = optarg;
      GetMaterialAZ(target, A, Z);
      break;
    case 's':
      PrintMode+=1;
      break;
    case 'p':
      PrintMode=0;
      break;
    case 'C':
      CMS = 1;
      break;
    case 'X':
      Example(argv);
      exit(0);
    case 'h':
    case '?':
    case '*':
      Usage(argv);
      exit(0);
    }
  }

  // Get target and residual nucleus mass
  Mt = M2 = MASS_PROTON; //MassFormula(A,Z);
  Hadron PhProd;

  // calculate kinematics in Laboratory Frame
  int err = PhProd.CalcgKinemaLAB(E0, M1, t1, M2, Q, omega, tq, p1);

  // calculate kinematics in CM Frame and Get CM angle
  //PhProd.CalcHKinemaCM(E0, t1, E0+Mt, Mt, M1, M2);

  // print out kinematics
  //PhProd.PrintgKinema(PrintMode, CMS, target, E0, A, Z, Mt, M1, t1, M2, Q, omega, tq, p1);

  // calculate kinematics in CM Frame
  PhProd.CalcHKinemaCM(E0, t1, E0+Mt, Mt, M1, M2);
  double tCM = PhProd.getHKinema(-1);

  return 0;

}

