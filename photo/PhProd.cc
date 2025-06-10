#include <iostream.h>
#include <fstream.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include "const.h"
#include "calc.h"
#include "FF.h"
#include "calcXS.h"
#include "PhProd.h"
#include "Phys.h"
#include "Math.h"

int
Usage(char *argv[]){

  cout << "\n Usage:" << argv[0] << " [-swRhXv] [-D <days>][-H <hours>][-d <angle>]" << endl;
  cout << "\n Description: calculate real photon kinematics" << endl;
  cout << "\n Options:" << endl;
  // cout << "\t -C        calculate in Center-of-Mass system [def]:Lab" << endl;
  cout << "\t -w \t     energy weighted average cross section [def]:off" << endl;
  cout << "\t -d <angle>  dtheta [deg] (def:0.005) " << endl;
  cout << "\t -D <days>   running days" << endl;
  cout << "\t -H <hours>  running hours" << endl;
  cout << "\t -i \t     Integrate up to angle theta_pi" << endl;
  cout << "\t -R \t     output rate [/s]  " << endl;
  //  cout << "\t -v \t     verbose mode      " << endl;
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


//
// Function name : GetInput
//
// Description   : Get Input in Extended Format
// Input         : 
// Return        : 
//
void  
GetInput(double &E0, double &t1, double &Q, double &omega, double &tq, double &p1, 
	 double &M1, double &A, int &Z, double &Mt, double &p1CM, double &t1CM_fwd, 
	 double &t1CM_bck, double &E0CM, double &qCM_fwd, double &qCM_bck){


  cin >> E0 >> t1 >> Q >> omega >> tq >> p1 >> M1 
      >> A >> Z >> Mt >> p1CM >> t1CM_fwd >> t1CM_bck
      >> E0CM >> qCM_fwd >> qCM_bck;

  return;

}


//
// Function name : GetInput
//
// Description   : Get Input in short Format
// Input         : 
// Return        : 
//
void  
GetInput(double &E0, double &t1, double &Q, double &omega, double &tq, double &p1, 
	 double &M1, double &A, int &Z, double &Mt, double &ELab, double &tLab, int &CMS){


  cin >> E0 >> t1 >> Q >> omega >> tq >> p1 >> M1 
      >> A >> Z >> Mt >> ELab >> tLab >> CMS;

  return;

}

 
void
PrintOut(int PrintMode, double t1, double Q, double Fnc, double Fc, 
	 double Finc, double Xnc, double Xc, double Xint, double Xinc, double Xtot){


  if (PrintMode>>3) { // show unit

    if (PrintMode&2) { // Counting Rate 

      printf("*---------------------------------------------------------------------------------------------\n");
      printf("* theta   Q     Fnc     Fc    Finc      Rnc         Rc         Rint        Rinc        Rtot\n");
      printf("* [deg] [GeV]                           \n");
      printf("*---------------------------------------------------------------------------------------------\n");

    } else {  // Cross Sections

      printf("*---------------------------------------------------------------------------------------------\n");
      printf("* theta   Q     Fnc     Fc    Finc      Xnc         Xc         Xint        Xinc        Xtot\n");
      printf("* [deg] [GeV]                         [mb/sr]     [mb/sr]     [mb/sr]     [mb/sr]     [mb/sr]\n");
      printf("*---------------------------------------------------------------------------------------------\n");
    }

  }

  double unit = PrintMode&2 ? 1 : mcb2mb ;

  printf("%6.3f", t1*r2d); // [deg] 
  printf("%7.3f", Q*M2G);  // [GeV/c]
  printf("%7.3f", Fnc);    
  printf("%7.3f", Fc);    
  printf("%7.3f", Finc);
  printf("%12.4e", Xnc*unit);
  printf("%12.4e", Xc*unit);
  printf("%12.4e", Xint*unit);
  printf("%12.4e", Xinc*unit);
  printf("%12.4e", Xtot*unit);
  printf("\n");

  return;

}


void
kinemaInCMorLab(int CMS, double &theta, double t1CM, double t1, double &p, double p1CM, double p1, 
		double &E, double E0CM, double E0, double &Q, double qCM, double qLab, double &tq,
		double tqLab){


  theta = CMS!=0 ? t1CM    : t1 ;
  p     = CMS!=0 ? p1CM    : p1 ;
  E     = CMS!=0 ? E0CM    : E0 ;
  Q     = CMS!=0 ? qCM     : qLab ; // 4-vector Q is invaliant
  //  q     = CMS!=0 ? qCM     : qLab ;
  tq    = CMS!=0 ? t1CM    : tqLab ; 

  return ;

}


//
// Function name : CalcCrossSection
//
// Description   : Calculate form factors and cross sections
// Input         : A, Z, Q, tq, theta, E, p, M1
// Return        : Fnc, Fc, Finc, Xnc, Xc, Xint, Xinc, Xtot
//
void
CalcCrossSection(double A, int Z, double Q, double tq, double theta, double E, double p, double M1,
		 double &Fnc, double &Fc, double &Finc, double &Xnc, double &Xc, double &Xint, 
		 double &Xinc, double &Xtot){

  //                                   //
  // Calculate Model Form Factors here //
  //                                   //

  // if Model is Cornell, three mom q has to be parsed, but 4vec Q is parsed for now.//
 
  // Note: Sevelal models are selected manualy here.
  // Currently Available options are
  // Nuclear Coherent   1. Cornell  and 2. Sergey 
  // Nuclear Incoherent 1. Cornell and 2. Cascade
  // Cornell parametrizations are selected by default. (Jan.29,2005 IN)


  // Nuclear Coherent Form Factor   
  Cornell Model_NC;
  //Sergey Model_NC;
  Fnc  = Model_NC.NuclearCoherentFF(A, Z, Q, tq); 

  // Coulomb Form Factor for Primakoff
  Cornell Model_C;
  Fc   = Model_C.CoulombFF(A, Z, Q, tq);

  // Nuclear Incoherent Model
  Cornell Model_NI;
  //  Cascade Model_NI;
  Finc = Model_NI.NuclearIncoherentFF(A, Z, theta);



  // Calculate Cross Sections from Form Factors
  Xnc  = NuclearCoherentXS(A, theta, E, Fnc);    
  Xc   = PrimakoffXS(p, M1, Z, E, Q, Fc, theta);
  Xint = InterferenceXS(Xc, Xnc);
  Xinc = NuclearIncoherentXS(Finc);
  Xtot = TotalXS(Xnc, Xc, Xint, Xinc);


  return;

}


//
// Function name : GetEfficiency
//
// Description   : Calculate Detection Efficiency from phase space database
// Input         : theta [rad]
// Return        : eff
//
double
GetEfficiency(double theta){

  // File Operation 
  // Currently hardcorded here //
  char * PhaseSpace_DB = "/home/itaru/kinema/db/PhaseSpace_5GeV.db";

  ifstream in_file;
  in_file.open(PhaseSpace_DB);
  if (in_file.fail())
    {
      cerr <<"GetEfficiency(): File opening Error. \n "<< PhaseSpace_DB
           <<" is not found."<<endl;
    exit(-1);
   }

  // Scan PhaseSpace database until finding interporation point of theta_pi angle
  int ch;
  int j=0;
  int OutOfRange=1;
  double x[MAXDATA],y[MAXDATA];
  for (int i=0; i<MAXDATA; i++) x[i]=y[i]=0; // Initialize Arrays

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
    //cerr << "GetEfficiency: Warning theta is out of database coverage. " << endl;
    //cerr << "               Detection Efficiency is 0 at theta_pi > " << x[j] << endl;
  }

  // Perform Linear Interporation 
  eff = LinearInterporate(x[j], y[j], x[j-1], y[j-1], theta*r2d);

  return eff;

}



//
// Function name : CalcRate
//
// Description   : Calculate Counting Rate from Cross section and Luminosity
// Input         : 
// Return        : 
//
void 
CalcRate(double Nt, double Qeq, double theta, double eff, double runtime,
	 double &Xnc, double &Xc, double &Xint, double &Xinc, double &Xtot){ 

  double Lumi   = Nt*Qeq/cm2mcb ; // [/s/mcb]
  double dOmega = sin(theta)*dtheta*dphi ; // [sr]

  // Get theta dependent detection efficiency map form the phase space database
  // Since GetEfficiency routine needs some work, this function call is disabled.
  // Jan.29, 2005 (IN)
  //eff = GetEfficiency(theta);

  Xnc  = Xnc  * Lumi * dOmega * eff * runtime ; // [mcb/sr] * [/s/mcb] * [sr] * [s]
  Xc   = Xc   * Lumi * dOmega * eff * runtime ; // [mcb/sr] * [/s/mcb] * [sr] * [s]
  Xint = Xint * Lumi * dOmega * eff * runtime ; // [mcb/sr] * [/s/mcb] * [sr] * [s]
  Xinc = Xinc * Lumi * dOmega * eff * runtime ; // [mcb/sr] * [/s/mcb] * [sr] * [s]
  Xtot = Xtot * Lumi * dOmega * eff * runtime ; // [mcb/sr] * [/s/mcb] * [sr] * [s]

  return ;

}



//
// Function name : recalcKinema
//
// Description   : reCalculate kinematics
// Input         : 
// Return        : 
//
void 
recalcKinema(double ELab, double tLab, double M1, double Mt, int CMS, 
	     double &Q, double &tq, double &theta, double &E, double &p){

  // Get target and residual nucleus mass
  double M2 = Mt ; 
  Hadron PhProd;

  double omega, p1;
  // calculate kinematics in Laboratory Frame
  int err = PhProd.CalcgKinemaLAB(ELab, M1, tLab, M2, Q, omega, tq, p1);

  // calculate kinematics in CM Frame and Get CM angle
  PhProd.CalcHKinemaCM(ELab, tLab, ELab+Mt, Mt, M1, M2);

  if (!CMS) { //Lab
    E = ELab;
    theta = tLab ;
    p = p1 ;
  }else{ // CMS
    PhProd.getCMSKinema(tq, theta, E, p);
  }

  return;
}


//
// Function name : calcEnergyWeightSum
//
// Description   : Recalculate Kinematics and Energy Weight Sum Cross Section
// Input         : ELab, tLab, M1, Mt, CMS, A, Z, E0
// Return        : &Fnc, &Fc, &Finc, &Xnc, &Xc, &Xint, &Xinc, &Xtot
//
void 
calcEnergyWeightSum(double ELab, double tLab, double M1, double Mt, int CMS,  
		    double A, int Z, double E0, 
		    double &Fnc,  double &Fc, double &Finc, 
		    double &Xnc,  double &Xc, double &Xint, 
		    double &Xinc, double &Xtot){

  RESET = 0;

  double Q, tq, theta, E, p;
  Q = tq = theta = E = p = 0;

  double F[nF], X[nX];
  for (int i=0; i<nF; i++) F[i] = 0;
  for (int i=0; i<nX; i++) X[i] = 0;

  double Egamma=0;
  for (Egamma=ELab*EMin; Egamma<=ELab*EMax; Egamma+=ELab*dE) {

    // recalculate kinematics 
    recalcKinema(Egamma, tLab, M1, Mt, CMS, Q, tq, theta, E, p);

    // calculate Cross Section
    CalcCrossSection(A, Z, Q, tq, theta, E, p, M1, F[0], F[1], F[2], X[0], X[1], X[2], X[3], X[4]);

    if (Egamma >= ELab*(EMax-dE)) RESET = 1; 
    // integrate over the energy range
    EnergyWeightSum(Egamma, E0*dE, F, X, Fnc, Fc, Finc, Xnc, Xc, Xint, Xinc, Xtot);

  }

  return ;

}


//
// Function name : Energy Weight Sum
//
// Description   : Calculate Energy Weighted Sum of Form Factors and Cross Sections
// Input         : kgamma, dk, F[], X[], &Fnc, &Fc, &Finc, &Xnc, &Xc, &Xint, &Xinc, &Xtot
// Return        : &Fnc, &Fc, &Finc, &Xnc, &Xc, &Xint, &Xinc, &Xtot
//
void 
EnergyWeightSum(double kgamma, double dk, double F[], double X[], 
		double &Fnc,  double &Fc, double &Finc, 
		double &Xnc,  double &Xc, double &Xint, 
		double &Xinc, double &Xtot){

  static double Fsum[nF]={0,0,0};
  static double Xsum[nX]={0,0,0,0,0};

  for (int i=0; i<nF; i++) Fsum[i]+=F[i]/kgamma*dk ;
  for (int i=0; i<nX; i++) Xsum[i]+=X[i]/kgamma*dk ;

  Fnc  = Fsum[0];
  Fc   = Fsum[1];
  Finc = Fsum[2];

  Xnc  = Xsum[0];
  Xc   = Xsum[1];
  Xint = Xsum[2];
  Xinc = Xsum[3];
  Xtot = Xsum[4];

  // Please let me know if there is smarter way to initialize static variables runtime.
  if (RESET) {
    for (int i=0; i<nF; i++) Fsum[i] -= Fsum[i];
    for (int i=0; i<nX; i++) Xsum[i] -= Xsum[i];
  }

  return;
}


//
// Function name : xs2rate
//
// Description   : Calculate Counting Rate from Cross Sections
// Input         : double A, int Z, double Xtarg, double Qeq, double theta, double eff, double runtime
// Return        : &Fnc, &Fc, &Finc, &Xnc, &Xc, &Xint, &Xinc, &Xtot
//
void 
xs2rate(double A, int Z, double Xtarg, double Qeq, double theta, double eff, double runtime,
	double &Fnc,  double &Fc, double &Finc, 
	double &Xnc,  double &Xc, double &Xint, 
	double &Xinc, double &Xtot){

  // Calculate Number of nucleus in the target
  double Xrad=RadiationLength(A,Z);
  double Nt=CalcNt(A, Xrad, Xtarg);

  // Calculate Counting Rate
  CalcRate(Nt, Qeq, theta, eff, runtime, Xnc, Xc, Xint, Xinc, Xtot);

  return ;

}


//
// Function name : IntegralXS
//
// Description   : Calculate cross section integral over the acceptance 
// Input         : double theta, double dtheta 
// Return        : &Fnc, &Fc, &Finc, &Xnc, &Xc, &Xint, &Xinc, &Xtot
//
void 
IntegralXS(double theta, double dtheta, double &Fnc, double &Fc, double &Finc, 
	 double &Xnc, double &Xc, double &Xint, double &Xinc, double &Xtot){

  double dOmega = sin(theta)*dtheta*dphi ; // [sr]

  Fnc  *= dOmega ;
  Fc   *= dOmega ;
  Finc *= dOmega ;

  Xnc  *= dOmega ;
  Xc   *= dOmega ;
  Xint *= dOmega ;
  Xinc *= dOmega ;
  Xtot *= dOmega ;

  return ;

}


//
// Function name : Integral
//
// Description   : Calculate Counting Rate from Cross Sections
// Input         : double A, int Z, double Xtarg, double Qeq, double theta, double eff, double runtime
// Return        : &Fnc, &Fc, &Finc, &Xnc, &Xc, &Xint, &Xinc, &Xtot
//
void 
Integral(double dtheta, double &Fnc, double &Fc, double &Finc, 
	 double &Xnc, double &Xc, double &Xint, double &Xinc, double &Xtot){

  static double Fsum[nF]={0,0,0};
  static double Xsum[nX]={0,0,0,0,0};

  // cumulative sum
  Fsum[0] += Fnc  ;
  Fsum[1] += Fc   ;
  Fsum[2] += Finc ;

  Xsum[0] += Xnc  ;
  Xsum[1] += Xc   ;
  Xsum[2] += Xint ;
  Xsum[3] += Xinc ;
  Xsum[4] += Xtot ;

  // redump sum into Fs and Xs
  Fnc  = Fsum[0];
  Fc   = Fsum[1];
  Finc = Fsum[2];

  Xnc  = Xsum[0];
  Xc   = Xsum[1];
  Xint = Xsum[2];
  Xinc = Xsum[3];
  Xtot = Xsum[4];

  return ;

}



int
main(int argc, char *argv[]) {

  int opt;

  // Defults //
  double day = 0;
  double hour = 0;
  int PrintMode = 0; // 
  int CMS = 0; // calc XS in Lab frame otherwise, CMS
  int EnergyWeight = 0;
  int RateMode = 0 ;
  int CalcIntegral = 0;

  while (EOF != (opt = getopt(argc, argv, "h?CXiswd:D:H:R"))) {
    switch (opt) {
    case 's':
      PrintMode += 8;
      break;
    case 'X':
      Example(argv);
      exit(0);
    case 'w':
      EnergyWeight = 1;
      break;
    case 'd':
      dtheta = atof(optarg)*d2r ;
      break;
    case 'D':
      day = atof(optarg) ;
      runtime *= day*24*3600 ;
      RateMode = 1;
      break;
    case 'i':
      CalcIntegral = 1;
      break;
    case 'H':
      hour = atof(optarg) ;
      runtime *= hour*3600 ;
      RateMode = 1;
      break;
    case 'R':
      RateMode = 1;
      break;
    case 'C':
      cerr << "Error: -C option is not supported any longer" << endl;
      exit(-1);
      /* yet ready
   case 'v':
      VERBOSE = 1;
      break;
      */
    case 'h':
    case '?':
    case '*':
      Usage(argv);
      exit(0);
    }
  }

  if (RateMode) PrintMode += 2;

  // Initialization of Variables
  double Fnc, Fc, Finc, Xnc, Xc, Xint, Xinc, Xtot;
  Fnc = Fc = Finc = Xnc = Xc = Xint = Xinc = Xtot = 0;

  double E0, t1, Q, omega, tq, p1, M1, A, Mt, p1CM, t1CM[2], E0CM, qCM[2];
  double p, theta, ELab, tLab, t0;
  int Z = 0 ;
  E0 = t1 = Q = omega = tq = p1 = M1 = A = Mt = p1CM = t1CM[0] = t1CM[1] = E0CM = qCM[0] = qCM[1] = 0;

  // Get Input Kinematics (theta can be CM or Lab)
  GetInput(E0, theta, Q, omega, tq, p, M1, A, Z, Mt, ELab, tLab, CMS);
  double Qeq = EquivPhoton(E0);



  // integral over theta_pi angle 
  t0 = CalcIntegral ? 0 : tLab;
  for (double t=t0; t<=tLab ; t+=dtheta) {

    // calculate cross sections
    if (EnergyWeight) {
      // this routine requires laboratory frame ELab and tLab
      calcEnergyWeightSum(ELab, t, M1, Mt, CMS, A, Z, E0, Fnc, Fc, Finc, Xnc, Xc, Xint, Xinc, Xtot);
    } else {
      if (CMS) cerr << "CMS integral is not implemented yet" << endl;
      CalcCrossSection(A, Z, Q, tq, theta, E0, p, M1, Fnc, Fc, Finc, Xnc, Xc, Xint, Xinc, Xtot);
    }

    // Calculate Counting Rate
    if (RateMode) {
      xs2rate(A, Z, Xtarg, Qeq, t, eff, runtime, Fnc, Fc, Finc, Xnc, Xc, Xint, Xinc, Xtot);
    } else if (CalcIntegral) { 
      IntegralXS(t, dtheta, Fnc, Fc, Finc, Xnc, Xc, Xint, Xinc, Xtot);
   }
    Integral(dtheta, Fnc, Fc, Finc, Xnc, Xc, Xint, Xinc, Xtot);

  }// end of t-loop


  // Print out cross sections
  PrintOut(PrintMode, theta, Q, Fnc, Fc, Finc, Xnc, Xc, Xint, Xinc, Xtot);

  return 0;

}

