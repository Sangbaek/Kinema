#ifndef PHPROD_H
#define PHPROD_H

//
// global constants
//

// Photon Energy range
const double dE = 0.01; // 0.01;
const double EMax=0.95;
const double EMin=0.85;

// Total Number of Photon per sec.
const double Ngamma=8e6;

// Equivalent Photons
double EquivPhoton(double E0){

  double Einv = 0;
  for (double E=EMin; E<=EMax; E+=dE) Einv += 1/E ;

  return Ngamma/Einv/dE ;

}

// Target Radiation Length 
const double Xtarg = 5e-2;   // 5[%]

// Solid Angle Acceptance;
double dphi = 2*M_PI ;
double dtheta = 0.05*d2r ;


// Detection Efficiency
double eff = 0.7 ;

// Running Time [s]
double runtime = 1 ;

// Number of form factors and cross sections
const int nF = 3;
const int nX = 5;

// Zero static variables
int RESET = 0;


int Example(char *argv[]);
int Usage(char *argv[]);

void GetInput(double &E0, double &t1, double &Q, double &omega, double &tq, double &p1, 
	       double &M1, double &A, int &Z, double &Mt, double &p1CM, double &t1CM_fwd, 
	       double &t1CM_bck, double &E0CM, double &qCM_fwd, double &qCM_bck);

void GetInput(double &E0, double &t1, double &Q, double &omega, double &tq, double &p1, 
	      double &M1, double &A, int &Z, double &Mt, double &ELab, double &tLab, int &CMS);

void recalcKinema(double ELab, double tLab, double M1, double Mt, int CMS, 
		  double &Q, double &tq, double &theta, double &E, double &p);

void calcEnergyWeightSum(double ELab, double tLab, double M1, double Mt, int CMS,  
		    double A, int Z, double E0, 
		    double &Fnc,  double &Fc, double &Finc, 
		    double &Xnc,  double &Xc, double &Xint, 
		    double &Xinc, double &Xtot);

void EnergyWeightSum(double kgamma, double dk, double F[], double X[], 
		     double &Fnc,  double &Fc, double &Finc, 
		     double &Xnc,  double &Xc, double &Xint, 
		     double &Xinc, double &Xtot);

 
void kinemaInCMorLab(int CMS, double &theta, double t1CM, double t1, double &p, double p1CM, double p1, 
		     double &E, double E0CM, double E0, double &q, double qCM, double qLab, double &tq,
		     double tqLab);

void CalcCrossSection(double A, int Z, double Q, double tq, double theta, double E, double p, double M1, 
		      double &Fnc, double &Fc, double &Finc, double &Xnc, double &Xc, double &Xint, 
		      double &Xinc, double &Xtot);

double GetEfficiency(double theta);

void xs2rate(double A, int Z, double Xtarg, double Qeq, double theta, double eff, double runtime,
	     double &Fnc,  double &Fc, double &Finc, 
	     double &Xnc,  double &Xc, double &Xint, 
	     double &Xinc, double &Xtot);


//
// Function name : CalcNt
//
// Description   : Calculate Number of Target Nucleus
// Input         : double A, double Xrad, double Xtarg
// Return        : double Nt
//
double CalcNt(double A, double Xrad, double Xtarg){ 

  return Xtarg*Xrad/A*AVOGADRO ;

};


void PrintOut(int PrintMode, double t1, double Q, double Fnc, double Fc, double Finc, 
	      double Xnc, double Xc, double Xint, double Xinc, double Xtot);



#endif /* PHPROD_H */


