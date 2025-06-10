//
// Fermi.cc
//
#include <iostream.h>
#include <fstream.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>

int
Usage(char *argv[]){

  cout << "\n Usage:" << argv[0] << "[-h] [-x]" << endl;
  cout << "\n Description: " << endl;
  cout << "\t calculate error propagation   " << endl;
  cout << "\n Options:" << endl;
  cout << "\t -Z <Z>  nucleus Z [def]=82 " << endl;
  cout << "\t -o \t output r,rho_r,integral every step" << endl;
  cout << "\t -H \t Harmonic Oscillator instead of Fermi model" << endl;
  cout << "\t -h \t show this help    " << endl;
  cout << "\t -x \t show example    " << endl;
  cout << endl;
  exit(0);

}

int 
Example(char *argv[]){

  cout << "\n Exapmle: " << endl;
  cout << "\t" << argv[0] << " -x" << endl << endl;;
  cout << "\t" << argv[0] << " -Z 6" << endl;

  cout << endl;
  exit(0);

}


double 
Fermi(double c, double z, double w, double r){

  return (1+w*r*r/c/c)/(1+exp((r-c)/z)) ;

}

double 
HarmonicOscillator(double a, double alpha_ho, double r){

  double a2       = a*a;
  double rho_r = (1 + alpha_ho*(r*r/a2))*exp(-r*r/a2);

  return rho_r ;

}



void 
target(int Z, double &c, double &z, double &w, double &a, double &alpha_ho){

  switch (Z) {
  case 82 : // 208Pb
    c     = 6.624;
    z     = 0.549;
    w     = 0;
    break;
  case 50: // 120Sn
    c     = 5.315;
    z     = 0.576;
    w     = 0;
    break;
  case 12: // Mg
    c     = 3.108;
    z     = 0.607;
    w     = -0.163;
    /* 2par
    c     = 2.942;
    z     = 0.538;
    */
    break;
  case 20: // Ca
    c     = 3.725;
    z     = 0.591;
    w     = 0-0.169;
    break;
  case 6 :  // 12C
    c     = 2.355;
    z     = 0.5224;
    w     = -0.149;
    a     = 1.687;   // For Harmonic Oscillator
    alpha_ho = 1.067;
  break;
  }

  if (c*z==0) {
    cerr << "target: Error, Z=" << Z << " is not in the database" << endl;
    exit(-1);
  }

  return;
}


int
main(int argc, char *argv[]) {

  int TextOut = 0;
  int calcHO  = 0;

  // Default 208Pb target
  int Z=82; 
  double c, z, w, a, alpha_ho;
  c=z=w=a=alpha_ho=0;

  // integral step //
  double dr = 0.01;

  int opt;
   while (EOF != (opt = getopt(argc, argv, "Hhxo?Z:"))) {
    switch (opt) {
    case 'Z':
      Z = atoi(optarg);
      break;
    case 'H':
      calcHO = 1;
      break;
    case 'x':
      Example(argv);
      break;
    case 'o':
      TextOut = 1;
      break;
    case 'h':
    case '?':
    case '*':
      Usage(argv);
    }

   }

   // get Fermi Model parameters for target Z
   target(Z, c, z, w, a, alpha_ho);

   double r   = 0;
   double rho_r = 0;
   double integral = 0;
   for (int i=0; i<3000; i++){
     rho_r = calcHO ? HarmonicOscillator(a, alpha_ho, r) : Fermi( c, z, w, r);
     //if (rho_r<0) break;

     integral += rho_r * 4*M_PI*r*r * dr;
     if (TextOut) cout << r << " " << rho_r << " " << integral << endl;

     r += dr;

   }

   if (!TextOut) {
     printf("r = %6.2f ",   r);
     printf("rho(r) = %8.3e ", rho_r);
     printf("integral = %6.3f ",integral);
     printf("\n");
     printf("rho_0=%6.3e\n",1/integral);
   }

   return 0;

}

