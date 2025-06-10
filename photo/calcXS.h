#ifndef calcXS_H
#define calcXS_H

/* Constants from Cornell Paper PRL33, 1400 (1974) */
const double CORNELL_CONST_bn  = 1.48;     // constants for nuclear coherent
const double CORNELL_CONST_bc  = 7.92e-6;  // constants for Coulomb pi0->2gamma decay width [MeV]
const double CORNELL_CONST_phi = 1.;       // nuclear coherent and Coulomb phase shift [rad]
const double CORNELL_CONST_bb  = 6.3;      // constants for nuclear incoherent

double NuclearCoherentXS(double A, double t1, double E0, double Fnc);
double PrimakoffXS(double p1, double m1, int Z, double E0, double q, double Fc, double t1);
double InterferenceXS(double Xc, double Xnc);
double NuclearIncoherentXS(double Finc);
double TotalXS(double Xnc, double Xc, double Xint, double Xinc);


#endif /* calcXS_H */




