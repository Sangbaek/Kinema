#ifndef RELATIVISTIC_H
#define RELATIVISTIC_H

double Energy2Mom(double E, double M);
double Mom2Energy(double p, double M);

double calcBeta(double p, double E);
double calcBeta(double p, double T, double M);

double calcgamma(double beta);

double Mandelstam_t(double Ma, double M1, double Ea, double E1, double theta);
double Mandelstam_s(double Ma, double Mb, double Ea, double Eb, double theta);


#endif /* RELATIVISTIC_H */

