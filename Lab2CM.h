#ifndef LAB2CM_H
#define LAB2CM_H

int VPLab2CM(double beta, double q, double omega, 
		             double &q_cm, double &omega_cm) ;

double HLab2CM(double W, double M1, double M2);

void printVPLab2CM(double q_cm, double omega_cm);

double HadronEnergyCM(double W, double M1, double M2);

double calcPhotonEquivalentEnergyCM(double W, double Mt);

double calcbetaCM(double q, double omega, double M);

double calcbetaCM(double q, double W);

double calcgammaCM(double omega, double W, double M);

double calcJacobianLab(double gammaCM, double pLab, double pCM, 
		       double betaCM,  double ELab, double tpqLab);

double calcJacobianCM(double gammaCM, double pLab, double pCM, 
		      double betaCM,  double ECM,  double tpqCM);


double calcJacobianLab2CM(double tCM, double betaCM, double beta1CM);

#endif /* LAB2CM */










