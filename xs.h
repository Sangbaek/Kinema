#ifndef XS_H
#define XS_H

class VP
{
 private:
  double me; 
  double k,kp;
  double kgamma, kgamma_cm, Ppikgamma; 
  double Qomega_cm,  Qomega_cm2, qCMomegaCM, QqCM; 
  double epslnL,  rhoTL,   rhoTT,   rhoTLp,   rhoTTp;

  double Gamma, Jacobian;
  double epsln, epslnTL, epslnTT, epslnTLp, epslnTTp;
 
 public:
  int calcVirtualPhoton(double Ee, double Eep, double q, double Q2, 
                        double te, double W, double Mt, double M1, double M2); 

  double calcPpikgamma(double W, double Mt, double M1, double M2);
  double calcQomegacm(double Q2, double  q, double  W, double Mt);
  double calcqcmomegacm(double Q2, double  q, double  W, double Mt);
  double calcQqcm(double Q2, double  q, double  W, double Mt);

  double calcVPhPolarization(double q, double Q2, double te);
  double calcVPhPolarization(double, double, double, double, double );

  double calcVPhKinVarTL(double  epsln);
  double calcVPhKinVarTLp(double epsln);
  double calcVPhKinVarTTp(double epsln);

  double calcPhotonEquivalentEnergy(double W, double Mt);	
  double calcPhotonEquivalentEnergy(double, double, double);	

  double calcVPhFlux(double Ee, double Eep, double kgamma, double Q2, 
		                                       double epsln);

  double calcJacobian(double Mt, double Ee, double Eep, double W);

  double getGamma();

  void   printVPh();
  void   printVPhFlux(double);
		  
  void   getVP(double &Gamma_, double &epsln_, double &epslnTT_, 
	       double &epslnTL_, double &epslnTLp_, double &Jacobian_);

};
 
#endif /* XS_H */



















