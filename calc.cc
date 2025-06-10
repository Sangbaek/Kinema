#include <math.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "const.h"
#include "calc.h"
#include "Math.h"
#include "Lab2CM.h"
#include "relativistic.h"
#include "Phys.h"

using namespace std;

//
// Function name : Calc'e'Kinema(double, double, double,  )
//
// Description   : calculate Elastic Electron Scattering kinematics
// Input         : Mt, Ee, te
// Return        : Rec, Eep, q2, Q2, t2, ERec, Prec
//
int
CalceKinema(double Mt, double Ee, double te, double &Rec, double &Eep, 
            double &q2, double &Q2, double &tq, double &ERec, double &PRec) {

  double cte = cos(te);
  Rec = Mt / (Mt + Ee*(1 - cte));
  Eep = Ee*Rec;
  q2  = Ee*Ee + Eep*Eep - 2*Ee*Eep*cte;
  tq  = acos( (Ee*Ee + q2 - Eep*Eep) / (2*Ee*sqrt(q2)) );
  ERec = Ee - Eep;
  PRec = sqrt(ERec*ERec + 2*ERec*Mt) ;
  Q2  = 4 * Ee * Eep * sin(te/2) * sin(te/2);

  return 0;
}// end of CalceKinema


//
// Function name : Print'e'Kinema(double, double, double,  )
//
// Description   : print Elastic Electron Scattering kinematics
// Input         : Mt, Ee, te, Rec, Eep, q2, t2, ERec, Prec
// Return        : 
//
int
PrinteKinema(double Mt, double Ee, double te, double Rec, double Eep, 
        	   double q2, double Q2, double tq, double ERec, double PRec) {

  cout << "  " << endl;
  cout << "----------------------------------------------" << endl;
  cout << " Incident Energy       = " << Ee << " [MeV]" <<endl;
  cout << " Scattering Angle      = " << te*rad2deg << " [deg.]" <<endl;
  cout << " Target Mass           = " << Mt << " [MeV]" <<endl;
  cout << " Final Eenergy         = " << Eep << " [MeV]" <<endl;
  cout << " Momentum Transfer q   = " << sqrt(q2/1e6) << " [GeV]\n" 
       << "                       = " << sqrt(q2/hc/hc) << " [fm^-1]" << endl;
  cout << " Momentum Transfer q^2 = " << q2/1e6 << " [GeV^2]\n" 
       << "                       = " << q2/hc/hc << " [fm^-2]" << endl;
  cout << " 4-Mom. Transfer Q^2   = " << Q2/1e6 << " [GeV^2]\n" 
       << "                       = " << Q2/hc/hc << " [fm^-2]" << endl;
  cout << " Theta_q               = " << tq*rad2deg   << " [deg.]" << endl;
  cout << " Recoil Factor         = " << Rec  << endl;
  cout << " Recoil Energy         = " << ERec << " [MeV]" << endl;
  cout << " Recoil Momentum       = " << PRec << "[MeV/c]" << endl;
  cout << "----------------------------------------------" << endl;
  cout << "  " << endl;

  return 0;
}// end of PrinteKinema

//
// Function name : Calc'p'Kinema(double, double, double,  )
//
// Description   : calculate Elastic Electron Scattering kinematics
//               : from thetaq to thetae
// Input         : Mt, Ee, tq
// Return        : Rec, Eep, q2, Q2, te, ERec, Prec
//
int
CalcpKinema(double Mt, double Ee, double tq, double &Rec, double &Eep, 
            double &q2, double &Q2, double &te, double &ERec, double &PRec) {

  double ctq, cte;
  ctq  = cos(tq);
  ERec = 2*Mt*ctq*ctq/((1+Mt/Ee)*(1+Mt/Ee)-ctq*ctq);
  Eep  = Ee - ERec;
  q2   = ERec*ERec + 2*ERec*Mt;
  cte  = (Ee*Ee + Eep*Eep - q2)/(2*Ee*Eep);
  te   = acos(cte);
  PRec = sqrt(ERec*ERec + 2*ERec*Mt) ;
  Rec = Mt / (Mt + Ee*(1 - cte));
  Q2  = 4 * Ee * Eep * sin(te/2) * sin(te/2);

  return 0;
}// end of CalcpKinema


//
// Function name : Calc'Pp'Kinema(double, double, double)
//
// Description   : calculate Elastic Electron Scattering kinematics
//               : from PRec
// Input         : Mt, Ee, PRec
// Return        : Rec, Eep, q2, Q2, te, ERec, tq 
//
int
CalcPpKinema(double Mt, double Ee, double PRec, double &Rec, double &Eep, 
            double &q2, double &Q2, double &te, double &ERec, double &tq) {

  double ctq, cte;
  ERec = sqrt(PRec*PRec + Mt*Mt) - Mt;
  ctq  = sqrt(ERec)*(1+Mt/Ee)/sqrt(ERec+2*Mt);
  tq   = acos(ctq);
  Eep  = Ee - ERec;
  q2   = ERec*ERec + 2*ERec*Mt;
  cte  = (Ee*Ee + Eep*Eep - q2)/(2*Ee*Eep);
  te   = acos(cte);
  Rec = Mt / (Mt + Ee*(1 - cte));
  Q2  = 4 * Ee * Eep * sin(te/2) * sin(te/2);

  return 0;
}// end of CalcPpKinema


//
// Function name : CalcExKinema(double, double, double,  )
//
// Description   : calculate Electron Scattering kinematics
//               : from given final energy
// Input         : Mt, Ee, tq
// Return        : Rec, Eep, q2, te, ERec, Prec
//
int
CalcExKinema(double Mt, double Ee, double te, double Eep, double &Rec, 
	     double &Ex, double &q2, double &tq, double &ERec, double &PRec) {

  double cte = cos(te);
  Rec = Mt / (Mt + Ee*(1 - cte));
  q2  = Ee*Ee + Eep*Eep - 2*Ee*Eep*cte;
  tq  = acos( (Ee*Ee + q2 - Eep*Eep) / (2*Ee*sqrt(q2)) );
  Ex  = Ee - Eep*Rec;

  ERec = Ee - Eep;
  PRec = sqrt(ERec*ERec + 2*ERec*Mt) ;

  return 0;
}// end of CalcpKinema



//
// Function name : CalcWQKinema(double *4, double& *5)
//
// Description   : calculate Electron Scattering kinematics
//               : without using electron mass; me->0 approximation
//               : for forward angle scattering
// Input         : Mt, Ee, 'Eep', te
// Return        : omega, W, Q2, q(3-mom.transfer), tq
//
int 
CalcWQKinema(double     Mt, double   Ee, double te, double   Eep,
	   double &omega, double &W, double &q, double &Q2, double &tq) {
  // calculate inelastic kinematics 
  double me = MASS_ELECTRON;
  double k  = Energy2Mom(Ee,  me); // incidenet electron momentum
  double kp = Energy2Mom(Eep, me); // scattered electron momentum

  // This program assumes proton as a target
  Q2    = -2*me*me + 2*Ee*Eep - 2*k*kp*cos(te);

  omega = Ee - Eep;
  q     = sqrt(omega*omega + Q2);
  W     = sqrt(2*omega*Mt - Q2 + Mt*Mt);  // this needs check
 
  double ste_2 =  Q2/4/k/kp;
  te    =   asin(sqrt(ste_2))*2;

  double   ctq = (k*k + q*q - kp*kp)/2/q/k;
  tq    = acos(ctq);

  return 0;
} // end of CalcWQKinema(inelastic, elastic)



//
// Function name : Calc'te'Kinema(double *4, double& *5)
//
// Description   : calculate Inelastic Electron Scattering kinematics
// Input         : Mt, Ee, 'te', W
// Return        : omega, Eep, q(3-mom.tranfer), 'Q2', tq
//
int 
CalcteKinema(double     Mt, double   Ee, double te, double   W,
	     double &omega, double &Eep, double &q, double &Q2, double &tq) {

  double ste_2 = sin(te/2);
  Eep   = (Mt*Mt + 2*Mt*Ee - W*W)/2/(Mt + 2*Ee*ste_2*ste_2);
  omega = Ee - Eep;
  
  double q4vec2 = -4*Ee*Eep*ste_2*ste_2;
  q     = sqrt(omega*omega - q4vec2) ;
  Q2    = -q4vec2;

  double   ctq = (Ee*Ee + q*q - Eep*Eep)/2/q/Ee;
  tq    = acos(ctq);

  return 0;

}//end of Calc'te'Kinema(inelastic)


//
// Function name : CalcKinema(double *4, double& *5)
//
// Description   : calculate Inelastic Electron Scattering kinematics
// Input         : Mt, Ee, 'Q2', W
// Return        : omega, Eep, q(3-mom.transfer), 'te', tq
//
int 
CalcKinema(double     Mt, double   Ee, double Q2, double   W,
	   double &omega, double &Eep, double &q, double &te, double &tq) {
  // calculate inelastic kinematics 

  // This program assumes proton as a target
  omega = (W*W - Mt*Mt + Q2) / 2 / Mt;
  Eep   =   Ee - omega;
  q     =   sqrt( omega*omega + Q2 );

  double ste_2 =  Q2/4/Ee/Eep;
  te    =   asin(sqrt(ste_2))*2;

  double   ctq = (Ee*Ee + q*q - Eep*Eep)/2/q/Ee;
  tq    = acos(ctq);

  return 0;
} // end of CalcKinema(inelastic)


//
// Function name : PrintKinema(double *9)
//
// Description   : print inElastic Electron Scattering kinematics
// Input         : Mt, Ee, Q2, W, omega, Eep, q, te, tq
// Return        : 
//
int
PrintKinema(double     Mt, double   Ee, double Q2, double   W,
	    double omega, double Eep, double q, double te, double tq) {

  cout << " " << endl;
  cout << " ------------------------------------------------ " << endl;
  cout << " " << endl;
  cout << "\t Incident Energy  = " << Ee << " [MeV]" << endl;
  cout << "\t Scattering Angle = " << rad2deg*te << " [deg.]" << endl;
  cout << "\t Final Energy     = " << Eep << " [MeV]" << endl;
  cout << " " << endl;
  cout << " w        = " << omega << " [MeV]" << endl;
  cout << " W        = " << W << " [MeV]" << endl;
  cout << " q        = " << q << " [MeV/c]     " << q/hc << "[fm^-1]" << endl;
  cout << " theta_q  = " << rad2deg*tq << " [deg.]" << endl;
  cout << " q^2(4vec)= " << -Q2*1e-6  << " [GeV/c]^2     " 
                   << -Q2/hc/hc << " [fm^-2]" << endl;
  cout << " Q^2      = " << -Q2*1e-6  << " [GeV/c]^2     " 
                   << Q2/hc/hc << " [fm^-2]" << endl;

  cout << " " << endl;

  return 0;
} // end of PrintKinema


//
// Function name : Calc'H'Kinema(double *5, double& *4 )
//
// Description   : calculate exclusive ([H]adron) kinematics
// Input         : q, omega, tpq(Lab), Mt, M1, M2
// Return        : 
//
int 
Hadron::CalcHKinema(double q, double omega, double tpq, double Mt,
		    double Mass1, double Mass2){

  M1 = Mass1 ;
  M2 = Mass2 ;
  double ctpq = cos(tpq);
  double stpq = sin(tpq);
  double AA  = omega + Mt;
  double BB  = M1*M1 - AA*AA - q*q - M2*M2;
  double A = 4*( q*q*ctpq*ctpq - AA*AA );
  double B = 4*q*ctpq*( BB + 2*AA*AA );
  double C = BB*BB - 4*AA*AA*( M2*M2 + q*q );

  //  p[0] = (-B - sqrt(B*B - 4*A*C))/4/C (High Momentum Branch) //
  p1[0] = QuadraticSolution(A, B, C, -1);
  //  p[1] = (-B + sqrt(B*B - 4*A*C))/4/C (Low Momentum Branch)  //
  p1[1] = QuadraticSolution(A, B, C,  1);

  for (int i=0; i<=1; i++) {

    E1[i]      = sqrt(p1[i]*p1[i] + M1*M1);

    p2[i]      = sqrt(p1[i]*p1[i] + q*q - 2*p1[i]*q*ctpq);
    E2[i]      = sqrt(p2[i]*p2[i] + M2*M2);

    tpq2Lab[i] = - asin( p1[i]/p2[i] * stpq ) ;
    // Additional 0.001 is necessary to test the sign. UGLY!!! //
    if ( p1[i]*ctpq + p2[i]*cos(tpq2Lab[i]) > q + 0.001 ) {
      tpq2Lab[i] = - M_PI - tpq2Lab[i] ;
    }

  }

  return 0;

} // end of calc'H'Kinema(exclusive)

//
// Function name : Calc'H'Kinema(double *5, double& *4 )
//
// Description   : calculate exclusive ([H]adron) kinematics
// Input         : W, q, omega, tpq(CM), Mt, M1, M2
// Return        : 
//
int 
Hadron::CalcHKinema(double W,  double q, double omega, double tpqCM, 
		    double Mt, double Mass1, double Mass2){

  M1 = Mass1 ;
  M2 = Mass2 ;
  double betaCM  = calcBeta(q, omega + Mt);
  double gammaCM = calcgamma(betaCM);

  for (int i=0; i<=1; i++) {

    tpq1CM[i] = tpqCM;
    tpq2CM[i] = M_PI - tpq1CM[i];

    E1CM[i] = (W*W + M1*M1 - M2*M2)/2/W ;
    E2CM[i] = (W*W + M2*M2 - M1*M1)/2/W ;

    p1CM[i] = Energy2Mom(E1CM[i], M1);
    p2CM[i] = Energy2Mom(E2CM[i], M2);

    //            Laboratory kinematics                   //
    E1[i] = gammaCM*( E1CM[i] + betaCM*p1CM[i]*cos(tpq1CM[i]) );
    E2[i] = gammaCM*( E2CM[i] + betaCM*p2CM[i]*cos(tpq2CM[i]) );

    p1[i] = Energy2Mom(E1[i], M1);
    p2[i] = Energy2Mom(E2[i], M2);

    tpq1Lab[i] =   p1CM[i]/p1[i] * sin(tpq1CM[i]);
    tpq2Lab[i] = - p2CM[i]/p2[i] * sin(tpq2CM[i]);

    Jacobian1[i] 
      = calcJacobianCM(gammaCM, p1[i], p1CM[i], betaCM, E1CM[i], tpq1CM[i]);
    Jacobian2[i] 
      = calcJacobianCM(gammaCM, p2[i], p2CM[i], betaCM, E2CM[i], tpq2CM[i]);

  }


  cerr << " This routine Hadron::CalcHKinema(double W, double q, double omega, double tpqCM, \n" 
       << "double Mt, double Mass1, double Mass2) requires check before being used" << endl;
  cerr << "Exit forced (I. Nakagawa, March 4, 2003)" << endl;
  exit(-1);

  return 0;

} // end of calc'H'Kinema(exclusive)


//
// Function name : Print'H'Kinema(double *5, double& *4 )
//
// Description   : Print exclusive kinematics
// Input         : Mode
// Return        :
//
int 
Hadron::PrintHKinema(int Mode){

  if (!Mode) {
   printf(" ---------------------------------------------------------------------------\n");
   printf("\t\t LAB \t\t\t |\t     CM \t     CM/Lab\n"                                    );  
   printf(" Branch  Mass\t   T  \t   p    theta_pq |    T       p   theta_pq  Jacobian \n");
   printf("  H/L    [MeV]   [MeV]  [MeV/c]  [deg.]  |  [MeV]  [MeV/c]  [deg.]          \n");
   printf(" ---------------------------------------------------------------------------\n");
  }

  char Branch='H' ;
  for (int i=0; i<=1; i++){

    if (i) Branch='L';
    printf("   %1c   %7.2f %7.2f %7.2f %10.5f %7.2f %7.2f %8.3f  %7.2f\n",
           Branch, M1, E1[i]-M1, p1[i], tpq1Lab[i]*rad2deg, 
           E1CM[i]-M1, p1CM[i], tpq1CM[i]*rad2deg, Jacobian1[i]);
    printf("   %1c   %7.2f %7.2f %7.2f %10.5f %7.2f %7.2f %8.3f  %7.2f\n",
           Branch, M2, E2[i]-M2, p2[i], tpq2Lab[i]*rad2deg, 
           E2CM[i]-M2, p2CM[i], tpq2CM[i]*rad2deg, Jacobian2[i]);

  }

  if (!Mode) printf("\n");

  return 0;
} // end of PrintKinema(exclusive)

//
// Function name : Calc'H'KinemaCM(double *5, double& *4)
//
// Description   : calculate exclusive kinematics in CM frame
// Input         : q, omega, tpq, W, M1, M2
// Return        : 
//
int 
Hadron::CalcHKinemaCM(double q, double omega, double tpq, double W, 
                      double Mt, double Mass1, double Mass2){ 

  M1 = Mass1;
  M2 = Mass2;
  tpq1Lab[0] = tpq1Lab[1] = tpq;

  double tantpq  = tan(tpq);
  double tantpq2 = tantpq * tantpq;

  //  double *ctpqCM = new double(2);
  double ctpqCM[2];
  betaCM  = calcbetaCM(q, omega, Mt);
  gammaCM = calcgammaCM(omega, W, Mt);

  double gammaCM2 = gammaCM * gammaCM;

  int i;
  double A, B, C;
  double betaFrac;
  for (i=0; i<=1; i++){ // loop for high- and low- momentum

    E1CM[i]    = HadronEnergyCM(W, M1, M2);
    p1CM[i]    = Energy2Mom(E1CM[i], M1);

    E2CM[i]    = HadronEnergyCM(W, M2, M1);
    p2CM[i]    = Energy2Mom(E2CM[i], M2);

    beta1CM[i] = calcBeta(p1CM[i], E1CM[i]);
    betaFrac   = betaCM/beta1CM[i];

    A = gammaCM2 * tantpq2 + 1;
    B = 2 * gammaCM2 * betaFrac * tantpq2 ;
    C = gammaCM2 * tantpq2 * betaFrac * betaFrac - 1;

    ctpqCM[i] = QuadraticSolution(A, B, C, -i);
    tpq1CM[i] = acos(ctpqCM[i]) ;
    tpq2CM[i] = M_PI - tpq1CM[i] ;

    Jacobian1[i] 
      = calcJacobianCM(gammaCM, p1[i], p1CM[i], betaCM, E1CM[i], tpq1CM[i]);
    Jacobian2[i] 
      = calcJacobianCM(gammaCM, p2[i], p2CM[i], betaCM, E2CM[i], tpq2CM[i]);

  }

  return 0;

}



//
// Function name : Calc'H'KinemaCM(double *5, double& *4)
//
// Description   : calculate exclusive kinematics in CM frame
// Input         : q, tpq, W, Mt, M1, M2
// Return        : 
//
int 
Hadron::CalcHKinemaCM(double q, double tpq, double W, 
                      double Mt, double Mass1, double Mass2){ 

  M1 = Mass1;
  M2 = Mass2;
  tpq1Lab[0] = tpq1Lab[1] = tpq;

  double tantpq  = tan(tpq);
  double tantpq2 = tantpq * tantpq;

  double ctpqCM[2];
  betaCM  = calcBeta(q, W);
  gammaCM = calcgamma(betaCM);
  double gammaCM2 = gammaCM * gammaCM;

  // Photon Energy in CMS
  E0CM = HadronEnergyCM(W, 0, Mt);
  p0CM = Energy2Mom(E0CM, 0);

  int i;
  double A, B, C;
  double betaFrac;
  for (i=0; i<=1; i++){ // loop for high- and low- momentum

    E1CM[i]    = HadronEnergyCM(W, M1, M2);
    p1CM[i]    = Energy2Mom(E1CM[i], M1);

    E2CM[i]    = HadronEnergyCM(W, M2, M1);
    p2CM[i]    = Energy2Mom(E2CM[i], M2);

    beta1CM[i] = calcBeta(p1CM[i], E1CM[i]);
    betaFrac   = betaCM/beta1CM[i];

    A = gammaCM2 * tantpq2 + 1;
    B = 2 * gammaCM2 * betaFrac * tantpq2 ;
    C = gammaCM2 * tantpq2 * betaFrac * betaFrac - 1;

    ctpqCM[i] = QuadraticSolution(A, B, C, -i);
    tpq1CM[i] = acos(ctpqCM[i]) ;
    tpq2CM[i] = M_PI - tpq1CM[i] ;

    // 3 momentum transfers in CMS
    qCM[i] = sqrt(p0CM*p0CM + p1CM[i]*p1CM[i] - 2*p0CM*p1CM[i]*ctpqCM[i]);

    Jacobian1[i] 
      = calcJacobianLab2CM(tpq1CM[i], betaCM, beta1CM[i]);

  }

  return 0;

}



//
// Function name : Calc'H'KinemaLab(double *5, double& *4)
//
// Description   : calculate exclusive kinematics in Lab frame
// Input         : q, omega, tpqCM, W, Mt, M1, M2
// Return        : 
//
int 
Hadron::CalcHKinemaLab(double q, double omega, double tpqCM, double W, 
                      double Mt, double Mass1, double Mass2){ 

  M1 = Mass1;
  M2 = Mass2;
  tpq1CM[0] = tpq1CM[1] = tpqCM;
  tpq2CM[0] = tpq2CM[1] = M_PI - tpq1CM[0] ;
  double ctpqCM[2];
  ctpqCM[0] = ctpqCM[1] = cos(tpqCM);
  double ttpqLab2[2]; // tan^2(theta_pq_lab)
  ttpqLab2[0] = ttpqLab2[1] = 0;

  betaCM  = calcbetaCM(q, omega, Mt);
  gammaCM = calcgammaCM(omega, W, Mt);
  double gammaCM2 = gammaCM * gammaCM;

  int i;
  double betaFrac;
  for (i=0; i<=1; i++){ // loop for high- and low- momentum

    E1CM[i]    = HadronEnergyCM(W, M1, M2);
    p1CM[i]    = Energy2Mom(E1CM[i], M1);
    E2CM[i]    = HadronEnergyCM(W, M2, M1);
    p2CM[i]    = Energy2Mom(E2CM[i], M2);

    beta1CM[i] = calcBeta(p1CM[i], E1CM[i]);
    betaFrac   = betaCM/beta1CM[i];

    ttpqLab2[i] = (1 - ctpqCM[i]*ctpqCM[i])/gammaCM2;
    ttpqLab2[i] = ttpqLab2[i]/(betaFrac + ctpqCM[i])/(betaFrac + ctpqCM[i]);

    tpq1Lab[i] = atan( sqrt(ttpqLab2[i]) );

    p1[i] = gammaCM*(p1CM[i]*ctpqCM[i] + betaCM*E1CM[i])/cos(tpq1Lab[i]);
    E1[i] = Mom2Energy(p1[i], M1);

    p2[i] = sqrt(p1[i]*p1[i] + q*q - 2*p1[i]*q*cos(tpq1Lab[i]));
    E2[i] = Mom2Energy(p2[i], M2);

    tpq2Lab[i] = - asin( p1[i]/p2[i] * sin(tpq1Lab[i]) ) ;
    // Additional 0.001 is necessary to test the sign. UGLY!!! //
    if ( p1[i]*cos(tpq1Lab[i]) + p2[i]*cos(tpq2Lab[i]) > q + 0.001 ) {
      tpq2Lab[i] = - M_PI - tpq2Lab[i] ;
    }

    Jacobian1[i] 
      = calcJacobianCM(gammaCM, p1[i], p1CM[i], betaCM, E1CM[i], tpq1Lab[i]);
    Jacobian2[i] 
      = calcJacobianCM(gammaCM, p2[i], p2CM[i], betaCM, E2CM[i], tpq2Lab[i]);

  }

  return 0;

}





//
// Function name : CalcgKinemaLAB(double *5, double& *4)
//
// Description   : calculate gamma kinematics in LAB frame
// Input         : E0, M1, t1, M2
// Return        : Q, omega, tq, p1
//
int 
Hadron::CalcgKinemaLAB(double E0, double M1, double t1, double M2, 
		     double &Q, double &omega, double &tq, double &p1){

  double ct1 = cos(t1) ;
  double st1 = sqrt(1 - cos(t1)*cos(t1)) ;

  double A = 2*( M2 + E0 );
  double B = (2*E0*E0 - M1*M1 )/A ;
  double C = 2*E0*ct1/A;

  double D = C*C - 1;
  double E = 2*C*( E0 - B );
  double F = E0*E0 - 2*E0*B + B*B - M1*M1; 

  p1 = QuadraticSolution(D, E, F, -1);
  double E1 = Mom2Energy(p1,M1);

  double q2 = E0*E0 + p1*p1 - 2*E0*p1*ct1;
  double q = sqrt(q2);

  omega = B - C*p1 ;  // omega = -(E0 - E1)

  Q = sqrt(q*q - omega*omega);

  double stq = -p1/q*st1 ;
  tq = asin(stq) ;

  return 0;

}

/* Disabled May 11, 2004
   B = 2*( E0*E0 - M1*M1 )/A  is wrong, It supposed to be B = (2*E0*E0 - M1*M1 )/A
//
// Function name : CalcgKinemaLAB(double *5, double& *4)
//
// Description   : calculate gamma kinematics in LAB frame
// Input         : E0, M1, t1, M2
// Return        : Q, omega, tq, p1
//
int 
Hadron::CalcgKinemaLAB(double E0, double M1, double t1, double M2, 
		       double &Q, double &omega, double &tq, double &p1){

  double ct1 = cos(t1) ;
  double st1 = sqrt(1 - cos(t1)*cos(t1)) ;

  double A = 2*( M2 - E0 );
  double B = 2*( E0*E0 - M1*M1 )/A ;
  double C = 2*E0*ct1/A;

  double D = C*C - 1;
  double E = -2*C*( E0 + B );
  double F = E0*E0 + 2*E0*B + B*B - M1*M1; 

  p1 = QuadraticSolution(D, E, F, -1);

  double q2 = E0*E0 + p1*p1 - 2*E0*p1*ct1;
  double q = sqrt(q2);

  omega = B - C*p1 ;  // omega = -(E0 - E1)

  Q = sqrt(q*q - omega*omega);

  // Redefine omega as standard definition "omega = E0 - E1".
  // See PrimEx Logbook P88 for details.
  omega = -omega;     

  double stq = -p1/q*st1 ;
  tq = asin(stq) ;

  return 0;

}
*/

//
// Function name : CalcgKinema3bodyLAB(double *7, double& *4)
//
// Description   : calculate kinematics (gamma,Mt) -> (M1,M2,M3) in LAB frame
//                 presently it assumes Mt=M3 and M1=M2
// Input         : E0, M3, M1, t1, E1, M2, t2
// Return        : Q, omega, tq, p1
//
int 
Hadron::CalcgKinema3bodyLAB(double E0, double M3, double M1, double t1, double E1, double M2, double t2,
				   double &Q, double &omega, double &tq, double &p2){

  double ct1  = cos(t1) ;
  double ct2  = cos(t2) ;
  double ct12 = cos(t1-t2);
  
  double st1 = sin(t1);
  double st2 = sin(t2);

  double p1 = Energy2Mom(E1,M1);

  // here M3=Mt, M1=M2 assumed
  double A = -(M1*M1 + M2*M2 - 2*(E0*E1 - E0*p1*ct1))/2/M3;

  double B = 1 + E0/M3 + E1/M3 ;
  double C = (p1*ct12 + E0*ct2)/M3 ;
  double D = E1 + A - E0 ;

  double AA = B*B - C*C ;
  double BB = 2*C*D ;
  double CC = B*B*M2*M2 - D*D ;

  p2 = QuadraticSolution(AA, BB, CC, -1);
  double E2 = Mom2Energy(p2, M2);

  //omega = A + (E0*E2 - E0*p2*ct2)/M3 + (E1*E2 - p1*p2*ct12)/M3 ; this is redundant
  omega = E0-E1-E2 ;
  double q = sqrt(omega*omega + 2*omega*M3) ;

  Q = sqrt(q*q - omega*omega) ;
  double stq = (p1*st1 + p2*st2)/q ;
  tq = asin(stq);

  cout << p2*M2G << " " << E2*M2G << " " << omega+E1+E2 << " " << E0-E1-E2 << " " << q*M2G << " " << Q*M2G << " " << tq*r2d << endl;

  return 0;

}



//
// Function name : PrintgKinema
//
// Description   : print gamma kinematics
// Input         : 
// Return        : 
//
void
Hadron::PrintgKinema(int PrintMode, int CMS, char target[], double E0, double A, int Z, 
		     double Mt, double M1, double t1, double M2, double Q, double omega, 
		     double tq, double p1){

  // Get Kinematic Energy for p1
  double T1 = Mom2Energy(p1,M1) - M1;

  if (PrintMode&1){
    printf("-----------------------------------------------------------------------------\n");
    printf("    E0   target   t1      Q       omega    tq      p1     T1     t1CM  p1CM\n");
    printf("  [GeV]         [deg.] [GeV/c]    [MeV]  [deg.] [GeV/c] [GeV]   [deg] [GeV/c]\n");
    printf("-----------------------------------------------------------------------------\n");
  }

  switch(PrintMode>>3){
  case 1:
    printf("%7.3f",E0*M2G);
    if (PrintMode&1) printf("%6s ",target);
    printf("%7.3f",t1*r2d);
    printf("%8.3f",Q*M2G);
    printf("%10.3f",omega);
    printf("%8.2f",tq*r2d);
    printf("%7.3f",p1*M2G);
    printf("%7.3f",T1*M2G);
    printf("%8.3f",tpq1CM[0]*r2d); //forward branch
    printf("%7.3f",p1CM[0]*M2G);
    printf("\n");
    break;
  case 3: // No longer in use
    printf("%7.1f",E0);
    printf("%10.7f",t1);
    printf("%13.5f",Q);
    printf("%13.5f",omega);
    printf("%10.7f",tq);
    printf("%13.5f",p1);
    printf("%12.5f",M1);
    printf("%9.3f",A);
    printf("%5d",Z);
    printf("%12.3f",Mt);
    printf("%10.3f",p1CM[0]);
    printf("%10.6f",tpq1CM[0]); //forward branch
    printf("%10.6f",tpq1CM[1]); //backward branch
    printf("%10.3f",E0CM);
    printf("%10.3f",qCM[0]);
    printf("%10.3f",qCM[1]);
    printf("\n");
    break;
  case 0:

    CMS ?  printf("%10.3f",E0CM)      : printf("%7.1f",E0);
    CMS ?  printf("%11.8f",tpq1CM[0]) : printf("%11.8f",t1);
    printf("%13.5f",Q);
    printf("%13.5f",omega);
    CMS ?  printf("%12.7f",tpq1CM[0]) : printf("%12.7f",tq);
    CMS ?  printf("%10.3f",p1CM[0])   :printf("%13.5f",p1);
    printf("%12.5f",M1);
    printf("%9.3f",A);
    printf("%5d",Z);
    printf("%12.3f",Mt);

    printf("%7.1f",E0);
    printf("%12.8f",t1);
    printf("%3d",CMS);

    printf("\n");
    break;
  }

  if (PrintMode&3) printf("\n");

  return ;




}



