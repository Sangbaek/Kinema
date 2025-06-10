#ifndef CALC_H
#define CALC_H

int CalceKinema(double Mt, double Ee, double te, 
		double &, double &, double &, double &, double &, double &, double &) ;
int PrinteKinema(double Mt, double Ee, double te, double Rec, double Eep, 
		 double q2, double Q2, double tq, double ERec, double PRec) ;

int CalcpKinema(double Mt, double Ee, double tq, 
		double &, double &, double &, double &, double &, double &, double &) ;

int CalcPpKinema(double Mt, double Ee, double PRec, 
		double &, double &, double &, double &, double &, double &, double &) ;

int CalcExKinema(double Mt, double Ee, double te, double Eep, double &Rec, 
		 double &Ex, double &q2, double &tq, double &ERec, double &PRec) ;

int CalcWQKinema(double Mt, double Ee, double te, double Eep,
		 double &omega, double &W, double &q, double &Q2, double &tq);

int CalcteKinema(double Mt, double Ee, double te, double W,
		 double &, double &, double &, double &, double &);

int CalcKinema(double Mt, double Ee, double Q2, double W,
	       double &, double &, double &, double &, double &);
int PrintKinema(double Mt, double Ee, double Q2, double W,
	        double omega, double Eep, double q, double te, double tq);

class Hadron
{
 private:
  double gammaCM, betaCM; // photon kinematics
  double qCM[2], p0CM, E0CM; // incident momentum and energy in CMS
  double M1, p1[2], E1[2], M2, p2[2], E2[2]; // laboratory frame
  double p1CM[2], E1CM[2], p2CM[2], E2CM[2], beta1CM[2], beta2CM[2];
  double tpq1Lab[2], tpq2Lab[2], tpq1CM[2], tpq2CM[2];
  double Jacobian1[2], Jacobian2[2];
 public:
  int CalcHKinema(double q, double omega, double tpq, double Mt, double M1, double M2);
  int CalcHKinema(double W,  double q, double omega, double tpqCM, 
		  double Mt, double Mass1, double Mass2);

  int PrintHKinema(int Mode);

  int CalcHKinemaCM(double q, double omega, double tpq, double W, 
		    double Mt, double M1, double M2); 

  int CalcHKinemaCM(double q, double tpq, double W, 
		    double Mt, double M1, double M2); 

  int CalcHKinemaLab(double q, double omega, double tpqCM, double W, 
		     double Mt, double M1, double M2); 

  int CalcgKinemaLAB(double E0, double M1, double t1, double M2, 
		     double &Q, double &omega, double &tq, double &p1);

  int CalcgKinema3bodyLAB(double E0, double M3, double M1, double t1, double E1, double M2, double t2,
			  double &Q, double &omega, double &tq, double &p2);

  void PrintgKinema(int PrintMode, double E0, double M1, double t1, double M2, 
		   double  q, double omega, double tq, double p1);

  void PrintgKinema(int PrintMode, int CMS, char target[], double E0, double A, int Z, double Mt, double M1, 
		    double t1, double M2, double Q, double omega, double tq, double p1);

  /* return th*_piq: (HL>0 high momentum brunch, HL<0 low momentum branch) */
  double getHKinema(int HL){ 
    return HL<0 ? tpq2CM[1]: tpq2CM[0] ;
  };
  /* return th_pq: (HL>0 high momentum brunch, HL<0 low momentum branch) */
  double getHKinemaLab(int HL){ 
    return HL<0 ? tpq1Lab[1]: tpq1Lab[0] ;
  };

  void getHadronKinema(double _p1[2], double _tpq1Lab[2], double _tpq1CM[2]){ 

    _p1[0]      = p1[0];
    _p1[1]      = p1[1];
    _tpq1Lab[0] = tpq1Lab[0];
    _tpq1Lab[1] = tpq1Lab[1];
    _tpq1CM[0]  = tpq1CM[0];
    _tpq1CM[1]  = tpq1CM[1];

    return;
  };

  void getCMSKinema(double &_tq, double &_theta, double &_E, double &_p){

      _tq    = -tpq1CM[0];
      _theta =  tpq1CM[0];
      _E     =  E0CM;
      _p     =  p1CM[0];

    return;
  };


};


#endif /* CALC_H */

