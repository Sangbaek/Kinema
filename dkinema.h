//
// dkinema.h
//

#ifndef DKINEMA_H
#define DKINEMA_H

void Usage();
void getParameters(double &p1, double &p2);

double CalcDtq(double Ee, double Eep, double te, double   q, double   tq,
	       double dte, double dEep) ;

void printDkinema(double te, double dte, double Eep, double dEep, 
		  double tq, double dtq);

#endif /* DKINEMA_H */

