#ifndef FF_H
#define FF_H

// Maximum row of model data 
const int MAXDATA=500;

class Cornell
{

public:
  double NuclearCoherentFF(double A, int Z, double q, double tq);
  double CoulombFF(double A, int Z, double q, double tq);
  double NuclearIncoherentFF(double A, int Z, double theta);

};


class Sergey
{

public:
  double NuclearCoherentFF(double A, int Z, double q, double tq);

};

class Cascade
{

public:

  double NuclearIncoherentFF(double A, int Z, double theta);

};



#endif /* FF_H */









