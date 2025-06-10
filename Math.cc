//
//   Math.cc
//

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include "Math.h"

using namespace std;

//
// Function name : QuadraticSolution(double A, double B, double C, double pm)
//
// Description   : calculate quadratic solusion: A*x*x + B*x + C = 0
// Input         : coefficiencies A, B, C , pm
// Return        : returns plus solution x if pm>=0, else minus solution x 
//
double
QuadraticSolution(double A, double B, double C, double pm){

  if ( B*B<4*A*C ) {
    cerr << "Math.cc: complex solusion in QuadraticSolution" << endl;
    exit(-1);
  }

  return pm>=0 ? (-B + sqrt(B*B - 4*A*C))/2/A : (-B - sqrt(B*B - 4*A*C))/2/A ;

}

//
// Function name : QuadraticSolution(double A, double B, double C)
//
// Description   : calculate quadratic solusion: A*x*x + B*x + C = 0
// Input         : coefficiencies A, B, C 
// Return        : x[0]=plus solution, x[1]=minus solution 
//
double *
QuadraticSolution(double A, double B, double C){

  if ( B*B<4*A*C ) {
    cerr << "Math.cc: complex solusion in QuadraticSolution" << endl;
    exit(-1);
  }

  double x[2];
  x[0] = (-B + sqrt(B*B - 4*A*C))/2/A; 
  x[1] = (-B - sqrt(B*B - 4*A*C))/2/A; 

  return x;
}



//
// Function name : det(double (*a)[2])
//
// Description   : calculate determinant A
// Input         : 2x2 Matrix a[2][2]
// Return        : detA
//
double 
det(double (*a)[2]){

  return a[0][0]*a[1][1] - a[0][1]*a[1][0] ;

} // End-of-det()

//
// Function name : det(double (*a)[3])
//
// Description   : calculate determinant A
// Input         : 3x3 Matrix a[3][3]
// Return        : detA
//
double 
det(double (*a)[3]){

  double A0[2][2] = {{a[1][1], a[1][2]}, {a[2][1], a[2][2]}};
  double A1[2][2] = {{a[1][0], a[2][0]}, {a[1][2], a[2][2]}};
  double A2[2][2] = {{a[1][0], a[1][1]}, {a[2][0], a[2][1]}};

  return a[0][0]*det(A0) - a[0][1]*det(A1) + a[0][2]*det(A2) ;

} // End-of-det()



//
// Function name : InvMatrix(double (*a)[2], double b[2], double x[2]){
//
// Description   : calculate Solutions for 2x2 Matrix
// Input         : 2x2 Matrix a[2][2], b[2]
// Return        : x[2]
//
void
InvMatrix(double (*a)[2], double b[2], double x[2]){

  double detA = det(a);

  if (!detA) {
    cerr << " Math: Error. Determinant 2x2 Matrix A equal zero" << endl;
    return ;
  } else {
    x[0] = (a[1][1]*b[0] - a[0][1]*b[1])/detA;
    x[1] = (a[0][0]*b[1] - a[1][0]*b[0])/detA;
  }

  return ;

}


//
// Function name : InvMatrix(double (*a)[3], double b[3], double x[3]){
//
// Description   : calculate Solutions for 3x3 Matrix A=a[3][3] -> Ax=b
// Input         : 3x3 Matrix a[3][3], b[3]
// Return        : x[3]
//
void
InvMatrix(double (*a)[3], double b[3], double x[3]){

  double detA = det(a);

  if (!detA) {
    cerr << " Math: Error. Determinant 3x3 Matrix A equal zero" << endl;
    return ;
  } else {
    double c[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};

    c[0][0] = a[1][1]*a[2][2] - a[2][1]*a[1][2];
    c[0][1] = a[2][1]*a[0][2] - a[0][1]*a[2][2];
    c[0][2] = a[0][1]*a[1][2] - a[1][1]*a[0][2];
    c[1][0] = a[2][0]*a[1][2] - a[1][0]*a[2][2];
    c[1][1] = a[0][0]*a[2][2] - a[2][0]*a[0][2];
    c[1][2] = a[1][0]*a[0][2] - a[0][0]*a[1][2];
    c[2][0] = a[1][0]*a[2][1] - a[1][1]*a[2][0];
    c[2][1] = a[0][1]*a[2][0] - a[0][0]*a[2][1];
    c[2][2] = a[0][0]*a[1][1] - a[0][1]*a[1][0];

    x[0] = (c[0][0]*b[0] + c[0][1]*b[1] + c[0][2]*b[2])/detA;
    x[1] = (c[1][0]*b[0] + c[1][1]*b[1] + c[1][2]*b[2])/detA;
    x[2] = (c[2][0]*b[0] + c[2][1]*b[1] + c[2][2]*b[2])/detA;

  }

  return ;

}


//
// Function name : InvMatrix(double (*a)[3], double b[3], double x[3]){
//
// Description   : calculate Solutions for 3x3 Matrix A=a[3][3] -> Ax=b
//               : Inverse Matrix C=A^-1. x=Cb
// Input         : a[3][3], b[3]
// Return        : c[3][3], x[3] 
//
void
InvMatrix(double (*a)[3], double b[3], double x[3], double (*c)[3]){

  double detA = det(a);

  if (!detA) {
    cerr << " Math: Error. Determinant 3x3 Matrix A equal zero" << endl;
    return ;
  } else {

    c[0][0] = (a[1][1]*a[2][2] - a[2][1]*a[1][2])/detA;
    c[0][1] = (a[2][1]*a[0][2] - a[0][1]*a[2][2])/detA;
    c[0][2] = (a[0][1]*a[1][2] - a[1][1]*a[0][2])/detA;
    c[1][0] = (a[2][0]*a[1][2] - a[1][0]*a[2][2])/detA;
    c[1][1] = (a[0][0]*a[2][2] - a[2][0]*a[0][2])/detA;
    c[1][2] = (a[1][0]*a[0][2] - a[0][0]*a[1][2])/detA;
    c[2][0] = (a[1][0]*a[2][1] - a[1][1]*a[2][0])/detA;
    c[2][1] = (a[0][1]*a[2][0] - a[0][0]*a[2][1])/detA;
    c[2][2] = (a[0][0]*a[1][1] - a[0][1]*a[1][0])/detA;

    x[0] = c[0][0]*b[0] + c[0][1]*b[1] + c[0][2]*b[2];
    x[1] = c[1][0]*b[0] + c[1][1]*b[1] + c[1][2]*b[2];
    x[2] = c[2][0]*b[0] + c[2][1]*b[1] + c[2][2]*b[2];


  }

  return ;

}

//
// Class name  : 
// Method name : LinearInterporate(double x1, y1,  double x2, y2, double x)
//
// Description : Linear Interporate between (x1,y1) and (x2,y2) for (x,y)
// Input       : x1,y1, x2,y2, x
// Return      : y
//

double
LinearInterporate(double x1, double y1, double x2, double y2, double x){

  if ( (x2 - x1) ) {
    return (y2 - y1)/(x2 - x1)*(x - x1) + y1;
  } else {
    cerr << " Error! LinearInterporate: denominator zero " << endl;
    exit (-1);
  }

} // end of LinearInterporate()




//
// Function name : Abs(double x)
//
// Description   : return absolute of x
// Input         : x
// Return        : abs(x)
//
double 
Abs(double X) { 
  return X>=0 ? X : -(X) ; 
}



