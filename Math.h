
#ifndef MATH_H
#define MATH_H

double * QuadraticSolution(double A, double B, double C);

double QuadraticSolution(double A, double B, double C, double pm);

double det(double (*a)[2]);
double det(double (*a)[3]);

void InvMatrix(double (*a)[2], double b[2], double x[2]);
void InvMatrix(double (*a)[3], double b[3], double x[3]);
void InvMatrix(double (*a)[3], double b[3], double x[3], double (*c)[3]);

double LinearInterporate(double x1, double y1, double x2, double y2, double x);

double Abs(double X);

#endif /* MATH_H */



