#include <iostream.h>

void
calcXXX(double y, double x[]){

  x[0] = y;
  x[1] = -y;

  return;

}

void 
main(){

  double y=1;
  double x[2];
  x[0] = x[2] = 0;

  calcXXX(y, x);

  cout << x[0] << " " << x[1] << endl;
}

