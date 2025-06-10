#include <string.h>

//
// Class  name  : 
// Method name : GetMaterialAZ(char target[], double &A, int &Z)
//
// Description : return A and Z of materical char "target"
//             : 
// Input       : target
// Return      : A, Z
//

void
GetMaterialAZ(char target[], double &A, int &Z){


  if ((!strncmp("Pb",  target,2))||(!strncmp("pb",  target,2))) {
    Z = 82;
    A = 207.19;
  } else if ((!strncmp("12C",  target,3))||(!strncmp("C", target,1))||(!strncmp("c", target,1))) {
    Z = 6;
    A = 12.01;
  } else if ((!strncmp("Sn",  target,2))||(!strncmp("sn", target,1))) {
    Z = 50;
    A = 118.69;
  } else if ((!strncmp("Be",  target,2))||(!strncmp("be", target,1))) {
    Z = 4;
    A = 9.01;
  } else if ((!strncmp("4He",  target,3))||(!strncmp("He",  target,2))||(!strncmp("he", target,1))) {
    Z = 2;
    A = 4.00;

  } else if ((!strncmp("H",  target,3))||(!strncmp("H",  target,2))||(!strncmp("h", target,1))) {
    Z = 1;
    A = 1.00794;
  }
    
  return ;

}
