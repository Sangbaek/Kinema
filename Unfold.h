#ifndef UNFOLD_H
#define UNFOLD_H
const int SIZE = 1000;
const int FSIZ =   32;

class Unfold
{
 private:
  float Ex_x[SIZE], Cnts_x[SIZE], dCnts_x[SIZE];
  float Ex_c[SIZE], Cnts_c[SIZE];
  float dEx_x, dEx_c;
  float LSsum;
  float LS[SIZE];
  int LSbin[2];
  int Nx, Nc;
  int at0, atThr;

 public:
  void Initialize();
  void Unfolding(char OutFile[FSIZ]);
  int getData(char DataFile[FSIZ]); 
  int getLineShape(char LSFile[FSIZ]);
  int SearchX1prime(float x);
  float getExpArea(int bin);
  float getLSArea(float x1, float x2);
  int FillLS(float dEx);

};

float LinearInterporate(float x1, float y1, float x2, float y2, float x);
float LinearIntegration(float x1, float y1, float x2, float y2);


//
// Class name  : Unfold
// Method name : Initialize
//
// Description : Initialize arrays
// Input       : 
// Return      : 
//

void
Unfold::Initialize(){
  int j;
  for (j=0; j<=SIZE; ++j) {
    Ex_x[j]    = 0;
    Cnts_x[j]  = 0;
    dCnts_x[j] = 0;
    Ex_c[j]    = 0;
    Cnts_c[j]  = 0;
    LS[j]      = 0;
  }
}

#endif /* UNFOLD_H */
