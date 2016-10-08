// Examples from Valassi and Chierici, 2014

#include "Blue.h"

void myTest2014(int theCase=1)
{
  TMatrixD *measAinp=NULL, *measBinp=NULL;
  TMatrixD *measTot=NULL;
  TMatrixD *Uinp=NULL, *covInp=NULL;
  int nMeas=2;
  TString caseName="unknown";
  Info_t info;
  info.add(Form("Case number in the macro theCase=%d",theCase));

  if (theCase==1) {
    caseName="2";
    nMeas=2;
    measAinp= new TMatrixD(1,1);
    measBinp= new TMatrixD(1,1);
    assignValues(*measAinp,0, "103.00");
    assignValues(*measBinp,0, "98.0");
    measTot= new TMatrixD(2,1);
    assignValuesT(*measTot,0, "103. 98.");
    covInp= new TMatrixD(2,2);
    assignValues(*covInp, 0, "15.00  0.  ");
    assignValues(*covInp, 1, " 0.   10.00");
    Uinp= new TMatrixD(2,1);
    assignValuesT(*Uinp, 0, "1 1");
    info.add("Eq.3");
  }
  else if (theCase==2) {
    caseName="2 (B is 2 measurements)";
    nMeas=2;
    measAinp= new TMatrixD(1,1);
    measBinp= new TMatrixD(2,1);
    assignValues(*measAinp,0, "103.00");
    assignValuesT(*measBinp,0, "99.0 101.");
    measTot= new TMatrixD(3,1);
    assignValuesT(*measTot,0, "103. 99. 101.");
    covInp= new TMatrixD(3,3);
    assignValues(*covInp, 0, "15.00  0.    0.  ");
    assignValues(*covInp, 1, " 0.   16.00 28.00");
    assignValues(*covInp, 2, " 0.   28.00 64.00");
    Uinp= new TMatrixD(3,1);
    assignValuesT(*Uinp, 0, "1 1 1");
    info.add("Eq.4");
  }

  BLUEResult_t blue;
  if (!blue.estimate(*measTot,*Uinp,*covInp)) return;

  print("lambda ",blue.getLambda());
}

