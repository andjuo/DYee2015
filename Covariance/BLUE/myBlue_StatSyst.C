#include "myBlue.C"

void myBlue_StatSyst(int theCase, int printInfo=0)
{
  TMatrixD *measAinp=NULL, *measBinp=NULL;
  TMatrixD *covInp=NULL, *covStat=NULL, *covSyst=NULL;
  int nMeas=2;
  TString caseName="unknown";
  Info_t info;
  info.add(Form("myBlue_StatSyst: Case number in the macro theCase=%d",theCase));

  if ((theCase==1) || (theCase==2) || (theCase==3)) {
    caseName="4.3.3. Breakdown of error contributions";
    if (theCase==1) caseName.Append(" (uncorrelated)");
    else if (theCase==2) caseName.Append(" (correlated)");
    else if (theCase==3) caseName.Append(" (anti-correlated)");
    nMeas=2;
    measAinp= new TMatrixD(2,1);
    measBinp= new TMatrixD(2,1);
    assignValuesT(*measAinp,0,"10.50 9.50");
    assignValuesT(*measBinp,0,"13.50 14.00");
    covInp= new TMatrixD(4,4);
    covStat= new TMatrixD(*covInp);
    covSyst= new TMatrixD(*covInp);
    info.add("Eq.(51) values:");
    info.add("B_A^e=10.50 +- 1.00(stat)");
    info.add("B_B^e=13.50 +- 0.21(stat) +- 2.99(syst)");
    info.add("B_A^t= 9.50 +- 3.00(stat)");
    info.add("B_B^t=14.00 +- 0.21(stat) +- 2.99(syst)");
    info.add("To get the previous case, stat^2+syst^2=9");
    assignValues(*covStat,0,"1.00 0.     0.   0.   ");
    assignValues(*covStat,1,"0.   0.0441 0.   0.   ");
    assignValues(*covStat,2,"0.   0.     9.00 0.   ");
    assignValues(*covStat,3,"0.   0.     0.   0.0441");
    assignValues(*covSyst,0,"0.   0.      0.   0.   ");
    assignValues(*covSyst,1,"0.   8.9559  0.   8.9559");
    assignValues(*covSyst,2,"0.   0.      0.   0.   ");
    assignValues(*covSyst,3,"0.   8.9559  0.   8.9559");
    if (0)
    for (int ir=0; ir<covStat->GetNrows(); ir++) {
      for (int ic=0; ic<covStat->GetNcols(); ic++) {
	(*covStat)(ir,ic)= 0.01*trunc(100.*(*covStat)(ir,ic));
	(*covSyst)(ir,ic)= 0.01*trunc(100.*(*covSyst)(ir,ic));
      }
    }
    if (theCase==1) {
      (*covSyst)(1,3)=0;
      (*covSyst)(3,1)=0;
    }
    if (theCase==3) {
      (*covSyst)(1,3)=-(*covSyst)(1,3);
      (*covSyst)(3,1)=-(*covSyst)(3,1);
    }
    (*covInp)=(*covStat)+(*covSyst);
  }
  else {
    std::cout << "not ready for case=" << theCase << "\n";
    return;
  }

  info.caseName=caseName;

  combine(*measAinp,*measBinp,*covInp,(printInfo) ? &info : NULL);
}
