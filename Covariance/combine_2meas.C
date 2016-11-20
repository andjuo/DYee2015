#include "inputs.h"
#include "Blue.h"

void combine_2meas(int theCase=1, int corrAll=0, double scale=1.)
{
  TMatrixD measEE(2,1), measMM(2,1);
  TMatrixD covEEtot(2,2), covMMtot(2,2), covEM(2,2);
  std::vector<TMatrixD*> uncEEv, uncMMv;
  std::vector<TMatrixD*> covEEv, covMMv;
  std::vector<TMatrixD*> covKindV;
  std::vector<TString> labelV;

  if (valueEquals(theCase,"1 -1 10 11")) {
    // TOP-13-008
    assignValuesT(measEE,0,"0.705 0.304");
    assignValuesT(measMM,0,"0.685 0.328");
    measEE *= scale;
    measMM *= scale;
    TMatrixD uncStat(2,2), uncSyst(2,2);
    uncStat(0,0)=0.013;
    uncStat(1,1)=0.009;
    uncSyst(0,0)=0.037;
    uncSyst(1,1)=0.020;
    if (theCase!=10) {
      uncEEv.push_back(new TMatrixD(uncStat));
      covEEv.push_back(new TMatrixD(uncStat,TMatrixD::kMult,uncStat));
      labelV.push_back("stat");
    }
    if (theCase!=-1) {
      uncEEv.push_back(new TMatrixD(uncSyst));
      covEEv.push_back(new TMatrixD(uncSyst,TMatrixD::kMult,uncSyst));
      labelV.push_back("syst");
    }
    double rho=-0.950;
    for (unsigned int i=0; i<covEEv.size(); i++) {
      (*uncEEv[i]) *= scale;
      (*covEEv[i]) *= scale*scale;
      TMatrixD *m= covEEv[i];
      if (corrAll || (labelV[i]=="stat")) {
	double v= rho* sqrt((*m)(0,0) * (*m)(1,1));
	(*m)(0,1)=v;
	(*m)(1,0)=v;
      }
      covEEtot += (*m);
    }
    uncStat(0,0)=0.013;
    uncStat(1,1)=0.009;
    uncSyst(0,0)=0.024;
    uncSyst(1,1)=0.014;
    if (theCase!=10) {
      uncMMv.push_back(new TMatrixD(uncStat));
      covMMv.push_back(new TMatrixD(uncStat,TMatrixD::kMult,uncStat));
    }
    if (theCase!=-1) {
      uncMMv.push_back(new TMatrixD(uncSyst));
      covMMv.push_back(new TMatrixD(uncSyst,TMatrixD::kMult,uncSyst));
    }
    rho=-0.957;
    for (unsigned int i=0; i<covMMv.size(); i++) {
      (*uncMMv[i]) *= scale;
      (*covMMv[i]) *= scale*scale;
      TMatrixD *m= covMMv[i];
      if (corrAll || (labelV[i]=="stat")) {
	double v= rho* sqrt((*m)(0,0) * (*m)(1,1));
	(*m)(0,1)=v;
	(*m)(1,0)=v;
      }
      covMMtot += (*m);
    }
  }
  else {
    std::cout << "not ready for case " << theCase << "\n";
    return;
  }


  BLUEResult_t blue;
  if (!blue.estimate(measEE,covEEtot,measMM,covMMtot,covEM)) {
    std::cout << "estimate failed\n";
    return;
  }

  TMatrixD measLL( *blue.getEst() );
  TMatrixD covLL( *blue.getCov() );

  for (int ir=0; ir<measLL.GetNrows(); ir++) {
    std::cout << "combined value: " << measLL(ir,0) << " +- " << sqrt(covLL(ir,ir)) << "\n";
  }
  printVec("included systematics : ",labelV);

  for (int ir=0; ir<measLL.GetNrows(); ir++) {
    std::cout << "combined value: " << measLL(ir,0);
    for (unsigned int iSrc=0; iSrc<labelV.size(); iSrc++) {
      TMatrixD covSrc(4,4);
      covSrc += blue.measCovMatrix( *covEEv[iSrc], 1 );
      covSrc += blue.measCovMatrix( *covMMv[iSrc], 0 );
      TMatrixD covFinSrc= blue.contributedCov(covSrc);
      std::cout << " +- " << sqrt(covFinSrc(ir,ir));
      std::cout << " (" << labelV[iSrc] << ")";
    }
    std::cout << "\n";
  }

  //TMatrixD covEEFinal( blue.contributedCov(blue.measCovMatrix(covEE,1)) );
  //TMatrixD covMMFinal( blue.contributedCov(blue.measCovMatrix(covMM,0)) );

}
