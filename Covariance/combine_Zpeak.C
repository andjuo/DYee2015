#include "inputs.h"
#include "Blue.h"

void combine_Zpeak(int theCase=1, int subVer=0)
{
  double sigmaZee=0, sigmaZmm=0;
  std::vector<double> dSigmaEE,dSigmaMM,dSigmaEM;
  std::vector<TString> labelV;
  double dLumiRel=2;

  if (valueEquals(theCase,"1 -1 10 11")) {
    // SMP-15-004-pas
    sigmaZee= 1920.;
    addToVector(dSigmaEE,2,20.,60.);
    sigmaZmm= 1900.;
    addToVector(dSigmaMM,2,10.,50.);
    addToVector(dSigmaEM,2,0.,0.8*sqrt(50*60));
    addToVector(labelV,"stat syst");
    if ((theCase==-1) || (theCase==10)) {
      dSigmaEE.clear(); dSigmaEE.push_back(20);
      dSigmaMM.clear(); dSigmaMM.push_back(10);
      dSigmaEM.clear(); dSigmaEM.push_back(0);
      labelV.clear(); labelV.push_back("stat");
    }
    if ((theCase==11) || (theCase==10)) {
      addToVector(dSigmaEE,1,90.);
      addToVector(dSigmaMM,1,90.);
      addToVector(dSigmaEM,1,90.);
      addToVector(labelV,"lumi");
    }
  }
  else if (valueEquals(theCase,"7 77 777 70 -7")) {
    // DY 7TeV
    sigmaZee= 984.6;
    addToVector(dSigmaEE,2,0.9,7.3);
    sigmaZmm= 989.5;
    addToVector(dSigmaMM,2,0.8,9.8);
    addToVector(dSigmaEM,2,0.,0.);
    addToVector(labelV,"stat expSyst");
    if (theCase==-7) {
      dSigmaEE.clear(); addToVector(dSigmaEE,1,0.9);
      dSigmaMM.clear(); addToVector(dSigmaMM,1,0.8);
      dSigmaEM.clear(); addToVector(dSigmaEM,1,0.);
      labelV.clear(); addToVector(labelV,"stat");
    }
    if ((theCase==70) || (theCase==777)) {
      addToVector(dSigmaEE,1,21.4);
      addToVector(dSigmaMM,1,21.9);
      addToVector(dSigmaEM,1,sqrt(21.4*21.9));
      addToVector(labelV,"th.syst");
    }
    if ((theCase==77) || (theCase==777)) {
      addToVector(dSigmaEE,1,21.7);
      addToVector(dSigmaMM,1,21.8);
      addToVector(dSigmaEM,1,sqrt(21.7*21.8));
      addToVector(labelV,"lumi");
    }
  }
  else if (valueEquals(theCase,"8 80 88 -8")) {
    // DY 8TeV
    sigmaZee= 1141;
    addToVector(dSigmaEE,2,11.,25.);
    sigmaZmm= 1135;
    addToVector(dSigmaMM,2,11.,25.);
    addToVector(labelV,"exp theo");
    addToVector(dSigmaEM,2,0.,25.);
    if ((theCase==-8) || (theCase==80)) {
      dSigmaEE.clear(); addToVector(dSigmaEE,1,11.);
      dSigmaMM.clear(); addToVector(dSigmaMM,1,11.);
      dSigmaEM.clear(); addToVector(dSigmaEM,1,0.);
      labelV.clear(); addToVector(labelV,"exp");
    }
    if ((theCase==88) || (theCase==80)) {
      addToVector(dSigmaEE,1,30.);
      addToVector(dSigmaMM,1,30.);
      addToVector(dSigmaEM,1,30.);
      addToVector(labelV,"lumi");
    }
  }
  else if ((subVer==0) && valueEquals(theCase,"13 1313 131313 130 -13")) {
    // DY 13TeV
    sigmaZee= 1861.0;
    addToVector(dSigmaEE,2,2.,1000.);
    sigmaZmm= 1860.0;
    addToVector(dSigmaMM,2,1.,14.);
    addToVector(dSigmaEM,2,0.,0.);
    addToVector(labelV,"stat expSyst");
    if (theCase==-13) {
      dSigmaEE.clear(); addToVector(dSigmaEE,1,2.);
      dSigmaMM.clear(); addToVector(dSigmaMM,1,1.);
      dSigmaEM.clear(); addToVector(dSigmaEM,1,0.);
      labelV.clear(); addToVector(labelV,"stat");
    }
    if ((theCase==130) || (theCase==131313)) {
      addToVector(dSigmaEE,1,25.);
      addToVector(dSigmaMM,1,25.);
      addToVector(dSigmaEM,1,sqrt(25.*25.));
      addToVector(labelV,"th.syst");
    }
    if ((theCase==1313) || (theCase==131313)) {
      addToVector(dSigmaEE,1,round(sigmaZee*0.027));
      addToVector(dSigmaMM,1,round(sigmaZmm*0.027));
      addToVector(dSigmaEM,1,sqrt(dSigmaEE.back()*dSigmaMM.back()));
      addToVector(labelV,"lumi");
    }
  }
  else if ((subVer==1) && valueEquals(theCase,"13 1313 131313 130 -13")) {
    // DY 13TeV
    // ee: 1949 +- 2 (stat) +- 28 (exp.syst) +- 32 (theo) +- 45 (lumi) pb
    // mm: 1907 +- 1 (Stat) +- 14 (exp.syst) +- 29 (theo) +- 44 (lumi) pb
    dLumiRel=0.023;
    sigmaZee= 1949.0;
    addToVector(dSigmaEE,2,2.,28.);
    sigmaZmm= 1907.0;
    addToVector(dSigmaMM,2,1.,14.);
    addToVector(dSigmaEM,2,0.,0.);
    addToVector(labelV,"stat expSyst");
    if (theCase==-13) {
      dSigmaEE.clear(); addToVector(dSigmaEE,1,2.);
      dSigmaMM.clear(); addToVector(dSigmaMM,1,1.);
      dSigmaEM.clear(); addToVector(dSigmaEM,1,0.);
      labelV.clear(); addToVector(labelV,"stat");
    }
    if ((theCase==130) || (theCase==131313)) {
      addToVector(dSigmaEE,1,32.);
      addToVector(dSigmaMM,1,29.);
      addToVector(dSigmaEM,1,sqrt(32.*29.));
      addToVector(labelV,"th.syst");
    }
    if ((theCase==1313) || (theCase==131313)) {
      addToVector(dSigmaEE,1,round(sigmaZee*dLumiRel));
      addToVector(dSigmaMM,1,round(sigmaZmm*dLumiRel));
      addToVector(dSigmaEM,1,sqrt(dSigmaEE.back()*dSigmaMM.back()));
      addToVector(labelV,"lumi");
    }
  }
  else if ((subVer==2) && valueEquals(theCase,"13 1313 131313 130 -13")) {
    // DY 13TeV
    // ee: 1949 +- 2 (stat) +- 28 (exp.syst) +- 32 (theo) +- 45 (lumi) pb
    // mm: 1907 +- 1 (Stat) +- 25 (exp.syst) +- 30 (theo) +- 44 (lumi) pb
    dLumiRel=0.023;
    sigmaZee= 1949.0;
    addToVector(dSigmaEE,2, 2.,28.);
    sigmaZmm= 1907.0;
    addToVector(dSigmaMM,2, 1.,25.);
    addToVector(dSigmaEM,2, 0.,0.);
    addToVector(labelV,"stat expSyst");
    if (theCase==-13) {
      dSigmaEE.clear(); addToVector(dSigmaEE,1, 2.);
      dSigmaMM.clear(); addToVector(dSigmaMM,1, 1.);
      dSigmaEM.clear(); addToVector(dSigmaEM,1, 0.);
      labelV.clear(); addToVector(labelV,"stat");
    }
    if ((theCase==130) || (theCase==131313)) {
      addToVector(dSigmaEE,1, 32.);
      addToVector(dSigmaMM,1, 30.);
      addToVector(dSigmaEM,1, sqrt(32.*30.));
      addToVector(labelV,"th.syst");
    }
    if ((theCase==1313) || (theCase==131313)) {
      addToVector(dSigmaEE,1,round(sigmaZee*dLumiRel));
      addToVector(dSigmaMM,1,round(sigmaZmm*dLumiRel));
      addToVector(dSigmaEM,1,sqrt(dSigmaEE.back()*dSigmaMM.back()));
      addToVector(labelV,"lumi");
    }
  }
  else {
    std::cout << "not ready for case " << theCase << "\n";
    return;
  }

  std::cout << "sigmaZee=" << sigmaZee; printVec(", dSigmaEE : ",dSigmaEE);
  std::cout << "sigmaZmm=" << sigmaZmm; printVec(", dSigmaMM : ",dSigmaMM);
  printVec("dSigmaEM : ",dSigmaEM);

  double totCovEE= totUnc(dSigmaEE,1);
  double totCovMM= totUnc(dSigmaMM,1);
  double totCovEM= totUnc(dSigmaEM,1);

  TMatrixD measEE(1,1), measMM(1,1);
  TMatrixD covEE(1,1), covMM(1,1), covEM(1,1);
  measEE(0,0) = sigmaZee;
  measMM(0,0) = sigmaZmm;
  covEE(0,0) = totCovEE;
  covMM(0,0) = totCovMM;
  covEM(0,0) = totCovEM;

  BLUEResult_t blue;
  if (!blue.estimate(measEE,covEE,measMM,covMM,covEM)) {
    std::cout << "estimate failed\n";
    return;
  }

  TMatrixD measLL( *blue.getEst() );
  TMatrixD covLL( *blue.getCov() );

  std::cout << "combined value: " << measLL(0,0) << " +- " << sqrt(covLL(0,0)) << "\n";
  printVec("included systematics : ",labelV);

  if ((theCase==1) || (theCase==11)) {
    std::vector<double> dSigmaLL_chk;
    addToVector(dSigmaLL_chk,2,10.,40.);
    if (theCase==11) addToVector(dSigmaLL_chk,1,90.);
    std::cout << " check result: " << 1910 << " +- " << totUnc(dSigmaLL_chk,0)
	      << "\n";
  }
  if (theCase==8) {
    std::vector<double> dSigmaLL_chk;
    addToVector(dSigmaLL_chk,2,8.,25.);
    std::cout << " check result : " << 1138 << " +- " << totUnc(dSigmaLL_chk,0)
	      << "\n";
  }

  std::stringstream ssNiceEE, ssNiceMM, ssNiceLL;
  ssNiceEE << "sigmaEE = " << measEE(0,0);
  ssNiceMM << "sigmaMM = " << measMM(0,0);
  ssNiceLL << "sigmaLL = " << measLL(0,0);
  for (unsigned int i=0; i<dSigmaEE.size(); i++) {
    TMatrixD covTotSrc(2,2);
    covTotSrc(0,0) = pow(dSigmaEE[i],2);
    covTotSrc(1,1) = pow(dSigmaMM[i],2);
    covTotSrc(0,1) = pow(dSigmaEM[i],2);
    covTotSrc(1,0) = pow(dSigmaEM[i],2);

    TMatrixD covFinSrc= blue.contributedCov(covTotSrc);
    std::cout << "source " << labelV[i] << " :\n";
    //covFinSrc.Print();
    std::cout << " sqrt( @(0,0) ) = " << sqrt(covFinSrc(0,0)) << "\n";
    ssNiceEE << " +- " << dSigmaEE[i] << " (" << labelV[i] << ")";
    ssNiceMM << " +- " << dSigmaMM[i] << " (" << labelV[i] << ")";
    ssNiceLL << " +- " << sqrt(covFinSrc(0,0)) << " (" << labelV[i] << ")";
  }
  std::cout << ssNiceEE.str() << "\n";
  std::cout << ssNiceMM.str() << "\n";
  std::cout << ssNiceLL.str() << "\n";


}
