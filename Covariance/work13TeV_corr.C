#include "inputs.h"
#include "crossSection.h"
#include "Blue.h"
#include "analyseBLUEResult.C"
#include "work13TeV.C"


// -------------------------------------------------------

int addLumiCov(TMatrixD &eeCov, TMatrixD &mmCov, TMatrixD &emCov,
	       double lumiRelUnc, const TH1D *h1cs=NULL);

// -------------------------------------------------------

void work13TeV_corr(int printCanvases=0, int corrCase=1, int includeLumiUnc=0)
{
  std::vector<int> eeCovIdx, mmCovIdx;
  if (1) {
    addToVector(eeCovIdx,7, _varYield, _varBkg, _varDetRes, _varEff,
		_varRhoFile, _varAcc, _varFSRRes);
    addToVector(mmCovIdx,7, _varYield, _varBkg, _varDetRes, _varEff,
		_varRhoFile, _varAcc, _varFSRRes);
  }
  else {
    addToVector(eeCovIdx,6,  _varBkg, _varDetRes, _varEff,
		_varRhoFile, _varAcc, _varFSRRes);
    addToVector(mmCovIdx,6,  _varBkg, _varDetRes, _varEff,
		_varRhoFile, _varAcc, _varFSRRes);
  }
   TString fileTag="_corr";
  TString plotTag=" (corr)";

  finfomap_t eeCovFNames,mmCovFNames;
  for (unsigned int i=0; i<eeCovIdx.size(); i++) {
    TVaried_t var= TVaried_t(eeCovIdx[i]);
    eeCovFNames[var] = "cov_ee_" + variedVarName(var) + TString("_1000.root");
  }
  for (unsigned int i=0; i<mmCovIdx.size(); i++) {
    TVaried_t var= TVaried_t(mmCovIdx[i]);
    mmCovFNames[var] = "cov_mumu_" + variedVarName(var) + TString("_1000.root");
  }

  std::vector<TH2D*> eeCovH2D, mmCovH2D;
  std::vector<TMatrixD> eeCovV, mmCovV;
  TMatrixD eeCovStat(43,43);
  TMatrixD eeCovSyst(eeCovStat);
  TMatrixD mmCovStat(eeCovStat);
  TMatrixD mmCovSyst(eeCovStat);

  if (!loadCovData(eeCovFNames,eeCovH2D,eeCovV,eeCovStat,eeCovSyst) ||
      !loadCovData(mmCovFNames,mmCovH2D,mmCovV,mmCovStat,mmCovSyst)) {
    std::cout << "failed to load covariances\n";
    return;
  }

  TMatrixD eeCovTot(eeCovStat,TMatrixD::kPlus,eeCovSyst);
  TMatrixD mmCovTot(mmCovStat,TMatrixD::kPlus,mmCovSyst);
  TMatrixD emCovTot(TMatrixD::kZero, eeCovTot);

  double lumiUnc=0.027;
  TString fileTagExtra,plotTagExtra;
  if (includeLumiUnc) {
    fileTagExtra="_wLumi"; plotTagExtra=" (wLumi)";
    if (includeLumiUnc==2) {
      eeCovTot.Zero(); mmCovTot.Zero(); emCovTot.Zero();
      lumiUnc=0.27;
      fileTagExtra="_wLumi10x"; plotTagExtra=" (wLumi10x)";
    }
    if (!addLumiCov(eeCovTot,mmCovTot,emCovTot, lumiUnc)) return;
  }

  if (corrCase==0) {
    eeCovTot = reduceCorrelations(eeCovTot, 1.);
    mmCovTot = reduceCorrelations(mmCovTot, 1.);
    if (includeLumiUnc!=2) emCovTot.Zero();
    fileTag="_uncorr";
    plotTag=" (uncorr)";
  }
  else {
    fileTag="_corr";
    plotTag=" (corr)";
  }
  if (fileTagExtra.Length()) fileTag+=fileTagExtra;
  if (plotTagExtra.Length()) plotTag+=plotTagExtra;

  BLUEResult_t *blue=
    combineData(&eeCovTot,&mmCovTot,&emCovTot,fileTag,plotTag,printCanvases);

  if (includeLumiUnc) {
    TMatrixD totFinalCov(*blue->getCov());
    TMatrixD ca(0,0), cb(0,0);
    TH1D *h1csLL= convert2histo1D(*blue->getEst(),totFinalCov,
				  "h1csLL_test","h1csLL test",NULL);
    if (!addLumiCov(ca,cb,totFinalCov,-lumiUnc,h1csLL)) return;
    plotCovCorr(totFinalCov,NULL,"h2finCovNoLumi","cFinCovNoLumi");
  }
}

// -------------------------------------------------------

int addLumiCov(TMatrixD &eeCov, TMatrixD &mmCov, TMatrixD &emCov,
	       double lumiRelUnc, const TH1D *h1cs)
{
  if (lumiRelUnc==0.) {
    std::cout << "addLumiCov: lumiUnc=0., no effect\n";
    return 1;
  }

  TString eeCSFName="cs_DYee_13TeV_El3.root";
  TString mmCSFName="cs_DYmm_13TeVMuApproved_cs.root";

  TH1D* h1csEE=loadHisto(eeCSFName, "h1PreFSRCS", "h1csEE", 1,h1dummy);
  TH1D* h1csMM=loadHisto(mmCSFName, "h1CS", "h1csMM", 1,h1dummy);

  if (!h1csEE || !h1csMM) return 0;

  double sign= (lumiRelUnc>0) ? 1 : -1;

  for (int ir=0; ir<eeCov.GetNrows(); ir++) {
    for (int ic=0; ic<eeCov.GetNcols(); ic++) {
      eeCov(ir,ic) +=
	h1csEE->GetBinContent(ir+1) * h1csEE->GetBinContent(ic+1)
	* pow( lumiRelUnc, 2) * sign;
    }
  }
  for (int ir=0; ir<mmCov.GetNrows(); ir++) {
    for (int ic=0; ic<mmCov.GetNcols(); ic++) {
      mmCov(ir,ic) +=
	h1csMM->GetBinContent(ir+1) * h1csMM->GetBinContent(ic+1)
	* pow( lumiRelUnc, 2) * sign;
    }
  }

  if (h1cs==NULL) {
    for (int ir=0; ir<eeCov.GetNrows(); ir++) {
      for (int ic=0; ic<eeCov.GetNcols(); ic++) {
	emCov(ir,ic) +=
	  h1csEE->GetBinContent(ir+1) * h1csMM->GetBinContent(ic+1)
	  * pow( lumiRelUnc, 2) * sign;
      }
    }
  }
  else {
    for (int ir=0; ir<emCov.GetNrows(); ir++) {
      for (int ic=0; ic<emCov.GetNcols(); ic++) {
	emCov(ir,ic) +=
	  h1cs->GetBinContent(ir+1) * h1cs->GetBinContent(ic+1)
	  * pow( lumiRelUnc, 2) * sign;
      }
    }
  }

  if (1) {
    double diagMultEE=1., diagMultMM=1.;
    for (int ir=0; ir<eeCov.GetNrows(); ir++) {
      diagMultEE *= eeCov(ir,ir);
      diagMultMM *= mmCov(ir,ir);
    }
    std::cout << "diagMultEE=" << diagMultEE << ", diagMultMM=" << diagMultMM << "\n";
    //return 0;
  }

  return 1;
}

// -------------------------------------------------------
