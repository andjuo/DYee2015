#include "inputs.h"
#include "crossSection.h"
#include "Blue.h"
#include "analyseBLUEResult.C"
#include "work13TeV.C"


// -------------------------------------------------------

int addLumiCov(TMatrixD &eeCov, TMatrixD &mmCov, TMatrixD &emCov,
	       double lumiRelUnc, const TH1D *h1cs=NULL);

int loadUncData(int isEE,
		const finfomap_t &fNames,
		const finfomap_t &h1names,
		const TH1D *h1binning,
		const TH1D *h1centralVal, // if the uncertainty is relative
		std::vector<TH1D*> &h1uncV,
		TString tag);

int adjustMMUnc(const finfomap_t &mmCovFNames, // needed for first value
		std::vector<TMatrixD> &mmCovV);

int adjustUnc(const finfomap_t &covFNames,
	      std::vector<TMatrixD> &covV,
	      const finfomap_t &uncFiles,
	      const finfomap_t &uncStat,
	      const finfomap_t &uncSyst,
	      const TH1D *h1binning,
	      const TH1D *h1central,
	      int plotCmp, // 1 - plot uncertainties, 2 - and covariances
	      int change_covV, // 1 - change values, 2 - verify change by plot
	      TString tag);

int compareUnc(const std::vector<TString> &labelV,
	       const std::vector<TH1D*> &h1V_target,
	       std::vector<TMatrixD> &covV,
	       int plotCmp,
	       int change_covV,
	       TString tag);



// -------------------------------------------------------

void work13TeV_corr(int printCanvases=0, int corrCase=1, int includeLumiUnc=0,
		    int changeNames=1)
{
  std::vector<int> eeCovIdx, mmCovIdx;
  std::vector<TString> eeOldFName, eeNewFName; // for file name changing
  std::vector<TString> mmOldFName, mmNewFName; // for file name changing
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

  if (changeNames) {
    eeOldFName.push_back("cov_ee_varYield_2000.root");
    eeNewFName.push_back("cov_ee_varYield_5000.root");
    eeOldFName.push_back("cov_ee_varRhoFile_2000.root");
    eeNewFName.push_back("cov_ee_varRhoFile_5000.root");
    mmOldFName.push_back("cov_mumu_varYield_2000.root");
    mmNewFName.push_back("dymm_cov_varYield5000_20170202-covOnly.root");
    mmOldFName.push_back("cov_mumu_varRhoFile_2000.root");
    mmNewFName.push_back("cov_mumu_varRhoFile_5000.root");
  }

  TString fileTag="_corr";
  TString plotTag=" (corr)";

  finfomap_t eeCovFNames,mmCovFNames; // covariance derived by randomization
  finfomap_t eeUncTarget_stat, eeUncTarget_syst; // uncertainties from other
  finfomap_t mmUncTarget_stat, mmUncTarget_syst; // ... groups

  for (unsigned int i=0; i<eeCovIdx.size(); i++) {
    TVaried_t var= TVaried_t(eeCovIdx[i]);
    TString fname = "cov_ee_" + variedVarName(var) + TString("_2000.root");
    if (eeOldFName.size()>0) {
      for (unsigned int i=0; i<eeOldFName.size(); i++) {
	if (fname==eeOldFName[i]) fname=eeNewFName[i];
      }
    }
    eeCovFNames[var] = fname;
  }
  for (unsigned int i=0; i<mmCovIdx.size(); i++) {
    TVaried_t var= TVaried_t(mmCovIdx[i]);
    TString fname = "cov_mumu_" + variedVarName(var) + TString("_2000.root");
    if (mmOldFName.size()>0) {
      for (unsigned int i=0; i<mmOldFName.size(); i++) {
	if (fname==mmOldFName[i]) fname=mmNewFName[i];
      }
    }
    mmCovFNames[var] = fname;
  }

  std::cout << "ee covariance names:\n";
  for (finfomap_t::const_iterator it= eeCovFNames.begin();
       it!=eeCovFNames.end(); it++) {
    TString fname= it->second;
    std::cout << "  -- " << variedVarName(it->first) << "\t\t"
	      << it->second << "\n";
  }
  std::cout << "mm covariance names:\n";
  for (finfomap_t::const_iterator it= mmCovFNames.begin();
       it!=mmCovFNames.end(); it++) {
    TString fname= it->second;
    std::cout << "  -- " << variedVarName(it->first) << "\t\t"
	      << it->second << "\n";
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

  adjustMMUnc(mmCovFNames,mmCovV);
  return;


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

  TH1D* h1csEE=loadHisto(eeCSFName, eeCSH1Name, "h1csEE", 1,h1dummy);
  TH1D* h1csMM=loadHisto(mmCSFName, mmCSH1Name, "h1csMM", 1,h1dummy);

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

int adjustMMUnc(const finfomap_t &mmCovFNames,
		std::vector<TMatrixD> &mmCovV)
{
  finfomap_t mmUncFiles, mmUncStat, mmUncSyst;

  TH1D* h1csEE=loadHisto(eeCSFName, eeCSH1Name, "h1csEE", 1,h1dummy);
  TH1D* h1csMM=loadHisto(mmCSFName, mmCSH1Name, "h1csMM", 1,h1dummy);
  if (!h1csEE || !h1csMM) return 0;

  TString inpPath="DYmm-uncertainties-2016Dec07/";
  mmUncFiles[_varYield]= "dymm_unc_fromNote.root";
  mmUncStat [_varYield]= "h1csMM_statUnc";
  mmUncSyst [_varYield]= "";
  mmUncFiles[_varBkg]= inpPath + "ROOTFile_RelUnc_Stat_Syst_Tot_BkgEst.root";
  mmUncStat [_varBkg]= "h_RelUnc_Stat";
  mmUncSyst [_varBkg]= "h_RelUnc_Syst";
  mmUncFiles[_varDetRes]= inpPath + "ROOTFile_RelUnc_Stat_Syst_Tot_DetRes.root";
  mmUncStat [_varDetRes]= "h_RelUnc_Stat";
  mmUncSyst [_varDetRes]= "h_RelUnc_Syst";
  mmUncFiles[_varEff]= inpPath + "ROOTFile_RelUnc_Stat_Syst_Tot_EffSF.root";
  mmUncStat [_varEff]= "h_RelUnc_Stat";
  mmUncSyst [_varEff]= "h_RelUnc_Syst";
  mmUncFiles[_varRhoFile]= "";
  mmUncStat [_varRhoFile]= "";
  mmUncSyst [_varRhoFile]= "";
  //mmUncFiles[_varRhoFile]= inpPath + "ROOTFile_RelUnc_Stat_Syst_Tot_EffSF.root";
  //mmUncStat [_varRhoFile]= "h_RelUnc_Stat";
  //mmUncSyst [_varRhoFile]= "h_RelUnc_Syst";
  mmUncFiles[_varFSRRes]= inpPath + "ROOTFile_RelUnc_Stat_Syst_Tot_FSR.root";
  mmUncStat [_varFSRRes]= "h_RelUnc_Stat";
  mmUncSyst [_varFSRRes]= "h_RelUnc_Syst";
  mmUncFiles[_varAcc]= "dymm_unc_fromNote.root";
  mmUncStat [_varAcc]= "";
  mmUncSyst [_varAcc]= "h1csMM_theoUnc";

  const TH1D *h1central= (0) ? h1csMM : NULL;
  int plotCmp=2;
  int change_covV=1;
  int res=adjustUnc( mmCovFNames, mmCovV,
		     mmUncFiles,mmUncStat,mmUncSyst,
		     h1csMM,h1central,
		     plotCmp,change_covV,
		     "_mm");
  if (!res) {
    std::cout << "error in adjustMMUnc\n";
    return 0;
  }
  return 1;
}
// -----------------------------------------------------------------------

int loadUncData(int isEE,
		const finfomap_t &fNames,
		const finfomap_t &h1names,
		const TH1D *h1binning,
		const TH1D *h1centralVal, // if the uncertainty is relative
		std::vector<TH1D*> &h1uncV,
		TString tag)
{
  TString lepTag=(isEE) ? "ee" : "mm";

  std::cout << "\n";
  if (h1centralVal) printHisto(h1centralVal);


  for (finfomap_t::const_iterator it= fNames.begin(); it!=fNames.end(); it++) {
    TString fname= it->second;
    TString h1newName= "h1unc_" + lepTag + "_" + variedVarName(it->first) + tag;
    TH1D *h1unc= cloneHisto(h1binning, h1newName,h1newName);
    h1unc->Reset();
    h1uncV.push_back(h1unc);
    if (fname.Length()==0) continue;

    finfomap_t::const_iterator itHName= h1names.find(it->first);
    if (itHName==h1names.end()) {
      std::cout << "could not find the histo name for "
		<< variedVarName(it->first) << "\n";
      return 0;
    }
    if (itHName->second.Length()==0) continue;

    TH1D *h1= loadHisto(fname,itHName->second,itHName->second + "_tmp",h1dummy);
    if (!h1) return 0;
    std::cout << "got " << h1newName << "\n";
    printHisto(h1);
    std::cout << "\n\n" << std::endl;
    if (!copyContents(h1unc, h1)) return 0;
    if ((h1centralVal!=NULL) && (fname.Index("fromNote")==-1)) {
      h1unc->Multiply(h1centralVal);
    }
  }
  return 1;
}

// -------------------------------------------------------

int adjustUnc(const finfomap_t &covFNames,
	      std::vector<TMatrixD> &covV,
	      const finfomap_t &uncFiles,
	      const finfomap_t &uncStat,
	      const finfomap_t &uncSyst,
	      const TH1D *h1binning,
	      const TH1D *h1central,
	      int plotCmp, // 1 - plot uncertainties, 2 - and covariances
	      int change_covV, // 1 - change values, 2 - verify change by plot
	      TString tag)
{
  if (uncFiles.size()!=covV.size()) {
    std::cout << "sizes are different!\n";
  }

  std::vector<TString> labelV;
  for (finfomap_t::const_iterator it= covFNames.begin();
       it!=covFNames.end(); it++) {
    TVaried_t var= it->first;
    if (uncFiles.find(var)==uncFiles.end()) {
      std::cout << "mmUncFiles does not contain " << variedVarName(var) << "\n";
      return 0;
    }
    labelV.push_back(variedVarName(var) + tag);
  }

  std::vector<TH1D*> h1uncStat, h1uncSyst;
  if (!loadUncData(0,uncFiles,uncStat,h1binning,h1central,h1uncStat,tag+"_stat")) {
    std::cout << "failed to load statistical component\n";
    return 0;
  }
  if (!loadUncData(0,uncFiles,uncStat,h1binning,h1central,h1uncSyst,tag+"_syst")) {
    std::cout << "failed to load systematic component\n";
    return 0;
  }

  int res=compareUnc(labelV,h1uncStat,covV,plotCmp,change_covV,tag);
  if (res && plotCmp && (change_covV==2))
    res=compareUnc(labelV,h1uncStat,covV,2,0,tag+"_adj");

  if (!res) {
    std::cout << "error in adjustUnc\n";
    return 0;
  }
  return 1;
}

// -------------------------------------------------------

int compareUnc(const std::vector<TString> &labelV,
	       const std::vector<TH1D*> &h1V_target,
	       std::vector<TMatrixD> &covV,
	       int plotCmp,
	       int change_covV,
	       TString tag)
{
  for (unsigned int i=0; i<labelV.size(); i++) {
    TString h2Name=h1V_target[i]->GetName() + TString("_cov");
    TH2D *h2cov= convert2histo(covV[i],h1V_target[i],h2Name,h2Name,NULL);
    TH1D *h1uncFromCov= uncFromCov(h2cov,NULL,0);
    if (plotCmp) {
      TString canvName="cCmpUnc_" + labelV[i] + "_" + tag;
      TH1D *h1= cloneHisto(h1V_target[i],
			   h1V_target[i]->GetName() + TString("_plotClone")+tag,
			   h1V_target[i]->GetTitle());
      hsRed.SetStyle(h1);
      plotHisto(h1,canvName,1,1,"LP",labelV[i] + " target");
      plotHistoSame(h1uncFromCov,canvName,"LP",labelV[i] + " cov");
    }
    if (plotCmp==2) {
      TString canvName2= "cCov_" + labelV[i] + "_" + tag;
      PlotCovCorrOpt_t opt;
      plotCovCorr(covV[i],h1V_target[1],"h2cov_"+labelV[i]+"_"+tag,
		  canvName2,opt);
    }

    if (change_covV) {
      if (h1V_target[i]->Integral()!=double(0)) {
	std::cout << "adjusting covV of " << labelV[i] << "\n";
	if (!changeCov(covV[i],h1V_target[i])) return 0;
      }
      else {
	std::cout << "no adjustment of covV for " << labelV[i] << "\n";
      }
    }
  }
  return 1;
}

// -------------------------------------------------------
