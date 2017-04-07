#include "inputs.h"
#include "crossSection.h"
#include "Blue.h"
#include "analyseBLUEResult.h"
#include "CovStruct.h"

// -------------------------------------------------------
// -------------------------------------------------------

int adjustMMUnc(const finfomap_t &mmCovFNames, // needed for first value
		std::vector<TMatrixD> &mmCovV,
		CovStruct_t &covM, int showChangedCov,
		int plotCmp);

int adjustEEUnc(const finfomap_t &eeCovFNames, // needed for first value
		std::vector<TMatrixD> &eeCovV,
		CovStruct_t &eeCovM, int showChangedCov,
		int plotCmp);

// -------------------------------------------------------

const int useUncorrYieldUnc=0;
TString eeCSFName="cs_DYee_13TeV_El3.root";
TString eeCSH1Name="h1PreFSRCS";
TString mmCSFName="cs_DYmm_13TeVMuApproved_cs.root";
TString mmCSH1Name="h1CS";

// -------------------------------------------------------

void work13TeV_corr(int printCanvases=0, int corrCase=1, int includeLumiUnc=0,
		    int plotChangedCov=0,
		    int excludeSyst=0, // 1 - all, 2 - Acc only
		    TString saveTag="")
{
  closeCanvases(5);
  const int changeNames=1;

  std::string showCanvs;
  //showCanvs+=" plotChannelInputCovOnFileAdj";
 // whether relative uncertainty in Acc in the channels is similar
  //showCanvs+=" plotChannelRelSystAccUnc";
  // compare randomized vs provided uncertainties
  //showCanvs+=" plotAdjUnc";
  //showCanvs+=" plotAdjCov";
  //showCanvs+=" plotInputContributedCovs"; // 4 canvases in each channel
  showCanvs+= " plotFinalCovByType"; // 4 canvases in the combined channel

  std::string showCombinationCanvases;
  if (0) {
    showCombinationCanvases="ALL";
  }
  else {
    showCombinationCanvases+= " showCombiCSCheck";
    showCombinationCanvases+= " showCombiCSUnc";
  }

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

  // xxCovStat is the covariance of measured yield
  // xxCovSyst equals to other covariances
  // xxCovV contains all covariance matrices
  if (!loadCovData("ee",eeCovFNames,eeCovH2D,eeCovV,eeCovStat,eeCovSyst) ||
      !loadCovData("mm",mmCovFNames,mmCovH2D,mmCovV,mmCovStat,mmCovSyst)) {
    std::cout << "failed to load covariances\n";
    return;
  }

  // Adjust components of the uncertainties
  CovStruct_t eeCovS(eeCovStat); // creates and nullifies
  CovStruct_t mmCovS(mmCovStat);

  int plotCmpUnc= (hasValue("plotAdjUnc",showCanvs)) ? 1 : 0;
  if (hasValue("plotAdjCov",showCanvs)) plotCmpUnc=2;

  adjustMMUnc(mmCovFNames,mmCovV,mmCovS,plotChangedCov,plotCmpUnc);
  adjustEEUnc(eeCovFNames,eeCovV,eeCovS,plotChangedCov,plotCmpUnc);

  if (excludeSyst) {
    if (excludeSyst==1) {
      mmCovS.editSystNonAcc().Zero();
      eeCovS.editSystNonAcc().Zero();
    }
    mmCovS.editSystAcc().Zero();
    eeCovS.editSystAcc().Zero();
  }

  TH1D* h1csEE_tmp=loadHisto(eeCSFName, eeCSH1Name, "h1csEE_tmp", 1,h1dummy);
  TH1D* h1csMM_tmp=loadHisto(mmCSFName, mmCSH1Name, "h1csMM_tmp", 1,h1dummy);
  if (!h1csEE_tmp || !h1csMM_tmp) return;

  // plot input covariance in ee and mm channels
  if (hasValue("plotChannelInputCovOnFileAdj",showCanvs)) {
    double ymin=0, ymax=0;
    ymin=1e-9; ymax=20.;
    mmCovS.Plot("mmCov_onFile_adj",h1csMM_tmp);
    mmCovS.PlotUnc("mmCovUnc_onFile_adj",h1csMM_tmp,NULL,ymin,ymax);
    eeCovS.Plot("eeCov_onFile_adj",h1csEE_tmp);
    eeCovS.PlotUnc("eeCovUnc_onFile_adj",h1csEE_tmp,NULL,ymin,ymax);
    //return;
  }

  // check whether mmCovS.Syst_Acc = eeCovS.Syst_Acc
  if (hasValue("plotChannelRelSystAccUnc",showCanvs)) {
    TH1D *h1mmAcc_rel= uncFromCov(mmCovS.covSyst_Acc,"h1mmAcc_rel",
				  h1csMM_tmp,h1csMM_tmp,0);
    TH1D *h1eeAcc_rel= uncFromCov(eeCovS.covSyst_Acc,"h1eeAcc_rel",
				  h1csEE_tmp,h1csEE_tmp,0);
    hsRed.SetStyle(h1eeAcc_rel);
    plotHisto(h1mmAcc_rel,"cAccUnc",1,1,"LP","#mu#mu");
    plotHistoSame(h1eeAcc_rel,"cAccUnc","LP","ee");
    //std::cout << "stopping\n"; return;
  }


  // create Acc theory-correlated uncertainty
  TMatrixD emCovTot(TMatrixD::kZero, eeCovS.covSyst_Acc);
  if (!excludeSyst) {
    double accCorrFactor=1.;
    emCovTot= constructEMAccCov(eeCovS,mmCovS,accCorrFactor);
  }

  if (corrCase==0) {
    double reductionFactor=1.;
    eeCovS.ReduceCorrelations(reductionFactor);
    mmCovS.ReduceCorrelations(reductionFactor);
    if (includeLumiUnc!=2) emCovTot.Zero();
  }

  if (0) {
    if (!eeCovS.PrintZRangeUnc("ee",h1csEE_tmp,0) ||
	!mmCovS.PrintZRangeUnc("mm",h1csMM_tmp,0)) {
      std::cout << "error\n";
    }
    return;
  }

  TMatrixD eeCovTot(eeCovS.Sum());
  TMatrixD mmCovTot(mmCovS.Sum());

  double lumiUnc=0.026;
  TString fileTagExtra,plotTagExtra;
  if (includeLumiUnc) {
    fileTagExtra="_wLumi"; plotTagExtra=" (wLumi)";
    if (includeLumiUnc==2) {
      eeCovTot.Zero(); mmCovTot.Zero(); emCovTot.Zero();
      lumiUnc=0.26;
      fileTagExtra="_wLumi10x"; plotTagExtra=" (wLumi10x)";
    }
    if (!addLumiCov(eeCovTot,mmCovTot,emCovTot, lumiUnc, h1csEE_tmp,h1csMM_tmp))
      return;
  }

  if (corrCase==0) {
    fileTag="_uncorr";
    plotTag=" (uncorr)";
  }
  else {
    fileTag="_corr";
    plotTag=" (corr)";
  }
  if (fileTagExtra.Length()) fileTag+=fileTagExtra;
  if (plotTagExtra.Length()) plotTag+=plotTagExtra;
  if (plotChangedCov) fileTag.Append("_plotChCov");

  if (!changeNames) fileTag.Append("_noChangeNames");
  if (useUncorrYieldUnc) fileTag.Append("_uncorrYieldUnc");
  if (excludeSyst) fileTag.Append(Form("_excludeSyst%d",excludeSyst));
  if (saveTag.Length()) fileTag.Append(saveTag);


  std::cout << "showCombinationCanvases=" << showCombinationCanvases << "\n";
  BLUEResult_t *blue=
    combineData(&eeCovTot,&mmCovTot,&emCovTot,fileTag,plotTag,printCanvases,
		showCombinationCanvases);
  if (!blue) { std::cout << "failed to get BLUE result\n"; return; }

  if (includeLumiUnc) {
    TMatrixD totFinalCov(*blue->getCov());
    TH1D *h1csLL_test= convert2histo1D(*blue->getEst(),totFinalCov,
				       "h1csLL_test","h1csLL test",NULL);
    if (!addLumiCov(totFinalCov,-lumiUnc,h1csLL_test)) return;
    plotCovCorr(totFinalCov,NULL,"h2finCovNoLumi","cFinCovNoLumi");
  }

  // make analysing plots
  if (!detailedCovPlots(blue,corrCase,
			h1csEE_tmp,h1csMM_tmp,
			eeCovS,mmCovS,showCanvs)) return;

  if (printCanvases || (saveTag.Length()>0)) {
    TString outFName="dyll-combi-" + fileTag + ".root";
    TFile fout(outFName,"recreate");
    blue->write(fout,fileTag);
    SaveCanvases("ALL","dyll-combi-"+fileTag,&fout);
    //h2cov_fromYield->Write();
    //h1_dCS_fromYield->Write();
    writeTimeTag(&fout);
    fout.Close();
    std::cout << "file <" << fout.GetName() << "> created\n";
  }
}

// -------------------------------------------------------

int adjustMMUnc(const finfomap_t &mmCovFNames,
		std::vector<TMatrixD> &mmCovV,
		CovStruct_t &mmCovS, int plotChangedCov,
		int plotCmp)
{
  finfomap_t mmUncFiles, mmUncStat, mmUncSyst;
  fvarweightmap_t relUncW;

  TH1D* h1csEE_loc=loadHisto(eeCSFName, eeCSH1Name, "h1csEE_loc", 1,h1dummy);
  TH1D* h1csMM_loc=loadHisto(mmCSFName, mmCSH1Name, "h1csMM_loc", 1,h1dummy);
  if (!h1csEE_loc || !h1csMM_loc) return 0;

  TString inpPath="./";
  TString inpFile=inpPath + "DYmm_ROOTFile_Input-20170207.root";

  TH1D *h1csMM= cloneHisto(h1csMM_loc,"h1csMM","h1csMM");
  if (0) {
    h1csMM= loadHisto(inpFile, "h_DiffXSec", "h1csMM_KL", 1,h1dummy);
    plotHisto(h1csMM_loc,"cCSmm",1,1,"LPE","old #mu#mu cs");
    plotHistoSame(h1csMM,"cCSmm","LPE","new #mu#mu cs");
    printRatio(h1csMM,h1csMM_loc);
    return 0;
  }

  mmUncFiles[_varYield]= inpFile;
  mmUncStat [_varYield]= "h_RelStatUnc";
  mmUncSyst [_varYield]= "";
  relUncW   [_varYield]= (useUncorrYieldUnc) ? 1. : 0.;
  mmUncFiles[_varBkg]= inpFile;
  mmUncStat [_varBkg]= "h_RelUnc_Stat_BkgEst";
  mmUncSyst [_varBkg]= "h_RelUnc_Syst_BkgEst";
  relUncW   [_varBkg]= 1;
  mmUncFiles[_varDetRes]= inpFile;
  mmUncStat [_varDetRes]= "h_RelUnc_Stat_DetRes";
  mmUncSyst [_varDetRes]= "h_RelUnc_Syst_DetRes";
  relUncW   [_varDetRes]= 1;
  mmUncFiles[_varEff]= inpFile;
  mmUncStat [_varEff]= "";
  mmUncSyst [_varEff]= "";
  relUncW   [_varEff]= 0;
  mmUncFiles[_varRhoFile]= inpFile;
  mmUncStat [_varRhoFile]= "h_RelUnc_Stat_EffSF";
  mmUncSyst [_varRhoFile]= "h_RelUnc_Syst_EffSF";
  relUncW   [_varRhoFile]= 1;
  mmUncFiles[_varFSRRes]= inpFile;
  mmUncStat [_varFSRRes]= "h_RelUnc_Stat_FSR";
  mmUncSyst [_varFSRRes]= "h_RelUnc_Syst_FSR";
  relUncW   [_varFSRRes]= 1;
  mmUncFiles[_varAcc]= inpFile;
  mmUncStat [_varAcc]= "";
  mmUncSyst [_varAcc]= "h_RelUnc_Syst_Acc";
  relUncW   [_varAcc]= 1;

  const TH1D *h1central= (1) ? h1csMM : NULL;
  int change_covV=1+(plotChangedCov!=0) ? 1:0;
  int res=adjustUnc( "mm",
		     mmCovFNames, mmCovV,
		     mmUncFiles,mmUncStat,mmUncSyst,relUncW,
		     h1csMM,h1central,
		     plotCmp,change_covV,
		     "_mm",
		     mmCovS,
		     "");
  if (!res) {
    std::cout << "error in adjustMMUnc\n";
    return 0;
  }

  return 1;
}
// -----------------------------------------------------------------------

int adjustEEUnc(const finfomap_t &eeCovFNames,
		std::vector<TMatrixD> &eeCovV,
		CovStruct_t &eeCovS, int plotChangedCov,
		int plotCmp)
{
  finfomap_t eeUncFiles, eeUncStat, eeUncSyst;
  fvarweightmap_t relUncW;

  TH1D* h1csEE_loc=loadHisto(eeCSFName, eeCSH1Name, "h1csEE_loc", 1,h1dummy);
  TH1D* h1csMM_loc=loadHisto(mmCSFName, mmCSH1Name, "h1csMM_loc", 1,h1dummy);
  if (!h1csEE_loc || !h1csMM_loc) return 0;

  TString inpPath="./";
  TString inpFile=inpPath + "DYee_ROOTFile_Input-20170208.root";

  TH1D *h1csEE= cloneHisto(h1csEE_loc,"h1csEE","h1csEE");
  if (0) {
    h1csEE= loadHisto(inpFile, "h_DiffXSec", "h1csEE_RC", 1,h1dummy);
    plotHisto(h1csEE_loc,"cCSee",1,1,"LPE","old #it{ee} cs");
    plotHistoSame(h1csEE,"cCSee","LPE","new #it{ee} cs");
    printRatio(h1csEE,h1csEE_loc);
    return 0;
  }

  eeUncFiles[_varYield]= inpFile;
  eeUncStat [_varYield]= "h_RelStatUnc";
  eeUncSyst [_varYield]= "";
  relUncW   [_varYield]= (useUncorrYieldUnc) ? 0.01 : 0.;
  eeUncFiles[_varBkg]= inpFile;
  eeUncStat [_varBkg]= "h_RelUnc_Stat_BkgEst";
  eeUncSyst [_varBkg]= "h_RelUnc_Syst_BkgEst";
  relUncW   [_varBkg]= 0.01;
  eeUncFiles[_varDetRes]= inpFile;
  eeUncStat [_varDetRes]= "h_RelUnc_Stat_DetRes"; // this is actually systematic
  eeUncSyst [_varDetRes]= "h_RelUnc_Syst_DetRes";
  relUncW   [_varDetRes]= 0.01;
  eeUncFiles[_varEff]= inpFile;
  eeUncStat [_varEff]= "";
  eeUncSyst [_varEff]= "";
  relUncW   [_varEff]= 0;
  eeUncFiles[_varRhoFile]= inpFile;
  eeUncStat [_varRhoFile]= "h_RelUnc_Stat_EffSF";
  eeUncSyst [_varRhoFile]= "h_RelUnc_Syst_EffSF";
  relUncW   [_varRhoFile]= 0.01;
  eeUncFiles[_varFSRRes]= inpFile;
  eeUncStat [_varFSRRes]= "";
  eeUncSyst [_varFSRRes]= "h_RelUnc_Syst_FSR";
  relUncW   [_varFSRRes]= 0.01;
  eeUncFiles[_varAcc]= inpFile;
  eeUncStat [_varAcc]= "";
  eeUncSyst [_varAcc]= "h_RelUnc_Syst_Acc";
  relUncW   [_varAcc]= 0.01;

  const TH1D *h1central= (1) ? h1csEE : NULL;
  int change_covV=1 + (plotChangedCov!=0) ? 1:0;
  int res=adjustUnc( "ee",
		     eeCovFNames, eeCovV,
		     eeUncFiles,eeUncStat,eeUncSyst,relUncW,
		     h1csEE,h1central,
		     plotCmp,change_covV,
		     "_ee",
		     eeCovS,
		     "dyee13TeV_detRes_correction");

  if (!res) {
    std::cout << "error in adjustEEUnc\n";
    return 0;
  }
  return 1;
}
// -----------------------------------------------------------------------
