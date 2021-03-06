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
		int plotCmp, int noUncAdjustment=0,
		TVersion_t inpVer=_verUndef,
		TSubVersion_t inpSubVer=_subverUndef);

int adjustEEUnc(const finfomap_t &eeCovFNames, // needed for first value
		std::vector<TMatrixD> &eeCovV,
		CovStruct_t &eeCovM, int showChangedCov,
		int plotCmp, int noUncAdjustment=0,
		TVersion_t inpVer=_verUndef,
		TSubVersion_t inpSubVer=_subverUndef);

// -------------------------------------------------------

const int useInputYieldUnc=1;
//TString eeCSFName="cs_DYee_13TeV_El3.root";
//TString eeCSH1Name="h1PreFSRCS";
//TString mmCSFName="cs_DYmm_13TeVMuApproved_cs.root";
//TString mmCSH1Name="h1CS";
const int modifyAccUnc=0;

// -------------------------------------------------------

void work13TeV_corr(int printCanvases=0, int corrCase=1, int includeLumiUnc=0,
		    int plotChangedCov=0,
		    int excludeSyst=0, // 1 - all, 2 - Acc only
		    TString saveTag="",
		    int noUncAdjustment=0,
		    int includeRhoSyst=0)
{
  closeCanvases(5);
  const int changeNames=1;

  TVersion_t inpVer=_verMuApproved;
  TSubVersion_t inpSubVer=_subverUndef;
  inpVer=_verMuMay2017; inpSubVer=_subverUndef;
  inpVer=_verMuMay2017; inpSubVer=_subver1;
  inpVer=_verMuMay2017; inpSubVer=_subver2;
  inpVer=_verMuMay2017; inpSubVer=_subver3;
  //inpVer=_verElMay2017false;
  inpVer=_verMuNov2017; inpSubVer=_subver1;

  eeCSHisto1DName="h1PreFSRCS";
  mmCSHisto1DName="h1CS";

  if (inpVer==_verMuApproved) {
    eeCSFName="cs_DYee_13TeV_El3.root";
    mmCSFName="cs_DYmm_13TeVMuApproved_cs.root";
  }
  else if (inpVer==_verMuMay2017) {
    eeCSFName="cs_DYee_13TeV_ElMay2017.root";
    mmCSFName="cs_DYmm_13TeVMuMay2017_cs.root";
  }
  else if (inpVer==_verElMay2017false) {
    eeCSFName="cs_DYee_13TeV_ElMay2017false.root";
    mmCSFName="cs_DYmm_13TeVMuMay2017_cs.root";
  }
  else if (inpVer==_verMuNov2017) {
    eeCSFName="cs_DYee_13TeV_ElNov2017false.root";
    mmCSFName="cs_DYmm_13TeVMuNov2017_cs.root";
  }

  std::string showCanvs;
  showCanvs+=" plotChannelInputCovOnFileAdj";
 // whether relative uncertainty in Acc in the channels is similar
  showCanvs+=" plotChannelRelSystAccUnc";
  // compare randomized vs provided uncertainties
  showCanvs+=" plotAdjUnc";
  //showCanvs+=" plotAdjCov";
  //showCanvs+=" plotInputContributedCovs"; // 4 canvases in each channel
  showCanvs+=" plotFinalCov";
  showCanvs+= " plotFinalCovByType"; // 4 canvases in the combined channel

  std::cout << "eeCSFName=" << eeCSFName << ", eeCSHisto1DName=" << eeCSHisto1DName << "\n";
  std::cout << "mmCSFName=" << mmCSFName << ", mmCSHisto1DName=" << mmCSHisto1DName << "\n";


  std::string showCombinationCanvases;
  if (0) {
    showCombinationCanvases="ALL";
  }
  else {
    showCombinationCanvases+= " showCombiCS";
    showCombinationCanvases+= " showCombiCSCheck";
    showCombinationCanvases+= " showCombiCSUnc";
  }

  std::vector<int> eeCovIdx, mmCovIdx;
  std::vector<TString> eeOldFName, eeNewFName; // for file name changing
  std::vector<TString> mmOldFName, mmNewFName; // for file name changing
  if (1) {
    if (!modifyAccUnc) { // all factors
      addToVector(eeCovIdx,7, _varYield, _varBkg, _varDetRes, _varEff,
		  _varRhoFile, _varAcc, _varFSRRes);
      addToVector(mmCovIdx,7, _varYield, _varBkg, _varDetRes, _varEff,
		  _varRhoFile, _varAcc, _varFSRRes);

    }
    else { // all factors and _varTheory
      addToVector(eeCovIdx,8, _varYield, _varBkg, _varDetRes, _varEff,
		  _varRhoFile, _varAcc, _varFSRRes, _varTheory);
      addToVector(mmCovIdx,8, _varYield, _varBkg, _varDetRes, _varEff,
		  _varRhoFile, _varAcc, _varFSRRes, _varTheory);
    }
    if (includeRhoSyst) { // symm case is managed by file names
      // if _varRhoSystFile is present, the uncertainty is handled
      // in adjustXXUnc functions
      addToVector(eeCovIdx,1,_varRhoSystFile);
      addToVector(mmCovIdx,1,_varRhoSystFile);
    }
  }
  else { // all factors, but no _varYield
    addToVector(eeCovIdx,6,  _varBkg, _varDetRes, _varEff,
		_varRhoFile, _varAcc, _varFSRRes);
    addToVector(mmCovIdx,6,  _varBkg, _varDetRes, _varEff,
		_varRhoFile, _varAcc, _varFSRRes);
  }

  if (changeNames && (inpVer==_verMuApproved)) {
    eeOldFName.push_back("cov_ee_varYield_2000.root");
    eeNewFName.push_back("cov_ee_varYield_5000.root");
    eeOldFName.push_back("cov_ee_varRhoFile_2000.root");
    eeNewFName.push_back("cov_ee_varRhoFile_5000.root");
    mmOldFName.push_back("cov_mumu_varYield_2000.root");
    mmNewFName.push_back("dymm_cov_varYield5000_20170202-covOnly.root");
    mmOldFName.push_back("cov_mumu_varRhoFile_2000.root");
    mmNewFName.push_back("cov_mumu_varRhoFile_5000.root");
  }
  else if (changeNames && (inpVer==_verMuMay2017)) {
    eeOldFName.push_back("cov_ee_ElMay2017_varYield_2000.root");
    eeNewFName.push_back("cov_ee_ElMay2017_varYieldPoisson_2000.root");
    mmOldFName.push_back("cov_mumu_varYield_2000.root");
    mmNewFName.push_back("cov_mumu_varYieldPoisson_2000.root");
    if (includeRhoSyst==1) {
      eeOldFName.push_back("cov_ee_ElMay2017_varRhoSystFile_2000.root");
      eeNewFName.push_back("cov_ee_ElMay2017false_varRhoSystFile_2000.root");
      //mmOldFName.push_back("cov_mumu_varRhoSystFile_2000.root");
      //mmNewFName.push_back("");
    }
    else if (includeRhoSyst==2) {
      eeOldFName.push_back("cov_ee_ElMay2017_varRhoSystFile_2000.root");
      eeNewFName.push_back("cov_ee_ElMay2017false_varRhoSystFileSymm_2000.root");
      //mmOldFName.push_back("cov_mumu_varRhoSystFile_2000.root");
      //mmNewFName.push_back("");
    }
  }
  else if (changeNames && (inpVer==_verElMay2017false)) {
    eeOldFName.push_back("cov_ee_varYield_2000.root");
    eeNewFName.push_back("cov_ee_ElMay2017false_varYieldPoisson_2000.root");
    //eeOldFName.push_back("cov_ee_varRhoFile_2000.root");
    //eeNewFName.push_back("cov_ee_varRhoFile_5000.root");
    mmOldFName.push_back("cov_mumu_varYield_2000.root");
    mmNewFName.push_back("cov_mumu_varYieldPoisson_2000.root");
    //mmOldFName.push_back("cov_mumu_varRhoFile_2000.root");
    //mmNewFName.push_back("cov_mumu_varRhoFile_5000.root");
  }
  else if (changeNames && (inpVer==_verMuNov2017)) {
    eeOldFName.push_back("cov_ee_ElMay2017_varYield_2000.root");
    eeNewFName.push_back("cov_ee_ElMay2017_varYieldPoisson_2000.root");
    eeOldFName.push_back("cov_ee_ElMay2017_varRhoSystFile_2000.root");
    eeNewFName.push_back("cov_ee_ElMay2017false_varRhoSystFile_2000.root");
    //eeOldFName.push_back("cov_ee_varRhoFile_2000.root");
    //eeNewFName.push_back("cov_ee_varRhoFile_5000.root");
    mmOldFName.push_back("cov_mumu_varYield_2000.root");
    mmNewFName.push_back("cov_mumu_varYieldPoisson_2000.root");
    //mmOldFName.push_back("cov_mumu_varRhoFile_2000.root");
    //mmNewFName.push_back("cov_mumu_varRhoFile_5000.root");
  }

  if (modifyAccUnc) {
    eeOldFName.push_back("cov_ee_ElMay2017_varTheory_2000.root");
    eeNewFName.push_back("");
    mmOldFName.push_back("cov_mumu_varTheory_2000.root");
    mmNewFName.push_back("");
  }

  TString fileTag="_corr";
  TString plotTag=" (corr)";

  finfomap_t eeCovFNames,mmCovFNames; // covariance derived by randomization
  finfomap_t eeUncTarget_stat, eeUncTarget_syst; // uncertainties from other
  finfomap_t mmUncTarget_stat, mmUncTarget_syst; // ... groups

  for (unsigned int i=0; i<eeCovIdx.size(); i++) {
    TVaried_t var= TVaried_t(eeCovIdx[i]);
    TString fname = "cov_ee_";
    if ((inpVer==_verMuMay2017) || (inpVer==_verMuNov2017)) {
      fname.Append(versionName(_verElMay2017) + "_");
    }
    else if (inpVer==_verElMay2017false) {
      fname.Append(versionName(inpVer) + "_");
    }
    fname+= variedVarName(var) + TString("_2000.root");
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

  PlotCovCorrOpt_t ccOpt;

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
  HERE("\nCovariances loaded");

  // Adjust components of the uncertainties
  CovStruct_t eeCovS(eeCovStat); // creates and nullifies
  CovStruct_t mmCovS(mmCovStat);

  int plotCmpUnc= (hasValue("plotAdjUnc",showCanvs)) ? 1 : 0;
  if (hasValue("plotAdjCov",showCanvs)) plotCmpUnc=2;

  HERE("call adjustMMUnc");
  if (!adjustMMUnc(mmCovFNames,mmCovV,mmCovS,plotChangedCov,plotCmpUnc,noUncAdjustment,inpVer,inpSubVer)) return;
  HERE("call adjustEEUnc");
  if (!adjustEEUnc(eeCovFNames,eeCovV,eeCovS,plotChangedCov,plotCmpUnc,noUncAdjustment,inpVer,inpSubVer)) return;
  HERE("unc adjusted\n");

  //return;

  if (excludeSyst) {
    if (excludeSyst==1) {
      mmCovS.editSystNonAcc().Zero();
      eeCovS.editSystNonAcc().Zero();
    }
    mmCovS.editSystAcc().Zero();
    eeCovS.editSystAcc().Zero();
  }

  TH1D* h1csEE_tmp=loadHisto(eeCSFName,eeCSHisto1DName,"h1csEE_tmp", 1,h1dummy);
  TH1D* h1csMM_tmp=loadHisto(mmCSFName,mmCSHisto1DName,"h1csMM_tmp", 1,h1dummy);
  if (!h1csEE_tmp || !h1csMM_tmp) return;

  // plot input covariance in ee and mm channels
  if (hasValue("plotChannelInputCovOnFileAdj",showCanvs)) {
    double ymin=0, ymax=0, yminRel=0, ymaxRel=0;
    ymin=1e-9; ymax=20.;
    yminRel=1e-3; ymaxRel=10;
    PlotCovCorrOpt_t ccOpt_loc(ccOpt);
    ccOpt_loc.yTitleOffset=1.2;
    PlotCovCorrOpt_t ccOpt_locRel(ccOpt_loc);
    //ccOpt_locRel.noLogScale();
    TString canvNameTag="_onFile";
    if (!noUncAdjustment) canvNameTag+= "_adj";
    mmCovS.Plot("mmCov"+canvNameTag,h1csMM_tmp,ccOpt);
    mmCovS.PlotUnc("mmCovUnc"+canvNameTag,h1csMM_tmp,NULL,ymin,ymax,ccOpt_loc,
		   "#delta#sigma_{#mu#mu}");
    TCanvas *cUncRelMM=mmCovS.PlotUnc("mmCovUncRel"+canvNameTag,
	      h1csMM_tmp,h1csMM_tmp, yminRel,ymaxRel,ccOpt_locRel,
				      "#delta#sigma_{#mu#mu}/#sigma_{#mu#mu}");
    eeCovS.Plot("eeCov"+canvNameTag,h1csEE_tmp,ccOpt);
    eeCovS.PlotUnc("eeCovUnc"+canvNameTag,h1csEE_tmp,NULL,ymin,ymax,ccOpt_loc,
		   "#delta#sigma_{ee}");
    TCanvas *cUncRelEE=eeCovS.PlotUnc("eeCovUncRel"+canvNameTag,
	      h1csEE_tmp,h1csEE_tmp, yminRel,ymaxRel,ccOpt_locRel,
				      "#delta#sigma_{ee}/#sigma_{ee}");
    double dyLeg=0.48;
    moveLegend(cUncRelMM,0,dyLeg);
    moveLegend(cUncRelEE,0,dyLeg);
    //return;
  }

  // check whether mmCovS.Syst_Acc = eeCovS.Syst_Acc
  if (hasValue("plotChannelRelSystAccUnc",showCanvs)) {
    TH1D *h1mmAcc_rel= uncFromCov(mmCovS.covSyst_Acc,"h1mmAcc_rel",
				  h1csMM_tmp,h1csMM_tmp,0);
    TH1D *h1eeAcc_rel= uncFromCov(eeCovS.covSyst_Acc,"h1eeAcc_rel",
				  h1csEE_tmp,h1csEE_tmp,0);
    setNiceMassAxisLabel(h1mmAcc_rel,2,0,"",3,"acc", 1+2);
    setNiceMassAxisLabel(h1eeAcc_rel,2,0,"",3,"acc", 1+2);
    hsRed.SetStyle(h1eeAcc_rel);
    plotHistoAuto(h1eeAcc_rel,"cAccUnc",1,1,"LP","ee");
    plotHistoAuto(h1mmAcc_rel,"cAccUnc",1,1,"LP","#mu#mu");
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

  TMatrixD eeCovTot_noLumi(eeCovTot);
  TMatrixD mmCovTot_noLumi(mmCovTot);
  TMatrixD emCovTot_noLumi(emCovTot);

  double lumiUnc=0.026;
  if ((inpVer==_verMuMay2017) || (inpVer==_verElMay2017false)
      || (inpVer==_verMuNov2017)) {
    lumiUnc=0.023;
  }

  TString fileTagExtra,plotTagExtra;
  TMatrixD lumiCov(2*eeCovTot.GetNrows(),2*eeCovTot.GetNcols());
  lumiCov.Zero();
  if (includeLumiUnc) {
    fileTagExtra="_wLumi"; plotTagExtra=" (wLumi)";
    if (includeLumiUnc==2) {
      eeCovTot.Zero(); mmCovTot.Zero(); emCovTot.Zero();
      fileTagExtra="_wLumi10x"; plotTagExtra=" (wLumi10x)";
    }
    lumiCov-= BLUEResult_t::measCovMatrix(eeCovTot,1);
    lumiCov-= BLUEResult_t::measCovMatrix(mmCovTot,0);
    lumiCov-= BLUEResult_t::measCovMatrix(emCovTot,2);
    if (!addLumiCov(eeCovTot,mmCovTot,emCovTot, lumiUnc, h1csEE_tmp,h1csMM_tmp))
      return;
    lumiCov+= BLUEResult_t::measCovMatrix(eeCovTot,1);
    lumiCov+= BLUEResult_t::measCovMatrix(mmCovTot,0);
    lumiCov+= BLUEResult_t::measCovMatrix(emCovTot,2);
    //lumiCov.Print();
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
  if (useInputYieldUnc) fileTag.Append("_inpYieldUnc");
  if (excludeSyst) fileTag.Append(Form("_excludeSyst%d",excludeSyst));
  if (saveTag.Length()) fileTag.Append(saveTag);


  std::cout << "showCombinationCanvases=" << showCombinationCanvases << "\n";

  double internalScale=1.;
  BLUEResult_t *blue=
    combineData(&eeCovTot,&mmCovTot,&emCovTot,fileTag,plotTag,printCanvases,
		showCombinationCanvases,internalScale,&ccOpt);
  if (!blue) { std::cout << "failed to get BLUE result\n"; return; }

  TH1D *h1csLL_test= convert2histo1D(*blue->getEst(),*blue->getCov(),
				     "h1csLL_test","h1csLL test",h1csEE_tmp);
  setNiceMassAxisLabel(h1csLL_test,2,0,"",1,"",1+0*2);

  if (includeLumiUnc) {
    TMatrixD totFinalCov(*blue->getCov());
    if (!addLumiCov(totFinalCov,-lumiUnc,h1csLL_test)) return;
    plotCovCorr(totFinalCov,NULL,"h2finCovNoLumi","cFinCovNoLumi",
		PlotCovCorrOpt_t(1,1,0));
    TMatrixD lumiContrCov(blue->contributedCov(lumiCov));
    TH2D *h2lumiContrCov=NULL;
    plotCovCorr(lumiContrCov,h1csLL_test,"h2lumiContrCov","cLumiContrCov",
		PlotCovCorrOpt_t(1,1,1),&h2lumiContrCov);
    TH1D *h1uncLumi= uncFromCov(h2lumiContrCov);
    setNiceMassAxisLabel(h1uncLumi,2,0,"",2,"",1+0*2);
    h1uncLumi->GetYaxis()->SetNoExponent(false);
    HistoStyle_t(6,2).SetStyle(h1uncLumi);
    plotHistoAuto(h1uncLumi,"canvContrUnc",1,1,"LP","lumi");
    if (0) {
      TH1D *h1uncLumiRel= uncFromCov(h2lumiContrCov,h1csLL_test);
      printHisto(h1uncLumiRel);
      return;
    }
  }

  // make analysing plots
  if (!detailedCovPlots(blue,corrCase,
			h1csEE_tmp,h1csMM_tmp,
			eeCovS,mmCovS,showCanvs)) return; //,&ccOpt)) return;

  if (printCanvases || (saveTag.Length()>0)) {
    TString outFName="dyll-combi-" + fileTag + ".root";
    TFile fout(outFName,"recreate");
    blue->write(fout,fileTag);
    SaveCanvases("ALL","dyll-combi-"+fileTag,&fout);
    //h2cov_fromYield->Write();
    //h1_dCS_fromYield->Write();
    if (1) {
      if (!writeHistosFromCanvases("cCombiCS2Theory canvContrUnc",1,&fout)) {
	std::cout << "failed to get histos from canvas\n";
      }
    }

    if (1) {
      // save inputs to BLUE method
      fout.mkdir("BLUE_inp");
      fout.cd("BLUE_inp");
      eeCovTot_noLumi.Write("covTotEE_noLumi");
      mmCovTot_noLumi.Write("covTotMM_noLumi");
      emCovTot_noLumi.Write("covTotEM_noLumi");
      lumiCov.Write("covTotLumi");
      h1csEE_tmp->Write("h1csEE");
      h1csMM_tmp->Write("h1csMM");
      h1csLL_test->Write("h1csLL");
      fout.cd();
    }

    writeTimeTag(&fout);
    fout.Close();
    std::cout << "file <" << fout.GetName() << "> created\n";
  }
}

// -------------------------------------------------------

int adjustMMUnc(const finfomap_t &mmCovFNames,
		std::vector<TMatrixD> &mmCovV,
		CovStruct_t &mmCovS, int plotChangedCov,
		int plotCmp, int noUncAdjustment,
		TVersion_t inpVer, TSubVersion_t inpSubVer)
{
  finfomap_t mmUncFiles, mmUncStat, mmUncSyst;
  fvarweightmap_t relUncW;

  //TH1D* h1csEE_loc=loadHisto(eeCSFName,eeCSHisto1DName,"h1csEE_loc", 1,h1dummy);
  TH1D* h1csMM_loc=loadHisto(mmCSFName,mmCSHisto1DName,"h1csMM_loc", 1,h1dummy);
  if (!h1csMM_loc) return 0;

  TString inpPath="./";
  TString inpFile=inpPath + "DYmm_ROOTFile_Input-20170207.root";
  if ((inpVer==_verMuMay2017) || (inpVer==_verElMay2017false)) {
    inpFile="./DYmm_ROOTFile_Input-20170504.root";
    if (inpSubVer==_subver2) {
      inpFile="./DYmm_ROOTFile_Input_Cov_v2-20170529.root";
    }
    else if (inpSubVer==_subver3) {
      inpFile="./DYmm_ROOTFile_Input_Cov_v3-20170614.root";
      std::cout << "version3\n";
    }
  }
  else if (inpVer==_verMuNov2017) {
    inpFile="./ROOTFile_Input_Cov_v4-KPLee20171124.root";
  }

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
  relUncW   [_varYield]= (useInputYieldUnc) ? 1. : 0.;
  mmUncFiles[_varBkg]= inpFile;
  mmUncStat [_varBkg]= "h_RelUnc_Stat_BkgEst";
  mmUncSyst [_varBkg]= "h_RelUnc_Syst_BkgEst";
  relUncW   [_varBkg]= 1;
  if (((inpVer==_verMuMay2017) || (inpVer==_verElMay2017false)) &&
      valueEquals(int(inpSubVer),Form("%d",_subver1))) {
    mmUncFiles[_varBkg]= "DYmm_ROOTFile_RelUnc_Stat_Syst_Tot_BkgEst-20170519.root";
    mmUncStat[_varBkg]= "h_RelUnc_Stat";
    mmUncSyst[_varBkg]= "h_RelUnc_Syst";
    relUncW  [_varBkg]= 0.01;
  }
  //if (inpVer==_verMuNov2017) relUncW[_varBkg]=100;
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
  if (modifyAccUnc) {
    mmUncFiles[_varTheory]= "";  // place holder
    mmUncStat [_varTheory]= "";
    mmUncSyst [_varTheory]= "";
    relUncW   [_varTheory]= 0;
  }

  if (findVaried(mmCovFNames,_varRhoSystFile,0)!=-1) {
    mmUncSyst [_varRhoFile]= "";
    mmUncFiles[_varRhoSystFile]= inpFile;
    mmUncStat [_varRhoSystFile]= (1) ? "h_RelUnc_Syst_EffSF" : "";
    mmUncSyst [_varRhoSystFile]= "";
    relUncW   [_varRhoSystFile]= 1;
  }

  const TH1D *h1central= (1) ? h1csMM : NULL;
  int change_covV=(noUncAdjustment) ? 0:1;
  if(change_covV && plotChangedCov) change_covV++;
  std::string specFlags="";
  if (modifyAccUnc==1) specFlags+=" dy13TeV_acc_correction";
  int res=adjustUnc( "mm",
		     mmCovFNames, mmCovV,
		     mmUncFiles,mmUncStat,mmUncSyst,relUncW,
		     h1csMM,h1central,
		     plotCmp,change_covV,
		     "_mm",
		     mmCovS,
		     specFlags);
  if (!res) {
    std::cout << "error in adjustMMUnc\n";
    return 0;
  }

  // study DetRes
  if (0) {
    int idx= findVaried(mmCovFNames,_varDetRes);
    if (idx==-1) {
      std::cout << "failed to find detRes (needed for the electron channel)\n";
      return 0;
    }
    std::cout << "idx(_varDetRes)=" << idx << "\n";
    TH1D *h1KL= loadHisto(mmUncFiles[_varDetRes],mmUncStat[_varDetRes],
			      "h1KL",h1dummy);
    h1KL->Multiply(h1csMM_loc);
    TH1D* h1uncFromCovDetRes= uncFromCov(mmCovV[idx],
					 "h1uncFromCovDetRes",
					 h1KL);
    printRatio(h1KL,h1uncFromCovDetRes);
    return 0;
  }


  return 1;
}
// -----------------------------------------------------------------------

int adjustEEUnc(const finfomap_t &eeCovFNames,
		std::vector<TMatrixD> &eeCovV,
		CovStruct_t &eeCovS, int plotChangedCov,
		int plotCmp, int noUncAdjustment,
		TVersion_t inpVer, TSubVersion_t inpSubVer)
{
  finfomap_t eeUncFiles, eeUncStat, eeUncSyst;
  fvarweightmap_t relUncW;

  TH1D* h1csEE_loc=loadHisto(eeCSFName, eeCSHisto1DName, "h1csEE_loc", 1,h1dummy);
  //TH1D* h1csMM_loc=loadHisto(mmCSFName, mmCSHisto1DName, "h1csMM_loc", 1,h1dummy);
  if (!h1csEE_loc) return 0;

  TString inpPath="./";
  TString inpFile=inpPath + "DYee_ROOTFile_Input-20170208.root";
  if ((inpVer==_verMuMay2017) || (inpVer==_verElMay2017false)) {
    //inpFile="DYee_ROOTFile_Input_v1-20170518.root";
    //inpFile="./DYee_ROOTFile_Input_v2-20170526.root";
    inpFile="./DYee_ROOTFile_Input_v3-20170529.root";
  }
  else if (inpVer==_verMuNov2017) {
    inpFile="./ROOTFile_Input_Electron_v4_RC20171128.root";
  }
  if (inpSubVer==_subverUndef) { ; } // subver has no effect

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
  relUncW   [_varYield]= (useInputYieldUnc) ? 0.01 : 0.;
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
  if (modifyAccUnc) {
    eeUncFiles[_varTheory]= ""; // place holder
    eeUncStat [_varTheory]= "";
    eeUncSyst [_varTheory]= "";
    relUncW   [_varTheory]= 0;
  }

  if (0 && ((inpVer==_verMuMay2017) || (inpVer==_verElMay2017false))) {
    eeUncFiles[_varAcc]= "DYmm_ROOTFile_Input-20170504.root";
    eeUncSyst [_varAcc]= "h_RelUnc_Syst_Acc";
    relUncW   [_varAcc]= 1.;
  }

  if (findVaried(eeCovFNames,_varRhoSystFile,0)!=-1) {
    eeUncSyst [_varRhoFile]= "";
    eeUncFiles[_varRhoSystFile]= inpFile;
    eeUncStat [_varRhoSystFile]= (1) ? "h_RelUnc_Syst_EffSF" : "";
    eeUncSyst [_varRhoSystFile]= "";
    relUncW   [_varRhoSystFile]= 0.01;
  }

  const TH1D *h1central= (1) ? h1csEE : NULL;
  int change_covV=(noUncAdjustment) ? 0:1;
  if(change_covV && plotChangedCov) change_covV++;
  std::string specFlags;
  specFlags="dyee13TeV_detRes_correction";
  if (modifyAccUnc==1) specFlags+=" dy13TeV_acc_correction";
  int res=adjustUnc( "ee",
		     eeCovFNames, eeCovV,
		     eeUncFiles,eeUncStat,eeUncSyst,relUncW,
		     h1csEE,h1central,
		     plotCmp,change_covV,
		     "_ee",
		     eeCovS,
		     specFlags);

  if (!res) {
    std::cout << "error in adjustEEUnc\n";
    return 0;
  }
  return 1;
}
// -----------------------------------------------------------------------
