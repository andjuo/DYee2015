#include "CovStruct.h"
#include <TPaveText.h>

// -------------------------------------------------------
// -------------------------------------------------------

CovStruct_t* detailedCovPlots(const BLUEResult_t *blue,
			      int corrCase,
			      const TH1D *h1csEE,
			      const TH1D *h1csMM,
			      const CovStruct_t &eeCovS,
			      const CovStruct_t &mmCovS,
			      std::string showCanvs,
			      const PlotCovCorrOpt_t *optCC_user)
{

  if (!blue || !h1csEE || !h1csMM) {
    std::cout << "detailedCovPlots: null ptr\n";
    return 0;
  }
  PlotCovCorrOpt_t optCC(1,1,1,1.8,0.15,0.15);
  if (optCC_user) optCC.assign(*optCC_user);

  if (hasValue("plotInputContributedCovs",showCanvs)) {
    for (int i=1; i<=4; i++) {
      if ((eeCovS.plotStatYieldOnly ||
	   mmCovS.plotStatYieldOnly) && (i>1)) break;
      TString tagName="ee_" + eeCovS.GetPartName(i);
      TMatrixD partM= eeCovS.GetPart(i);
      //plotCovCorr(partM,h1csEE_tmp,"h2_"+tagName,"c2_"+tagName,optCC);
      TMatrixD corr = blue->contributedCorrPartial(partM,1);
      TString cName= "cContrCov_" + tagName;
      TH2D *h2= convert2histo(corr,h1csEE,"h2"+tagName,"h2"+tagName);
      logAxis(h2);
      plotHisto(h2,cName,optCC);
    }
    for (int i=1; i<=4; i++) {
      TString tagName="mm_" + mmCovS.GetPartName(i);
      TMatrixD partM= mmCovS.GetPart(i);
      //plotCovCorr(partM,h1csMM_tmp,"h2_"+tagName,"c2_"+tagName,optCC);
      TMatrixD corr = blue->contributedCorrPartial(partM,0);
      TString cName= "cContrCov_" + tagName;
      TH2D *h2= convert2histo(corr,h1csMM,"h2"+tagName,"h2"+tagName);
      logAxis(h2);
      plotHisto(h2,cName,optCC);
    }
  }

  TH1D* h1csLL_tmp= cloneHisto(h1csEE, "h1csLL_tmp", "h1csLL_tmp");
  h1csLL_tmp->GetXaxis()->SetTitle( niceMassAxisLabel(2,"",0) );
  h1csLL_tmp->GetYaxis()->SetTitle( niceMassAxisLabel(2,"",1) );

  TH1D *h1csLL_badBinning= convert2histo1D(*blue->getEst(),*blue->getCov(),
			   "h1csLL_badBinning","h1csLL_badBinning",NULL);
  if (!h1csLL_badBinning) return 0;
  if (!copyContents(h1csLL_tmp, h1csLL_badBinning)) return 0;

  TH1D* h1deltacsLL_tmp= cloneHisto(h1csLL_tmp,
				    "h1deltacsLL_tmp", "h1deltacsLL_tmp");
  h1deltacsLL_tmp->GetXaxis()->SetTitle( niceMassAxisLabel(2,"",0) );
  h1deltacsLL_tmp->GetYaxis()->SetTitle( "\\delta" + niceMassAxisLabel(2,"",1) );

  PlotCovCorrOpt_t optCCNoLog;
  optCCNoLog.setLogScale(0,0);

  if (0) { // check function constructMeasCov
    TMatrixD covM_ee_yield= blue->measCovMatrix(eeCovS.getStatYield(),1);
    TMatrixD covM_mm_yield= blue->measCovMatrix(mmCovS.getStatYield(),0);
    //covM_ee_yield.Print();
    //covM_mm_yield.Print();
    TMatrixD covM_yield(covM_ee_yield,TMatrixD::kPlus,covM_mm_yield);
    plotCovCorr(covM_yield,NULL,"h2covMeasYield_chk","cCovMeasYield_chk");
    TMatrixD covFin_yield= blue->contributedCov(covM_yield);
    TH2D* h2cov_fromYield=NULL;
    plotCovCorr(covFin_yield,h1csLL_tmp,
		"h2covFin_fromYield_chk","cCovFin_fromYield_chk",
		optCC,&h2cov_fromYield);
    TH1D *h1_dCS_fromYield= uncFromCov(covFin_yield,
				       "h1_dCS_fromYield_chk",
				       h1deltacsLL_tmp,NULL,0);
    printHisto(h2cov_fromYield);
    printHisto(h1_dCS_fromYield);
  }

  /*
  TMatrixD sumContributedCov(TMatrixD::kZero,eeCovS.getStatYield());
  CovStruct_t llCovS(eeCovS.getStatYield()); // creates and nullifies
  for (int iFlag=1; iFlag<=4; iFlag++) {
    if ((eeCovS.plotStatYieldOnly ||
	 mmCovS.plotStatYieldOnly) && (iFlag>1)) break;
    CovStruct_t::TCovPart_t part=CovStruct_t::TCovPart_t(iFlag);
    double accCorrFlag= (corrCase==1) ? 1. : 0.;
    TMatrixD measCov= constructMeasCov(eeCovS,mmCovS,part,accCorrFlag,blue);
    TString label= eeCovS.GetPartName(iFlag);
    if (0) plotCovCorr(measCov,NULL,"h2covMeas"+label,"cCovMeas_"+label,
		       optCCNoLog,NULL);
    TMatrixD covFin_contrib= blue->contributedCov(measCov);
    llCovS.SetPart(part,covFin_contrib);
    sumContributedCov+=covFin_contrib;
    if (hasValue("plotFinalCovByType",showCanvs)) {
      TH2D* h2cov_contrib=NULL;
      plotCovCorr(covFin_contrib,h1csLL_tmp,
		  "h2covFin_from"+label,"cCovFin_from_"+label,
		  PlotCovCorrOpt_t(),&h2cov_contrib);
      TH1D *h1_dCS_contrib= uncFromCov(covFin_contrib,
				       "h1_dCS_from_"+label,
				       h1deltacsLL_tmp,NULL,0);
      if (!h2cov_contrib || !h1_dCS_contrib) return 0;
      hsVec[iFlag-1].SetStyle(h1_dCS_contrib);
      plotHistoAuto(h1_dCS_contrib,"canvContrUnc",1,1,"LPE",label);
    }
  }
  */

  if (hasValue("plotFinalCov",showCanvs)) {
    plotCovCorr(*blue->getCov(),h1csLL_tmp,"h2covFin","cCovFin",optCC);
  }

  //TMatrixD sumContributedCov(TMatrixD::kZero,eeCovS.getStatYield());
  int res=1;
  CovStruct_t llCovS= constructCombCov(eeCovS,mmCovS,blue,corrCase,res);
  if (!res) {
    std::cout << "failed to construct CombCov\n";
    return NULL;
  }
  TMatrixD sumContributedCov= llCovS.Sum();
  if (hasValue("plotFinalCovByType",showCanvs)) {
    for (int iFlag=1; iFlag<=4; iFlag++) {
      if (llCovS.plotStatYieldOnly && (iFlag>1)) break;
      //CovStruct_t::TCovPart_t part=CovStruct_t::TCovPart_t(iFlag);
      //double accCorrFlag= (corrCase==1) ? 1. : 0.;
      TString label= llCovS.GetPartName(iFlag);
      TMatrixD covFin_contrib= llCovS.GetPart(iFlag);
      TH2D* h2cov_contrib=NULL;
      plotCovCorr(covFin_contrib,h1csLL_tmp,
		  "h2covFin_from"+label,"cCovFin_from_"+label,
		  optCC,&h2cov_contrib);
      TH1D *h1_dCS_contrib= uncFromCov(covFin_contrib,
				       "h1_dCS_from_"+label,
				       h1deltacsLL_tmp,NULL,0);
      if (!h2cov_contrib || !h1_dCS_contrib) return 0;
      hsVec[iFlag-1].SetStyle(h1_dCS_contrib);
      plotHistoAuto(h1_dCS_contrib,"canvContrUnc",optCC.logScaleX,1,"LPE",label);
    }
  }

  if (1) {
    TString eeStr,mmStr,llStr;
    if (!eeCovS.PrintZRangeUnc("ee",h1csEE,0,&eeStr) ||
	!mmCovS.PrintZRangeUnc("mm",h1csMM,0,&mmStr) ||
	!llCovS.PrintZRangeUnc("ll",h1csLL_tmp,0,&llStr)) {
      std::cout << "error\n";
    }
    if (1) {
      TCanvas *cx= new TCanvas("cIntCS","cIntCS",800,300);
      if (!cx) return 0;
      TPaveText *pt= new TPaveText(0.05,0.1,0.95,0.8);
      pt->AddText("Integrated cross section (60-120 GeV)");
      pt->AddText(eeStr);
      pt->AddText(mmStr);
      pt->AddText(llStr);
      pt->Draw();
    }

    if (0) {
      std::cout << "\n EE    MM    Combi-XSect\n";
      const TH1D *h1def=h1csEE;
      for (int ibin=1; ibin<=h1def->GetNbinsX(); ibin++) {
	std::cout << "ibin=" << ibin << " " << h1def->GetBinLowEdge(ibin)
		  << " -- "
		  << (h1def->GetBinLowEdge(ibin)+h1def->GetBinWidth(ibin))
		  << "  "
		  << h1csEE->GetBinContent(ibin) << "   "
		  << h1csMM->GetBinContent(ibin) << "   "
		  << h1csLL_tmp->GetBinContent(ibin) << "\n";
      }
    }
  }

  //TMatrixD sumContrChk( sumContributedCov, TMatrixD::kMinus, *blue->getCov() );
  //sumContrChk.Print();

  return new CovStruct_t(llCovS);
}


// -----------------------------------------------------------------------

int aggregateUnc(const finfomap_t &fnames,
		 const std::vector<TMatrixD> &covStatV,
		 const std::vector<TH1D*> &h1SystV,
		 CovStruct_t &covS)
{

  covS.Zero();

  int idx= findVaried(fnames,_varYield);
  if (idx==-1) {
    std::cout << "failed to find varYield\n";
    return 0;
  }
  covS.covStat_Yield= covStatV[idx];

  idx=0;
  for (finfomap_t::const_iterator it=fnames.begin(); it!=fnames.end();
       it++, idx++) {
    if (it->first == _varYield) continue;
    covS.covStat_nonYield += covStatV[idx];
  }

  if (fnames.size() != h1SystV.size()) {
    std::cout << "aggregateUnc: fnames.size=" << fnames.size() << ", "
	      << h1SystV.size() << "\n";
    return 0;
  }

  idx=0;
  for (finfomap_t::const_iterator it=fnames.begin(); it!=fnames.end();
       it++, idx++) {
    if (h1SystV[idx]->Integral()==0) {
      std::cout << "aggregate systematic error: skipping "
		<< variedVarName(it->first) << " (no input)\n";
      continue;
    }
    if (it->first == _varAcc) {
      covS.covSyst_Acc = convert2cov1D(h1SystV[idx],0);
    }
    else {
      TMatrixD tmpM= convert2cov1D(h1SystV[idx],0);
      std::cout << "h1SystV[idx].size=" << h1SystV[idx]->GetNbinsX()
		<< ", tmpM.size=" << tmpM.GetNrows() << " x " << tmpM.GetNcols()
		<< ", covS.covSyst_nonAcc.size="
		<< covS.covSyst_nonAcc.GetNrows() << " x "
		<< covS.covSyst_nonAcc.GetNcols() << "\n";
      covS.covSyst_nonAcc += convert2cov1D(h1SystV[idx],0);
    }
  }

  return 1;
}

// -------------------------------------------------------

int adjustUnc(TString lepTag,
	      const finfomap_t &covFNames,
	      std::vector<TMatrixD> &covV,
	      const finfomap_t &uncFiles,
	      const finfomap_t &uncStat,
	      const finfomap_t &uncSyst,
	      const fvarweightmap_t &relUncW,
	      const TH1D *h1binning,
	      const TH1D *h1central,
	      int plotCmp, // 1 - plot uncertainties, 2 - and covariances
	      int change_covV, // 1 - change values, 2 - verify change by plot
	      TString tag,
	      CovStruct_t &covS,
	      std::string specialArgs)
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

  std::vector<TH1D*> h1uncStatV, h1uncSystV;
  if (!loadUncData(lepTag,uncFiles,uncStat,relUncW,
		   h1binning,h1central,h1uncStatV,tag+"_stat")) {
    std::cout << "failed to load statistical component\n";
    return 0;
  }
  if (!loadUncData(lepTag,uncFiles,uncSyst,relUncW,
		   h1binning,h1central,h1uncSystV,tag+"_syst")) {
    std::cout << "failed to load systematic component\n";
    return 0;
  }

  // correct for the input: Stat_DetRes in the electron channel is a part of
  // Syst_DetRes
  if (hasValue("dyee13TeV_detRes_correction",specialArgs)) {
    int idx= findVaried(covFNames,_varDetRes);
    if (idx==-1) {
      std::cout << "failed to find detRes (needed for the electron channel)\n";
      return 0;
    }
    std::cout << "idx(_varDetRes)=" << idx << "\n";
    //printHisto(h1uncStatV[idx]);
    //printHisto(h1uncSystV[idx]);
    if (!addInQuadrature(h1uncSystV[idx],h1uncStatV[idx],1,1.)) {
      return 0;
    }
    h1uncStatV[idx]->Reset();
    std::cout << "cleared h1detRes_ee statistical error from RC\n";
    //return 0;
  }

  // Special correction for acceptance uncertainty
  if (hasValue("dy13TeV_acc_correction",specialArgs)) {
    int idxAcc= findVaried(covFNames,_varAcc);
    int idxTh= findVaried(covFNames,_varTheory);
    if ((idxAcc==-1) || (idxTh==-1)) {
      std::cout << "failed to find _varAcc or _varTheory the channel "
		<< lepTag << ": " << idxAcc << ", " << idxTh << "\n";
      return 0;
    }
    std::cout << "idxAcc, idxTh=" << idxAcc << ", " << idxTh << "\n";

    TH1D *h1th= loadHisto("DYmm_ROOTFile_Input-20170504.root",
			  "h_RelUnc_Syst_Acc","h1acc_th"+lepTag,1,h1dummy);
    TH1D *h1ee= loadHisto("DYee_ROOTFile_Input-20170208.root",
			  "h_RelUnc_Syst_Acc","h1acc_ee"+lepTag,1,h1dummy);

    if (1) {
      std::cout << "lepTag=" << lepTag << "loaded unc\n";
      printHisto(h1th);
      printHisto(h1ee);
    }

    h1th->Scale(0.5);
    h1th->Add(h1ee,0.5*0.01);
    if (h1central) h1th->Multiply(h1central);

    if (1) {
      std::cout << "lepTag=" << lepTag << "  before change\n";
      printHisto(h1th);
      printHisto(h1uncSystV[idxTh]);
      printHisto(h1uncSystV[idxAcc]);
    }

    h1uncSystV[idxTh]->Add(h1uncSystV[idxAcc]);
    int resTmp= addInQuadrature(h1uncSystV[idxTh],h1th,1,-1.) ? 1:0;
    copyContents(h1uncSystV[idxAcc],h1th);

    delete h1th;
    delete h1ee;

    if (1) {
      std::cout << "lepTag=" << lepTag << "  after change\n";
      printHisto(h1uncSystV[idxTh]);
      printHisto(h1uncSystV[idxAcc]);
    }
    if (!resTmp) {
      std::cout << "  removing negative entries\n";
      removeNegatives(h1uncSystV[idxTh]);
    }
    //if (lepTag=="ee") return 0;
  }

  // eliminate eff syst. uncertainty
  if (1 && change_covV) {
    int idx1= findVaried(covFNames,_varEff);
    int idx2= findVaried(covFNames,_varAcc);
    if ((idx1==-1) || (idx2==-1)) {
      std::cout << "failed to find eff or acc\n";
      return 0;
    }
    std::cout << "idx(_varEff)=" << idx1 << "\n";
    std::cout << "idx(_varAcc)=" << idx2 << "\n";
    if (change_covV) {
      covV[idx1].Zero();
      covV[idx2].Zero();
    }
    /*
    if (idx==0) { std::cout << "  --- idx=0, fix the code\n"; return 0; }
    if (!h1uncStatV[idx]) {
      h1uncStatV[idx]= cloneHisto(h1uncStatV[0],Form("h1eff_%d",isEE),"h1Eff");
      h1uncStatV[idx]->Reset();
    }
    */
  }

  int res=compareUncAndScale(labelV,h1uncStatV,covV,plotCmp,change_covV,tag);
  if (res && plotCmp && (change_covV==2))
    res=compareUncAndScale(labelV,h1uncStatV,covV,2,0,tag+"_adj");

  if (res) {
    res= aggregateUnc(covFNames,covV,h1uncSystV,covS);
  }

  if (!res) {
    std::cout << "error in adjustUnc\n";
    return 0;
  }

  return 1;
}

// -------------------------------------------------------
// -------------------------------------------------------

TMatrixD constructEMAccCov(const CovStruct_t &eeCovS, const CovStruct_t &mmCovS,
			   double corrFactor)
{
  TMatrixD emCovTot(TMatrixD::kZero, eeCovS.covSyst_Acc);
  for (int ir=0; ir<eeCovS.covSyst_Acc.GetNrows(); ir++) {
    for (int ic=0; ic<mmCovS.covSyst_Acc.GetNcols(); ic++) {
      emCovTot(ir,ic) = corrFactor *
	sqrt(eeCovS.covSyst_Acc(ir,ic)) *
	sqrt(mmCovS.covSyst_Acc(ir,ic));
    }
  }
  return emCovTot;
}

// -------------------------------------------------------

TMatrixD constructMeasCov(const CovStruct_t &eeCovS, const CovStruct_t &mmCovS,
			  CovStruct_t::TCovPart_t covFlag, double accCorrFactor,
			  const BLUEResult_t *blue)
{
  TMatrixD measCov_ee(blue->measCovMatrix(eeCovS.GetPart(covFlag),1));
  TMatrixD measCov_mm(blue->measCovMatrix(mmCovS.GetPart(covFlag),0));
  TMatrixD measCov( measCov_ee, TMatrixD::kPlus, measCov_mm );
  if (covFlag==CovStruct_t::_varCov_systAcc) {
    TMatrixD emCov( constructEMAccCov(eeCovS,mmCovS,accCorrFactor) );
    TMatrixD measEMCov( blue->measCovMatrix(emCov,2) );
    measCov += measEMCov;
  }
  return measCov;
}

// -------------------------------------------------------

CovStruct_t constructCombCov(const CovStruct_t &eeCovS,
			     const CovStruct_t &mmCovS,
			     const BLUEResult_t *blue,
			     int corrAccCase,
			     int &flag)
{
  CovStruct_t llCovS(eeCovS.getStatYield()); // creates and nullifies
  llCovS.plotStatYieldOnly=
    (eeCovS.plotStatYieldOnly || mmCovS.plotStatYieldOnly) ? 1:0;
  flag=0;
  if (!blue) return llCovS;

  for (int iFlag=1; iFlag<=4; iFlag++) {
    CovStruct_t::TCovPart_t part=CovStruct_t::TCovPart_t(iFlag);
    TMatrixD measCov= constructMeasCov(eeCovS,mmCovS,part,corrAccCase,blue);
    TString label= eeCovS.GetPartName(iFlag);
    if (0) {
      PlotCovCorrOpt_t optCCNoLog;
      optCCNoLog.setLogScale(0,0);
      plotCovCorr(measCov,NULL,"h2covMeas"+label,"cCovMeas_"+label,
		  optCCNoLog,NULL);
    }
    TMatrixD covFin_contrib= blue->contributedCov(measCov);
    llCovS.SetPart(part,covFin_contrib);
  }
  flag=1;
  return llCovS;
}


// -------------------------------------------------------

