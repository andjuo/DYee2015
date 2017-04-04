#include "inputs.h"
#include "crossSection.h"
#include "Blue.h"
#include "analyseBLUEResult.h"
#include "CovStruct.h"
#include "helper_8TeV.h"

// -------------------------------------------------------
// -------------------------------------------------------

const int useUncorrYieldUnc=0;
TString inpFName="out_dy8TeV_1D.root";

// -------------------------------------------------------

void work8TeV_1D(int printCanvases=0, int corrCase=1, int includeLumiUnc=0,
		 int plotChangedCov=0,
		 int excludeSyst=0, // 1 - all, 2 - Acc only
		 TString saveTag="",
		 double use_scale_user=-1)
{
  closeCanvases(5);
  const int showOnlyCovStatYield=1;

  eeCSFName=inpFName;
  mmCSFName=inpFName;
  theoryCSFName=inpFName;
  eeCSHisto1DName="h1ee_orig";
  mmCSHisto1DName="h1mm_orig";
  theoryCSHisto1DName="h1ll_orig";
  theoryCSHisto1D_perMBW=1;

  std::string showCanvs;
  showCanvs+=" plotChannelInputCov";
 // whether relative uncertainty in Acc in the channels is similar
  //showCanvs+=" plotChannelRelSystAccUnc";
  // compare randomized vs provided uncertainties
  //showCanvs+=" plotAdjUnc";
  //showCanvs+=" plotAdjCov";
  showCanvs+=" plotInputContributedCovs"; // 4 canvases in each channel
  showCanvs+= " plotFinalCovByType"; // 4 canvases in the combined channel

  std::string showCombinationCanvases;
  showCombinationCanvases="ALL";
  //
  //showCanvs+= " showCombiCSCheck";
  //showCanvs+= " showCombiCSUnc";

  TString fileTag="_corr";
  TString plotTag=" (corr)";

  std::vector<TH2D*> eeCovH2D, mmCovH2D;
  std::vector<TMatrixD> eeCovV, mmCovV;
  TMatrixD eeCovStat(41,41);
  TMatrixD eeCovSyst(eeCovStat);
  TMatrixD mmCovStat(eeCovStat);
  TMatrixD mmCovSyst(eeCovStat);

  
  TFile fin(inpFName);
  TH1D* h1csEE_tmp=loadHisto(fin, eeCSHisto1DName, "h1csEE_tmp", 1,h1dummy);
  TH1D* h1csMM_tmp=loadHisto(fin, mmCSHisto1DName, "h1csMM_tmp", 1,h1dummy);
  TH1D* h1csLL_orig=loadHisto(fin,"h1ll_orig","h1csLL_orig",1,h1dummy);
  TH2D *h2covEE= loadHisto(fin,"h2eeCov_orig","h2eeCov_orig",1,h2dummy);
  TH2D *h2covMM= loadHisto(fin,"h2mmCov_orig","h2mmCov_orig",1,h2dummy);
  TH2D *h2covLL_orig=loadHisto(fin,"h2llCov_orig","h2llCov_orig",1,h2dummy);
  fin.Close();

  if (!h1csEE_tmp || !h1csMM_tmp || !h1csLL_orig ||
      !h2covEE || !h2covMM || !h2covLL_orig) return;

  // xxCovStat is the covariance of measured yield
  // xxCovSyst equals to other covariances
  // xxCovV contains all covariance matrices
  //if (!loadCovData("ee",eeCovFNames,eeCovH2D,eeCovV,eeCovStat,eeCovSyst) ||
  //    !loadCovData("mm",mmCovFNames,mmCovH2D,mmCovV,mmCovStat,mmCovSyst)) {
  //  std::cout << "failed to load covariances\n";
  //  return;
  //}
  
  eeCovH2D.push_back(h2covEE);
  mmCovH2D.push_back(h2covMM);
  eeCovStat= convert2mat(h2covEE);
  mmCovStat= convert2mat(h2covMM);
  eeCovV.push_back(eeCovStat);
  mmCovV.push_back(mmCovStat);

  plotCovCorr(h2covEE,"cCovEE");

  // Adjust components of the uncertainties
  CovStruct_t eeCovS(eeCovStat); // creates and nullifies
  CovStruct_t mmCovS(mmCovStat);
  eeCovS.plotStatYieldOnly= showOnlyCovStatYield;
  mmCovS.plotStatYieldOnly= showOnlyCovStatYield;
  eeCovS.covStat_Yield=eeCovStat;
  mmCovS.covStat_Yield=mmCovStat;

  int plotCmpUnc= (hasValue("plotAdjUnc",showCanvs)) ? 1 : 0;
  if (hasValue("plotAdjCov",showCanvs)) plotCmpUnc=2;

  if (excludeSyst) {
    if (excludeSyst==1) {
      mmCovS.editSystNonAcc().Zero();
      eeCovS.editSystNonAcc().Zero();
    }
    mmCovS.editSystAcc().Zero();
    eeCovS.editSystAcc().Zero();
  }

  // plot input covariance in ee and mm channels
  if (hasValue("plotChannelInputCov",showCanvs)) {
    double ymin=0, ymax=0;
    ymin=1e-9; ymax=20.;
    mmCovS.Plot("mmCov",h1csMM_tmp);
    mmCovS.PlotUnc("mmCovUnc_",h1csMM_tmp,NULL,ymin,ymax);
    eeCovS.Plot("eeCov",h1csEE_tmp);
    eeCovS.PlotUnc("eeCovUnc_",h1csEE_tmp,NULL,ymin,ymax);
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

  if (0) {
    plotCovCorr(eeCovTot,h1csEE_tmp,"h2eeCov_tmp","cCovEEtot_chk");
    plotCovCorr(mmCovTot,h1csMM_tmp,"h2mmCov_tmp","cCovMMtot_chk");
    return;
  }

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

  if (useUncorrYieldUnc) fileTag.Append("_uncorrYieldUnc");
  if (excludeSyst) fileTag.Append(Form("_excludeSyst%d",excludeSyst));
  if (saveTag.Length()) fileTag.Append(saveTag);

  double use_scale=1e1;
  if (use_scale_user>0) use_scale=use_scale_user;
  BLUEResult_t *blue=
    combineData(&eeCovTot,&mmCovTot,&emCovTot,fileTag,plotTag,printCanvases,
		showCombinationCanvases,use_scale);
  if (!blue) { std::cout << "failed to get BLUE result\n"; return; }

  //return;

  if (includeLumiUnc) {
    TMatrixD totFinalCov(*blue->getCov());
    TH1D *h1csLL_test= convert2histo1D(*blue->getEst(),totFinalCov,
				       "h1csLL_test","h1csLL test",NULL);
    if (!addLumiCov(totFinalCov,-lumiUnc,h1csLL_test)) return;
    plotCovCorr(totFinalCov,NULL,"h2finCovNoLumi","cFinCovNoLumi");
  }

  //return;

  // make analysing plots
  if (!detailedCovPlots(blue,corrCase,
			h1csEE_tmp,h1csMM_tmp,
			eeCovS,mmCovS,showCanvs)) return;

  // check the values and the uncertainties
  if (1) {
    TString massStr= niceMassAxisLabel(2,"");
    TString sigmaStr= niceMassAxisLabel(2,"",1);
    TString axisStr= ";" + massStr + ";" + sigmaStr;
    TH1D *h1comb_loc= convert2histo1D(*blue->getEst(),*blue->getCov(),
				      "h1Combi_loc",
				      "combined measurement " + axisStr,
				      h1csEE_tmp,0);
    hsBlack.SetStyle(h1comb_loc);
    logAxis(h1comb_loc);
    TH1D *h1comb_diff= cloneHisto(h1csLL_orig,"h1comb_diff","h1comb_diff");
    h1comb_diff->Add(h1comb_loc,-1);
    removeError(h1comb_diff);
    printHisto(h1comb_diff,0);
    TH1D *h1combUnc_orig= errorAsCentral(h1csLL_orig);
    TH1D *h1combUnc= errorAsCentral(h1comb_loc);
    TH1D *h1combUnc_diff= cloneHisto(h1combUnc_orig,
				     "h1combUnc_diff","h1combUnc_diff");
    h1combUnc_diff->Add(h1combUnc,-1);
    printHisto(h1combUnc_diff,0);
 
    TH1D *h1comb_ratio= cloneHisto(h1csLL_orig,
			   "h1comb_ratio;"+massStr+";orig/new","h1comb_ratio");
    h1comb_ratio->Divide(h1comb_loc);
    removeError(h1comb_ratio);
    plotHisto(h1comb_ratio,"cCombRatio",1,0,"LP");

    plotHisto(h1combUnc_orig,"cCombUnc",1,1,"LP","orig");
    plotHistoSame(h1combUnc,"cCombUnc","LP","new comb");
    if (use_scale_user>0)
      std::cout << "use_scale_user=" << use_scale_user << "\n";
    return;
  }

  if (printCanvases || (saveTag.Length()>0)) {
    TString outFName="dyll8TeV-combi-" + fileTag + ".root";
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
