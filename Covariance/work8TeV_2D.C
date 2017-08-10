#include "inputs.h"
#include "crossSection.h"
#include "Blue.h"
#include "analyseBLUEResult.h"
#include "CovStruct.h"
#include "helper_8TeV.h"

// -------------------------------------------------------
// -------------------------------------------------------

const int useUncorrYieldUnc=0;
TString inpFName="out_dy8TeV_2D_mdf.root";

// -------------------------------------------------------

void work8TeV_2D(int printCanvases=0, int corrCase=1, int includeLumiUnc=0,
		 int plotChangedCov=0,
		 int excludeSyst=0, // 1 - all, 2 - Acc only
		 TString saveTag="",
		 double use_scale_user=-1,
		 int specSave=0)
{
  closeCanvases(5);
  const int showOnlyCovStatYield=1;
  const int logX=0;

  eeCSFName=inpFName;
  mmCSFName=inpFName;
  theoryCSFName=inpFName;
  eeCSHisto1DName="h1ee_orig";
  mmCSHisto1DName="h1mm_orig";
  theoryCSHisto1DName="h1ll_orig";
  theoryCSHisto1D_perMBW=1;
  TString covEE_H2NameOnFile="h2eeCov_orig";
  TString covMM_H2NameOnFile="h2mmCov_orig";
  TString h1combiOrig_FileName=inpFName;
  TString h1combiOrig_HNameOnFile="h1ll_orig";
  TString covCombiOrig_H2NameOnFile="h2llCov_orig";
  int scaleCov=0;
  int correctEEcov=0;
  TString specSaveTag="";

  if (1) {
    // this file contains Alexeys inputs
    h1combiOrig_FileName="dy8TeV-check-20170810.root";
    h1combiOrig_HNameOnFile="h1combiFlat";
    covCombiOrig_H2NameOnFile="h2covHepData";
  }

  if (0) {
    inpFName="dy8TeV-check-20170810.root";
    eeCSFName=inpFName;
    mmCSFName=inpFName;
    theoryCSFName=inpFName;
    eeCSHisto1DName="h1eeFlat";
    mmCSHisto1DName="h1mmFlat";
    theoryCSHisto1DName="h1combiFlat";
    theoryCSHisto1D_perMBW=1;
    h1combiOrig_FileName=inpFName;
    h1combiOrig_HNameOnFile="h1combiFlat";
    covEE_H2NameOnFile="h2covEE_materials";
    covMM_H2NameOnFile="h2covMM_materials";
    covCombiOrig_H2NameOnFile="h2covHepData";
    scaleCov=1;
    correctEEcov=0; // is not right
    specSaveTag="-materials";
  }
  else if (1) {
    if (0) {
      inpFName="resultCombiner2D_bug.root";
      specSaveTag="-bugVersion";
    }
    else {
      inpFName="resultCombiner2D_aj.root";
      specSaveTag="-ajVersion";
    }
    eeCSFName=inpFName;
    mmCSFName=inpFName;
    theoryCSFName=inpFName;
    eeCSHisto1DName="h1ee_blueInp";
    mmCSHisto1DName="h1mm_blueInp";
    theoryCSHisto1DName="h1ll_blueRes";
    theoryCSHisto1D_perMBW=1;
    covEE_H2NameOnFile="h2eeCov_blueInp";
    covMM_H2NameOnFile="h2mmCov_blueInp";
    if (0) {
      h1combiOrig_FileName=inpFName;
      h1combiOrig_HNameOnFile="h1ll_blueRes";
      covCombiOrig_H2NameOnFile="h2llCov_blueRes";
    }
    scaleCov=0;
    correctEEcov=0;
  }
  else if (0) {
    inpFName="resultCombiner2D_aj.root";
    specSaveTag="-ajVersion";
    eeCSFName=inpFName;
    mmCSFName=inpFName;
    theoryCSFName=inpFName;
    eeCSHisto1DName="h1ee_blueInp";
    mmCSHisto1DName="h1mm_blueInp";
    theoryCSHisto1DName="h1ll_blueRes";
    theoryCSHisto1D_perMBW=1;
    covEE_H2NameOnFile="h2eeCov_blueInp";
    covMM_H2NameOnFile="h2mmCov_blueInp";
    if (0) {
      eeCSHisto1DName=mmCSHisto1DName;
      covEE_H2NameOnFile=covMM_H2NameOnFile;
      specSaveTag+="-mmOnly";
    }
    else {
      mmCSHisto1DName=eeCSHisto1DName;
      covMM_H2NameOnFile=covEE_H2NameOnFile;
      specSaveTag+="-eeOnly";
    }
    if (0) {
      h1combiOrig_FileName=inpFName;
      h1combiOrig_HNameOnFile="h1ll_blueRes";
      covCombiOrig_H2NameOnFile="h2llCov_blueRes";
    }
    scaleCov=0;
    correctEEcov=0;
  }

  std::string showCanvs;
  showCanvs+=" plotChannelInputCovOnFile";
 // whether relative uncertainty in Acc in the channels is similar
  //showCanvs+=" plotChannelRelSystAccUnc";
  // compare randomized vs provided uncertainties
  //showCanvs+=" plotAdjUnc";
  showCanvs+=" plotAdjCov";
  showCanvs+=" plotInputContributedCovs"; // 4 canvases in each channel
  showCanvs+= " plotFinalCovByType"; // 4 canvases in the combined channel

  std::string showCombinationCanvases;
  //showCombinationCanvases="ALL";
  //
  //showCombinationCanvases+= " showCombiCS";
  showCombinationCanvases+= " showCombiCSInRanges2D132";
  showCombinationCanvases+= " showCombiCSCheck";
  showCombinationCanvases+= " showCombiCSUnc";

  TString fileTag, plotTag;

  std::vector<TH2D*> eeCovH2D, mmCovH2D;
  std::vector<TMatrixD> eeCovV, mmCovV;
  TMatrixD eeCovStat(132,132);
  TMatrixD eeCovSyst(eeCovStat);
  TMatrixD mmCovStat(eeCovStat);
  TMatrixD mmCovSyst(eeCovStat);

  TH1D *h1csLL_orig=NULL;
  TH2D *h2covLL_orig=NULL;
  
  TFile fin(inpFName);
  TH1D* h1csEE_tmp=loadHisto(fin, eeCSHisto1DName, "h1csEE_tmp", 1,h1dummy);
  TH1D* h1csMM_tmp=loadHisto(fin, mmCSHisto1DName, "h1csMM_tmp", 1,h1dummy);
  TH2D *h2covEE= loadHisto(fin,covEE_H2NameOnFile,"h2eeCov_orig",1,h2dummy);
  TH2D *h2covMM= loadHisto(fin,covMM_H2NameOnFile,"h2mmCov_orig",1,h2dummy);

  if (0) {
    //const int nPairs=3;
    //const int fix[nPairs][2] = { {94,118}, {95,119}, {96,120} };
    const int nPairs=1;
    const int fix[nPairs][2] = { { 96,120} };

    for (int ip=0; ip<nPairs; ip++) {
      int ibin=fix[ip][0];
      int jbin=fix[ip][1];
      double v= h2covMM->GetBinContent(ibin,jbin);
      double newV= 0.5*v;
      h2covMM->SetBinContent(ibin,jbin, newV);
      h2covMM->SetBinContent(jbin,ibin, newV);
    }
    specSaveTag+="-mdfCovMM";

    for (int ibin=1; ibin<=h2covMM->GetNbinsX(); ibin++) {
      for (int jbin=ibin; jbin<=h2covMM->GetNbinsY(); jbin++) {
	if ((h2covMM->GetBinContent(ibin,jbin)<0) &&
	    (fabs(h2covMM->GetBinContent(ibin,jbin))>1e-3)) {
	  std::cout << "ibin,jbin=" << ibin << "," << jbin << ", "
		    << h2covMM->GetBinContent(ibin,jbin) << "\n";
	}
      }
    }
    std::cout << "specSaveTag="<< specSaveTag << "\n";
    //return;
  }

  if (inpFName==h1combiOrig_FileName) {
    h1csLL_orig=loadHisto(fin,h1combiOrig_HNameOnFile,"h1csLL_orig",1,h1dummy);
    h2covLL_orig=loadHisto(fin,covCombiOrig_H2NameOnFile,"h2llCov_orig",1,h2dummy);
  }
  fin.Close();

  if (inpFName!=h1combiOrig_FileName) {
    TFile fin2(h1combiOrig_FileName);
    h1csLL_orig=loadHisto(fin2,h1combiOrig_HNameOnFile,"h1csLL_orig",1,h1dummy);
    h2covLL_orig=loadHisto(fin2,covCombiOrig_H2NameOnFile,"h2llCov_orig",1,h2dummy);
    fin2.Close();
  }

  if (!h1csEE_tmp || !h1csMM_tmp || !h1csLL_orig ||
      !h2covEE || !h2covMM || !h2covLL_orig) return;

  if (correctEEcov) {
    for (int ibin=132-12; ibin<=132; ibin++) {
      for (int jbin=1; jbin<=132; jbin++) {
	double factor=2.;
	if (jbin>120) factor*=2;
	h2covEE->SetBinContent(ibin,jbin,
			       h2covEE->GetBinContent(ibin,jbin)/factor);
      }
    }
  }

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

  PlotCovCorrOpt_t ccOpt;
  ccOpt.setLogScale(logX,logX);
  plotCovCorr(h2covEE,"cCovEE",ccOpt);

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
  if (hasValue("plotChannelInputCovOnFile",showCanvs)) {
    double ymin=0, ymax=0;
    ymin=1e-9; ymax=20.;
    mmCovS.Plot("mmCov_onFile",h1csMM_tmp,ccOpt);
    mmCovS.PlotUnc("mmCovUnc_onFile",h1csMM_tmp,NULL,ymin,ymax,ccOpt.logScaleX);
    eeCovS.Plot("eeCov_onFile",h1csEE_tmp,ccOpt);
    eeCovS.PlotUnc("eeCovUnc_onFile",h1csEE_tmp,NULL,ymin,ymax,ccOpt.logScaleX);
    //return;
  }

  // check whether mmCovS.Syst_Acc = eeCovS.Syst_Acc
  if (hasValue("plotChannelRelSystAccUnc",showCanvs)) {
    TH1D *h1mmAcc_rel= uncFromCov(mmCovS.covSyst_Acc,"h1mmAcc_rel",
				  h1csMM_tmp,h1csMM_tmp,0);
    TH1D *h1eeAcc_rel= uncFromCov(eeCovS.covSyst_Acc,"h1eeAcc_rel",
				  h1csEE_tmp,h1csEE_tmp,0);
    hsRed.SetStyle(h1eeAcc_rel);
    plotHisto(h1mmAcc_rel,"cAccUnc",logX,1,"LP","#mu#mu");
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
    plotCovCorr(eeCovTot,h1csEE_tmp,"h2eeCov_tmp","cCovEEtot_chk",ccOpt);
    plotCovCorr(mmCovTot,h1csMM_tmp,"h2mmCov_tmp","cCovMMtot_chk",ccOpt);
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
		showCombinationCanvases,use_scale,&ccOpt);
  if (!blue) { std::cout << "failed to get BLUE result\n"; return; }

  //return;

  if (includeLumiUnc) {
    TMatrixD totFinalCov(*blue->getCov());
    TH1D *h1csLL_test= convert2histo1D(*blue->getEst(),totFinalCov,
				       "h1csLL_test","h1csLL test",NULL);
    if (!addLumiCov(totFinalCov,-lumiUnc,h1csLL_test)) return;
    plotCovCorr(totFinalCov,NULL,"h2finCovNoLumi","cFinCovNoLumi",ccOpt);
  }

  //return;

  // make analysing plots
  CovStruct_t *llCovS_ptr=detailedCovPlots(blue,corrCase,
					   h1csEE_tmp,h1csMM_tmp,
					   eeCovS,mmCovS,showCanvs,&ccOpt);
  if (!llCovS_ptr) {
    std::cout << "detailedCovPlots did not return llCovS_ptr\n";
    return;
  }

  if (1) {
    std::cout << "\n\nmanual Z-range cross section calculation\n";
    int idxMin= 3*24; // skip 20-30,30-45,45-60 GeV ranges
    int idxMax_p1= idxMin+24+1;
    TH1D *h1rapBinW= new TH1D("h1rapBinW","h1rapBinW",
			      _nFlatBins2D,_flatBins2D_rapBinW);
    TH1D *h1ee_rapBinW= cloneHisto(h1rapBinW,"h1ee_rapBinW","h1ee_rapBinW");
    TH1D *h1mm_rapBinW= cloneHisto(h1rapBinW,"h1mm_rapBinW","h1mm_rapBinW");
    TH1D *h1ll_rapBinW= cloneHisto(h1rapBinW,"h1ll_rapBinW","h1ll_rapBinW");
    copyContents(h1ee_rapBinW,h1csEE_tmp);
    copyContents(h1mm_rapBinW,h1csMM_tmp);
    copyContents(h1ll_rapBinW,h1csLL_orig);
    //printHisto(perMassBinWidth(h1ee_rapBinW,1),0); // check binning
    eeCovS.PrintZRangeUnc("ee",idxMin,idxMax_p1,h1ee_rapBinW,0);
    mmCovS.PrintZRangeUnc("mm",idxMin,idxMax_p1,h1mm_rapBinW,0);
    llCovS_ptr->PrintZRangeUnc("ll",idxMin,idxMax_p1,h1ll_rapBinW,0);
    std::cout << "\n\n";
    //return;
  }


  // prepare combined result for plotting
  TString massStr= "idx"; //niceMassAxisLabel(2,"");
  TString sigmaStr= niceMassAxisLabel(2,"",1);
  TString axisStr= ";" + massStr + ";" + sigmaStr;
  TH1D *h1comb_loc= convert2histo1D(*blue->getEst(),*blue->getCov(),
				    "h1Combi_loc",
				    "combined measurement " + axisStr,
				    h1csEE_tmp,0);
  hsBlack.SetStyle(h1comb_loc);
  logAxis(h1comb_loc);
  TH1D *h1combUnc= errorAsCentral(h1comb_loc);

  // compare to theory uncertainties
  if (1) {
    TH1D *h1theory_onFile= loadHisto(theoryCSFName,theoryCSHisto1DName,
				     "h1cs_theory_loc_onFile",1,h1dummy);
    if (!h1theory_onFile) return;
    TH1D *h1theory= (theoryCSHisto1D_perMBW) ?
      h1theory_onFile : perMassBinWidth(h1theory_onFile);
    HistoStyle_t(kRed,7,1).SetStyle(h1theory);
    TH1D *h1theoryErr= errorAsCentral(h1theory);

    TString theoryIsOrig=1;
    TString thLabel= (theoryIsOrig) ? "orig.BLUE code" : "theory";
    plotHisto(h1comb_loc, "cTheoryVsCombi",0,1,"LPE1","BLUE result");
    plotHistoAuto(h1theory, "cTheoryVsCombi",0,1,"LPE1",thLabel);
    plotHisto(h1combUnc, "cTheoryVsCombiErr",0,1,"LP","BLUE result unc.");
    plotHistoAuto(h1theoryErr,"cTheoryVsCombiErr",0,1,"LP",thLabel+" unc.");
    printRatio(h1theory,h1comb_loc);
    printRatio(h1theoryErr,h1combUnc);
  }

  // check the values and the uncertainties
  if (1) {
    h1csLL_orig->SetLineColor(kRed);
    h1csLL_orig->SetMarkerColor(kRed);
    h1csLL_orig->SetStats(0);

    TH1D *h1comb_diff= cloneHisto(h1csLL_orig,"h1comb_diff","h1comb_diff");
    h1comb_diff->Add(h1comb_loc,-1);
    removeError(h1comb_diff);
    printHisto(h1comb_diff,0);
    TH1D *h1combUnc_orig= errorAsCentral(h1csLL_orig);
    TH1D *h1combUnc_diff= cloneHisto(h1combUnc_orig,
				     "h1combUnc_diff","h1combUnc_diff");
    h1combUnc_diff->Add(h1combUnc,-1);
    printHisto(h1combUnc_diff,0);
 
    TH1D *h1comb_ratio= cloneHisto(h1csLL_orig, "h1comb_ratio",
				   "h1comb_ratio;"+massStr+";orig/new");
    h1comb_ratio->Divide(h1csLL_orig,h1comb_loc);
    removeError(h1comb_ratio);
    plotHisto(h1comb_ratio,"cCombRatio",logX,0,"LP");

    plotHisto(h1combUnc,"cCombUnc",logX,1,"LP","new comb");
    plotHistoAuto(h1combUnc_orig,"cCombUnc",logX,1,"LP","orig");

    TH1D *h1combUnc_ratio= cloneHisto(h1combUnc_orig,"h1combUnc_ratio",
			      "h1combUnc_ratio;"+massStr+";origUnc/newUnc");
    h1combUnc_ratio->Divide(h1combUnc_orig,h1combUnc);
    removeError(h1combUnc_ratio);
    plotHisto(h1combUnc_ratio,"cCombUncRatio",logX,0,"LP");

    //printHisto(h1comb_ratio);
    printRatio(h1csLL_orig,h1comb_loc);
    printRatio(h1combUnc_orig,h1combUnc);

    if (use_scale_user>0)
      std::cout << "use_scale_user=" << use_scale_user << "\n";
    //return;
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

  if (specSave) {
    TString canvList;
    canvList+=" cErr cTheoryVsCombi cTheoryVsCombiErr";
    canvList+=" cCombRatio cCombUnc cCombUncRatio";
    canvList+=" c2CovStatYield_mmCov_onFile c2CovStatYield_eeCov_onFile";
    canvList+=" cCombiRange_0_24 cCombiRange_24_48 cCombiRange_48_72 ";
    canvList+=" cCombiRange_72_96 cCombiRange_96_120 cCombiRange_120_132";
    SaveCanvases(canvList,"dir-8TeV-spec"+specSaveTag);
  }

}

// -------------------------------------------------------
