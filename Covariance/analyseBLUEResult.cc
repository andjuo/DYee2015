#include "analyseBLUEResult.h"

// -------------------------------------------------------
// -----------------------------------------------------------------------

TString eeCSFName="cs_DYee_13TeV_El3.root";
TString mmCSFName="cs_DYmm_13TeVMuApproved_cs.root";
TString theoryCSFName="theory13TeVmm.root";
TString eeCSHisto1DName="h1PreFSRCS";
TString mmCSHisto1DName="h1CS";
TString theoryCSHisto1DName="h1cs_theory";
int theoryCSHisto1D_perMBW=0;

// -------------------------------------------------------
// -----------------------------------------------------------------------

int loadUncData(TString lepTag,
		const finfomap_t &fNames,
		const finfomap_t &h1names,
		const fvarweightmap_t &relUncW,
		const TH1D *h1binning,
		const TH1D *h1centralVal, // if the uncertainty is relative
		std::vector<TH1D*> &h1uncV,
		TString tag)
{
  //std::cout << "Entered loadUncData\n";
  //if (h1centralVal) printHisto(h1centralVal);


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

    fvarweightmap_t::const_iterator itRelW= relUncW.find(it->first);
    if (itRelW==relUncW.end()) {
      std::cout << "coult not find the relUncW for "
		<< variedVarName(it->first) << "\n";
      return 0;
    }

    TH1D *h1= loadHisto(fname,itHName->second,itHName->second + "_tmp",h1dummy);
    if (!h1) return 0;
    std::cout << "got " << h1newName << "\n";
    //printHisto(h1);
    //std::cout << "\n\n" << std::endl;
    if (itRelW->second!=0.) {
      if (!copyContents(h1unc, h1)) return 0;
    }
    if ((h1centralVal!=NULL) && (itRelW->second!=0.)) {
      std::cout << "scaling uncertainty for " << variedVarName(it->first) << "\n";
      h1unc->Multiply(h1centralVal);
      h1unc->Scale(itRelW->second);
    }
  }
  return 1;
}

// -------------------------------------------------------
// -------------------------------------------------------

int loadCovData(TString lepTag,
		const finfomap_t &covFNames,
		std::vector<TH2D*> &covAsH2D, std::vector<TMatrixD> &covV,
		TMatrixD &eeCovStat, TMatrixD &eeCovSyst)
{
  for (finfomap_t::const_iterator it= covFNames.begin();
       it!=covFNames.end(); it++) {
    TString fname= it->second;
    TString h2newName = "h2Cov_" + lepTag + "_" + variedVarName(it->first);
    TH2D *h2= loadHisto(fname,"h2Cov",h2newName,h2dummy);
    if (!h2) return 0;
    std::cout << "got " << h2newName << "\n";
    covAsH2D.push_back(h2);
    covV.push_back(convert2mat(h2));
    if (it->first==_varYield) eeCovStat+= covV.back();
    else eeCovSyst+= covV.back();
  }
  return 1;
}

// -------------------------------------------------------

int addLumiCov(TMatrixD &cov, double lumiRelUnc, const TH1D *h1cs)
{
  if (lumiRelUnc==0.) {
    std::cout << "addLumiCov(cov,lumiRelUnc,h1cs): lumiUnc=0., no effect\n";
    return 1;
  }

  if (!h1cs) {
    std::cout << "addLumiCov(cov,lumiRelUnc,h1cs): null histo\n";
    return 0;
  }

  double sign= (lumiRelUnc>0) ? 1 : -1;

  for (int ir=0; ir<cov.GetNrows(); ir++) {
    for (int ic=0; ic<cov.GetNcols(); ic++) {
      cov(ir,ic) +=
	h1cs->GetBinContent(ir+1) * h1cs->GetBinContent(ic+1)
	* pow( lumiRelUnc, 2) * sign;
    }
  }

  return 1;
}

// -------------------------------------------------------

int addLumiCov(TMatrixD &eeCov, TMatrixD &mmCov, TMatrixD &emCov,
	       double lumiRelUnc, const TH1D *h1csEE, const TH1D *h1csMM)
{
  if (lumiRelUnc==0.) {
    std::cout << "addLumiCov: lumiUnc=0., no effect\n";
    return 1;
  }
  if (!h1csEE || !h1csMM) {
    std::cout << "addLumiCov: null cs histo\n";
    return 0;
  }

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

  for (int ir=0; ir<eeCov.GetNrows(); ir++) {
    for (int ic=0; ic<eeCov.GetNcols(); ic++) {
      emCov(ir,ic) +=
	h1csEE->GetBinContent(ir+1) * h1csMM->GetBinContent(ic+1)
	* pow( lumiRelUnc, 2) * sign;
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

int compareUncAndScale(const std::vector<TString> &labelV,
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
      h1->GetYaxis()->SetTitleOffset(1.5);
      logAxis(h1,1+8,"M [GeV]","","");
      TCanvas *cx=plotHisto(h1,canvName,1,1,"LP",labelV[i] + " target");
      setLeftMargin(cx,0.15);
      moveLegend(cx,0.05,0.);
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

BLUEResult_t* combineData(const TMatrixD *covEE_inp,
			  const TMatrixD *covMM_inp,
			  const TMatrixD *covEM_inp,
			  TString outputFileTag, TString plotTag,
			  int printCanvases, std::string showCanvases,
			  double scale,
			  const PlotCovCorrOpt_t *ccOpt_user)
{

  if (showCanvases=="ALL") {
    showCanvases="";
    showCanvases+= " showCSCheck"; // plot 1D xs in channels vs theory
    //showCanvases+= " showCSCheckInRanges"; // plot 1D xs in channels vs theory
    //showCanvases+= " showCombiCS"; // plot 1D xs in the combined channel vs inp
    // plot 1D xs in the combined channel vs theory
    //showCanvases+= " showCombiCSCheckInputInRanges";
    // plot 1D xs combi vs theory (abs values and the ratio plot)
    //showCanvases+= " showCombiCSCheck";
    // plot 1D xs uncertainties combi vs channels
    //showCanvases+= " showCombiCSUnc";
    // plot input covariances
    showCanvases+= " showInputCov";
    // plot final covariance
    //showCanvases+= " showFinalCov";
    // plot contributed covariances from the channels
    //showCanvases+= " showContributedCov";
  }

  TString plotTagSafe=plotTag;
  plotTagSafe.ReplaceAll(" ","_");
  plotTagSafe.ReplaceAll("(","_");
  plotTagSafe.ReplaceAll(")","_");
  plotTagSafe.ReplaceAll(".","");
  plotTagSafe.ReplaceAll(";","_");

  //TString eeCSFName="cs_DYee_13TeV_El3.root";
  //TString mmCSFName="cs_DYmm_13TeVMuApproved_cs.root";

  TH1D* h1csEE=loadHisto(eeCSFName, eeCSHisto1DName, "h1csEE", 1,h1dummy);
  TH1D* h1csMM=loadHisto(mmCSFName, mmCSHisto1DName, "h1csMM", 1,h1dummy);
  TH1D* h1csTheory_onFile=loadHisto(theoryCSFName,
				    theoryCSHisto1DName, "h1cs_theory",
				    1,h1dummy);
  TH1D *h1csTheory= (theoryCSHisto1D_perMBW) ?
    h1csTheory_onFile : perMassBinWidth(h1csTheory_onFile);
  TH1D *h1csTheoryRed=(TH1D*)h1csTheory->Clone("h1csTheoryRed");

  TString massStr= niceMassAxisLabel(2,"");
  TString eemassStr= niceMassAxisLabel(0,"");
  TString mmmassStr= niceMassAxisLabel(1,"");
  TString sigmaStr= niceMassAxisLabel(2,"",1);
  TString axisStr= ";" + massStr + ";" + sigmaStr;
  TString DYeeStr="\\text{DY}\\!\\rightarrow\\!ee";
  TString DYmmStr="\\text{DY}\\!\\rightarrow\\!\\mu\\mu";
  TString DYllStr="\\text{DY}\\!\\rightarrow\\!\\ell\\!\\ell";

  HistoStyle_t hsEE(kGreen+1, 24, 1, 0.8);
  HistoStyle_t hsMM(46, 20, 1, 0.8);
  HistoStyle_t hsTheory(9, 1);
  HistoStyle_t hsTheoryRed(kRed, 1, 1, 0.8);
  HistoStyle_t hsCombi(kBlue, 5, 1, 1.);
  hsEE.SetStyle(h1csEE);
  hsMM.SetStyle(h1csMM);
  hsTheory.SetStyle(h1csTheory);
  hsTheoryRed.SetStyle(h1csTheoryRed);
  h1csTheory->SetStats(0);
  h1csTheoryRed->SetStats(0);

  h1csEE->GetXaxis()->SetTitle(eemassStr);
  h1csMM->GetXaxis()->SetTitle(mmmassStr);

  TMatrixD measEE = convert2mat1D(h1csEE);
  TMatrixD covEE  = (covEE_inp) ? (*covEE_inp) : convert2cov1D(h1csEE);
  TMatrixD measMM = convert2mat1D(h1csMM);
  TMatrixD covMM  = (covMM_inp) ? (*covMM_inp) : convert2cov1D(h1csMM);
  TMatrixD covEM = (covEM_inp) ? (*covEM_inp) : TMatrixD(TMatrixD::kZero,covEE);
  if (covEE_inp) { setError(h1csEE,*covEE_inp); }
  if (covMM_inp) { setError(h1csMM,*covMM_inp); }

  if (0) {
    std::cout << "\nMatrix sizes:  ";
    printMDim("measEE",measEE); printMDim(", covEE",covEE);
    printMDim(", measMM",measMM); printMDim(", covMM",covMM);
    printMDim(", covEM",covEM);
    std::cout << "\n";
  }


  //printHisto(h1csEE);
  //printHisto(h1csMM);

  PlotCovCorrOpt_t optCC(1,1,1,1.8,0.15,0.15);
  if (ccOpt_user) optCC.assign(*ccOpt_user);

  if (hasValue("showCSCheck",showCanvases)) {
    TH1D *h1frame= cloneHisto(h1csEE,"h1frame_inpCS","frame");
    h1frame->Reset();
    logAxis(h1frame,1+0*2,massStr,sigmaStr);
    h1frame->SetTitle("input cross section");
    h1frame->GetYaxis()->SetTitleOffset(1.5);
    h1frame->GetYaxis()->SetRangeUser(1e-8,1e3);
    TCanvas *cx=plotHisto(h1frame,"cCScheck",optCC.logScaleX,1, "");
    plotHistoSame(h1csEE,"cCScheck","LPE1", DYeeStr);
    setLeftMargin(cx,0.15); moveLegend(cx,0.05,0.);
    plotHistoSame(h1csMM,"cCScheck", "LPE1", DYmmStr);
    plotHistoSame(h1csTheory,"cCScheck","LP", "theory");
    //return NULL;
  }

  TGraphErrors *grCSEE= createGraph(h1csEE,"grEE", -1);
  TGraphErrors *grCSMM= createGraph(h1csMM,"grMM",  1);
  hsEE.SetStyle(grCSEE);
  hsMM.SetStyle(grCSMM);

  if (hasValue("showCSCheckInRanges",showCanvases)) {
    std::vector<TCanvas*> canvasCSV;
    for (int i=0; i<5; i++) {
      TString canvName="";
      TH2D *h2frame=NULL;
      TCanvas *cx= createMassFrame(i,"cCSCheckRange_",
				   "input cross section" + axisStr,
				   _massFrameCS,&canvName,&h2frame);
      canvasCSV.push_back(cx);
      //plotHistoSame(h1csEE,canvName,"LPE1", DYeeStr);
      //plotHistoSame(h1csMM,canvName,"LPE1", DYmmStr);
      plotGraphSame(grCSEE,canvName,"PE1", DYeeStr);
      plotGraphSame(grCSMM,canvName,"PE1", DYmmStr);
      plotHistoSame(h1csTheory,canvName,"LP", "theory");
      if (i==1) moveLegend(cx,0.,0.55);
      //if (i==4) moveLegend(cx,0.05,0.);
    }

    if (printCanvases) {
      TCanvas *cCSCheck= findCanvas("cCScheck");
      if (cCSCheck) canvasCSV.push_back(cCSCheck);

      TFile fout("foutCanvas_DYCSInp" + outputFileTag + ".root","RECREATE");
      SaveCanvases(canvasCSV,"dir-plot-DYCSInp" + outputFileTag,&fout);
      writeTimeTag();
      fout.Close();
      std::cout << "file <" << fout.GetName() << "> created\n";
    }
    //return NULL;
  }

  BLUEResult_t *blue= new BLUEResult_t();
  if (!blue->estimate(measEE,covEE,measMM,covMM,covEM,scale)) {
    std::cout << "combineData failed\n";
    return NULL;
  }

  TH1D *h1comb= convert2histo1D(*blue->getEst(),*blue->getCov(),
				"h1Combi",
				"combined measurement " + plotTag + axisStr,
				h1csEE,0);
  printHisto(h1comb);
  hsCombi.SetStyle(h1comb);
  logAxis(h1comb);

  // plot combined cross section
  std::vector<TCanvas*> canvasCombiCSV;
  if (1) {
    if (hasValue("showCombiCS",showCanvases)) {
      TH1D *h1frame= cloneHisto(h1csEE,"h1frame_combCS","frame");
      h1frame->Reset();
      logAxis(h1frame,1+0*2,massStr,sigmaStr);
      h1frame->SetTitle("cross section");
      h1frame->GetYaxis()->SetTitleOffset(1.5);
      h1frame->GetYaxis()->SetRangeUser(1e-8,1e3);
      TCanvas *cx=plotHisto(h1frame,"cCombiCS",optCC.logScaleX,1, "");
      plotHistoSame(h1csEE,"cCombiCS","LPE1", DYeeStr);
      setLeftMargin(cx,0.15); moveLegend(cx,0.05,0.);
      plotHistoSame(h1csMM,"cCombiCS", "LPE1", DYmmStr);
      plotHistoSame(h1comb,"cCombiCS", "LPE1", DYllStr);
      //return NULL;
    }

    if (hasValue("showCombiCSCheckInputInRanges",showCanvases)) {
      TH1D *h1theory_red= cloneHisto(h1csTheory,
				     "h1theory_red",h1csTheory->GetTitle());
      h1theory_red->SetLineColor(kRed);
      h1theory_red->SetMarkerColor(kRed);

      for (int i=0; i<5; i++) {
	TString canvName="";
	TH2D *h2frame=NULL;
	TCanvas *cx= createMassFrame(i,"cCombiRange_",
				     "cross section " + plotTag + axisStr,
				     _massFrameCS,&canvName,&h2frame);
	canvasCombiCSV.push_back(cx);
	plotGraphSame(grCSEE,canvName,"PE1", DYeeStr);
	plotGraphSame(grCSMM,canvName,"PE1", DYmmStr);
	plotHistoSame(h1theory_red,canvName,"LP", "theory");
	plotHistoSame(h1comb,canvName,"LPE1", DYllStr);
	if (i==1) moveLegend(cx,0.,0.55);
	//if (i==4) moveLegend(cx,0.05,0.);
      }
    }

    // compare combination to theory
    if (hasValue("showCombiCSCheck",showCanvases)) {
      h1csTheoryRed->GetYaxis()->SetRangeUser(1e-7,1e3);
      h1csTheoryRed->GetYaxis()->SetTitleOffset(1.4);
      plotHisto(h1csTheoryRed,"cCombiCS2Theory",optCC.logScaleX,1,"LPE","theory");
      plotHistoSame(h1comb,"cCombiCS2Theory","LPE1","combined " + plotTag);

      TH1D *h1ratio=(TH1D*)h1comb->Clone("h1ratio");
      h1ratio->Divide(h1csTheory);
      h1ratio->GetXaxis()->SetTitle(massStr);
      h1ratio->GetYaxis()->SetRangeUser(-0.5,1.5);
      h1ratio->GetYaxis()->SetTitle("combined/theory");
      plotRatio(h1ratio,"cCombiCS2TheoryRatio",optCC.logScaleX,0,"LPE1");
    }

    // compare the uncertainties
    if (hasValue("showCombiCSUnc",showCanvases)) {
      TH1D *h1EEErr= errorAsCentral(h1csEE,0);
      logAxis(h1EEErr,1,massStr,"\\delta"+sigmaStr);
      h1EEErr->GetYaxis()->SetRangeUser(1e-7,20);
      TH1D *h1MMErr= errorAsCentral(h1csMM,0);
      TH1D *h1LLErr= errorAsCentral(h1comb,0);

      TCanvas *cErr=plotHisto(h1EEErr,"cErr",optCC.logScaleX,1,"LPE1",DYeeStr);
      setLeftMargin(cErr,0.15); moveLegend(cErr,0.05,0.);
      plotHistoSame(h1MMErr,"cErr","LPE1",DYmmStr);
      plotHistoSame(h1LLErr,"cErr","LPE1",DYllStr);
    }

    if (printCanvases) {
      TFile fout("foutCanvas_DYCSCombi" + outputFileTag + ".root","RECREATE");
      SaveCanvases("cCombiCS cCombiCS2Theory cCombiCS2TheoryRatio cErr",
		   "dir-plot-DYCSCombi" + outputFileTag, &fout);
      SaveCanvases(canvasCombiCSV,"dir-plot-DYCSCombi" + outputFileTag,&fout);
      h1comb->Write();
      writeTimeTag();
      fout.Close();
      std::cout << "file <" << fout.GetName() << "> created\n";
    }
  }

  // plot the covariances
  if (hasValue("showInputCov",showCanvases)) {
    if (1 && (covEE_inp || covMM_inp || covEM_inp)) {
      std::cout << "plot input covariances\n";
      plotCovCorr(covEE,h1csEE,"h2eeCov","cEECov",optCC);
      plotCovCorr(covMM,h1csMM,"h2mmCov","cMMCov",optCC);
      std::cout << "covEM_inp is " << ((covEM_inp==NULL) ? "" : "NON ")
		<< "null";
      if (covEM_inp) std::cout << "; " << covEM_inp->NonZeros() << " non zeros";
      std::cout << "\n";
      if (covEM_inp) { // && (covEM_inp->NonZeros()>0)) {
	TMatrixD inpCov(*blue->getInpCov());
	TH1D *h1binDef= new TH1D("h1binDef","h1binDef;binIdx;count",
				 inpCov.GetNrows(),1,inpCov.GetNrows()+1);
	PlotCovCorrOpt_t optCCNoLog(optCC);
	optCCNoLog.logScaleX=0;
	optCCNoLog.logScaleY=0;
	optCCNoLog.autoZRangeCorr=0;
	plotCovCorr(inpCov,h1binDef,"h2TotInpCov","cTotInpCov",optCCNoLog);
      }
    }


    TMatrixD covEEFinal( blue->contributedCov(blue->measCovMatrix(covEE,1)) );
    TMatrixD covMMFinal( blue->contributedCov(blue->measCovMatrix(covMM,0)) );
    TMatrixD covEMFinal( blue->contributedCov(blue->measCovMatrix(covEM,2)) );
    TMatrixD finalCov(*blue->getCov());

    //std::cout << "h1comb->GetXaxis()->GetTitle=" << h1comb->GetXaxis()->GetTitle() << "\n";
    //std::cout << "h1csEE title " << h1csEE->GetXaxis()->GetTitle() << "\n";
    //std::cout << "h1csMM title " << h1csMM->GetXaxis()->GetTitle() << "\n";

    if (hasValue("showFinalCov",showCanvases)) {
      plotCovCorr(finalCov,h1comb,"h2finCov","cFinCov",optCC);
    }

    if (hasValue("showContributedCov",showCanvases)) {
      TH2D* h2eeCorrPart= convert2histo(covToCorrPartial(covEEFinal,finalCov),
					h1csEE,"h2eeCorrPart","h2eeCorrPart");
      logAxis(h2eeCorrPart);
      TCanvas *cCorrPartEE= plotHisto(h2eeCorrPart,"cCorrPartEE",optCC);
      cCorrPartEE->SetGrid(1,1);

      TH2D* h2mmCorrPart= convert2histo(covToCorrPartial(covMMFinal,finalCov),
					h1csMM,"h2mmCorrPart","h2mmCorrPart");
      logAxis(h2mmCorrPart);
      TCanvas *cCorrPartMM= plotHisto(h2mmCorrPart,"cCorrPartMM",optCC);
      cCorrPartMM->SetGrid(1,1);
    }

    if (printCanvases) {
      TFile fout("foutCanvas_DYCSCov" + outputFileTag + ".root","RECREATE");
      int count=SaveCanvases("cEECov cMMCov cTotInpCov " +
			     TString("cFinCov cCorrPartEE cCorrPartMM"),
			     "dir-plot-DYCSCov" + outputFileTag,&fout);
      h1comb->Write();
      std::cout << "saved " << count << " canvases\n";
      writeTimeTag();
      fout.Close();
      std::cout << "file <" << fout.GetName() << "> created\n";
    }
  }

  return blue;
}

// -------------------------------------------------------
// -------------------------------------------------------
