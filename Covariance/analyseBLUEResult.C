#include "inputs.h"
#include "crossSection.h"
#include "Blue.h"

// -------------------------------------------------------

BLUEResult_t* combineData(const TMatrixD *covEE_inp,
			  const TMatrixD *covMM_inp,
			  const TMatrixD *covEM_inp,
			  TString outputFileTag, TString plotTag,
			  int printCanvases)
{
  double scale=1; //100000.; // 1000000

  TString plotTagSafe=plotTag;
  plotTagSafe.ReplaceAll(" ","_");
  plotTagSafe.ReplaceAll("(","_");
  plotTagSafe.ReplaceAll(")","_");
  plotTagSafe.ReplaceAll(".","");
  plotTagSafe.ReplaceAll(";","_");

  TString eeCSFName="cs_DYee_13TeV_El3.root";
  TString mmCSFName="cs_DYmm_13TeVMuApproved_cs.root";

  TH1D* h1csEE=loadHisto(eeCSFName, "h1PreFSRCS", "h1csEE", 1,h1dummy);
  TH1D* h1csMM=loadHisto(mmCSFName, "h1CS", "h1csMM", 1,h1dummy);
  TH1D* h1csTheory=perMassBinWidth(loadHisto("theory13TeVmm.root",
					     "h1cs_theory", "h1cs_theory",
					     1,h1dummy));
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

  printHisto(h1csEE);
  printHisto(h1csMM);

  PlotCovCorrOpt_t optCC(1,1,1,1.8,0.15,0.15);

  if (1) {
    TH1D *h1frame= cloneHisto(h1csEE,"h1frame_inpCS","frame");
    h1frame->Reset();
    logAxis(h1frame,1+0*2,massStr,sigmaStr);
    h1frame->SetTitle("input cross section");
    h1frame->GetYaxis()->SetTitleOffset(1.5);
    h1frame->GetYaxis()->SetRangeUser(1e-8,1e3);
    TCanvas *cx=plotHisto(h1frame,"cCScheck",1,1, "");
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

  if (1) {
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
  if (scale!=1.) {
    measEE *= scale;
    covEE *= scale*scale;
    measMM *= scale;
    covMM *= scale*scale;
    covEM *= scale*scale;
  }
  if (!blue->estimate(measEE,covEE,measMM,covMM,covEM)) {
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
    if (1) {
      TH1D *h1frame= cloneHisto(h1csEE,"h1frame_combCS","frame");
      h1frame->Reset();
      logAxis(h1frame,1+0*2,massStr,sigmaStr);
      h1frame->SetTitle("cross section");
      h1frame->GetYaxis()->SetTitleOffset(1.5);
      h1frame->GetYaxis()->SetRangeUser(1e-8,1e3);
      TCanvas *cx=plotHisto(h1frame,"cCombiCS",1,1, "");
      plotHistoSame(h1csEE,"cCombiCS","LPE1", DYeeStr);
      setLeftMargin(cx,0.15); moveLegend(cx,0.05,0.);
      plotHistoSame(h1csMM,"cCombiCS", "LPE1", DYmmStr);
      plotHistoSame(h1comb,"cCombiCS", "LPE1", DYllStr);
      //return NULL;
    }

    if (1) {
      for (int i=0; i<5; i++) {
	TString canvName="";
	TH2D *h2frame=NULL;
	TCanvas *cx= createMassFrame(i,"cCombiRange_",
				     "cross section " + plotTag + axisStr,
				     _massFrameCS,&canvName,&h2frame);
	canvasCombiCSV.push_back(cx);
	plotGraphSame(grCSEE,canvName,"PE1", DYeeStr);
	plotGraphSame(grCSMM,canvName,"PE1", DYmmStr);
	//plotHistoSame(h1csTheory,canvName,"LP", "theory");
	plotHistoSame(h1comb,canvName,"LPE1", DYllStr);
	if (i==1) moveLegend(cx,0.,0.55);
	//if (i==4) moveLegend(cx,0.05,0.);
      }
    }

    // compare combination to theory
    if (1) {
      h1csTheoryRed->GetYaxis()->SetRangeUser(1e-7,1e3);
      h1csTheoryRed->GetYaxis()->SetTitleOffset(1.4);
      plotHisto(h1csTheoryRed,"cCombiCS2Theory",1,1,"LPE","theory");
      plotHistoSame(h1comb,"cCombiCS2Theory","LPE1","combined " + plotTag);

      TH1D *h1ratio=(TH1D*)h1comb->Clone("h1ratio");
      h1ratio->Divide(h1csTheory);
      h1ratio->GetXaxis()->SetTitle(massStr);
      h1ratio->GetYaxis()->SetRangeUser(-0.5,1.5);
      h1ratio->GetYaxis()->SetTitle("combined/theory");
      plotRatio(h1ratio,"cCombiCS2TheoryRatio",1,0,"LPE1");
    }

    // compare the uncertainties
    if (1) {
      TH1D *h1EEErr= errorAsCentral(h1csEE,0);
      logAxis(h1EEErr,1,massStr,"\\delta"+sigmaStr);
      h1EEErr->GetYaxis()->SetRangeUser(1e-7,20);
      TH1D *h1MMErr= errorAsCentral(h1csMM,0);
      TH1D *h1LLErr= errorAsCentral(h1comb,0);

      TCanvas *cErr=plotHisto(h1EEErr,"cErr",1,1,"LPE1",DYeeStr);
      setLeftMargin(cErr,0.15); moveLegend(cErr,0.05,0.);
      plotHistoSame(h1MMErr,"cErr","LPE1",DYmmStr);
      plotHistoSame(h1LLErr,"cErr","LPE1",DYllStr);
    }

    if (printCanvases) {
      TFile fout("foutCanvas_DYCSCombi" + outputFileTag + ".root","RECREATE");
      SaveCanvases("cCombiCS cCombiCS2Theory cCombiCS2TheoryRatio cErr",
		   "dir-plot-DYCSCombi" + outputFileTag, &fout);
      SaveCanvases(canvasCombiCSV,"dir-plot-DYCSCombi" + outputFileTag,&fout);
      writeTimeTag();
      fout.Close();
      std::cout << "file <" << fout.GetName() << "> created\n";
    }
  }

  // plot the covariances
  if (1) {
    if (1 && (covEE_inp || covMM_inp || covEM_inp)) {
      std::cout << "plot input covariances\n";
      plotCovCorr(covEE,h1csEE,"h2eeCov","cEECov",optCC);
      plotCovCorr(covMM,h1csMM,"h2mmCov","cMMCov",optCC);
      std::cout << "covEM_inp is " << ((covEM_inp==NULL) ? "" : "NON ")
		<< "null";
      if (covEM_inp) std::cout << "; " << covEM_inp->NonZeros() << " non zeros";
      std::cout << "\n";
      if (covEM_inp && (covEM_inp->NonZeros()>0)) {
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

    plotCovCorr(finalCov,h1comb,"h2finCov","cFinCov",optCC);

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

    if (printCanvases) {
      TFile fout("foutCanvas_DYCSCov" + outputFileTag + ".root","RECREATE");
      int count=SaveCanvases("cEECov cMMCov cTotInpCov " +
			     TString("cFinCov cCorrPartEE cCorrPartMM"),
			     "dir-plot-DYCSCov" + outputFileTag,&fout);
      std::cout << "saved " << count << " canvases\n";
      writeTimeTag();
      fout.Close();
      std::cout << "file <" << fout.GetName() << "> created\n";
    }
  }

  return blue;
}
