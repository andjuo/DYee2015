#include "inputs.h"
#include "crossSection.h"
#include "Blue.h"

#include <map>
//#include <pair>
#include <sstream>
#include <fstream>
#include <TCanvas.h>
#include <TFile.h>

//typedef std::pair<TVaried_t,TString> finfopair_t;
typedef std::map<TVaried_t,TString> finfomap_t;

int loadCovData(const finfomap_t &covIdx, std::vector<TH2D*> &covAsH2D,
		std::vector<TMatrixD> &covV,
		TMatrixD &totCovStat, TMatrixD &totCovSyst);

int combineData(const TMatrixD &measEE, const TMatrixD &covEE,
		const TMatrixD &measMM, const TMatrixD &covMM,
		const TH1D *h1binning);

void checkSplit(const BLUEResult_t &blue,
		const TMatrixD &covEE, const TMatrixD &covMM);

// -----------------------------------------------------------------------
// -----------------------------------------------------------------------

void work13TeV()
{

  TString eeCSFName="cs_DYee_13TeV_El3.root";
  TString mmCSFName="cs_DYmm_13TeVMuApproved_cs.root";

  std::vector<int> eeCovIdx, mmCovIdx;
  addToVector(eeCovIdx,7, _varYield, _varBkg, _varDetRes, _varEff,
	      _varRhoFile, _varAcc, _varFSRRes);
  addToVector(mmCovIdx,7, _varYield, _varBkg, _varDetRes, _varEff,
	      _varRhoFile, _varAcc, _varFSRRes);

  TH1D* h1csEE=loadHisto(eeCSFName, "h1PreFSRCS", "h1csEE", 1,h1dummy);
  TH1D* h1csMM=loadHisto(mmCSFName, "h1CS", "h1csMM", 1,h1dummy);
  TH1D* h1csTheory=perMassBinWidth(loadHisto("theory13TeVmm.root",
					     "h1cs_theory", "h1cs_theory",
					     1,h1dummy));

  if (0) {
    histoStyle(h1csEE, kGreen+1, 24, 1);
    histoStyle(h1csMM, 46, 20, 1);
    plotHisto(h1csEE,"cCScheck",1,1, "LPE1", "DY#rightarrowee");
    plotHistoSame(h1csMM,"cCScheck", "LPE1", "DY#rightarrow#mu#mu");
    plotHistoSame(h1csTheory,"cCScheck","LP", "theory");
    return;
  }

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
  TMatrixD eeCovStat(h1csEE->GetNbinsX(),h1csEE->GetNbinsX());
  TMatrixD eeCovSyst(eeCovStat);
  TMatrixD mmCovStat(eeCovStat);
  TMatrixD mmCovSyst(eeCovStat);

  if (!loadCovData(eeCovFNames,eeCovH2D,eeCovV,eeCovStat,eeCovSyst) ||
      !loadCovData(mmCovFNames,mmCovH2D,mmCovV,mmCovStat,mmCovSyst)) {
    std::cout << "failed to load covariances\n";
    return;
  }

  TMatrixD measEE(convert2mat1D(h1csEE));
  TMatrixD measMM(convert2mat1D(h1csMM));
  TMatrixD eeCovTot(eeCovStat,TMatrixD::kPlus,eeCovSyst);
  TMatrixD mmCovTot(mmCovStat,TMatrixD::kPlus,mmCovSyst);

  if (0) {
    TMatrixD cov100 = reduceCorrelations(eeCovTot,-2.);
    TMatrixD cov066 = reduceCorrelations(cov100,0.666);
    TMatrixD cov000 = reduceCorrelations(cov100,1.0);
    plotCovCorr(eeCovTot,h1csEE,"h2covOrig","cCov orig");
    plotCovCorr(cov100,h1csEE,"h2cov100","cCov100corr fake tripled");
    plotCovCorr(cov066,h1csEE,"h2cov066","cCov066corr fake one third");
    plotCovCorr(cov000,h1csEE,"h2cov000","cCov000corr fake no corr");
    return;
  }

  if (!combineData(measEE,eeCovStat,measMM,mmCovStat, h1csEE)) {
    std::cout << "failed to combined data. stopping\n";
    return;
  }


  return;
}

// -----------------------------------------------------------------------
// -----------------------------------------------------------------------

int loadCovData(const finfomap_t &covFNames,
		std::vector<TH2D*> &covAsH2D, std::vector<TMatrixD> &covV,
		TMatrixD &eeCovStat, TMatrixD &eeCovSyst)
{
  int isEE= (covFNames.begin()->second.Index("mumu")==-1) ? 1 : 0;
  TString lepTag=(isEE) ? "ee" : "mm";

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


// -----------------------------------------------------------------------

int combineData(const TMatrixD &measEE_inp, const TMatrixD &covEE_inp,
		const TMatrixD &measMM_inp, const TMatrixD &covMM_inp,
		const TH1D* h1binning)
{
  const int nXbins= h1binning->GetNbinsX();
  const TArrayD *xb= h1binning->GetXaxis()->GetXbins();
  const Double_t *xBins= xb->GetArray();

  TMatrixD measEE(measEE_inp), covEE(covEE_inp);
  TMatrixD measMM(measMM_inp), covMM(covMM_inp);

  //if (0) { measEE=measMM; covEE=covMM; }
  //else if (0) { measMM=measEE; covMM=covEE; }

  if (1) { // uncorrelated case
    covEE= reduceCorrelations(covEE_inp,1.);
    covMM= reduceCorrelations(covMM_inp,1.);
  }

  TMatrixD covEM(TMatrixD::kZero,covEE);
  BLUEResult_t blue;
  if (!blue.estimate(measEE,covEE, measMM,covMM, covEM)) {
    std::cout << "combineData failed to estimate\n";
    return 0;
  }

  if (0) {
    TString fCovFName="fCov.dat";
    std::ofstream fCov(fCovFName.Data());
    for (int ir=0; ir< blue.getCov()->GetNrows(); ir++) {
      for (int ic=0; ic< blue.getCov()->GetNcols(); ic++) {
	fCov << Form("    %6d    %6d   %13.6e",ir+1,ic+1,blue.getCov(ir,ic))
	     << "\n";
      }
    }
    fCov.close();
    std::cout << "file <" << fCovFName << "> saved\n";
  }

  TString massStr= niceMassAxisLabel(2,"");
  TString eemassStr= niceMassAxisLabel(0,"");
  TString mmmassStr= niceMassAxisLabel(1,"");
  TString sigmaStr= niceMassAxisLabel(2,"",1);
  TString axisStr= ";" + massStr + ";" + sigmaStr;

  if (1) {
    TH1D *h1EEErr= createErrHisto("h1EEErr","EE unc",covEE,nXbins,xBins);
    TH1D *h1MMErr= createErrHisto("h1MMErr","MM unc",covMM,nXbins,xBins);
    TH1D *h1LLErr= createErrHisto("h1LLErr","Comb unc",*blue.getCov(),nXbins,xBins);
    histoStyle(h1EEErr,kBlue,24,1,0.7);
    histoStyle(h1LLErr,kRed,7,2,1.);
    histoStyle(h1MMErr,kBlack,20,1,0.7);
    h1EEErr->SetTitle("uncertainties;"+massStr+";d"+sigmaStr);
    h1EEErr->GetYaxis()->SetTitleOffset(1.5);
    int setLogX= (h1EEErr->GetNbinsX()<50) ? 1:0;
    logAxis(h1EEErr);
    h1EEErr->SetStats(0);
    TCanvas *cErr=plotHisto(h1EEErr,"cErr",setLogX,1,"LP","DY#rightarrowee");
    cErr->SetLeftMargin(0.15); gPad->Modified(); cErr->Update();
    moveLegend(cErr,0.05,0.);
    cErr->SetGrid(1,1);
    plotHistoSame(h1MMErr,"cErr","LP","DY#rightarrow#mu#mu");
    plotHistoSame(h1LLErr,"cErr","LP","combined");
    //return 1;
  }

  if (1) {
    TGraphErrors *grEE= createGraph("EE measurement", measEE,covEE,
				    nXbins,xBins,-1);
    TGraphErrors *grMM= createGraph("MM measurement", measMM,covMM,
				    nXbins,xBins, 1);
    TGraphErrors *grLL= createGraph("combined measurement",
				    *blue.getEst(),*blue.getCov(),
				    nXbins,xBins, 0);
    grLL->SetTitle("combined measurement" + axisStr);
    grLL->GetYaxis()->SetTitleOffset(1.3);
    grEE->SetLineColor(kBlue);
    grEE->SetMarkerColor(kBlue);
    grLL->SetLineColor(kRed);
    grLL->SetMarkerColor(kRed);

    graphStyle(grEE,kBlue, 20, 1,0.6);
    graphStyle(grMM,kBlack,24, 1,0.6);
    graphStyle(grLL,kRed, 5, 1, 1.);
    //grEE->SetMarkerStyle(7);
    //grMM->SetMarkerStyle(24);

    //grLL->GetXaxis()->SetRangeUser(10,50);
    //grLL->GetYaxis()->SetRangeUser(20,250);

    TCanvas *cMeas= new TCanvas("cMeas","cMeas", 600,600);
    cMeas->SetGrid(1,1);
    int setLogXY= (grLL->GetN()<50) ? 1:0;
    if (setLogXY) {
      cMeas->SetLogx();
      cMeas->SetLogy();
    }
    logAxis(grLL,1+0*2);
    grLL->Draw("ALPE1");
    grEE->Draw("PE1 same");
    grMM->Draw("PE1 same");
    cMeas->Update();

    //grEE->Print("range");
    //grMM->Print("range");
    //grLL->Print("range");

    if (1) {
      TGraphErrors *grLL_line= new TGraphErrors(*grLL);
      grLL_line->SetLineStyle(3);

      const int nDivs=5;
      const double divX[nDivs][2] = { { 15., 60.},
				      { 60., 102.},
				      { 100., 200.},
				      { 200., 600.},
				      { 500., 3000.} };
      const double divY[nDivs][2] = { { 3., 300. }, // 10-60
				      { 3., 300. }, // 60-100
				      { 2e-2, 10.}, // 100-200
				      { 1e-4, 0.1}, // 200-600
				      { 1e-9,1e-3} }; // >600
      for (int iD=0; iD<nDivs; iD++) {
	TString rangeStr=Form("%d_%d",int(divX[iD][0]),int(divX[iD][1]));
	TString titleStr="combined measurement: " + rangeStr + axisStr;
	TH2D *h2frame= new TH2D("h2frame" + rangeStr,titleStr,
				100,divX[iD][0],divX[iD][1],
				100,divY[iD][0],divY[iD][1]);
	h2frame->SetDirectory(0);
	h2frame->SetStats(0);
	TString cName="cMeas_" + rangeStr;
	TCanvas *cMeasD= new TCanvas(cName,cName,600,600);
	cMeasD->SetGrid(1,1);
	cMeasD->SetLogy();
	logAxis(h2frame);
	if (divY[iD][0]<1e-5) h2frame->GetYaxis()->SetNoExponent(0);
	if (divX[iD][0]>490) {
	  cMeasD->SetLeftMargin(0.15);
	  h2frame->GetYaxis()->SetTitleOffset(2.1);
	}
	else {
	  cMeasD->SetLeftMargin(0.11);
	  h2frame->GetYaxis()->SetTitleOffset(1.4);
	}
	h2frame->Draw();
	grLL_line->Draw("L");
	grLL->Draw("PE1 same");
	grEE->Draw("PE1 same");
	grMM->Draw("PE1 same");
	cMeasD->Update();
      }
    }


    //return 1;
  }

  if (0) {
    // check split
    checkSplit(blue,covEE,covMM);
    return 1;
  }

  TMatrixD covEEFinal( blue.contributedCov(blue.measCovMatrix(covEE,1)) );
  TMatrixD covMMFinal( blue.contributedCov(blue.measCovMatrix(covMM,0)) );
  TMatrixD covEMFinal( blue.contributedCov(blue.measCovMatrix(covEM,2)) );
  TMatrixD finalCov(*blue.getCov());

  h1binning->GetXaxis()->SetTitle(massStr);
  TCanvas *cFinCov= plotCovCorr(finalCov,h1binning,"h2finCov","cFinCov",
				NULL,NULL,1,1.8);
  cFinCov->SetGrid(1,1);

  h1binning->GetXaxis()->SetTitle(eemassStr);
  TH2D* h2eeCorrPart= convert2histo(covToCorrPartial(covEEFinal,finalCov),
				    h1binning,"h2eeCorrPart","h2eeCorrPart");
  logAxis(h2eeCorrPart);
  TCanvas *cCorrPartEE= plotHisto(h2eeCorrPart,"cCorrPartEE",1,1,1.8,1,1.);
  cCorrPartEE->SetGrid(1,1);

  h1binning->GetXaxis()->SetTitle(mmmassStr);
  TH2D* h2mmCorrPart= convert2histo(covToCorrPartial(covMMFinal,finalCov),
				    h1binning,"h2mmCorrPart","h2mmCorrPart");
  logAxis(h2mmCorrPart);
  TCanvas *cCorrPartMM= plotHisto(h2mmCorrPart,"cCorrPartMM",1,1,1.8,1,1.);
  cCorrPartMM->SetGrid(1,1);

  return 1;
}

// -----------------------------------------------------------------------

// check BLUEResult object

void checkSplit(const BLUEResult_t &blue,
		const TMatrixD &covEE, const TMatrixD &covMM)
{
  TMatrixD measCovEE(blue.measCovMatrix(covEE,1));
  TMatrixD measCovMM(blue.measCovMatrix(covMM,0));
  if (0) {
    TMatrixD measCovTot(measCovEE, TMatrixD::kPlus, measCovMM);
    TMatrixD chkDiff(measCovTot, TMatrixD::kMinus, *blue.getInpCov());
      std::cout << "check difference: "; chkDiff.Print();
  }

  TMatrixD measCorrEE( covToCorr(measCovEE) );
  TMatrixD measCorrMM( covToCorr(measCovMM) );

  TCanvas *cSplit= new TCanvas("cSplit","cSplit",1200,600);
  cSplit->Divide(2,1);
  cSplit->cd(1);
  measCorrEE.Draw("COLZ");
  cSplit->cd(2);
  measCorrMM.Draw("COLZ");
  cSplit->Update();

  if (1) {
    TMatrixD corrEE( covToCorr(covEE) );
    TMatrixD measCovEM(blue.measCovMatrix(corrEE,2));
    if (0) { // to have something on a diagonal
      measCovEM += removeNaNs(measCorrEE);
      measCovEM += removeNaNs(measCorrMM);
    }
    TCanvas *cEMCovChk= new TCanvas("cEMCovChk","cEMCovChk",600,600);
    measCovEM.Draw("COLZ");
    cEMCovChk->Update();
  }

}

// -----------------------------------------------------------------------
