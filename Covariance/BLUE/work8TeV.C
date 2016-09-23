#include "Blue.h"
#include <map>
//#include <pair>
#include <sstream>
#include <fstream>
#include <TCanvas.h>
#include <TFile.h>

typedef std::pair<TString,TString> finfo_t;

int loadData(int is2D, const std::map<TString,finfo_t> &fnames,
	     TMatrixD &measEE, TMatrixD &covEE,
	     TMatrixD &measMM, TMatrixD &covMM);

int loadData(const finfo_t &fnames, TMatrixD &meas, TMatrixD &cov,
	     int is2D, int isEleChannel);

int combineData(const TMatrixD &measEE, const TMatrixD &covEE,
		const TMatrixD &measMM, const TMatrixD &covMM,
		int nXbins, const double *xBins);

const int nXbins_1D=41;
const double xBins_1D[nXbins_1D+1] = { 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 81, 86, 91, 96, 101, 106, 110, 115, 120, 126, 133, 141, 150, 160, 171, 185, 200, 220, 243, 273, 320, 380, 440, 510, 600, 1000, 1500, 2000 };

const int nXbins_2D=132;
const double xBins_2D[nXbins_2D+1] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132 };


// -----------------------------------------------------------------------
// -----------------------------------------------------------------------

void work8TeV(int is2D=0)
{
  int nBins=(is2D) ? 132 : 41;
  TString path="dir-8TeV-data/";

  std::map<TString,finfo_t> fNames;
  fNames["ee1D"] = finfo_t(path+"ee_results_2014_abs_adhoc_PI_Bayes.dat",
			   path+"covForBlue-1D-absolute_EE_Bayes.dat");
  fNames["mm1D"] = finfo_t(path+"mumu_results_2014_abs_PI_Bayes.dat",
			   path+"covForBlue-1D-absolute_MuMu.dat");
  fNames["ee2D"] = finfo_t(path+"ee_results_2D_2014_abs_adhoc_PI_Bayes.dat",
			   path+"covForBlue-2D-absolute_EE_Bayes.dat");
  fNames["mm2D"] = finfo_t(path+"mumu_results_2D_2014_abs_PI_Bayes.dat",
			   path+"covForBlue-2D-absolute_MuMu.dat");

  TMatrixD measEE(nBins,1), covEE(nBins,nBins);
  TMatrixD measMM(nBins,1), covMM(nBins,nBins);

  if (!loadData(is2D,fNames,measEE,covEE,measMM,covMM)) {
    std::cout << "failed to load data. stopping\n";
    return;
  }

  int nXBins=(is2D) ? nXbins_2D : nXbins_1D;
  const double *xBins= (is2D) ? xBins_2D : xBins_1D;

  if (!combineData(measEE,covEE,measMM,covMM, nXBins, xBins)) {
    std::cout << "failed to combined data. stopping\n";
    return;
  }

  return;
}

// -----------------------------------------------------------------------
// -----------------------------------------------------------------------

int loadData(int is2D, const std::map<TString,finfo_t> &fnames,
	     TMatrixD &measEE, TMatrixD &covEE,
	     TMatrixD &measMM, TMatrixD &covMM)
{
  TString eeKey= (is2D) ? "ee2D" : "ee1D";
  TString mmKey= (is2D) ? "mm2D" : "mm1D";
  std::map<TString,finfo_t>::const_iterator itEE=fnames.find(eeKey);
  std::map<TString,finfo_t>::const_iterator itMM=fnames.find(mmKey);
  int res=
    loadData(itEE->second,measEE,covEE,is2D,1) &&
    loadData(itMM->second,measMM,covMM,is2D,0);
  if (res && 0) {
    print("measEE",measEE);
    print("measMM",measMM);
  }
  if (!res) { std::cout << "failed to load data\n"; return 0; }
  return res;
}


// -----------------------------------------------------------------------

int loadData(const finfo_t &fnames, TMatrixD &meas, TMatrixD &cov,
	     int is2D, int isEleChannel)
{

  TMatrixD binW(meas);
  TMatrixD measErr(meas);

  // load the cross section
  std::ifstream finCS(fnames.first.Data());
  if (!finCS.is_open()) {
    std::cout << "failed to open the file <" << fnames.first << ">\n";
    return 0;
  }

  std::string s;
  if (finCS.peek()=='#') {
    std::getline(finCS,s);
    //std::cout << "skipping <" << s << ">\n";
  }
  int ibin, mMin, mMax;
  double cs,csErr;

  while (!finCS.eof()) {
    finCS >> ibin >> mMin >> mMax >> cs >> csErr;
    std::cout << ", " << mMin;
    if (finCS.rdstate()) break;
    binW(ibin-1,0)= mMax-mMin;
    meas(ibin-1,0)= cs;
    measErr(ibin-1,0)= csErr;
  }
  finCS.close();

  // load the covariance
  std::ifstream finCov(fnames.second.Data());
  if (!finCov.is_open()) {
    std::cout << "failed to open the file <" << fnames.second << ">\n";
    return 0;
  }

  int jbin;
  double var;
  while (!finCov.eof()) {
    finCov >> ibin >> jbin >> var;
    if (finCov.rdstate()) break;
    if (isEleChannel) {
      if (is2D==0) var /= (binW(ibin-1,0) * binW(jbin-1,0));
      else {
	double dx= ((ibin<=120) ? 0.1 : 0.2) * ((jbin<=120) ? 0.1 : 0.2);
	var /= dx;
      }
    }
    cov(ibin-1,jbin-1) = var;
  }
  finCov.close();

  if (0) {
    std::cout << "measurement : \n";
    for (int ir=0; ir<meas.GetNrows(); ir++) {
      std::cout << "ir=" << ir << ", " << meas(ir,0) << " +- " << measErr(ir,0)
		<< " covErr=" << sqrt(cov(ir,ir)) << "\n";
    }
  }

  return 1;
}

// -----------------------------------------------------------------------

int combineData(const TMatrixD &measEE, const TMatrixD &covEE,
		const TMatrixD &measMM, const TMatrixD &covMM,
		int nXbins, const double *xBins)
{
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

  if (1) {
    TH1D *h1EEErr= createErrHisto("h1EEErr","EE unc",covEE,nXbins,xBins);
    TH1D *h1MMErr= createErrHisto("h1MMErr","MM unc",covMM,nXbins,xBins);
    TH1D *h1LLErr= createErrHisto("h1LLErr","Comb unc",*blue.getCov(),nXbins,xBins);
    h1EEErr->SetLineColor(kBlue);
    h1EEErr->SetMarkerColor(kBlue);
    h1EEErr->SetMarkerStyle(24);
    h1EEErr->SetMarkerSize(0.5);
    h1LLErr->SetLineColor(kRed);
    h1LLErr->SetMarkerColor(kRed);
    h1MMErr->SetMarkerStyle(20);
    h1MMErr->SetMarkerSize(0.5);
    TCanvas *cErr= new TCanvas("cErr","cErr",600,600);
    if (h1EEErr->GetNbinsX()<50) cErr->SetLogx();
    cErr->SetLogy();
    h1EEErr->Draw("LPE1");
    h1LLErr->Draw("LPE1 same");
    h1MMErr->Draw("LPE1 same");
    cErr->Update();
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

    grEE->SetLineColor(kBlue);
    grEE->SetMarkerColor(kBlue);
    grLL->SetLineColor(kRed);
    grLL->SetMarkerColor(kRed);

    grEE->SetMarkerStyle(7);
    grMM->SetMarkerStyle(7);

    //grLL->GetXaxis()->SetRangeUser(10,50);
    //grLL->GetYaxis()->SetRangeUser(20,250);

    TCanvas *cMeas= new TCanvas("cMeas","cMeas", 600,600);
    if (grLL->GetN()<50) {
      cMeas->SetLogx();
      cMeas->SetLogy();
    }
    grLL->Draw("ALPE1");
    grEE->Draw("PE1 same");
    grMM->Draw("PE1 same");
    cMeas->Update();
    //return 1;
  }

  if (0) {
    // check split
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

    return 1;
  }

  TMatrixD covEEFinal( blue.contributedCov(blue.measCovMatrix(covEE,1)) );
  TMatrixD covMMFinal( blue.contributedCov(blue.measCovMatrix(covMM,0)) );
  TMatrixD covEMFinal( blue.contributedCov(blue.measCovMatrix(covEM,2)) );

  TCanvas *cFinCov= new TCanvas("cFinCov","cFinCov",600,600);
  TMatrixD finalCov(*blue.getCov());
  finalCov.Draw("COLZ");
  cFinCov->Update();

  TCanvas *cFinCorr=new TCanvas("cFinCorr","cFinCorr",600,600);
  TMatrixD finalCorr(*blue.getCorr());
  finalCorr.Draw("COLZ");
  cFinCorr->Update();

  TCanvas *cCorrPartEE= new TCanvas("cCorrPartEE","cCorrPartEE",600,600);
  TMatrixD eeCorrPart(covToCorrPartial(covEEFinal,finalCov));
  eeCorrPart.Draw("COLZ");
  cCorrPartEE->Update();

  TCanvas *cCorrPartMM= new TCanvas("cCorrPartMM","cCorrPartMM",600,600);
  TMatrixD mmCorrPart(covToCorrPartial(covMMFinal,finalCov));
  mmCorrPart.Draw("COLZ");
  cCorrPartMM->Update();

  return 1;
}

// -----------------------------------------------------------------------
