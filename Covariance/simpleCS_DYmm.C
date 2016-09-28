#include <TFile.h>
#include <TROOT.h>
#include <TString.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include "../RooUnfold/src/RooUnfoldResponse.h"
#include "../RooUnfold/src/RooUnfoldBayes.h"
#include "../RooUnfold/src/RooUnfoldInvert.h"
#include <iostream>

// ---------------------------------------------------------

void printHisto(const TH1D *h1, int extraRange=0) {
  std::cout << "\nhisto " << h1->GetName() << " " << h1->GetTitle() << "\n";
  int d=(extraRange) ? 1 : 0;
  for (int ibin=1-d; ibin<=h1->GetNbinsX()+d; ibin++) {
    std::cout << "ibin=" << ibin << " " << h1->GetBinLowEdge(ibin)
	      << " " << (h1->GetBinLowEdge(ibin)+h1->GetBinWidth(ibin))
	      << "  " << h1->GetBinContent(ibin) << " +- "
	      << h1->GetBinError(ibin) << "\n";
  }
}

// ---------------------------------------------------------


// ----------------------------------------------------

void compareHistos(const TH1D *h1a, const TH1D *h1b);
TH1D* loadGraph(TFile &fin, TString graphNameOnFile, TString histoName,
		int posNegErrs, int plotIt=0);

// ----------------------------------------------------

void simpleCS_DYmm(int checkBackgrounds=0, int loadTheory=0)
{
  TString inpPath="/mnt/sdb/andriusj/v20160915_CovInput_ApprovedResults/";
  //inpPath="/mnt/sdb/andriusj/v20160527_1st_CovarianceMatrixInputs_76X/";
  TString dataInputFName=inpPath + "ROOTFile_Histograms_Data.root";
  const TString ending[2]= { "HLTv4p2", "HLTv4p3" };
  TString h1MeasNameBase="h_data_";
  TString h1YieldNameBase="h_yield_";

  TString bkgInputFName=inpPath + "ROOTFile_Histograms_Bkg.root";
  const int nBkgs=8;
  const TString bkgNames[nBkgs] = { "h_ZZ", "h_WZ", "h_WW", "h_ttbar",
				    "h_DYtautau", "h_tW", "h_WJets", "h_QCD" };
  int scaleMCbkg=1;

  TString corrInputFName=inpPath + "ROOTFile_Input6_CrossCheck.root";
  TString gAccEffName="g_AccEff";
  TString gEffSFNameBase="g_EffSF_";
  TString h1FinalXSName="h_DiffXsec_Data";
  TString DetResName="UnfoldRes_DetectorResol";
  TString FSRResName="UnfoldRes_FSR";
  int nDetResIters=17;
  const double lumiTot= 2832.673;
  const double lumi[2]= { 865.919, lumiTot-865.919 };

  TString h1AccEffName="h_AccEff";
  TString h1EffSFName="h_EffSF";

  TString theoryFName=(loadTheory) ? "theory13TeVmm.root" : "";
  TString theoryHistoName="h1cs_theory";

  TFile finData(dataInputFName);
  if (!finData.IsOpen()) {
    std::cout << "failed to open <" << finData.GetName() << ">\n";
    return;
  }
  TH1D* h1meas[2];
  TH1D* h1yield[2];
  h1meas[0]=(TH1D*)finData.Get(h1MeasNameBase + ending[0]);
  h1meas[1]=(TH1D*)finData.Get(h1MeasNameBase + ending[1]);
  h1yield[0]=(TH1D*)finData.Get(h1YieldNameBase + ending[0]);
  h1yield[1]=(TH1D*)finData.Get(h1YieldNameBase + ending[1]);
  for (int i=0; i<2; i++) {
    if (!h1meas[i] || !h1yield[i]) {
      std::cout << "failed to load histos from <"
		<< finData.GetName() << ">\n";
      return;
    }
    else {
      h1meas[i]->SetDirectory(0);
      h1yield[i]->SetDirectory(0);
    }
  }
  finData.Close();

  TH1D *h1bkgTot= (TH1D*)h1meas[0]->Clone("h1bkgTot");
  h1bkgTot->Reset();
  h1bkgTot->SetDirectory(0);
  h1bkgTot->Sumw2();
  TH1D *h1bkg_fromMC= (TH1D*)h1bkgTot->Clone("h1bkg_fromMC");
  TH1D *h1bkg_fromData=(TH1D*)h1bkgTot->Clone("h1bkg_fromData");
  TFile finBkg(bkgInputFName);
  if (!finBkg.IsOpen()) {
    std::cout << "failed to open the file <" << finBkg.GetName() << ">\n";
    return;
  }
  for (int iBkg=0; iBkg<nBkgs; iBkg++) {
    TH1D *h1= (TH1D*)finBkg.Get(bkgNames[iBkg]);
    if (!h1) {
      std::cout << "Failed to get <" << bkgNames[iBkg] << ">\n";
      return;
    }
    if ((bkgNames[iBkg]=="h_ZZ") || (bkgNames[iBkg]=="h_WZ")) {
      h1bkg_fromMC->Add(h1);
    }
    else h1bkg_fromData->Add(h1);
    h1bkgTot->Add(h1);
  }
  finBkg.Close();

  // background check requires efficiency correction scale factor

  // load other corrections and the results

  TFile finCorr(corrInputFName);
  if (!finCorr.IsOpen()) {
    std::cout << "failed to open <" << finCorr.GetName() << ">\n";
    return;
  }
  TH1D *h1AccEff=loadGraph(finCorr,gAccEffName,h1AccEffName,0);
  TH1D* h1EffSF[2];
  h1EffSF[0]=loadGraph(finCorr,gEffSFNameBase+ending[0],h1EffSFName+ending[0],0);
  h1EffSF[1]=loadGraph(finCorr,gEffSFNameBase+ending[1],h1EffSFName+ending[1],0);
  RooUnfoldResponse *detResResp=(RooUnfoldResponse*)finCorr.Get(DetResName);
  RooUnfoldResponse *detFSRRes=(RooUnfoldResponse*)finCorr.Get(FSRResName);
  TH1D *h1FinalXS=(TH1D*)finCorr.Get(h1FinalXSName);
  if (!h1AccEff || !h1EffSF[0] || !h1EffSF[1] || !h1FinalXSName) {
    std::cout << "failed to get either " << h1AccEffName << ", or "
	      << h1EffSFName << ", or " << h1FinalXSName << " from file <"
	      << finCorr.GetName() << ">\n";
    return;
  }
  if (!detResResp || !detFSRRes) {
    std::cout << "failed to get either " << DetResName << " or "
	      << FSRResName << "from file <"
	      << finCorr.GetName() << ">\n";
    return;
  }
  h1AccEff->SetDirectory(0);
  h1EffSF[0]->SetDirectory(0);
  h1EffSF[1]->SetDirectory(0);
  h1FinalXS->SetDirectory(0);
  finCorr.Close();

  if (checkBackgrounds) {
    std::cout << "\n\n Check background\n";
    TH1D *h1Bkg[2];
    for (int i=0; i<2; i++) {
      h1Bkg[i] = (TH1D*)h1meas[i]->Clone("h1bkg_from_input_" + ending[i]);
      h1Bkg[i]->Add(h1yield[i],-1);
      h1Bkg[i]->SetDirectory(0);
      h1Bkg[i]->Scale(lumiTot/lumi[i]);
      //compareHistos(h1bkgTot,h1Bkg[i]);
      if (scaleMCbkg) {
	TH1D *h1bkg_check=(TH1D*)h1bkg_fromMC->Clone(Form("h1bkg_check_%d",i));
	h1bkg_check->Multiply(h1EffSF[i]);
	//h1bkg_check->Divide(h1EffSF[i]);
	h1bkg_check->Add(h1bkg_fromData);
	compareHistos(h1Bkg[i],h1bkg_check);
      }
      else compareHistos(h1Bkg[i],h1bkgTot);
    }
    compareHistos(h1Bkg[0],h1Bkg[1]);
    std::cout << "\n\n";
    return;
  }


  if (1) {
    printHisto(h1yield[0],1);
    printHisto(h1yield[1],1);
    printHisto((TH1D*)detResResp->Hmeasured(),1);
  }

  detResResp->UseOverflow(true);

  // correct for detRes
  RooUnfoldBayes detResBayes4p2( detResResp, h1yield[0], nDetResIters );
  RooUnfoldBayes detResBayes4p3( detResResp, h1yield[1], nDetResIters );
  TH1D *h1Unf[2];
  h1Unf[0]=(TH1D*)detResBayes4p2.Hreco();
  h1Unf[1]=(TH1D*)detResBayes4p3.Hreco();

  // correct for the scale factor
  // correct for efficiency and acceptance
  TH1D *h1postFSR_4p2=(TH1D*)h1Unf[0]->Clone("h1postFSR_" + ending[0]);
  TH1D *h1postFSR_4p3=(TH1D*)h1Unf[1]->Clone("h1postFSR_" + ending[1]);
  h1postFSR_4p2->Divide(h1EffSF[0]);
  h1postFSR_4p3->Divide(h1EffSF[1]);

  h1postFSR_4p2->Divide(h1AccEff);
  h1postFSR_4p3->Divide(h1AccEff);

  // add contributions
  TH1D* h1postFSR=(TH1D*)h1postFSR_4p2->Clone("h1postFSR_");
  h1postFSR->Add(h1postFSR_4p3);

  if (1) {
    printHisto(h1yield[0]);
    printHisto(h1yield[1]);
    printHisto(h1postFSR_4p2);
    printHisto(h1postFSR_4p3);
  }

  // FSR unfolding
  RooUnfoldInvert fsrInvert( detFSRRes, h1postFSR, "rooInv");
  TH1D *h1preFSR=(TH1D*)fsrInvert.Hreco();

  // normalize
  h1preFSR->Scale(1/lumiTot);

  // final cross section
  TH1D *h1finalCalc=(TH1D*)h1preFSR->Clone("h1finalCalc");
  for (int ibin=1; ibin<=h1finalCalc->GetNbinsX(); ibin++) {
    double w= h1finalCalc->GetBinWidth(ibin);
    h1finalCalc->SetBinContent( ibin, h1finalCalc->GetBinContent(ibin)/w );
    h1finalCalc->SetBinError  ( ibin, 0 ); //h1finalCalc->GetBinError(ibin)/w );
  }

  // theory
  TH1D *h1theory=NULL;
  if (theoryFName.Length()) {
    TFile finTh(theoryFName);
    if (!finTh.IsOpen()) {
      std::cout << "Theory file <" << finTh.GetName() << "> not found\n";
      return;
    }
    TH1D *h1th_inp=(TH1D*)finTh.Get(theoryHistoName);
    if (!h1th_inp) {
      std::cout << "failed to get " << theoryHistoName << " from file <"
		<< finTh.GetName() << ">\n";
      return;
    }
    h1th_inp->SetDirectory(0);
    finTh.Close();
    h1theory=(TH1D*)h1th_inp->Clone("h1theory");
    for (int ibin=1; ibin<=h1th_inp->GetNbinsX(); ibin++) {
      double w= h1th_inp->GetBinWidth(ibin);
      h1theory->SetBinContent(ibin, h1th_inp->GetBinContent(ibin)/w );
      h1theory->SetBinError  (ibin, h1th_inp->GetBinError(ibin)/w );
    }
    delete h1th_inp;
    h1theory->SetStats(0);
  }

  // plot

  h1FinalXS->SetLineColor(kBlue);
  h1FinalXS->SetMarkerColor(kBlue);
  h1finalCalc->SetMarkerStyle(24);
  if (h1theory) {
    h1theory->SetLineColor(kRed);
    h1theory->SetMarkerColor(kRed);
    h1theory->SetMarkerStyle(7);
  }

  TCanvas *cx= new TCanvas("cs","cs",600,600);
  cx->SetLogx();
  cx->SetLogy();
  h1FinalXS->Draw("hist");
  //h1FinalXS->Draw("LPE1 same");
  //h1finalCalc->Draw("LPE1 same");
  h1finalCalc->Draw("LPsame");
  if (h1theory) h1theory->Draw("LPE1 same");
  cx->Update();

  // print ratio
  if (1) compareHistos(h1FinalXS,h1finalCalc);
  if (1 && h1theory) compareHistos(h1theory,h1finalCalc);
}

// ----------------------------------------------------
// ----------------------------------------------------

TH1D* loadGraph(TFile &fin, TString graphNameOnFile, TString histoName,
		int posNegErrs, int plotIt)
{
  TGraphAsymmErrors *gr= (TGraphAsymmErrors*)fin.Get(graphNameOnFile);
  if (!gr) {
    std::cout << "failed to get " << graphNameOnFile << " from <"
	      << fin.GetName() << ">\n";
    return NULL;
  }

  Double_t *xLow= gr->GetEXlow();
  Double_t *xHigh= gr->GetEXhigh();
  Double_t *yLow= gr->GetEYlow();
  Double_t *yHigh= gr->GetEYhigh();
  int n=gr->GetN();
  Double_t *x= gr->GetX();
  Double_t *y= gr->GetY();

  Double_t *xb= new Double_t[n+1];

  //std::cout << "n=" << n << "\n";
  double factor=1.;
  //factor=1.05; // to see the difference
  for (int i=0; i<n; i++) {
    xb[i]=x[i]-xLow[i]*factor;
    //std::cout << "i=" << i << ", xLow=" << xLow[i] << ", xHigh=" << xHigh[i] << ", x=" << x[i] << "\n";
  }
  xb[n]= x[n-1]+xHigh[n-1];
  //std::cout << "last point " << xb[n] << "\n";


  TH1D *h1= new TH1D(histoName, histoName, n, xb);
  h1->SetDirectory(0);
  for (int ibin=1; ibin<=h1->GetNbinsX(); ibin++) {
    h1->SetBinContent(ibin, y[ibin-1]);
    double w=0.5;
    if (posNegErrs==-1) w=1.;
    else if (posNegErrs==1) w=0;
    h1->SetBinError  (ibin, w*yLow[ibin-1] + (1.-w)*yHigh[ibin-1]);
  }

  if (plotIt) {
    TString cName= "c" + histoName;
    TCanvas *c= new TCanvas(cName,cName,600,600);
    c->SetLogx();
    h1->SetLineColor(kRed);
    h1->Draw("LPE1");
    gr->Draw("LPE same");
    c->Update();
  }

  return h1;
}

// ----------------------------------------------------
// ----------------------------------------------------

void compareHistos(const TH1D *h1a, const TH1D *h1b)
{
  std::cout << "compare " << h1a->GetName() << " to "
	    << h1b->GetName() << "\n";
  for (int ibin=1; ibin<=h1a->GetNbinsX(); ibin++) {
    if ( (h1a->GetBinLowEdge(ibin) != h1b->GetBinLowEdge(ibin)) ||
	 (h1a->GetBinWidth(ibin) != h1b->GetBinWidth(ibin)) ) {
      std::cout << "bining mismatch at ibin=" << ibin << "\n";
	return;
    }
    std::cout << "ibin=" << ibin << " " << h1a->GetBinLowEdge(ibin)
	      << " " << (h1a->GetBinLowEdge(ibin)+h1a->GetBinWidth(ibin))
	      << "  " << h1a->GetBinContent(ibin)
	      << "  " << h1b->GetBinContent(ibin)
	      << "  " << h1a->GetBinContent(ibin)/h1b->GetBinContent(ibin);
    if (fabs(h1a->GetBinContent(ibin)/h1b->GetBinContent(ibin) - 1.)>1e-3) {
      std::cout << " *";
    }
    std::cout << "\n";
  }
}

// ----------------------------------------------------
// ----------------------------------------------------

