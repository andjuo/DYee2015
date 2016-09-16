#include <TFile.h>
#include <TROOT.h>
#include <TString.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include "../RooUnfold/src/RooUnfoldResponse.h"
#include "../RooUnfold/src/RooUnfoldBayes.h"
#include <iostream>

// ----------------------------------------------------

void compareHistos(const TH1D *h1a, const TH1D *h1b);

// ----------------------------------------------------

void simpleCS(int loadTheory=0, int loadSet2=0)
{
  TString dataInputFName="/mnt/sdb/andriusj/v3_09092016_CovarianceMatrixInputs/ROOTFile_Input1_Histograms_Data.root";
  TString h1DataYieldName="h_yield";

  TString corrInputFName="/mnt/sdb/andriusj/v3_09092016_CovarianceMatrixInputs/ROOTFile_Input6_CrossCheck.root";
  TString h1AccEffName="h_AccEff";
  TString h1EffSFName="h_EffSF";
  TString h1FinalXSName="h_diffXsec_Meas";
  TString DetResName="Unfold_DetectorRes";
  TString FSRResName="Unfold_FSRCorr";
  int nDetResIters=15;
  int nFSRResIters=15;
  double lumi=2316.969;

  TString theoryFName=(loadTheory) ? "theory13TeVmm.root" : "";
  TString theoryHistoName="h1cs_theory";

  if (loadSet2) {
    corrInputFName="dyee_test_dressed_El2skim2.root";
    h1AccEffName="h1EffPUAcc";
    h1EffSFName="h1rho_inPostFsrAcc";
    h1FinalXSName="h1_preFsr_Mweighted";
    DetResName="rooUnf_detResRespPU";
    FSRResName="rooUnf_fsrResp";
  }

  TFile finData(dataInputFName);
  if (!finData.IsOpen()) {
    std::cout << "failed to open <" << finData.GetName() << ">\n";
    return;
  }
  TH1D *h1yield=(TH1D*)finData.Get(h1DataYieldName);
  if (!h1yield) {
    std::cout << "failed to get " << h1DataYieldName << " from <"
	      << finData.GetName() << ">\n";
    return;
  }
  h1yield->SetDirectory(0);
  finData.Close();

  TFile finCorr(corrInputFName);
  if (!finCorr.IsOpen()) {
    std::cout << "failed to open <" << finCorr.GetName() << ">\n";
    return;
  }
  TH1D *h1AccEff=(TH1D*)finCorr.Get(h1AccEffName);
  TH1D *h1EffSF=(TH1D*)finCorr.Get(h1EffSFName);
  TH1D *h1FinalXS=(TH1D*)finCorr.Get(h1FinalXSName);
  RooUnfoldResponse *detResResp=(RooUnfoldResponse*)finCorr.Get(DetResName);
  RooUnfoldResponse *detFSRRes=(RooUnfoldResponse*)finCorr.Get(FSRResName);
  if (!h1AccEff || !h1EffSF || !h1FinalXSName) {
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
  h1EffSF->SetDirectory(0);
  h1FinalXS->SetDirectory(0);
  finCorr.Close();

  // correct for detRes
  RooUnfoldBayes detResBayes( detResResp, h1yield, nDetResIters );
  TH1D *h1Unf=(TH1D*)detResBayes.Hreco();

  // correct for efficiency and acceptance
  TH1D *h1postFSR=(TH1D*)h1Unf->Clone("h1postFSR");
  h1postFSR->Divide(h1EffSF);
  h1postFSR->Divide(h1AccEff);

  // correct for FSR
  h1postFSR->Scale(1/lumi);
  std::cout << "first FSR mig bin: " << ((TH2D*)(detFSRRes->Hresponse()))->GetBinContent(1,1) << "\n";
  RooUnfoldBayes fsrBayes( detFSRRes, h1postFSR, nFSRResIters );
  TH1D *h1preFSR=(TH1D*)fsrBayes.Hreco();

  // final cross section
  TH1D *h1final=(TH1D*)h1preFSR->Clone("h1final");
  for (int ibin=1; ibin<=h1final->GetNbinsX(); ibin++) {
    double w= h1final->GetBinWidth(ibin);
    h1final->SetBinContent( ibin, h1final->GetBinContent(ibin)/w );
    h1final->SetBinError  ( ibin, 0 ); //h1final->GetBinError(ibin)/w );
    if (loadSet2) {
      h1FinalXS->SetBinContent( ibin, h1FinalXS->GetBinContent(ibin)/w );
      h1FinalXS->SetBinError  ( ibin, h1FinalXS->GetBinError(ibin)/w );
    }
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
  h1final->SetMarkerStyle(24);
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
  //h1final->Draw("LPE1 same");
  h1final->Draw("LPsame");
  if (h1theory) h1theory->Draw("LPE1 same");
  cx->Update();

  // print ratio
  if (1) compareHistos(h1FinalXS,h1final);
  if (1 && h1theory) compareHistos(h1theory,h1final);
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
    std::cout << "\n";
  }
}
