#include <TFile.h>
#include <TROOT.h>
#include <TString.h>
#include <TH1D.h>
#include <TCanvas.h>
#include "../RooUnfold/src/RooUnfoldResponse.h"
#include "../RooUnfold/src/RooUnfoldBayes.h"
#include <iostream>

// ----------------------------------------------------

void compareHistos(const TH1D *h1a, const TH1D *h1b);

// ----------------------------------------------------

void simpleCS()
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
  RooUnfoldBayes fsrBayes( detFSRRes, h1postFSR, nFSRResIters );
  TH1D *h1preFSR=(TH1D*)fsrBayes.Hreco();

  // final cross section
  TH1D *h1final=(TH1D*)h1preFSR->Clone("h1final");
  for (int ibin=1; ibin<=h1final->GetNbinsX(); ibin++) {
    double w= h1final->GetBinWidth(ibin);
    h1final->SetBinContent( ibin, h1final->GetBinContent(ibin)/w );
    h1final->SetBinError  ( ibin, 0 ); //h1final->GetBinError(ibin)/w );
  }

  // plot

  h1FinalXS->SetLineColor(kBlue);
  h1FinalXS->SetMarkerColor(kBlue);
  h1final->SetMarkerStyle(24);

  TCanvas *cx= new TCanvas("cs","cs",600,600);
  cx->SetLogx();
  cx->SetLogy();
  h1FinalXS->Draw("hist");
  //h1FinalXS->Draw("LPE1 same");
  //h1final->Draw("LPE1 same");
  h1final->Draw("LPsame");
  cx->Update();

  // print ratio
  if (1) compareHistos(h1FinalXS,h1final);
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
