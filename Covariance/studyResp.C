// Small macro for testing purposes

#include "inputs.h"
//#include "../RooUnfold/src/RooUnfoldResponse.h"
//#include "../RooUnfold/src/RooUnfoldBayes.h"
#include "crossSection.h"

void studyResp()
{
  TString inpFile="dymm_test_RECO_Mu76X_withOverflow.root";
  //inpFile="dymm_test_RECO_Mu76X.root";
  std::cout << "inpFile=<" << inpFile << ">\n";
  RooUnfoldResponse *detResResp= loadRooUnfoldResponse(inpFile,
						       "rooUnf_detResRespPU",
						       "DRR");
  if (!detResResp) return;
  detResResp->UseOverflow(kFALSE); std::cout << "useOverflow(kFALSE)\n";
  //detResResp->UseOverflow(kTRUE); std::cout << "useOverflow(kTRUE)\n";

  // independent histograms
  TH1D* h1measInd= loadHisto(inpFile, "h1recoSel_MWPU", "h1measInd",h1dummy);
  TH1D* h1trueInd= loadHisto(inpFile, "h1_postFsrInAccSel_MWPU", "h1trueInd",h1dummy);
  if (!h1measInd || !h1trueInd) {
    std::cout << "failed to get h1measInd or h1trueInd\n";
    return;
  }

  // dependent histograms
  TH1D* h1meas= cloneHisto(detResResp->Hmeasured(),"h1meas","h1meas_orig",h1dummy);
  TH1D* h1true= cloneHisto(detResResp->Htruth(), "h1true","h1true_orig",h1dummy);
  TH2D* h2r= cloneHisto(detResResp->Hresponse(),"h2r","h2r_orig",h2dummy);

  plotHisto(h1meas, "cMeas", 1,1, "LPE1","measFromResp");
  plotHisto(h1true, "cTrue", 1,1, "LPE1","trueFromResp");

  //
  // 1. Check whether the migration matrix corresponds to meas and true
  // histograms in detResResp by calculating the profiles
  //
  if (1) {
    TH1D* h1countX= cloneHisto(h1meas,"h1countX","h1countX");
    h1countX->Reset();
    TH1D* h1countY= cloneHisto(h1true,"h1countY","h1countY");
    h1countY->Reset();

    for (int ibin=1; ibin<= h2r->GetNbinsX(); ibin++) {
      double sum=0;
      double sumErr2=0;
      for (int jbin=1; jbin<= h2r->GetNbinsY(); jbin++) {
	sum+= h2r->GetBinContent(ibin,jbin);
	sumErr2+= pow(h2r->GetBinError(ibin,jbin),2);
      }
      h1countX->SetBinContent(ibin, sum);
      h1countX->SetBinError  (ibin, sqrt(sumErr2));
    }

    for (int jbin=1; jbin<= h2r->GetNbinsY(); jbin++) {
      double sum=0;
      double sumErr2=0;
      for (int ibin=1; ibin<= h2r->GetNbinsX(); ibin++) {
	sum+= h2r->GetBinContent(ibin,jbin);
	sumErr2+= pow(h2r->GetBinError(ibin,jbin),2);
      }
      h1countY->SetBinContent(jbin, sum);
      h1countY->SetBinError  (jbin, sqrt(sumErr2));
    }

    histoStyle(h1countX, kRed,24);
    histoStyle(h1countY, kRed,24);

    plotHistoSame(h1countX, "cMeas", "LPE", "countX");
    printRatio(h1meas,h1countX);

    plotHistoSame(h1countY, "cTrue", "LPE", "countY");
    printRatio(h1true,h1countY);
  }

  if (1) {
    histoStyle(h1measInd, kGreen+1, 5);
    histoStyle(h1trueInd, kGreen+1, 5);
    plotHistoSame(h1measInd, "cMeas", "PE", "measInd");
    plotHistoSame(h1trueInd, "cTrue", "PE", "trueInd");

    printRatio(h1meas,h1measInd,1);
    printRatio(h1true,h1trueInd,1);
  }

  if (1) {
    int nIters=4;
    TH1D *h1measClone=cloneHisto(h1meas, "h1measClone","h1measClone");
    TH1D *h1trueClone=cloneHisto(h1true, "h1trueClone","h1trueClone");
    //h1measClone->SetBinContent(0,0);
    //h1measClone->SetBinError(0,0);
    RooUnfoldBayes bayesDetRes( detResResp, h1measClone, nIters );
    TH1D *h1Unf= (TH1D*) bayesDetRes.Hreco()->Clone("h1Unf");
    RooUnfoldBayes bayesDetResInd( detResResp, h1measInd, nIters );
    TH1D *h1UnfInd= (TH1D*) bayesDetResInd.Hreco()->Clone("h1UnfInd");
    printRatio(h1Unf, h1true,1);
    printRatio(h1UnfInd, h1true,1);
  }
}
