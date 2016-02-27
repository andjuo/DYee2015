// Small macro for testing purposes

#include "inputs.h"
//#include "../RooUnfold/src/RooUnfoldResponse.h"
//#include "../RooUnfold/src/RooUnfoldBayes.h"
#include "crossSection.h"

void studyResp()
{
  RooUnfoldResponse *detResResp= loadRooUnfoldResponse("dymm_test_RECO.root",
						       "rooUnf_detResResp",
						       "DRR");
  if (!detResResp) return;
  TH1D* h1meas= cloneHisto(detResResp->Hmeasured(),"h1meas","h1meas_orig",h1dummy);
  TH1D* h1true= cloneHisto(detResResp->Htruth(), "h1true","h1true_orig",h1dummy);
  TH2D* h2r= cloneHisto(detResResp->Hresponse(),"h2r","h2r_orig",h2dummy);

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

  plotHisto(h1meas, "cMeas", 1,1, "LPE1");
  plotHistoSame(h1countX, "cMeas", "LPE");
  printRatio(h1meas,h1countX);

  plotHisto(h1true, "cTrue", 1,1, "LPE1");
  plotHistoSame(h1countY, "cTrue", "LPE");
  printRatio(h1true,h1countY);
}
