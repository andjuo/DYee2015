#include "inputs.h"
#include "crossSection.h"

void checkRooUnfoldResponse()
{
  TString inpFName="dyee_test_dressed_El3.root";
  TString rooUnfRespName= "rooUnf_detResRespPU";
  TString h1missName= "h1detRes_Miss";
  TString h1fakeName= "h1detRes_Fake";
  TString h1measName= "h1detRes_Meas";
  TString h1trueName= "h1detRes_True";
  TString h2migName = "h1detRes_Mig";

  TFile fin(inpFName);
  if (!fin.IsOpen()) {
    std::cout << "failed to open the file <" << fin.GetName() << ">\n";
    return;
  }
  RooUnfoldResponse *resp= loadRooUnfoldResponse(fin,rooUnfRespName,
						 rooUnfRespName,1);
  TH1D *h1miss= loadHisto(fin, h1missName, "h1miss","h1miss",1,h1dummy);
  TH1D *h1fake= loadHisto(fin, h1fakeName, "h1fake","h1fake",1,h1dummy);
  TH1D *h1meas= loadHisto(fin, h1measName, "h1meas","h1meas",1,h1dummy);
  TH1D *h1true= loadHisto(fin, h1trueName, "h1true","h1true",1,h1dummy);
  TH2D *h2mig = loadHisto(fin, h2migName, "h2mig", "h2mig", 1, h2dummy);

  fin.Close();
  if (!resp) {
    std::cout << "failed to get the response\n";
    return;
  }
  if (!h1miss || !h1fake || !h1meas || !h1true || !h2mig) {
    std::cout << "failed to load distributions\n";
    return;
  }

  TH1D *h1_fullMeas=cloneHisto(h1meas,"h1_fullMeas","h1_fullMeas");
  h1_fullMeas->Add(h1fake);
  TH1D *h1_fullTrue=cloneHisto(h1true,"h1_fullTrue","h1_fullTrue");
  h1_fullTrue->Add(h1miss);

  // compare histograms to rooUnfResponse
  //printRatio(h2mig, (TH2D*)resp->Hresponse(),1,1,1e-3);

  //printHisto(h1miss);
  //printHisto(h1fake);
  //printRatio(h1fake, (TH1D*)resp->Hfakes(),1,1,1e-3);
  //printRatio(h1_fullMeas, (TH1D*)resp->Hmeasured(),1,1,1e-3);
  //printRatio(h1_fullTrue, (TH1D*)resp->Htruth(),1,1,1e-3);

  // try to construct a new rooUnfResponse and compare histograms
  RooUnfoldResponse *respNew=
    new RooUnfoldResponse(h1_fullMeas, h1_fullTrue, h2mig,
			  "respNew","respNew");

  //printRatio(h2mig, (TH2D*)respNew->Hresponse(),1,1,1e-3);
  //printRatio(h1_fullMeas, (TH1D*)respNew->Hmeasured(),1,1,1e-3);
  //printRatio(h1_fullTrue, (TH1D*)respNew->Htruth(),1,1,1e-3);
  //printRatio(h1fake, (TH1D*)respNew->Hfakes(),1,1,1e-3);


  if (1) {
    std::cout << "try to unfold data\n";
    TH1D *h1data= loadHisto("cs_DYee_13TeV_El3.root","h1Signal","h1Signal",1,h1dummy);
    h1data->SetTitle("DYee signal;M_{ee} [GeV];signal yield");
    int nIters=17;
    RooUnfoldBayes bayes( resp, h1data, nIters, false );
    RooUnfoldBayes bayesNew( respNew, h1data, nIters, false );
    printRatio( (TH1D*)bayes.Hreco(), (TH1D*)bayesNew.Hreco(), 1,1,1e-3 );

    TMatrixD cov( bayes.Ereco(RooUnfold::kCovariance) );
    TMatrixD covNew( bayesNew.Ereco(RooUnfold::kCovariance) );
    TH2D *h2cov= convert2histo(cov,h1true,"h2cov","h2cov");
    TH2D *h2covNew= convert2histo(covNew,h1true,"h2covNew","h2covNew");
    printRatio( h2cov, h2covNew, 1,0,1e-3 );
  }
}
