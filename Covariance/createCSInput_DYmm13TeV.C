#include "crossSection.h"
#include <TGraphAsymmErrors.h>

// ---------------------------------------------------------

TH1D* convert(TGraphAsymmErrors *gr, TString hName, TString hTitle, int plotIt=0);

// ---------------------------------------------------------
// ---------------------------------------------------------

void createCSInput_DYmm13TeV()
{
  TString srcPath="/media/ssd/v20160214_1st_CovarianceMatrixInputs/";

  TH1D* h1Dummy=NULL;
  //TH2D* h2Dummy=NULL;

  double lumi42= 848.104;
  double lumi43= 2765.229-lumi42;
  double lumiTot= lumi42 + lumi43;

  TString fname1= srcPath + TString("Input1/ROOTFile_Histograms_Data.root");
  TFile fin1(fname1);
  if (!fin1.IsOpen()) {
    std::cout << "failed to open the file <" << fin1.GetName() << ">\n";
    return;
  }
  TH1D* h1Yield_42= loadHisto(fin1, "h_data_HLTv4p2", "h1Yield_42",1,h1Dummy);
  TH1D* h1Yield_43= loadHisto(fin1, "h_data_HLTv4p3", "h1Yield_43",1,h1Dummy);
  if (!h1Yield_42) return;
  if (!h1Yield_43) return;
  //plotHisto(h1Yield_42,"cYield42",1,1);
  //plotHisto(h1Yield_43,"cYield43",1,1);

  TH1D* h1Signal_42= loadHisto(fin1,"h_yield_HLTv4p2","h1Signal_42",1,h1Dummy);
  TH1D* h1Signal_43= loadHisto(fin1,"h_yield_HLTv4p3","h1Signal_43",1,h1Dummy);
  if (!h1Signal_42) return;
  if (!h1Signal_43) return;
  //plotHisto(h1Signal_42,"cSignal42",1,1);
  //plotHisto(h1Signal_43,"cSignal43",1,1);
  fin1.Close();


  TString fname2= srcPath + TString("Input2/ROOTFile_Histograms_Bkg.root");
  TFile fin2(fname2);
  if (!fin2.IsOpen()) {
    std::cout << "failed to open the file <" << fin2.GetName() << ">\n";
    return;
  }
  TH1D* h1BkgTot= cloneHisto(h1Yield_42, "h1BkgTot","h1BkgTot");
  h1BkgTot->Reset();
  //printHisto(h1BkgTot);
  std::vector<TH1D*> hBkgV;
  std::vector<double> bkgWeights;
  hBkgV.reserve(_bkgLast);
  bkgWeights.reserve(_bkgLast);
  for (TBkg_t b=_bkgZZ; b<_bkgLast; next(b)) {
    std::cout << "b=" << bkgName(b) << "\n";
    TH1D *h1= loadHisto(fin2, "h_"+bkgName(b), "h_"+bkgName(b),1,h1Dummy);
    if (!h1) return;
    double w=1.;
    if (1) {
      if (b==_bkgZZ) { w= 15.4/996944.; }
      else if (b==_bkgWZ) { w= 66.1/978512.; }
      else if (b==_bkgWW) { w= 118.7/993640.; }
      if (w!=double(1.)) h1->Scale( w * lumiTot );
    }
    bkgWeights.push_back(1);
    hBkgV.push_back(h1);
    h1BkgTot->Add(h1);
  }
  fin2.Close();
  plotHisto(h1BkgTot, "cBkgTot",1,1);


  if (0) {
    TH1D* h1Bkg42_chk= cloneHisto(h1Yield_42, "h1Bkg42_chk", "h1Bkg42_chk");
    h1Bkg42_chk->Add(h1Signal_42, -1);
    TH1D* h1Bkg43_chk= cloneHisto(h1Yield_43, "h1Bkg43_chk", "h1Bkg43_chk");
    h1Bkg43_chk->Add(h1Signal_43, -1);
    plotHistoSame(h1Bkg42_chk, "cBkgTot","hist");
    plotHistoSame(h1Bkg43_chk, "cBkgTot","hist");

    TH1D *h1ratio42= cloneHisto(h1Bkg42_chk, "h1ratio42","ratio42");
    TH1D *h1ratio43= cloneHisto(h1Bkg43_chk, "h1ratio43","ratio43");
    h1ratio42->Divide(h1BkgTot);
    h1ratio43->Divide(h1BkgTot);
    h1ratio43->SetLineColor(kGreen);
    if (1) {
      h1ratio42->Scale(lumiTot/lumi42);
      h1ratio43->Scale(lumiTot/lumi43);
    }
    h1ratio42->GetYaxis()->SetRangeUser(0,2.);
    plotHisto(h1ratio42, "cRatio42",1);
    plotHistoSame(h1ratio43, "cRatio42", "hist");
    printHisto(h1ratio42);
    TH1D *h1ratio43_to_42= cloneHisto(h1ratio43, "h1ratio43_to_42","ratio43_to_42");
    h1ratio43_to_42->Divide(h1ratio42);
    printHisto(h1ratio43_to_42);

    if (1) {
      TH1D *h1BkgSum= cloneHisto(h1Bkg42_chk, "h1BkgSum", "BkgSum");
      h1BkgSum->Add(h1Bkg43_chk);
      TH1D *h1Bkg_ChkDiff= cloneHisto(h1BkgSum, "h1Bkg_ChkDiff", "Bkg_ChkDiff");
      h1Bkg_ChkDiff->Add(h1BkgTot, -1);
      plotHisto(h1Bkg_ChkDiff, "cBkg_ChkDiff",1,0);
    }
  }

  TString fname6= srcPath + TString("Input6/ROOTFile_Input6_CrossCheck.root");
  TFile fin6(fname6);
  if (!fin6.IsOpen()) {
    std::cout << "failed to open the file <" << fin6.GetName() << ">\n";
    return;
  }
  TGraphAsymmErrors *grAcc=(TGraphAsymmErrors*)fin6.Get("g_Acc");
  TGraphAsymmErrors *grEff=(TGraphAsymmErrors*)fin6.Get("g_Eff");
  TGraphAsymmErrors *grAccEff=(TGraphAsymmErrors*)fin6.Get("g_AccEff");
  TGraphAsymmErrors *grEffSF_HLTv4p2=(TGraphAsymmErrors*)fin6.Get("g_EffSF_HLTv4p2");
  TGraphAsymmErrors *grEffSF_HLTv4p3=(TGraphAsymmErrors*)fin6.Get("g_EffSF_HLTv4p3");
  RooUnfoldResponse *rooUnfDetRes= loadRooUnfoldResponse(fin6,"UnfoldRes_DetectorResol", "rooUnfDetRes");
  RooUnfoldResponse *rooFSRRes= loadRooUnfoldResponse(fin6,"UnfoldRes_FSR", "rooUnfFSR");
  TH1D *h1_DiffXSec_data= loadHisto(fin6,"h_DiffXsec_Data", "h1KL_xsec", 1, h1Dummy);
  fin6.Close();
  if (!grAcc || !grEff || !grAccEff || !grEffSF_HLTv4p2 || !grEffSF_HLTv4p3 ||
      !rooUnfDetRes || !rooFSRRes || !h1_DiffXSec_data) return;
  //plotHisto(*rooUnfDetRes,"cDetRes",1,1);
  //plotHisto(*rooFSRRes,"cFSRRes",1,1);

  int plotAccEff=0;
  TH1D *h1Acc= convert(grAcc, "h1Acc","Acceptance;M_{#mu#mu} [GeV];A", plotAccEff);
  TH1D *h1Eff= convert(grEff ,"h1Eff","Efficiency;M_{#mu#mu} [GeV];#epsilon", plotAccEff);
  TH1D *h1EffAcc= convert(grAccEff, "h1EffAcc","Acc #times Eff;M_{#mu#mu} [GeV];A #times #epsilon", plotAccEff);
  if (!h1Acc || !h1Eff || !h1EffAcc) return;
  int plotIt_sf=0;
  TH1D *h1SF42= convert(grEffSF_HLTv4p2, "h1SF42","ScaleFactors 4.2;M_{#mu#mu} [GeV];eff.s.f. (v.4.2)", plotIt_sf);
  if (!h1SF42) return;
  TH1D *h1SF43= convert(grEffSF_HLTv4p3, "h1SF43","ScaleFactors 4.3;M_{#mu#mu} [GeV];eff.s.f. (v.4.3)", plotIt_sf);
  if (!h1SF43) return;
  

  MuonCrossSection_t muCS("muCS","",lumi42,lumi43);
  muCS.editCSa().h1Yield( h1Yield_42 );
  muCS.editCSa().h1Eff( h1Eff );
  muCS.editCSa().h1Acc( h1Acc );
  muCS.editCSa().h1EffAcc( h1EffAcc );
  muCS.editCSa().h1Rho( h1SF42 );
  muCS.editCSa().detRes( *rooUnfDetRes );
  muCS.editCSa().fsrRes( *rooFSRRes );

  muCS.editCSb().h1Yield( h1Yield_43 );
  muCS.editCSb().h1Eff( h1Eff );
  muCS.editCSb().h1Acc( h1Acc );
  muCS.editCSb().h1EffAcc( h1EffAcc );
  muCS.editCSb().h1Rho( h1SF43 );
  muCS.editCSb().detRes( *rooUnfDetRes );
  muCS.editCSb().fsrRes( *rooFSRRes );

  //muCS.h1Bkg( h1BkgTot );
  muCS.setBkgV(hBkgV,bkgWeights);
  muCS.h1Theory(h1_DiffXSec_data);

  if (0) {
    TH1D *h1preFSRcs= muCS.calcCrossSection();
    if (!h1preFSRcs) return;
    plotHisto(h1preFSRcs, "cPreFsrCS",1,1);
    h1_DiffXSec_data->SetMarkerStyle(24);
    h1_DiffXSec_data->SetLineColor(kGreen);
    h1_DiffXSec_data->SetMarkerColor(kGreen);
    plotHistoSame(h1_DiffXSec_data, "cPreFsrCS", "LPE1");
    printRatio(h1preFSRcs, h1_DiffXSec_data);
  }
  else {
    TCanvas *c= muCS.plotCrossSection();
    if (!c) return;
  }

  if (0) {
    h1Signal_42->SetMarkerStyle(24);
    plotHisto(h1Signal_42,"cSignal42",1,1);
    plotHistoSame(muCS.csA().copy(muCS.csA().h1Signal()),"cSignal42","LPE");

    h1Signal_43->SetMarkerStyle(24);
    plotHisto(h1Signal_43,"cSignal43",1,1);
    plotHistoSame(muCS.csB().copy(muCS.csB().h1Signal()),"cSignal43","LPE");
  }


  if (!muCS.save("cs_DYmm_13TeV_v1")) {
    std::cout << "saving failed\n";
    return;
  }
}

// ---------------------------------------------------------
// ---------------------------------------------------------

TH1D* convert(TGraphAsymmErrors *gr, TString hName, TString hTitle,
	      int plotIt)
{
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
  factor=1.05; // to see the difference
  for (int i=0; i<n; i++) {
    xb[i]=x[i]-xLow[i]*factor;
    //std::cout << "i=" << i << ", xLow=" << xLow[i] << ", xHigh=" << xHigh[i] << ", x=" << x[i] << "\n";
  }
  xb[n]= x[n-1]+xHigh[n-1];
  //std::cout << "last point " << xb[n] << "\n";


  TH1D *h1= new TH1D(hName, hTitle, n, xb);
  h1->SetDirectory(0);
  for (int ibin=1; ibin<=h1->GetNbinsX(); ibin++) {
    h1->SetBinContent(ibin, y[ibin-1]);
    h1->SetBinError  (ibin, 0.5*(yLow[ibin-1] + yHigh[ibin-1]));
  }

  if (plotIt) {
    TString cName= "c" + hName;
    TCanvas *c= new TCanvas(cName,cName,600,600);
    c->SetLogx();
    h1->SetLineColor(kRed);
    h1->Draw("LPE1");
    gr->Draw("LPE same");
    c->Update();
  }

  return h1;
}

// ---------------------------------------------------------