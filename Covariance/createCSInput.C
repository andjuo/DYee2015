#include "crossSection.h"

// ---------------------------------------------------------

// ---------------------------------------------------------
// ---------------------------------------------------------

void createCSInput(int analysisIs2D=0)
{
  TString anStr(Form("%dD",(analysisIs2D)?2:1));

  //TH1D* h1Dummy=NULL;
  TH2D* h2Dummy=NULL;

  TH2D* h2Yield_raw= loadHisto("dir-DYee-8TeV/bg-subtracted_yield_1D_peak20140220.root","Input/observedYield","h2ObsYield_raw",h2Dummy);
  if (!h2Yield_raw) return;
  //plotHisto(h2Yield_raw,"cRawH2");
  TH1D* h1Yield= flattenHisto(h2Yield_raw,"h1Yield");
  if (!h1Yield) return;
  //printHisto(h1Yield);
  //plotHisto(h1Yield,"cRawH1");

  TH2D *h2Bkg_raw= loadHisto("dir-DYee-8TeV/bg-subtracted_yield_1D_peak20140220.root","Input/ddBkgTotal","h2BkgYield_raw",h2Dummy);
  if (!h2Bkg_raw) return;
  //plotHisto(h2Bkg_raw,"cBkgRawH2");
  TH1D* h1Bkg= flattenHisto(h2Bkg_raw,"h1Bkg");
  if (!h1Bkg) return;
  //printHisto(h1Bkg);
  //plotHisto(h1Bkg,"cBkgH1");

  RooUnfoldResponse *rooUnfDetRes= loadRooUnfoldResponse("dir-DYee-8TeV/rooUnfDetRes_1D.root","rooUnfDetRes","rooUnfDetRes");
  if (!rooUnfDetRes) return;
  //plotHisto(*rooUnfDetRes,"cDetRes");

  RooUnfoldResponse *rooFSRRes= loadRooUnfoldResponse("dir-DYee-8TeV/rooUnfFSR_1D.root","rooUnfFSR","rooUnfFSR");
  if (!rooFSRRes) return;
  //plotHisto(*rooFSRRes,"cFSRRes");
  //printHisto(rooFSRRes->Hfakes());

  TH2D *h2Eff_raw= loadHisto("dir-DYee-8TeV/efficiency_1D.root","hEfficiency","h2Eff_raw",h2Dummy);
  if (!h2Eff_raw) return;
  TH1D *h1Eff= flattenHisto(h2Eff_raw,"h1Eff");
  if (!h1Eff) return;
  //plotHisto(h1Eff,"cEffH1");

  TH2D *h2EffPass_raw= loadHisto("dir-DYee-8TeV/efficiency_1D.root","hSumPass","h2EffPass_raw",h2Dummy);
  if (!h2EffPass_raw) return;
  TH1D *h1EffPass= flattenHisto(h2EffPass_raw,"h1EffPass");
  if (!h1EffPass) return;
  //plotHisto(h1EffPass,"cEffPassH1");

  TH2D *h2EffTot_raw= loadHisto("dir-DYee-8TeV/efficiency_1D.root","hSumTotal","h2EffTot_raw",h2Dummy);
  if (!h2EffTot_raw) return;
  TH1D *h1EffTot= flattenHisto(h2EffTot_raw,"h1EffTot");
  if (!h1EffTot) return;
  //plotHisto(h1EffTot,"cEffTotH1");

  //TH1D *h1EffTest=(TH1D*)h1EffPass->Clone("h1EffTest");
  //h1EffTest->Divide(h1EffTot);
  //h1EffTest->SetMarkerStyle(24);
  //plotHistoSame(h1EffTest,"cEffH1","LPE");

  TH2D *h2Acc_raw= loadHisto("dir-DYee-8TeV/acceptance_1D.root","hAcceptance","h2Acc_raw",h2Dummy);
  if (!h2Acc_raw) return;
  TH1D *h1Acc= flattenHisto(h2Acc_raw,"h1Acc");
  if (!h1Acc) return;
  //plotHisto(h1Acc,"cAccH1");

  TH2D *h2AccPass_raw= loadHisto("dir-DYee-8TeV/acceptance_1D.root","hSumPass","h2AccPass_raw",h2Dummy);
  if (!h2AccPass_raw) return;
  TH1D *h1AccPass= flattenHisto(h2AccPass_raw,"h1AccPass");
  if (!h1AccPass) return;
  //plotHisto(h1AccPass,"cAccPassH1");

  TH2D *h2AccTot_raw= loadHisto("dir-DYee-8TeV/acceptance_1D.root","hSumTotal","h2AccTot_raw",h2Dummy);
  if (!h2AccTot_raw) return;
  TH1D *h1AccTot= flattenHisto(h2AccTot_raw,"h1AccTot");
  if (!h1AccTot) return;
  //plotHisto(h1AccTot,"cAccTotH1");

  //TH1D *h1AccTest=(TH1D*)h1AccPass->Clone("h1AccTest");
  //h1AccTest->Divide(h1AccTot);
  //h1AccTest->SetMarkerStyle(24);
  //plotHistoSame(h1AccTest,"cAccH1","LPE");

  //TH1D *h1AccPassClone=(TH1D*)h1AccPass->Clone("h1AccPassClone");
  //h1AccPassClone->SetMarkerStyle(24);
  //plotHistoSame(h1AccPassClone,"cEffTotH1","LPE");

  TH1D *h1EffAcc=(TH1D*)h1EffPass->Clone("h1EffAcc");
  h1EffAcc->Divide(h1EffPass,h1AccTot,1,1,"B");
  //plotHisto(h1EffAcc,"cEffAccH1");

  //TH1D* h1EffAccTest=(TH1D*)h1Eff->Clone("h1EffAccTest");
  //h1EffAccTest->Multiply(h1Acc);
  //h1EffAccTest->SetMarkerStyle(24);
  //plotHistoSame(h1EffAccTest,"cEffAccH1","LPE");


  TH1D *h1Rho= loadVectorD("dir-DYee-8TeV/covRhoFileSF_nMB41_asymHLT_Unregressed_energy-allSyst_100.root","scaleFactorFlatIdxArray","scaleFactorErrFlatIdxArray","h1Rho", "Rho;M_{ee} [GeV];#rho", h1Yield);
  if (!h1Rho) return;
  //plotHisto(h1Rho, "cRhoH1");

  TH1D* h1Theory= loadHisto("dir-DYee-8TeV/1Dabsxsec_NNLO_CTEQ12NNLO_41.root",
			    "invm_FEWZ_41","h1Theory", h1Rho);
  if (!h1Theory) return;
  //plotHisto(h1Theory, "cTheoryH1");

  CrossSection_t cs("xs","8TeV_1D");
  cs.lumi(19712.);
  cs.h1Yield(h1Yield);
  cs.h1Bkg(h1Bkg);
  cs.h1Eff(h1Eff);
  cs.h1Rho(h1Rho);
  cs.h1Acc(h1Acc);
  cs.detRes(*rooUnfDetRes);
  cs.fsrRes(*rooFSRRes);
  cs.h1Theory(perMassBinWidth(h1Theory));

  cs.plotCrossSection();
  //return;

  std::cout << "cs.checkPtrs()=" << cs.checkPtrs() << "\n";
  if (!cs.save("cs_DYee_8TeV.root")) {
    std::cout << "saving failed\n";
    return;
  }

}

// ---------------------------------------------------------
// ---------------------------------------------------------


// ---------------------------------------------------------
