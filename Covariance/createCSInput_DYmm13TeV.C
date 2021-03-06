#include "crossSection.h"
#include <TGraphAsymmErrors.h>

// ---------------------------------------------------------

// moved to inputs.h
//TH1D* convert(TGraphAsymmErrors *gr, TString hName, TString hTitle, int plotIt=0);

// ---------------------------------------------------------
// ---------------------------------------------------------

void createCSInput_DYmm13TeV(int doSave=0, int removeNegativeSignal=0,
			     int checkStep=-1)
{
  TString srcPath="/media/ssd/v20160214_1st_CovarianceMatrixInputs/";
  srcPath="/mnt/sdb/andriusj/v20160214_1st_CovarianceMatrixInputs/";
  TVersion_t inpVer=_verMu1;
  TString inpVerTag="_v1";
  double lumi42= 848.104;
  double lumiTot= 2765.229;

  TString fname1Base="Input1/ROOTFile_Histograms_Data.root";
  TString fname2Base="Input2/ROOTFile_Histograms_Bkg.root";
  TString fname6Base="Input6/ROOTFile_Input6_CrossCheck.root";

  if (0) {
    srcPath="/mnt/sdb/andriusj/v20160527_1st_CovarianceMatrixInputs_76X/";
    inpVer=_verMu76X;
    inpVerTag="_76X";
    lumi42=865.919;
    lumiTot=2832.673;
  }
  else if (0) {
    srcPath="/mnt/sdb/andriusj/v20160915_CovInput_ApprovedResults/";
    inpVer=_verMuApproved;
    inpVerTag=versionName(inpVer);
    lumi42=865.919;
    lumiTot=2832.673;

    fname1Base="ROOTFile_Histograms_Data.root";
    fname2Base="ROOTFile_Histograms_Bkg.root";
    fname6Base="ROOTFile_Input6_CrossCheck.root";
  }
  else if (1) {
    srcPath="/media/sf_CMSData/DY13TeV-CovInputs/v20170504_Input_Cov_mm/";
    inpVer=_verMuMay2017;
    inpVer=_verMuNov2017;
    inpVerTag=versionName(inpVer);
    lumi42=843.404;
    lumiTot=2759.017;

    fname1Base="ROOTFile_Histograms_Data.root";
    fname2Base="ROOTFile_Histograms_Bkg.root";
    fname6Base="ROOTFile_Input6_CrossCheck.root";
  }

  int checkBackgrounds=0;
  if ((checkStep==3) || (checkStep==4)) checkBackgrounds=1;

  TH1D* h1Dummy=NULL;
  //TH2D* h2Dummy=NULL;

  double lumi43= lumiTot-lumi42;
  std::cout << "lumi42=" << Form("%7.3lf",lumi42) << ", lumi43=" << Form("%7.3lf",lumi43) << "\n";



  // --------- Load corrections

  TString fname6= srcPath + fname6Base;
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

  if (0 || (checkStep==1)) {
    TH1D *h1effAcc_check= cloneHisto(h1Eff,"h1effAcc_check","h1effAcc_check");
    for (int ibin=1; ibin<=h1Acc->GetNbinsX(); ibin++) {
      double a= h1Eff->GetBinContent(ibin);
      double da= h1Eff->GetBinError(ibin);
      double b= h1Acc->GetBinContent(ibin);
      double db= h1Acc->GetBinError(ibin);
      double val= a * b;
      double unc= sqrt( b*b * da*da + a*a * db*db );
      h1effAcc_check->SetBinContent( ibin, val );
      h1effAcc_check->SetBinError( ibin, unc );

      if (1) {
	std::cout << "ibin=" << ibin << ", eff=" << a << " +- " << da
		  << ", acc=" << b << " +- " << db
		  << ", (eff x acc)=" << val << " +- " << unc << ", KL value ="
		  << h1EffAcc->GetBinContent(ibin) << " +- " << h1EffAcc->GetBinError(ibin)
		  << "\n";
	if (1) {
	  std::cout << "  relEff=" << da/a << ", relAcc=" << db/b << ", rel effxAcc=" << unc/val << "; rel KL=" << h1EffAcc->GetBinError(ibin)/h1EffAcc->GetBinContent(ibin) << "\n";
	}
	if (1) {
	  std::cout << "  eff =" << grEff->GetY()[ibin-1] << " +" << grEff->GetEYhigh()[ibin-1] << " -" << grEff->GetEYlow()[ibin-1] << ", acc=" << grAcc->GetY()[ibin-1] << " +" << grAcc->GetEYhigh()[ibin-1] << " -" << grAcc->GetEYlow()[ibin-1] << ", eff x Acc = "<< grAccEff->GetY()[ibin-1] << " +" << grAccEff->GetEYhigh()[ibin-1] << " -" << grAccEff->GetEYlow()[ibin-1] << "\n";
	}
      }
    }
    plotHisto(h1EffAcc, "cEffAcc", 1,1, "LPE1", "origEffAcc");
    plotHistoSame(h1effAcc_check, "cEffAcc", "LPE", "effAcc_check");
    printRatio(h1EffAcc,h1effAcc_check);
    std::cout << "checkStep induced stop\n";
    return;
  }

  // ---------------- load the observed and background-subracted (signal) yield
  TString fname1= srcPath + fname1Base;
  TFile fin1(fname1);
  if (!fin1.IsOpen()) {
    std::cout << "failed to open the file <" << fin1.GetName() << ">\n";
    return;
  }
  TH1D* h1Yield_42= loadHisto(fin1, "h_data_HLTv4p2", "h1Yield_42",1,h1Dummy);
  TH1D* h1Yield_43= loadHisto(fin1, "h_data_HLTv4p3", "h1Yield_43",1,h1Dummy);
  if (!h1Yield_42) return;
  if (!h1Yield_43) return;
  //printHisto(h1Yield_42);
  //plotHisto(h1Yield_42,"cYield42",1,1);
  //plotHisto(h1Yield_43,"cYield43",1,1);

  TH1D* h1Signal_42= loadHisto(fin1,"h_yield_HLTv4p2","h1Signal_42",1,h1Dummy);
  TH1D* h1Signal_43= loadHisto(fin1,"h_yield_HLTv4p3","h1Signal_43",1,h1Dummy);
  if (!h1Signal_42) return;
  if (!h1Signal_43) return;
  //plotHisto(h1Signal_42,"cSignal42",1,1);
  //plotHisto(h1Signal_43,"cSignal43",1,1);
  fin1.Close();

  // --------- Load backgrounds
  TString fname2= srcPath + fname2Base;
  TFile fin2(fname2);
  if (!fin2.IsOpen()) {
    std::cout << "failed to open the file <" << fin2.GetName() << ">\n";
    return;
  }
  TH1D* h1BkgTot= cloneHisto(h1Yield_42, "h1BkgTot","h1BkgTot");
  h1BkgTot->Reset();
  //printHisto(h1BkgTot);
  std::vector<TH1D*> hBkgV;
  std::vector<double> bkgWeightUnc;
  hBkgV.reserve(_bkgLast);
  bkgWeightUnc.reserve(_bkgLast);
  for (TBkg_t b=_bkgZZ; b<_bkgLast; next(b)) {
    std::cout << "b=" << bkgName(b) << "\n";
    TString histoBkgName="h_"+bkgName(b);
    TH1D *h1=NULL;
    if (((inpVer==_verMuMay2017) || (inpVer==_verMuNov2017))
	&& ((b==_bkgZZ) || (b==_bkgWZ))) {
      TString tmpName="h_"+bkgName(b) + "_HLTv4p2";
      TH1D *h1a= loadHisto(fin2, tmpName, tmpName,1,h1Dummy);
      tmpName.ReplaceAll("4p2","4p3");
      TH1D *h1b= loadHisto(fin2, tmpName, tmpName,1,h1Dummy);
      if (!h1a || !h1b) return;
      if (0) {
	printHisto(h1a);
	printHisto(h1b);
	return;
      }
      TH1D *h1sf42= cloneHisto(h1SF42,"h1sf42"+bkgName(b),"sf42");
      removeError(h1sf42);
      TH1D *h1sf43= cloneHisto(h1SF43,"h1sf43"+bkgName(b),"sf43");
      removeError(h1sf43);
      h1a->Divide(h1sf42);
      h1a->Scale(lumiTot/lumi42);
      h1b->Divide(h1sf43);
      h1b->Scale(lumiTot/lumi43);
      h1=cloneHisto(h1b,"h_"+bkgName(b),"h_"+bkgName(b));
      if (1)
      for (int ibin=1; ibin<=h1->GetNbinsX(); ibin++) {
	h1->SetBinContent(ibin,
		  0.5*(h1a->GetBinContent(ibin)+h1b->GetBinContent(ibin)));
      }
      if (checkStep==2) {
	hsBlue.SetStyle(h1a);
	hsGreen.SetStyle(h1b);
	hsRed.SetStyle(h1);
	TString cNameChk="cBkg_"+bkgName(b);
	plotHisto(h1a,cNameChk,1,1,"LPE1","HLTv4p2");
	plotHistoSame(h1b,cNameChk,"LPE1","HLTv4p3");
	plotHistoSame(h1,cNameChk,"LPE1","avg");
	std::cout << "remember that HLTv4p2 and 4p3 are (un)scaled by SF\n";
      }
    }
    else {
      h1= loadHisto(fin2, "h_"+bkgName(b), "h_"+bkgName(b),1,h1Dummy);
    }
    if (!h1) return;

    int hasNegative=0;
    for (int ibin=1; ibin<=h1->GetNbinsX(); ibin++) {
      if (h1->GetBinContent(ibin)<0) {
	std::cout << "histo " << h1->GetName() << " has negative entries\n";
	break;
      }
    }
    if (hasNegative) {
      printHisto(h1);
    }

    double w=1., dw=0.;
    if (inpVer==_verMu1) {
      if (b==_bkgZZ) { w= 15.4/996944.; dw=0.11; }
      else if (b==_bkgWZ) { w= 66.1/978512.; dw= 0.40; }
      else if (b==_bkgWW) { w= 118.7/993640.; dw= 0.08; }
      if (w!=double(1.)) h1->Scale( w * lumiTot );
    }
    else if (inpVer==_verMu76X) {
      switch(b) {
      case _bkgZZ: dw=0.11; break;
      case _bkgWZ: dw=0.40; break;
      case _bkgWW: dw=0.08; break;
      default: dw=0.;
      }
    }
    else {
      std::cout << "bkg cross section uncertainty is not defined for gInpVer="
		<< int(inpVer) << " (" << versionName(inpVer) << ")\n";
    }
    bkgWeightUnc.push_back(dw);
    hBkgV.push_back(h1);
    printBin(h1BkgTot,1,0); std::cout << " + "; printBin(h1,1,1);
    h1BkgTot->Add(h1);
    std::cout << "  = "; printBin(h1BkgTot,1,1);
  }
  fin2.Close();
  plotHisto(h1BkgTot, "cBkgTot",1,1,"histo","h1BkgTot");

  if (checkStep==2) {
    if (1) {
      TH1D *h1_bkg_4p2= cloneHisto(h1BkgTot,"h1_bkg_4p2","h1_bkg_4p2");
      h1_bkg_4p2->Scale(lumi42/lumiTot);
      TH1D *h1_bkg_4p3= cloneHisto(h1BkgTot,"h1_bkg_4p3","h1_bkg_4p3");
      h1_bkg_4p3->Scale(lumi43/lumiTot);
      printRatio(h1_bkg_4p2,h1_bkg_4p3);
    }

    std::cout << "checkStep induced stop\n";
    return;
  }


  if (checkBackgrounds) {
    TH1D* h1Bkg42_chk= cloneHisto(h1Yield_42, "h1Bkg42_chk", "h1Bkg42_chk");
    h1Bkg42_chk->Add(h1Signal_42, -1);
    TH1D* h1Bkg43_chk= cloneHisto(h1Yield_43, "h1Bkg43_chk", "h1Bkg43_chk");
    h1Bkg43_chk->Add(h1Signal_43, -1);
    hsBlue.SetStyle(h1Bkg42_chk);
    hsGreen.SetStyle(h1Bkg43_chk);
    plotHistoSame(h1Bkg42_chk, "cBkgTot","hist","h1Bkg42_chk");
    TCanvas *ctmp=plotHistoSame(h1Bkg43_chk, "cBkgTot","hist","h1Bkg43_chk");
    ctmp->SetTitle("cBkgTot - compare shapes of bkgs");

    std::cout << "compare the backgrounds\n";
    printRatio(h1Bkg42_chk,h1Bkg43_chk);

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
    h1ratio42->SetTitle("ratio42 = (Nyield-Nsig)/NtotBkg");
    plotHisto(h1ratio42, "cRatio",1,0,"LPE","h1ratio42");
    ctmp=plotHistoSame(h1ratio43, "cRatio", "hist","h1ratio43");
    ctmp->SetTitle("cRatio - ratio of backgrounds - expectation is 1");
    //printHisto(h1ratio42);
    //printHisto(h1ratio43);
    //TH1D *h1ratio43_to_42= cloneHisto(h1ratio43, "h1ratio43_to_42","ratio43_to_42");
    //h1ratio43_to_42->Divide(h1ratio42);
    //printHisto(h1ratio43_to_42);
    printRatio(h1ratio43,h1ratio42);
    std::cout << "remember that some HLTv4p2/4p3 bkgs are (un)scaled by SF\n";

    if (0 || (checkStep==4)) {
      TH1D *h1BkgSum= cloneHisto(h1Bkg42_chk, "h1BkgSum", "BkgSum");
      h1BkgSum->Add(h1Bkg43_chk);
      TH1D *h1Bkg_ChkDiff= cloneHisto(h1BkgSum, "h1Bkg_ChkDiff", "Bkg_ChkDiff");
      h1Bkg_ChkDiff->Add(h1BkgTot, -1);
      removeError(h1Bkg_ChkDiff);
      h1Bkg_ChkDiff->GetYaxis()->SetRangeUser(-100,100);
      ctmp=plotHisto(h1Bkg_ChkDiff, "cBkg_ChkDiff",1,0);
      ctmp->SetTitle("cBkg_ChkDiff - difference in backgrounds - expectation is 0");
    }

    std::cout << "checkStep induced stop\n";
    return;
  }

  if (inpVer==_verMuApproved) {
    std::cout << "\n\n hacking the backgrounds to match inputs\n\n";
    HERE("");
    for (int i=0; i<2; i++) {
      TBkg_t iBkg=(i==0) ? _bkgZZ : _bkgWZ;
      if (!hBkgV[iBkg]) HERE("null bkg ptr");
      scaleBin(hBkgV[iBkg],15, 1214.56/1213.38);
      scaleBin(hBkgV[iBkg],16, 0.5*(2545.03/2551.31 + 2554.9/2563.36) );
      scaleBin(hBkgV[iBkg],17, 0.5*(2417.55/2423.37 + 2426.63/2434.19) );
      scaleBin(hBkgV[iBkg],43, 0.264957/0.2640006);
    }
  }

  // Bundle to one object


  MuonCrossSection_t muCS("muCS",inpVerTag,lumi42,lumi43,inpVer);
  muCS.editCSa().h1Yield( h1Yield_42 );
  muCS.editCSa().h1Eff( h1Eff );
  muCS.editCSa().h1Acc( h1Acc );
  muCS.editCSa().h1EffAcc( h1EffAcc );
  muCS.editCSa().h1Rho( h1SF42 );
  //rooUnfDetRes->UseOverflow(false);
  //rooFSRRes->UseOverflow(false);
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
  muCS.setBkgV(hBkgV, bkgWeightUnc);
  muCS.h1Theory(h1_DiffXSec_data);

  if (1) {
    std::cout << "Number of iterations: detRes " << muCS.csA().nItersDetRes()
	      << ", FSR " << muCS.csA().nItersFSR() << "\n";
    TH1D *h1preFSRcs= muCS.calcCrossSection(removeNegativeSignal);
    if (!h1preFSRcs) return;
    plotHisto(h1preFSRcs, "cPreFsrCS",1,1);
    h1_DiffXSec_data->SetMarkerStyle(24);
    h1_DiffXSec_data->SetLineColor(kGreen);
    h1_DiffXSec_data->SetMarkerColor(kGreen);
    plotHistoSame(h1_DiffXSec_data, "cPreFsrCS", "LPE1");
    printRatio(h1preFSRcs, h1_DiffXSec_data);
  }
  else {
    TCanvas *c= muCS.plotCrossSection("cs",0,removeNegativeSignal);
    if (!c) return;
  }

  //if (1) {
  //  printHisto(muCS.csA().h1UnfRhoEffAccCorr());
  //  printHisto(muCS.csB().h1UnfRhoEffAccCorr());
  //}

  if (0) {
    h1Signal_42->SetMarkerStyle(24);
    plotHisto(h1Signal_42,"cSignal42",1,1);
    plotHistoSame(muCS.csA().copy(muCS.csA().h1Signal()),"cSignal42","LPE");
    printRatio(h1Signal_42, muCS.csA().h1Signal());

    h1Signal_43->SetMarkerStyle(24);
    plotHisto(h1Signal_43,"cSignal43",1,1);
    plotHistoSame(muCS.csB().copy(muCS.csB().h1Signal()),"cSignal43","LPE");
    printRatio(h1Signal_43, muCS.csB().h1Signal());
  }

  if (0) {
    //printHisto(muCS.csB().h1UnfRhoEffAccCorr());
    //printHisto(muCS.csB().h1UnfRhoEffAccCorr());
    //printHisto((TH1D*)rooFSRRes->Hmeasured());

    TH1D *postFSR_a=cloneHisto(muCS.csA().h1UnfRhoEffAccCorr(), "postFSR_a","postFSR_A");
    postFSR_a->Scale(muCS.csA().lumi());
    TH1D *postFSR_b=cloneHisto(muCS.csB().h1UnfRhoEffAccCorr(), "postFSR_b","postFSR_B");
    postFSR_b->Scale(muCS.csB().lumi());
    printHisto(postFSR_a);
    printHisto(postFSR_b);
  }

  if (!doSave) { HERE("not saving"); return; }
  HERE("calling muCS.save()");
  if (!muCS.save("cs_DYmm_13TeV" + inpVerTag)) {
    std::cout << "saving failed\n";
    return;
  }
  HERE("saving done");
}

// ---------------------------------------------------------
// ---------------------------------------------------------

/*
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
  //factor=1.05; // to see the difference
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
*/
// ---------------------------------------------------------
