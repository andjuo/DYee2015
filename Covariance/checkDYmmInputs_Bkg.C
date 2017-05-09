#include "inputs.h"
#include "crossSection.h"

void checkDYmmInputs_Bkg(int iBkg=-1)
{
  TString fname="/media/sf_CMSData/DY13TeV-CovInputs/v20170504_Input_Cov/ROOTFile_Histograms_Bkg.root";
  TString hNameBase="h_ZZ_HLTv4p2";
  if (iBkg==1) hNameBase.ReplaceAll("ZZ","WZ");

  double lumi42=843.404;
  double lumiTot=2759.017;
  double lumi43= lumiTot-lumi42;

  TH1D *h1bkg4p2=NULL, *h1bkg4p3=NULL;

  TFile fin(fname);
  TH1D *h1a= loadHisto(fin,hNameBase,hNameBase,1,h1dummy);
  hNameBase.ReplaceAll("v4p2","v4p3");
  TH1D *h1b= loadHisto(fin,hNameBase,hNameBase,1,h1dummy);

  if (iBkg==111) {
    h1bkg4p2= cloneHisto(h1a,"h1bkg4p2","h1bkg4p2");
    h1bkg4p2->Reset();
    h1bkg4p3= cloneHisto(h1a,"h1bkg4p3","h1bkg4p3");
    h1bkg4p3->Reset();

    for (TBkg_t b=_bkgZZ; b<_bkgLast; next(b)) {
      std::cout << "b=" << bkgName(b) << "\n";
      TString histoBkgName="h_"+bkgName(b);
      TH1D *h1a=NULL, *h1b=NULL;
      TString tmpName="h_"+bkgName(b);
      if ((b==_bkgZZ) || (b==_bkgWZ)) {
	tmpName+="_HLTv4p2";
	h1a= loadHisto(fin, tmpName, tmpName,1,h1dummy);
	tmpName.ReplaceAll("4p2","4p3");
	h1b= loadHisto(fin, tmpName, tmpName,1,h1dummy);
      }
      else {
	h1a= loadHisto(fin,tmpName,tmpName+"_HLTv4p2",1,h1dummy);
	if (!h1a) return;
	h1b= cloneHisto(h1a,tmpName+"_HLTv4p3",tmpName+"_HLTv4p3");
	h1a->Scale(lumi42/lumiTot);
	h1b->Scale(lumi43/lumiTot);
      }
      if (!h1a || !h1b) return;
      h1bkg4p2->Add(h1a);
      h1bkg4p3->Add(h1b);
    }
  }

  fin.Close();

  TString fname6="/media/sf_CMSData/DY13TeV-CovInputs/v20170504_Input_Cov/ROOTFile_Input6_CrossCheck.root";
  TFile fin6(fname6);
  TGraphAsymmErrors *grEffSF_HLTv4p2=(TGraphAsymmErrors*)fin6.Get("g_EffSF_HLTv4p2");
  TGraphAsymmErrors *grEffSF_HLTv4p3=(TGraphAsymmErrors*)fin6.Get("g_EffSF_HLTv4p3");
  fin6.Close();
  int plotIt_sf=1;
  if (!grEffSF_HLTv4p2 || !grEffSF_HLTv4p3) return;
  TH1D *h1SF42= convert(grEffSF_HLTv4p2, "h1SF42","ScaleFactors 4.2;M_{#mu#mu} [GeV];eff.s.f. (v.4.2)", plotIt_sf);
  if (!h1SF42) return;
  TH1D *h1SF43= convert(grEffSF_HLTv4p3, "h1SF43","ScaleFactors 4.3;M_{#mu#mu} [GeV];eff.s.f. (v.4.3)", plotIt_sf);
  if (!h1SF43) return;
 
  TH1D *h1ax=cloneHisto(h1a,"h1a_scaled","h1a_scaled");
  h1ax->Divide(h1SF42);
  TH1D *h1bx=cloneHisto(h1b,"h1b_scaled","h1b_scaled");
  h1bx->Divide(h1SF43);
  h1bx->Scale(lumi42/lumi43);

  printRatio(h1a,h1b);
  printRatio(h1ax,h1bx);
  std::cout << "lumi42/lumi43=" << lumi42/lumi43 << "\n";


  plotHisto(h1a,"cCmp",1,1,"LPE1","v4p2");
  plotHistoSame(h1b,"cCmp","LPE1","v4p3");
  plotHisto(h1ax,"cCmpX",1,1,"LPE1","down-scaled by v4p2");
  plotHistoSame(h1bx,"cCmpX","LPE1","down-scaled by v4p3");

  std::cout << "\n\ninfo: hNameBase=" << hNameBase << "\n";

if (iBkg!=111) return;


  // ---------------- load the observed and background-subracted (signal) yield
  TString fname1= "/media/sf_CMSData/DY13TeV-CovInputs/v20170504_Input_Cov/ROOTFile_Histograms_Data.root";
  TFile fin1(fname1);
  if (!fin1.IsOpen()) {
    std::cout << "failed to open the file <" << fin1.GetName() << ">\n";
    return;
  }
  TH1D* h1Yield_42= loadHisto(fin1, "h_data_HLTv4p2", "h1Yield_42",1,h1dummy);
  TH1D* h1Yield_43= loadHisto(fin1, "h_data_HLTv4p3", "h1Yield_43",1,h1dummy);
  if (!h1Yield_42) return;
  if (!h1Yield_43) return;
  //printHisto(h1Yield_42);
  //plotHisto(h1Yield_42,"cYield42",1,1);
  //plotHisto(h1Yield_43,"cYield43",1,1);

  TH1D* h1Signal_42= loadHisto(fin1,"h_yield_HLTv4p2","h1Signal_42",1,h1dummy);
  TH1D* h1Signal_43= loadHisto(fin1,"h_yield_HLTv4p3","h1Signal_43",1,h1dummy);
  if (!h1Signal_42) return;
  if (!h1Signal_43) return;
  //plotHisto(h1Signal_42,"cSignal42",1,1);
  //plotHisto(h1Signal_43,"cSignal43",1,1);
  fin1.Close();

  // compare sum of backgrounds

  TH1D *h1totBkg4p2= cloneHisto(h1Yield_42,"h1totBkg4p2","h1totBkg4p2");
  h1totBkg4p2->Add(h1Signal_42,-1);
  TH1D *h1totBkg4p3= cloneHisto(h1Yield_43,"h1totBkg4p3","h1totBkg4p3");
  h1totBkg4p3->Add(h1Signal_43,-1);

  printRatio(h1bkg4p2,h1totBkg4p2);
  printRatio(h1bkg4p3,h1totBkg4p3);
}
