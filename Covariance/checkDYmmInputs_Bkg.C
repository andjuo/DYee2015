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

  TFile fin(fname);
  TH1D *h1a= loadHisto(fin,hNameBase,hNameBase,1,h1dummy);
  hNameBase.ReplaceAll("v4p2","v4p3");
  TH1D *h1b= loadHisto(fin,hNameBase,hNameBase,1,h1dummy);
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
  
}
