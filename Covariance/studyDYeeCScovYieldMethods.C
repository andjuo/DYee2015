#include "inputs.h"
#include "crossSection.h"

void studyDYeeCScovYieldMethods(int nSample=2000)
{
  const int nMethods=4;
  const TString methodName[nMethods] = { "_allRnd", "_keepYieldWithPosSignal",
					 "_nulifyNegSig",
					 "_keepYieldWithPosSignal2-bkgX1p3"};

  std::vector<TH2D*> h2covV;
  std::vector<TH1D*> h1uncV;
  TString fnameBase=Form("cov_ee_varYield_%d",nSample);
  if (nSample==2000) fnameBase.ReplaceAll("2000","2000_slim");
  for (int i=0; i<nMethods; i++) {
    TString fname= fnameBase + methodName[i] + ".root";
    TString h2NameFile= "h2cov_" + methodName[i];
    TString h2NameMem= "h2cov_" + methodName[i];
    if (methodName[i]=="_keepYieldWithPosSignal2-bkgX1p3") {
      h2NameFile.ReplaceAll("-bkgX1p3","");
    }
    TH2D *h2= loadHisto(fname,h2NameFile,h2NameMem,h2dummy);
    if (!h2) return;
    h2covV.push_back(h2);
    TH1D *h1= uncFromCov(h2,NULL,NULL);
    hsVec[i].SetStyle(h1);
    h1uncV.push_back(h1);
  }

  TH1D *h1unc_note= loadHisto("dyee_unc_fromNote.root","h1csEE_statUnc","h1csEE_statUnc",h1dummy);
  if (!h1unc_note) return;
  h1unc_note->SetStats(0);
  printHisto(h1unc_note);

  hsColor46.SetStyle(h1unc_note);
  hsGreen.SetStyle(h1uncV[2]);
  if (h1uncV.size()>3) hsViolet.SetStyle(h1uncV[3]);

  TString canvName="cCovYieldUnc";
  TString csLabel= "(\\delta" + niceMassAxisLabel(0,"",1);
  csLabel.ReplaceAll("ee}","ee})_{exp.unc.}");
  
  logAxis(h1unc_note,1+2,niceMassAxisLabel(0,"",0),csLabel,"");
  h1unc_note->GetYaxis()->SetNoExponent(false);
  h1unc_note->GetYaxis()->SetTitleOffset(1.5);
  h1unc_note->GetYaxis()->SetRangeUser(1e-7,10.);
  TCanvas *c=plotHisto(h1unc_note,canvName,1,1,"LP","from AN Note");
  setLeftRightMargins(c,0.15,0.15);
  moveLegend(c,0.02,0.);
  increaseLegend(c,0.13,-0.02);
  for (unsigned int i=0; i<h1uncV.size(); i++) {
    plotHistoSame(h1uncV[i],canvName,"LP",methodName[i]);
  }

  std::vector<TH1D*> h1ratioV;
  for (unsigned int i=0; i<h1uncV.size(); i++) {
    TH1D *h1= cloneHisto(h1uncV[i],h1uncV[i]->GetName() + TString("_div_AN"),
			 h1uncV[i]->GetTitle());
    h1->Divide(h1,h1unc_note);
    removeError(h1);
    h1ratioV.push_back(h1);
  }

  logAxis(h1ratioV[0],1,niceMassAxisLabel(0,"",0),"ratio","");
  //if (nSample==5000) 
    h1ratioV[0]->GetYaxis()->SetRangeUser(0.2,6);
  TCanvas *cr=plotRatio(h1ratioV[0],"cUncRatio",1,0,"LP",methodName[0]);
  moveLegend(cr,0.0,0.48);
  increaseLegend(cr,0.13,0.02);
  for (unsigned int i=1; i<h1ratioV.size(); i++) {
    plotHistoSame(h1ratioV[i],"cUncRatio","LP",methodName[i]);
  }
}
