#include "inputs.h"
#include "crossSection.h"

void studyDYmmCScovYieldMethods(int nSample=2000)
{
  const int nMethods=3;
  const TString methodName[nMethods] = { "_allRnd", "_keepYieldWithPosSignal",
					 "_nulifyNegSig" };

  std::vector<TH2D*> h2covV;
  std::vector<TH1D*> h1uncV;
  TString fnameBase=Form("cov_mumu_varYield_%d",nSample);
  if (nSample==5000) fnameBase.Append("_slim");
  for (int i=0; i<nMethods; i++) {
    TString fname= fnameBase + methodName[i] + ".root";
    TString h2NameFile= "h2cov_" + methodName[i];
    TH2D *h2= loadHisto(fname,h2NameFile,h2NameFile,h2dummy);
    if (!h2) return;
    h2covV.push_back(h2);
    TH1D *h1= uncFromCov(h2,NULL,NULL);
    hsVec[i].SetStyle(h1);
    h1uncV.push_back(h1);
  }

  TH1D *h1unc_note= loadHisto("dymm_unc_fromNote.root","h1csMM_statUnc","h1csMM_statUnc",h1dummy);
  if (!h1unc_note) return;
  h1unc_note->SetStats(0);

  hsColor46.SetStyle(h1unc_note);
  hsGreen.SetStyle(h1uncV[2]);

  TString canvName="cCovYieldUnc";
  TString csLabel= "(\\delta" + niceMassAxisLabel(1,"",1);
  csLabel.ReplaceAll("\\mu}","\\mu})_{exp.unc.}");

  logAxis(h1unc_note,1+2,niceMassAxisLabel(1,"",0),csLabel,
	  h1unc_note->GetTitle() + TString(Form(", nSample=%d",nSample)));
  h1unc_note->GetYaxis()->SetNoExponent(false);
  h1unc_note->GetYaxis()->SetTitleOffset(1.5);
  h1unc_note->GetYaxis()->SetRangeUser(1e-7,5.);
  TCanvas *c=plotHisto(h1unc_note,canvName,1,1,"LP","from AN Note");
  setLeftRightMargins(c,0.15,0.15);
  moveLegend(c,0.02,0.);
  increaseLegend(c,0.04,0.);
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

  h1ratioV[0]->GetYaxis()->SetRangeUser(0.7,2.6);
  logAxis(h1ratioV[0],1,niceMassAxisLabel(1,"",0),"ratio",
	  h1ratioV[0]->GetTitle() + TString(Form(", nSample=%d",nSample)));
  TCanvas *cr=plotRatio(h1ratioV[0],"cUncRatio",1,0,"LP",methodName[0]);
  moveLegend(cr,0.0,0.56);
  increaseLegend(cr,0.04,0.03);
  for (unsigned int i=1; i<h1ratioV.size(); i++) {
    plotHistoSame(h1ratioV[i],"cUncRatio","LP",methodName[i]);
  }
}
