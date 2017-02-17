#include "inputs.h"

void compareDYllRatios()
{
  TCanvas *cCorr= loadCanvas("dyll-combi-_corr_plotChCov.root","cCombiCS2TheoryRatio","cCombiCS2TheoryRatio_corr");
  //TCanvas *cUnCorr= loadCanvas("dyll-combi-_uncorr_plotChCov.root","cCombiCS2TheoryRatio","cCombiCS2TheoryRatio_uncorr");
  //TCanvas *cUnCorr= loadCanvas("dyll-combi-_corr_plotChCov_uncorrYieldUnc.root","cCombiCS2TheoryRatio","cCombiCS2TheoryRatio_uncorr");
  //TCanvas *cUnCorr= loadCanvas("dyll-combi-_corr_plotChCov-myRhoStat.root","cCombiCS2TheoryRatio","cCombiCS2TheoryRatio_corr-myRhoStat");
  TCanvas *cUnCorr= loadCanvas("dyll-combi-_corr_wLumi_plotChCov.root","cCombiCS2TheoryRatio","cCombiCS2TheoryRatio_wLumi_corr");
  if (!cCorr || !cUnCorr) return;

  cCorr->Draw();
  cUnCorr->Draw();

  std::vector<TH1D*> h1Vcorr, h1Vuncorr;

  if (!getHistosFromCanvas(cCorr,&h1Vcorr,NULL) ||
      !getHistosFromCanvas(cUnCorr,&h1Vuncorr,NULL)) return;

  std::cout << "got " << h1Vcorr.size() << " (corr) and "
	    << h1Vuncorr.size() << " (uncorr) histos\n";
  if ((h1Vcorr.size()!=1) || (h1Vuncorr.size()!=1)) return;


  TH1D *h1corr= cloneHisto(h1Vcorr[0],"h1corr","h1corr");
  TH1D *h1uncorr= h1Vuncorr[0];

  hsRed.SetStyle(h1corr);
  h1corr->SetTitle("");
  h1corr->GetXaxis()->SetTitleSize(0.04);
  h1corr->GetXaxis()->SetLabelSize(0.04);
  h1corr->GetYaxis()->SetTitleSize(0.04);
  h1corr->GetYaxis()->SetLabelSize(0.04);
  h1corr->GetYaxis()->SetTitleOffset(1.5);

  logAxis(h1corr,1);
  TCanvas *c=plotHisto(h1corr,"cCombiRatioCmp",1,0,"LPE1","correlated");
  plotHistoSame(h1uncorr,"cCombiRatioCmp","LPE1","uncorrelated");
  increaseLegend(c,0.,0.02);
  moveLegend(c,0.05,0.);
  setLeftMargin(c,0.15);
  
  return;
}
