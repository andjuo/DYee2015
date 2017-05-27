#include "inputs.h"


void studyDYmm_KL_varYield(int iMax=-1, int save=0)
{
  TString fname="DYmm-ROOTFile_RelStatUnc-KL-20170525.root";
  TString histoBaseNameA="h_DiffXSec_Smeared_";
  TString histoBaseNameB="h_FpoF_DiffXSec_Smeared_";
  const int nHistos=1000;
  std::vector<TH1D*> h1aV, h1bV;

  TH1D *h1relStatUncA=NULL, *h1relStatUncB=NULL;

  TFile fin(fname);
  if (!fin.IsOpen()) {
    std::cout << "failed to open the file\n";
    return;
  }

  h1relStatUncA= loadHisto(fin, "h_RelStatUnc", "h1relStatUncA",1,h1dummy);
  h1relStatUncB= loadHisto(fin, "h_FpoF_RelStatUnc", "h1relStatUncB",1,h1dummy);
  if (!h1relStatUncA || !h1relStatUncB) return;

  for (int i=0; i<nHistos; i++) {
    TString hNameA= histoBaseNameA + makeNumberStr(i,3);
    TString hNameB= histoBaseNameB + makeNumberStr(i,3);
    TH1D *h1A= loadHisto(fin, hNameA, hNameA, 1, h1dummy);
    TH1D *h1B= loadHisto(fin, hNameB, hNameB, 1, h1dummy);
    if (!h1A || !h1B) return;
    h1A->SetLineColor(i%9+1);
    h1B->SetLineColor(i%9+1);
    h1aV.push_back(h1A);
    h1bV.push_back(h1B);

    if (i==0) {
      printHisto(h1A);
      printHisto(h1B);
    }

    if (1 && (i<20)) {
      plotHistoAuto(h1A,"c"+histoBaseNameA,1,1,"hist","",3);
      plotHistoAuto(h1B,"c"+histoBaseNameB,1,1,"hist","",3);
    }
    if ((iMax>0) && (i>=iMax)) break;
  }
  fin.Close();

  TH1D *h1avgA=NULL, *h1avgB=NULL;
  TH2D *h2covA=NULL, *h2covB=NULL;
  TH2D *h2corrA=NULL, *h2corrB=NULL;
  deriveCovariance(h1aV, histoBaseNameA, "cov" + histoBaseNameA,
		   &h1avgA,&h2covA);
  deriveCovariance(h1bV, histoBaseNameB, "cov" + histoBaseNameB,
		   &h1avgB,&h2covB);

  PlotCovCorrOpt_t opt;
  plotCovCorr(h2covA,"cCov"+histoBaseNameA,opt,&h2corrA);
  plotCovCorr(h2covB,"cCov"+histoBaseNameB,opt,&h2corrB);

  TH1D *h1relUncA= uncFromCov(h2covA, h1avgA);
  TH1D *h1relUncB= uncFromCov(h2covB, h1avgB);
  h1relStatUncA->SetStats(0);
  h1relUncA->SetStats(0);

  hsBlue.SetStyle(h1relStatUncA);
  hsBlue2.SetStyle(h1relStatUncB);
  h1relStatUncB->SetMarkerStyle(2);
  hsGreen.SetStyle(h1relUncA);
  hsRed.SetStyle(h1relUncB);

  TCanvas *cx=plotHisto(h1relStatUncA,"cRelUnc",1,1,"LP","KL A");
  plotHistoAuto(h1relUncA,"cRelUnc",1,1,"LP","calc A");
  plotHistoAuto(h1relStatUncB,"cRelUnc",1,1,"LP","KL B");
  plotHistoAuto(h1relUncB,"cRelUnc",1,1,"LP","calc B");
  moveLegend(cx,0.45,0);

  if (1) {
    TString myFName="cov_mumu_varYieldPoisson_1000_veryslim.root";
    TH1D *h1myXS= loadHisto(myFName,"h1xsecAvg","h1xsecAvg",1,h1dummy);
    TH2D *h2myCov= loadHisto(myFName,"h2Cov","h2covVarYield",1,h2dummy);
    if (!h1myXS || !h2myCov) return;
    plotCovCorr(h2myCov,"cCov_varYieldPoisson",opt);
    TH1D *h1myRelUnc= uncFromCov(h2myCov, h1myXS);
    HistoStyle_t(kViolet,5).SetStyle(h1myRelUnc);
    plotHistoAuto(h1myRelUnc,"cRelUnc",1,1,"LP","varPoisson");
  }

  if (save) {
    TFile fout("cmp_KL_varYield.root","RECREATE");
    SaveCanvases("ALL","",&fout);
    writeTimeTag(&fout);
    fout.Close();
    std::cout << "file <" << fout.GetName() << "> created\n";
  }
}

