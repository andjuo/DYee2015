#include "inputs.h"


void studyDYmm_KL_varRho(int iMax=-1, int save=0)
{
  TString fname="/media/sf_CMSData/DY13TeV-CovInputs/v20170504_Input_Cov_mm/ROOTFile_Outputs_SysUncTool_EffCorr-20170602.root";
  TString histoBaseNameA="h_DiffXsec_Smeared_";
  //TString histoBaseNameB="h_FpoF_DiffXSec_Smeared_";
  const int nHistos=500;
  std::vector<TH1D*> h1aV;

  TH1D *h1xsCV=NULL;

  TFile fin(fname);
  if (!fin.IsOpen()) {
    std::cout << "failed to open the file\n";
    return;
  }

  h1xsCV= loadHisto(fin, "h_DiffXsec_CV", "h1xsCV",1,h1dummy);
  if (!h1xsCV) return;
  hsRed.SetStyle(h1xsCV);
  plotHistoAuto(h1xsCV,"c" + histoBaseNameA,1,1,"hist","",3);

  for (int i=0; i<nHistos; i++) {
    TString hNameA= histoBaseNameA + Form("%d",i);
    TH1D *h1A= loadHisto(fin, hNameA, hNameA, 1, h1dummy);
    if (!h1A) return;
    h1A->SetLineColor(i%9+1);
    h1aV.push_back(h1A);

    if (i==0) {
      printHisto(h1A);
    }

    if (1 && (i<20)) {
      plotHistoAuto(h1A,"c"+histoBaseNameA,1,1,"hist","",3);
    }
    if ((iMax>0) && (i>=iMax)) break;
  }
  fin.Close();


  TH1D *h1avgCS=NULL;
  TH2D *h2covCS=NULL;
  TH2D *h2corrCS=NULL;
  deriveCovariance(h1aV, histoBaseNameA, "cov" + histoBaseNameA,
		   &h1avgCS,&h2covCS);

  PlotCovCorrOpt_t opt;
  plotCovCorr(h2covCS,"cCov"+histoBaseNameA,opt,&h2corrCS);

  plotHisto(h1xsCV,"cXS",1,1,"LPE","KL");
  plotHistoSame(h1avgCS,"cXS","LPE","my");
  printRatio(h1xsCV,h1avgCS);

  TH1D *h1relUncA= uncFromCov(h2covCS, h1avgCS);
  TH1D *h1uncA= uncFromCov(h2covCS, NULL);
  h1relUncA->SetStats(0);
  h1uncA->SetStats(0);

  TH1D *h1uncKL= loadHisto("/media/sf_CMSData/DY13TeV-CovInputs/v20170504_Input_Cov_mm/ROOTFile_SysUnc_EffSF_TagProbe_UnbinnedFit-20170602.root","h_RMS","h1uncKL",1,h1dummy);
  if (!h1uncKL) return;

  hsBlue.SetStyle(h1uncKL);
  hsGreen.SetStyle(h1relUncA);
  hsRed.SetStyle(h1uncA);

  TCanvas *cxRel=plotHisto(h1relUncA,"cRelUnc",1,1,"LP","my");
  plotHistoAuto(h1uncKL,"cRelUnc",1,1,"LP","KL");
  moveLegend(cxRel,0.45,0);

  TCanvas *cx=plotHisto(h1uncA,"cUnc",1,1,"LP","my");
  plotHistoAuto(h1uncKL,"cUnc",1,1,"LP","KL");
  moveLegend(cx,0.45,0);
 
  if (1) {
    int checkKL=1;
    TString myFName="cov_mumu_varRhoFile_500.root";
    if (checkKL) myFName="cov_mumu_varRhoFile_500_KL.root";
    TH1D *h1myXS= loadHisto(myFName,"h1xsecAvg","h1xsecAvg",1,h1dummy);
    TH2D *h2myCov= loadHisto(myFName,"h2Cov","h2covVarRho",1,h2dummy);
    if (!h1myXS || !h2myCov) return;
    plotCovCorr(h2myCov,"cCov_varRho",opt);
    TH1D *h1myRelUnc= uncFromCov(h2myCov, h1myXS);
    HistoStyle_t(kViolet,5).SetStyle(h1myRelUnc);
    plotHistoAuto(h1myRelUnc,"cRelUnc",1,1,"LP",
		  (checkKL) ? "varRho500KL" : "varRho2000");
  }

  /*
  if (save) {
    TFile fout("cmp_KL_varYield.root","RECREATE");
    SaveCanvases("ALL","",&fout);
    writeTimeTag(&fout);
    fout.Close();
    std::cout << "file <" << fout.GetName() << "> created\n";
  }
  */
}

