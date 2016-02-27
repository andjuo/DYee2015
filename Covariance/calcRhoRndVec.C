#include "inputs.h"
#include "DYbinning.h"
#include "DYmm13TeV_eff.h"

void calcRhoRndVec(int nSamples=100)
{

  //TH2D* h2PostrFsrEventSpace=NULL;
  DYTnPEff_t tnpEff;
  TH1D* h1postFsrInAccSel_MWPU=NULL; // will provide mass binning

  EventSpace_t esPostFsr;

  // Load data
  TString fname="dymm_test_RECO.root";
  TFile fin(fname);
  if (!fin.IsOpen()) {
    std::cout << "failed to open the file <" << fin.GetName() << ">\n";
    return;
  }
  h1postFsrInAccSel_MWPU=loadHisto(fin, "h1_postFsrInAccSel_MWPU","h1postFsrInAccSel_MW",1, h1dummy);
  if (!esPostFsr.load(fin, "mainES")) return;
  if (!tnpEff.load(fin,"DYTnPEff")) {
    std::cout << "failed to load efficiencies\n";
    return;
  }
  fin.Close();

  // check the event space
  if (0) {
    TH1D* h1postFsrCheck= cloneHisto(h1postFsrInAccSel_MWPU,"h1postFsrCheck","h1postFsrCheck");
    h1postFsrCheck->Reset();

    for (int im=0; im<DYtools::nMassBins; im++) {
      const TH2D* h2= esPostFsr.h2ES(im);
      double sum=0, sumErr=0;
      for (int ibin=1; ibin<=h2->GetNbinsX(); ibin++) {
	for (int jbin=1; jbin<=h2->GetNbinsY(); jbin++) {
	  sum+= h2->GetBinContent(ibin,jbin);
	  sumErr+= pow(h2->GetBinError(ibin,jbin),2);
	}
      }
      sumErr=sqrt(sumErr);
      h1postFsrCheck->SetBinContent(im+1, sum);
      h1postFsrCheck->SetBinError  (im+1, sumErr);
    }

    plotHisto(h1postFsrInAccSel_MWPU,"cSel",1,1,"LPE");
    histoStyle(h1postFsrCheck,kBlue,24);
    plotHistoSame(h1postFsrCheck,"cSel","LPE");
    printRatio(h1postFsrInAccSel_MWPU,h1postFsrCheck);
  }

  // calculate event scale factors
  std::cout << "calculate h1rho_4p2\n";
  TH1D* h1rho_4p2= esPostFsr.calculateScaleFactor(tnpEff,0,"h1rho_4p2","scale factor 4p2");
  std::cout << "calculate h1rho_4p3\n";
  TH1D* h1rho_4p3= esPostFsr.calculateScaleFactor(tnpEff,1,"h1rho_4p3","scale factor 4p3");

  plotHisto(h1rho_4p2,"c4p2",1,0,"LPE");
  plotHisto(h1rho_4p3,"c4p3",1,0,"LPE");

  if (1) {
    TString fnameChk="../../v20160214_1st_CovarianceMatrixInputs/Input6/ROOTFile_Input6_CrossCheck.root";
    TFile finChk(fnameChk);
    if (!finChk.IsOpen()) {
      std::cout << "file <" << finChk.GetName() << "> could not be found\n";
      return;
    }
    TGraphAsymmErrors *gr4p2= (TGraphAsymmErrors*)finChk.Get("g_EffSF_HLTv4p2");
    TGraphAsymmErrors *gr4p3= (TGraphAsymmErrors*)finChk.Get("g_EffSF_HLTv4p3");
    finChk.Close();
    if (!gr4p2 || !gr4p3) {
      std::cout << "failed to load graphs\n";
      return;
    }
    TH1D *h1KP4p2= convert(gr4p2,"h1KP4p2","HLTv4p2 sf",0);
    TH1D *h1KP4p3= convert(gr4p3,"h1KP4p3","HLTv4p3 sf",0);
    histoStyle(h1KP4p2,kBlue,24);
    histoStyle(h1KP4p3,kBlue,24);
    plotHistoSame(h1KP4p2,"c4p2","LPE1");
    plotHistoSame(h1KP4p3,"c4p3","LPE1");
  }
}
