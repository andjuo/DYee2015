#include "crossSection.h"

// --------------------------------------------------------------

void work(MuonCrossSection_t &muCS, TVaried_t var, int nSample);

// --------------------------------------------------------------

void studyDYmmCS(TVaried_t var= _varNone, int nSample=10)
{
  MuonCrossSection_t muCS("muCS","v1");
  if (!muCS.load("cs_DYmm_13TeV_v1","v1")) {
    std::cout << "loading failed\n";
    return;
  }

  //TCanvas *ca=muCS.editCSa().plotCrossSection("csa");
  //TCanvas *cb=muCS.editCSb().plotCrossSection("csb");

  if (var==_varNone) {
    TCanvas *c= muCS.plotCrossSection("cs");
    if (!c) return;

    plotHisto(cloneHisto(muCS.h1Bkg(),"h1tmp1","tmp1"),"ctmp1",1,1,"hist");

    std::vector<double> bkgW;
    bkgW.reserve(muCS.bkgW().GetNoElements());
    for (int i=0; i<muCS.bkgW().GetNoElements(); i++) {
      bkgW.push_back(2.);
    }
    if (!muCS.recalcBkg(bkgW)) return;
    plotHistoSame(cloneHisto(muCS.h1Bkg(),"h1tmp2","tmp2"),"ctmp1","LPE");
 
    TCanvas *c2= muCS.plotCrossSection("cs2",1);
    if (!c2) return;
  }
  else {
    std::cout << "perform study\n";
    work(muCS,var,nSample);
  }
 
}


// --------------------------------------------------------------
// --------------------------------------------------------------

void work(MuonCrossSection_t &muCS, TVaried_t var, int nSample)
{
  std::vector<TH1D*> rndCSVec;
  int res=muCS.sampleRndVec(var,nSample,rndCSVec);
  if (!res) return;
  TH1D* h1Avg=NULL;
  TH2D* h2Cov=NULL;
  if (!muCS.deriveCov(rndCSVec,var,&h1Avg,&h2Cov)) return;
  
  TCanvas *c= muCS.plotCrossSection("cs");
  if (!c) return;
  h1Avg->SetLineColor(kGreen+1);
  h1Avg->SetMarkerColor(kGreen+1);
  plotHistoSame(h1Avg,"cs","LPE");

  TH2D* h2Corr=NULL;
  TCanvas *cx= plotCovCorr(h2Cov,"ccov",&h2Corr);
}

// --------------------------------------------------------------


