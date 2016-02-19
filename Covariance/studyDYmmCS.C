#include "crossSection.h"

void studyDYmmCS()
{
  MuonCrossSection_t muCS("muCS","v1");
  if (!muCS.load("cs_DYmm_13TeV_v1","v1")) {
    std::cout << "loading failed\n";
    return;
  }

  //TCanvas *ca=muCS.editCSa().plotCrossSection("csa");
  //TCanvas *cb=muCS.editCSb().plotCrossSection("csb");

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

