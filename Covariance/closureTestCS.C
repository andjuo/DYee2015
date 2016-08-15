#include "crossSection.h"

void closureTestCS()
{
  TVersion_t inpVer=_verEl2skim;
  TString inpVerTag=versionName(inpVer);

  TString csFName="cs_DYee_13TeV_El2.root";
  TString csOutFName="csClosure_DYee_13TeV_El2.root";
  CrossSection_t cs("cs",inpVerTag,_csPreFsrFullSp,inpVer);
  if (!cs.load(csFName,inpVerTag)) {
    std::cout << "loading failed\n";
    return;
  }

  // remove bkg and replace scale factors by 1.
  cs.removeBkg();
  cs.removeRho();

  // replace signal
  TH1D *h1reco_MC=cloneHisto((const TH1D*)cs.detRes()->Hmeasured(),"h1reco_MC","h1reco_MC");
  cs.h1Yield(h1reco_MC);

  // replace theory
  TH1D* h1preFSR_MC=cloneHisto((const TH1D*)cs.fsrRes()->Htruth(),"h1preFSR_MC","h1preFSR_MC");
  TH1D* h1preFSRCS_MC=cloneHisto(h1preFSR_MC,"h1preFSRCS_MC","h1preFSRCS_MC");
  h1preFSRCS_MC->Scale(1/cs.lumi());
  TH1D *h1preFSRCS_MC_perMBW= perMassBinWidth(h1preFSRCS_MC);
  cs.h1Theory(h1preFSRCS_MC_perMBW);

  if (0) {
    // replace Acc x Eff by my values
    TH1D* h1effAcc=loadHisto("dyee_test_El2skim_excludeGap.root","h1EffPUAcc",
			     "h1effAcc_my",1,h1dummy);
    cs.h1EffAcc(h1effAcc);
    csOutFName.ReplaceAll(".root","_AJeffAcc.root");
  }

  TH1D* h1cs= cs.calcCrossSection();
  printRatio(h1cs,cs.h1Theory());
  cs.plotCrossSection();

  cs.save(csOutFName);
}
