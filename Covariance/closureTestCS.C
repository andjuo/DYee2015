#include "crossSection.h"

void closureTestCS()
{
  TVersion_t inpVer=_verEl2;
  inpVer=_verEl3;
  TString inpVerTag=versionName(inpVer);

  TString csFName="cs_DYee_13TeV_" + inpVerTag + ".root";
  TString csOutFName="csClosure_DYee_13TeV_" + inpVerTag + ".root";
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
    TH1D* h1effAcc=loadHisto("dyee_test_dressed_El2skim3.root","h1EffPUAcc",
			     "h1effAcc_my",1,h1dummy);
    cs.h1EffAcc(h1effAcc);
    csOutFName.ReplaceAll(".root","_AJeffAcc.root");
  }
  if (1) {
    // replace Acc x Eff by modified code values
    TString fname="/home/andriusj/DY13TeV/DYanalysis-20160817/ElectronNtupler/test/Analysis_Codes/AccEff/dyee_preFSR_forAccEff_v2.root";
    TH1D *h1effAcc= loadHisto(fname,"h1effAcc","h1effAcc_v2",1,h1dummy);
    cs.h1EffAcc(h1effAcc);
    csOutFName.ReplaceAll(".root","_effAccV2.root");
  }
  if (0) {
    // replace Acc x Eff by my values
    TH1D* h1effAcc=loadHisto("dyee_test_El2skim2_excludeGap.root","h1EffPUAcc",
			     "h1effAcc_my",1,h1dummy);
    cs.h1EffAcc(h1effAcc);
    csOutFName.ReplaceAll(".root","_AJeffAcc.root");
  }
  if (0) {
    // replace detRes by my values
    RooUnfoldResponse *detRes=loadRooUnfoldResponse("dyee_test_RECO_El2skim.root","rooUnf_detResRespPU","detResRespPU");
    cs.detRes(*detRes);
  }
  if (0) {
    // replace fsrRes by my values
    RooUnfoldResponse *fsrRes=loadRooUnfoldResponse("dyee_test_El2skim2_excludeGap.root","rooUnf_fsrResp","fsrResp");
    cs.fsrRes(*fsrRes);
  }

  TH1D* h1cs= cs.calcCrossSection();
  printRatio(h1cs,cs.h1Theory());
  cs.plotCrossSection();

  cs.save(csOutFName);
}
