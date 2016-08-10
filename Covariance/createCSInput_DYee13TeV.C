#include "crossSection.h"

// ---------------------------------------------------------

// ---------------------------------------------------------
// ---------------------------------------------------------

void createCSInput_DYee13TeV(int doSave=0)
{
  TString srcPath="/mnt/sdb/andriusj/v1_31052016_CovarianceMatrixInputs/";
  TVersion_t inpVer=_verEl1;
  TString inpVerTag= versionName(inpVer);
  double lumiTot= 2300.; // approximate
  int setNiters_detRes=4;
  int setNiters_FSR=4;

  int checkBackgrounds=1;

  // ---------------- load the observed and background-subracted (signal) yield
  TString fname1=srcPath + "ROOTFile_Input1_Histograms_Data.root";
  TH1D* h1Yield= loadHisto(fname1, "h_data", "h1Yield", 1,h1dummy);
  TH1D* h1Signal= loadHisto(fname1, "h_yield", "h1Signal", 1,h1dummy);
  if (!h1Yield || !h1Signal) return;
  if (1) {
  //printHisto(h1Yield);
    histoStyle(h1Signal,kBlue,5,2);
    plotHisto(h1Yield,"cData",1,1,"LPE","Observed");
    plotHistoSame(h1Signal,"cData","LPE","Signal");
    printRatio(h1Yield,h1Signal);
  }

  // --------- Load backgrounds
  TString fname2=srcPath + "ROOTFile_Input2_Histograms_Bkg.root";
  //a) h_EMu: ttbar, diboson, dy->tautau and single top backgrounds estimated from data.
  //b) h_Dijet: QCD background estimated from data.
  //c) h_Wjet: W+jets background estimated from data.
  TH1D* h1Emu= loadHisto(fname2, "h_EMu", "h1Emu", 1,h1dummy);
  TH1D* h1Dijet= loadHisto(fname2, "h_Dijet", "h1QCD", 1,h1dummy);
  TH1D* h1Wjet= loadHisto(fname2, "h_Wjet", "h1Wjet", 1,h1dummy);
  if (!h1Emu || !h1Dijet || !h1Wjet) return;

  TString massAxisLabel= niceMassAxisLabel(0,"");
  TH1D* h1BkgTot= cloneHisto(h1Emu, "h1BkgTot","h1BkgTot;" + massAxisLabel + TString(";bkg.yield"));
  h1BkgTot->Add(h1Dijet);
  h1BkgTot->Add(h1Wjet);

  if (1) {
    histoStyle(h1Emu,kBlue,5,2);
    histoStyle(h1Dijet,kRed,24,2);
    histoStyle(h1Wjet,kGreen+1,27,2);
    plotHisto(h1BkgTot,"cBkg",1,1,"LPE","bkgTot");
    plotHistoSame(h1Emu,"cBkg","LPE","e-#mu");
    plotHistoSame(h1Dijet,"cBkg","LPE","QCD");
    plotHistoSame(h1Wjet,"cBkg","LPE","W+Jets");
  }

  if (checkBackgrounds) {
    TH1D *h1BkgFromSignal= cloneHisto(h1Yield, "h1BkgFromSignal","bkg from signal");
    h1BkgFromSignal->Add(h1Signal,-1);
    histoStyle(h1BkgFromSignal,kRed,24,2);
    plotHisto(h1BkgTot,"cBkgCheck",1,1,"LPE","bkgTot");
    plotHistoSame(h1BkgFromSignal,"cBkgCheck","LPE","bkg from signal");
    printRatio(h1BkgTot,h1BkgFromSignal);
  }

  // --------- Load corrections

  //"ROOTFile_Input6_CrossCheck.root"
  //a) UnfoldRes_DetectorRes: RooUnfoldResponse for unfolding correction for detector resolution.
  //b) UnfoldRes_FSR: RooUnfoldResponse for FSR correction.
  //c) h_[Acc/Eff/AccEff]: acceptance, efficiency and acceptance * efficiency for each mass bin.
  //d) h_EffSF: efficiency scale factor for data for each mass bin.
  //e) h_diffXsec_Meas: differential cross section distribution with all corrections applied.

  TString fname6= srcPath + "ROOTFile_Input6_CrossCheck.root";
  TFile fin6(fname6);
  if (!fin6.IsOpen()) {
    std::cout << "failed to open the file <" << fin6.GetName() << ">\n";
    return;
  }
  RooUnfoldResponse *rooUnfDetRes= loadRooUnfoldResponse(fin6,"UnfoldRes_DetectorRes","rooUnfDetRes");
  RooUnfoldResponse *rooUnfFSRRes= loadRooUnfoldResponse(fin6,"UnfoldRes_FSR","rooUnfFSR");
  TH1D *h1Acc= loadHisto(fin6,"h_Acc","h1Acc", "Acceptance;" + massAxisLabel + TString(";A"), 1,h1dummy);
  TH1D *h1Eff= loadHisto(fin6,"h_Eff","h1Eff", "Efficiency;" + massAxisLabel + TString(";#epsilon"), 1,h1dummy);
  TH1D *h1EffAcc= loadHisto(fin6,"h_AccEff","h1_EffAcc", "Acc #times Eff;" + massAxisLabel + TString(";A #times #epsilon"), 1,h1dummy);
  TH1D *h1Rho= loadHisto(fin6,"h_EffSF","h1Rho","Scale factor;" + massAxisLabel+ TString(";#rho"), 1,h1dummy);
  TH1D *h1_DiffXSec_data= loadHisto(fin6,"h_diffXsec_Meas","h1RC_xsec", 1,h1dummy);
  fin6.Close();
  if (!rooUnfDetRes || !rooUnfFSRRes) return;
  if (!h1Acc || !h1Eff || !h1EffAcc || !h1Rho || !h1_DiffXSec_data) return;


  if (1) {
    TH1D *h1effAcc_check= cloneHisto(h1Eff,"h1effAcc_check","h1effAcc_check");
    for (int ibin=1; ibin<=h1Acc->GetNbinsX(); ibin++) {
      double a= h1Eff->GetBinContent(ibin);
      double da= h1Eff->GetBinError(ibin);
      double b= h1Acc->GetBinContent(ibin);
      double db= h1Acc->GetBinError(ibin);
      double val= a * b;
      double unc= sqrt( b*b * da*da + a*a * db*db );
      h1effAcc_check->SetBinContent( ibin, val );
      h1effAcc_check->SetBinError( ibin, unc );


      const int detailed_checkEffAcc=0;
      if (detailed_checkEffAcc) {
	std::cout << "ibin=" << ibin << ", eff=" << a << " +- " << da
		  << ", acc=" << b << " +- " << db
		  << ", (eff x acc)=" << val << " +- " << unc << ", RC value ="
		  << h1EffAcc->GetBinContent(ibin) << " +- " << h1EffAcc->GetBinError(ibin)
		  << "\n";
	if (1) {
	  std::cout << "  relEff=" << da/a << ", relAcc=" << db/b << ", rel effxAcc=" << unc/val << "; rel RC=" << h1EffAcc->GetBinError(ibin)/h1EffAcc->GetBinContent(ibin) << "\n";
	}
      }
    }
    plotHisto(h1EffAcc, "cEffAcc", 1,1, "LPE1", "origEffAcc");
    plotHistoSame(h1effAcc_check, "cEffAcc", "LPE", "effAcc_check");
    printRatio(h1EffAcc,h1effAcc_check);
  }

  CrossSection_t cs("elCS",inpVerTag,_csPreFsrFullSp);
  cs.lumi(lumiTot);
  cs.h1Yield(h1Yield);
  cs.h1Bkg(h1BkgTot);
  cs.h1Eff(h1Eff);
  cs.h1Rho(h1Rho);
  cs.h1Acc(h1Acc);
  cs.h1EffAcc(h1EffAcc);
  rooUnfDetRes->UseOverflow(false);
  rooUnfFSRRes->UseOverflow(false);
  cs.detRes(*rooUnfDetRes);
  cs.fsrRes(*rooUnfFSRRes);
  cs.h1Theory(h1_DiffXSec_data);
  cs.nItersDetRes(setNiters_detRes);
  cs.nItersFSR(setNiters_FSR);

  if (1) {
    TH1D* h1cs= cs.calcCrossSection();
    printRatio(h1cs, h1_DiffXSec_data);
  }

  cs.plotCrossSection();

  //return;

  std::cout << "cs.checkPtrs()=" << cs.checkPtrs() << "\n";

  if (!doSave) { std::cout << "not saving" << std::endl; return; }
  std::cout << "calling cs.save()" << std::endl;
  if (!cs.save("cs_DYee_13TeV" + inpVerTag)) {
    std::cout << "saving failed\n";
    return;
  }

}

// ---------------------------------------------------------
// ---------------------------------------------------------


// ---------------------------------------------------------
