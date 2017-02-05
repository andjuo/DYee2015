#include "crossSection.h"
#include "DYbinning.h"
#include <TText.h>

// ---------------------------------------------------------
// ---------------------------------------------------------

void createCSInput_DYee13TeV_mb4X(int doSave=0)
{
  closeCanvases();

  TString srcPath="/mnt/sdb/andriusj/v1_31052016_CovarianceMatrixInputs/";
  double lumiTot= 2316.97; // approximate
  int setNiters_detRes=4;
  int setNiters_FSR=4;

  TVersion_t inpVer= DYtools::activeVer;
  std::cout << "\n\n\t** active version to " << versionName(inpVer)<<"\n\n";

  srcPath="/mnt/sdb/andriusj/v3_09092016_CovarianceMatrixInputs/";
  lumiTot=2316.969; // actually used
  setNiters_detRes=21;
  setNiters_FSR=21;

  TString inpVerTag= versionName(inpVer);
  std::cout << "inpVerTag=" << inpVerTag << "\n";

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
    std::cout << "\nNote about the next check: "
	      << "yield and signal does not necesarily match\n";
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
  // d) version 2 has also h_WZ and h_ZZ
  TH1D *h1WZ=NULL, *h1ZZ=NULL;
  if ((inpVer==_verEl3mb41) || (inpVer==_verEl3mb42)) {
    h1WZ= loadHisto(fname2, "h_WZ", "h1WZ", 1,h1dummy);
    h1ZZ= loadHisto(fname2, "h_ZZ", "h1ZZ", 1,h1dummy);
    if (!h1WZ || !h1ZZ) {
      std::cout << "failed to obtain h_WZ,h_ZZ\n";
      return;
    }
  }
  else {
    std::cout << "the code is not ready for this inpVer (h1WZ,h1ZZ)\n";
    return;
  }

  TString massAxisLabel= niceMassAxisLabel(0,"");
  TH1D* h1BkgTot= cloneHisto(h1Emu, "h1BkgTot","h1BkgTot;" + massAxisLabel + TString(";bkg.yield"));
  h1BkgTot->Add(h1Dijet);
  h1BkgTot->Add(h1Wjet);
  h1BkgTot->Add(h1WZ);
  h1BkgTot->Add(h1ZZ);

  if (1) {
    histoStyle(h1Emu,kBlue,5,2);
    histoStyle(h1Dijet,kRed,24,2);
    histoStyle(h1Wjet,kGreen+1,27,2);
    histoStyle(h1WZ,6,32,3);
    histoStyle(h1ZZ,46,26,3);
    plotHisto(h1BkgTot,"cBkg",1,1,"LPE","bkgTot");
    plotHistoSame(h1Emu,"cBkg","LPE","e-#mu");
    plotHistoSame(h1Dijet,"cBkg","LPE","QCD");
    plotHistoSame(h1Wjet,"cBkg","LPE","W+Jets");
    plotHistoSame(h1WZ,"cBkg","LPE","WZ");
    plotHistoSame(h1ZZ,"cBkg","LPE","ZZ");
  }

  if (checkBackgrounds) {
    TH1D *h1BkgFromSignal= cloneHisto(h1Yield, "h1BkgFromSignal","bkg from signal");
    h1BkgFromSignal->Add(h1Signal,-1);
    histoStyle(h1BkgFromSignal,kRed,24,2);
    plotHisto(h1BkgTot,"cBkgCheck",1,1,"LPE","bkgTot");
    plotHistoSame(h1BkgFromSignal,"cBkgCheck","LPE","bkg from signal");
    printRatio(h1BkgTot,h1BkgFromSignal);
    std::cout << "background check end\n\n";
  }

  // --------- Load corrections

  //"ROOTFile_Input6_CrossCheck.root"
  //a) UnfoldRes_DetectorRes: RooUnfoldResponse for unfolding correction for detector resolution.
  //b) UnfoldRes_FSR: RooUnfoldResponse for FSR correction.
  //c) h_[Acc/Eff/AccEff]: acceptance, efficiency and acceptance * efficiency for each mass bin.
  //d) h_EffSF: efficiency scale factor for data for each mass bin.
  //e) h_diffXsec_Meas: differential cross section distribution with all corrections applied.

  std::map<TVaried_t,TString> hNameMap;
  if ((inpVer==_verEl3mb41) || (inpVer==_verEl3mb42)) {
    hNameMap[_varDetRes]= "rooUnf_detResResp";
    hNameMap[_varFSRRes]= "rooUnf_fsrResp";
    hNameMap[_varEff]= "h1EffPU";
    hNameMap[_varRho]= "h1rho";
    hNameMap[_varAcc]= "h1Acc";
    hNameMap[_varEffAcc]= "h1EffPUAcc";
    hNameMap[_varLast]= "h1preFsrUnf";
  }
  else {
    std::cout << "code (hNameMap) is not ready for this inpVer\n";
    return;
  }

  TString fname6= srcPath + "ROOTFile_Input6_CrossCheck.root";
  if (inpVer==_verEl3mb41) {
    fname6="dyee_test_dressed_El3mb41.root";
  }
  else if (inpVer==_verEl3mb42) {
    fname6="dyee_test_dressed_El3mb42.root";
  }
  TFile fin6(fname6);
  if (!fin6.IsOpen()) {
    std::cout << "failed to open the file <" << fin6.GetName() << ">\n";
    return;
  }
  RooUnfoldResponse *rooUnfDetRes= loadRooUnfoldResponse(fin6,hNameMap[_varDetRes],"rooUnfDetRes");
  RooUnfoldResponse *rooUnfFSRRes= loadRooUnfoldResponse(fin6,hNameMap[_varFSRRes],"rooUnfFSR");
  TH1D *h1Acc= loadHisto(fin6,hNameMap[_varAcc],"h1Acc", "Acceptance;" + massAxisLabel + TString(";A"), 1,h1dummy);
  TH1D *h1Eff= loadHisto(fin6,hNameMap[_varEff],"h1Eff", "Efficiency;" + massAxisLabel + TString(";#epsilon"), 1,h1dummy);
  TH1D *h1EffAcc= loadHisto(fin6,hNameMap[_varEffAcc],"h1_EffAcc", "Acc #times Eff;" + massAxisLabel + TString(";A #times #epsilon"), 1,h1dummy);
  TH1D *h1Rho= loadHisto(fin6,hNameMap[_varRho],"h1Rho","Scale factor;" + massAxisLabel+ TString(";#rho"), 1,h1dummy);
  TH1D *h1_DiffXSec_data= loadHisto(fin6,hNameMap[_varLast],"h1my_xsec", 1,h1dummy);
  fin6.Close();
  if (!rooUnfDetRes || !rooUnfFSRRes) return;
  if (!h1Acc || !h1Eff || !h1EffAcc || !h1Rho || !h1_DiffXSec_data) return;

  TH1D* h1_DiffXSec_data_raw= h1_DiffXSec_data;
  h1_DiffXSec_data= perMassBinWidth(h1_DiffXSec_data_raw);

  //if (rooUnfFSRRes && (inpVer==_verEl3)) {
  //  std::cout << "\n\tscaling rooUnfFSRRes by luminosity\n";
  //  rooUnfFSRRes->Scale(lumiTot);
  //}

  if ((inpVer==_verEl3mb41) || (inpVer==_verEl3mb42)) {
    std::cout << "\n\nRebinning Ridhi's input\n";
  }


  TH1D *h1Yield4X=DYtools::rebin_43toLess(h1Yield);
  TH1D *h1Bkg4X=DYtools::rebin_43toLess(h1BkgTot);

  CrossSection_t cs("elCS",inpVerTag,_csPreFsrFullSp);
  cs.lumi(lumiTot);
  cs.h1Yield(h1Yield4X);
  cs.h1Bkg(h1Bkg4X);
  cs.h1Eff(h1Eff);
  cs.h1Rho(h1Rho);
  cs.h1Acc(h1Acc);
  cs.h1EffAcc(h1EffAcc);
  rooUnfDetRes->UseOverflow(false);
  rooUnfFSRRes->UseOverflow(false);
  cs.detRes(*rooUnfDetRes);
  cs.fsrRes(*rooUnfFSRRes);
  //cs.h1Theory(h1_DiffXSec_data);
  cs.h1Theory(NULL);
  cs.nItersDetRes(setNiters_detRes);
  cs.nItersFSR(setNiters_FSR);

  TH1D* h1cs= cs.calcCrossSection();
  if (h1cs) printRatio(h1cs, h1_DiffXSec_data);

  TCanvas *canv=cs.plotCrossSection();
  TText *info= new TText();
  info->DrawTextNDC(0.45,0.93,inpVerTag);
  canv->Update();

  TH1D *h1theoryRaw= loadHisto("theory13TeVmm.root","h1cs_theory","h1NNLOtheory",1,h1dummy);

  TString nMBStr= Form(" (%d)",DYtools::nMassBins);

  TH1D *h1theory4X= perMassBinWidth(DYtools::rebin_43toLess(h1theoryRaw));
  histoStyle(h1theory4X,kGreen+2,5,1);
  plotHistoSame(h1theory4X,"cs","LP","theory"+nMBStr);
  printRatio(h1cs,h1theory4X);

  TH1D *h1theory= perMassBinWidth(h1theoryRaw,0);
  histoStyle(h1theory,kRed,7,1);
  plotHistoSame(h1theory,"cs","LP","theory (43)");

  //cs.plotCrossSection_StepByStep("cCalcTest");

  //return;

  TH1D *h1Signal4X=cloneHisto(h1Yield4X,"h1Signal4X","h1Signal4X");
  h1Signal4X->Add(h1Bkg4X,-1);

  histoStyle(h1Yield4X,kBlue,5,1);
  histoStyle(h1Signal4X,kGreen+1,7,1);
  h1Yield4X->SetStats(0);
  h1Yield4X->SetTitle(Form("dielectron measured yield (%d mass bins)",DYtools::nMassBins));
  h1Signal4X->SetStats(0);
  TString cName= Form("cSignal%d",DYtools::nMassBins);
  plotHisto(h1Yield4X,cName,1,1,"LPE","yield"+nMBStr);
  plotHistoSame(h1Signal4X,cName,"LPE","signal"+nMBStr);

  std::cout << "cs.checkPtrs()=" << cs.checkPtrs() << "\n";

  if (!doSave) { std::cout << "not saving" << std::endl; return; }
  std::cout << "calling cs.save()" << std::endl;
  if (!cs.save("cs_DYee_13TeV_" + inpVerTag + ".root")) {
    std::cout << "saving failed\n";
    return;
  }

}

// ---------------------------------------------------------
// ---------------------------------------------------------


// ---------------------------------------------------------
