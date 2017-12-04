#include "crossSection.h"
#include <TText.h>

// ---------------------------------------------------------

// ---------------------------------------------------------
// ---------------------------------------------------------

void createCSInput_DYee13TeV(int doSave=0, int printHistosAtTheEnd=0)
{
  TString srcPath="/mnt/sdb/andriusj/v1_31052016_CovarianceMatrixInputs/";
  double lumiTot= 2316.97; // approximate
  int setNiters_detRes=4;
  int setNiters_FSR=4;
  bool forceOverflow=false;

  TVersion_t inpVer=_verEl1;
  TSubVersion_t inpSubVer=_subver1;
  inpVer=_verEl2;
  inpVer=_verEl3;
  inpVer=_verElMay2017;
  inpVer=_verElMay2017false;
  inpVer=_verElNov2017false;
  inpSubVer=_subver2; // take the final cs from another file

  if (inpVer==_verEl2) {
    srcPath="/mnt/sdb/andriusj/v2_08082016_CovarianceMatrixInputs/";
    lumiTot=2316.97; // approximate
    setNiters_detRes=15;
    setNiters_FSR=15;
  }
  else if (inpVer==_verEl3) {
    srcPath="/mnt/sdb/andriusj/v3_09092016_CovarianceMatrixInputs/";
    lumiTot=2316.969; // actually used
    setNiters_detRes=21;
    setNiters_FSR=21;
    forceOverflow=false;
  }
  else if ((inpVer==_verElMay2017) || (inpVer==_verElMay2017false)) {
    srcPath="/media/sf_CMSData/DY13TeV-CovInputs/v20170518_Input_Cov_ee/";
    lumiTot=2258.066; // new value
    setNiters_detRes=21;
    setNiters_FSR=21;
    forceOverflow=true;
    if (inpVer==_verElMay2017false) forceOverflow=false;
  }
  else if ((inpVer==_verElNov2017) || (inpVer==_verElNov2017false)) {
    srcPath="/media/sf_CMSData/DY13TeV-CovInputs/v20170518_Input_Cov_ee/";
    lumiTot=2258.066; // new value
    setNiters_detRes=21;
    setNiters_FSR=21;
    forceOverflow=true;
    if (inpVer==_verElNov2017false) forceOverflow=false;
  }

  TString inpVerTag= versionName(inpVer);
  std::cout << "inpVerTag=" << inpVerTag << "\n";

  int checkBackgrounds=1;

  // ---------------- load the observed and background-subracted (signal) yield
  TString fname1=srcPath + "ROOTFile_Input1_Histograms_Data.root";
  TH1D* h1Yield= loadHisto(fname1, "h_data", "h1Yield", 1,h1dummy);
  //printHisto(h1Yield,1); return;
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
  if (inpVer==_verEl1) {
    h1WZ= cloneHisto(h1Emu,"h1WZ","h1WZ");
    h1ZZ= cloneHisto(h1Emu,"h1ZZ","h1ZZ");
    h1WZ->Reset();
    h1ZZ->Reset();
  }
  else if ((inpVer==_verEl2) || (inpVer==_verEl3) || (inpVer==_verElMay2017) ||
	   (inpVer==_verElMay2017false) ||
	   (inpVer==_verElNov2017) || (inpVer==_verElNov2017false)) {
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
  if (inpVer==_verEl1) {
    hNameMap[_varDetRes]= "UnfoldRes_DetectorRes";
    hNameMap[_varFSRRes]= "UnfoldRes_FSR";
    hNameMap[_varEff]= "h_Eff";
    hNameMap[_varRho]= "h_EffSF";
    hNameMap[_varAcc]= "h_Acc";
    hNameMap[_varEffAcc]= "h_AccEff";
    hNameMap[_varLast]= "h_diffXsec_Meas";
  }
  else if (inpVer==_verEl2) {
    hNameMap[_varDetRes]= "UnfoldRes_DetectorRes";
    hNameMap[_varFSRRes]= "Unfold_FSRCorr";
    hNameMap[_varEff]= "h_EffCorr";
    hNameMap[_varRho]= "h_EfficiencySF";
    hNameMap[_varAcc]= "h_AccCorr";
    hNameMap[_varEffAcc]= "h_AccEff";
    hNameMap[_varLast]= "h_diffXsec_Meas";
  }
  else if ((inpVer==_verEl3) || (inpVer==_verElMay2017) ||
	   (inpVer==_verElMay2017false)) {
    hNameMap[_varDetRes]= "Unfold_DetectorRes";
    hNameMap[_varFSRRes]= "Unfold_FSRCorr";
    hNameMap[_varEff]= "h_Eff";
    hNameMap[_varRho]= "h_EffSF";
    hNameMap[_varAcc]= "h_Acc";
    hNameMap[_varEffAcc]= "h_AccEff";
    hNameMap[_varLast]= "h_diffXsec_Meas";
  }
  else if ((inpVer==_verElNov2017) || (inpVer==_verElNov2017false)) {
    hNameMap[_varDetRes]= "Unfold_DetectorRes";
    hNameMap[_varFSRRes]= "Unfold_FSRCorr";
    hNameMap[_varEff]= "h_Eff";
    hNameMap[_varRho]= "h_EffSF";
    hNameMap[_varAcc]= "h_Acc";
    hNameMap[_varEffAcc]= "h_AccEff";
    hNameMap[_varLast]= "h_diffXsec_Meas";
  }
  else {
    std::cout << "code (hNameMap) is not ready for this inpVer\n";
    return;
  }

  TString fname6= srcPath + "ROOTFile_Input6_CrossCheck.root";
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
  TH1D *h1_DiffXSec_data= loadHisto(fin6,hNameMap[_varLast],"h1RC_xsec", 1,h1dummy);
  fin6.Close();
  if (!rooUnfDetRes || !rooUnfFSRRes) return;
  if (!h1Acc || !h1Eff || !h1EffAcc || !h1Rho || !h1_DiffXSec_data) return;

  if (((inpVer==_verElNov2017) || (inpVer==_verElNov2017false)) &&
      (inpSubVer==_subver2)) {
    std::cout << "\n\tSpecial case for " << versionName(inpVer)
	      << subVersionName(inpSubVer) << "\n";
    if (0) h1_DiffXSec_data= loadHisto("/media/sf_CMSData/DY13TeV-CovInputs/v20170518_Input_Cov_ee/DiffXsec_Electron_v7_RC20171128.root","h_DiffXSec","h1RC_xsec_spec", 1,h1dummy);
    else  h1_DiffXSec_data= loadHisto("/media/sf_CMSData/DY13TeV-CovInputs/v20170518_Input_Cov_ee/ROOTFile_Input_Electron_v4_RC20171128.root","h_DiffXSec","h1RC_xsec_spec", 1,h1dummy);
    if (!h1_DiffXSec_data) return;
  }

  //if (rooUnfFSRRes && (inpVer==_verEl3)) {
  //  std::cout << "\n\tscaling rooUnfFSRRes by luminosity\n";
  //  rooUnfFSRRes->Scale(lumiTot);
  //}

  if (0 &&
      ((inpVer==_verEl3) ||
       (inpVer==_verElMay2017) || (inpVer==_verElMay2017false) ||
       (inpVer==_verElNov2017) || (inpVer==_verElNov2017false)
       )) {
    std::cout << "\n\nReplacing h1EffAcc with my version\n";
    TFile finAJ("dyee_test_dressed_" + versionName(inpVer) + TString(".root"));
    if (!finAJ.IsOpen()) {
      std::cout << "failed to open the file <" << finAJ.GetName() << ">\n";
      return;
    }
    TH1D *h1EffPU_aj= loadHisto(finAJ,"h1EffPU","h1EffPU_aj",1,h1dummy);
    TH1D *h1EffPUAcc_aj= loadHisto(finAJ,"h1EffPUAcc","h1EffPUAcc_aj",1,h1dummy);
    TH1D *h1rho_aj= loadHisto(finAJ, "h1rho","h1rho_aj",1,h1dummy);
    finAJ.Close();
    CrossSection_t csTmp("elCS_tmp",inpVerTag+"_tmp",_csPreFsrFullSp);
    csTmp.lumi(lumiTot);
    csTmp.h1Yield(h1Yield);
    csTmp.h1Bkg(h1BkgTot);
    csTmp.h1Eff(h1EffPU_aj);
    if (1) csTmp.h1Rho(h1Rho);
    //else csTmp.h1Rho(h1rho_aj);
    csTmp.h1Acc(h1Acc);
    csTmp.h1EffAcc(h1EffPUAcc_aj);
    rooUnfDetRes->UseOverflow(forceOverflow);
    rooUnfFSRRes->UseOverflow(forceOverflow);
    csTmp.detRes(*rooUnfDetRes);
    csTmp.fsrRes(*rooUnfFSRRes);
    csTmp.h1Theory(h1_DiffXSec_data);
    csTmp.nItersDetRes(setNiters_detRes);
    csTmp.nItersFSR(setNiters_FSR);
    TH1D *h1cs_aj= csTmp.calcCrossSection();
    TCanvas *canv_aj= csTmp.plotCrossSection();
    if (!canv_aj) { std::cout << "null canv_aj\n"; }
    printRatio(h1cs_aj,h1_DiffXSec_data);
    // copying
    copyContents(h1Eff, h1EffPU_aj);
    copyContents(h1EffAcc, h1EffPUAcc_aj);
    copyContents(h1_DiffXSec_data, h1cs_aj);
    copyContents(h1Rho, h1rho_aj);
    //return;
  }
  else if (1) {
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


  if ((inpVer==_verElMay2017) || (inpVer==_verElMay2017false) ||
      (inpVer==_verElNov2017) || (inpVer==_verElNov2017false)) {
    nullifyOverflow(h1Yield);
    nullifyOverflow(h1BkgTot);
    printHisto(h1Yield,1);
    printHisto(h1BkgTot,1);
  }

  CrossSection_t cs("elCS",inpVerTag,_csPreFsrFullSp);
  cs.useOverflow(forceOverflow);
  cs.lumi(lumiTot);
  cs.h1Yield(h1Yield);
  cs.h1Bkg(h1BkgTot);
  cs.h1Eff(h1Eff);
  cs.h1Rho(h1Rho);
  cs.h1Acc(h1Acc);
  cs.h1EffAcc(h1EffAcc);
  rooUnfDetRes->UseOverflow(forceOverflow);
  rooUnfFSRRes->UseOverflow(forceOverflow);
  cs.detRes(*rooUnfDetRes);
  cs.fsrRes(*rooUnfFSRRes);
  cs.h1Theory(h1_DiffXSec_data);
  cs.nItersDetRes(setNiters_detRes);
  cs.nItersFSR(setNiters_FSR);

  TH1D* h1cs= cs.calcCrossSection();
  if (h1cs) printRatio(h1cs, h1_DiffXSec_data);

  TCanvas *canv=cs.plotCrossSection();
  changeLegendEntry(canv,cs.h1Theory(),"Ridhi calc.");
  TText *info= new TText();
  info->DrawTextNDC(0.45,0.93,inpVerTag);
  canv->Update();


  TH1D *h1theoryRaw= loadHisto("theory13TeVmm.root","h1cs_theory","h1NNLOtheory",1,h1dummy);
  TH1D *h1theory= perMassBinWidth(h1theoryRaw,0);
  histoStyle(h1theory,kRed,7,1);
  plotHistoSame(h1theory,"cs","LP","theory");
  printRatio(h1cs,h1theory);

  if (printHistosAtTheEnd) {
    printHisto(h1Yield,1);
    printHisto(cs.h1Yield(),1);
    printHisto(cs.h1Signal(),1);
    printHisto(cs.h1Unf(),1);
    printHisto(cs.h1UnfRhoCorr(),1);
    printHisto(cs.h1UnfRhoEffAccCorr());
    TH1D *h1postFSRcs= cloneHisto(cs.h1UnfRhoEffAccCorr(),
				  "h1postFSRcs","h1postFSRcs");
    h1postFSRcs->Scale(1/cs.lumi());
    printHisto(h1postFSRcs,1);
    TH1D *h1preFSRcs= cloneHisto(cs.h1PreFsr(),
				 "h1preFSRcs","h1preFSRcs");
    h1preFSRcs->Scale(1/cs.lumi());
    printHisto(h1preFSRcs,1);
    printHisto(cs.h1PreFsr());
    printHisto(cs.h1PreFsrCS());
  }


  //cs.plotCrossSection_StepByStep("cCalcTest");

  //return;

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
