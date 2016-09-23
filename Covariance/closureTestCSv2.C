#include "crossSection.h"
#include "DYbinning.h"
#include <map>

void closureTestCSv2(int dataUnfold=0, int compareToTheory=0)
{
  TVersion_t inpVer=_verEl2skim3;
  TString inpVerTag=versionName(inpVer);
  TString csOutFName="csClosure_DYee_13TeV_" + inpVerTag + ".root";

  const TString tag_UseData="data";
  const TString tag_UseTheory="theory";
  const TString tag_UseUnity="unity";
  const TString tag_UseZero="zero";
  const TString tag_Load="combine:";
  std::vector<TString> dataInput, theoryInput;
  addToVector(dataInput, "cs_DYee_13TeV_El2.root  h1Yield  h1Bkg");
  addToVector(theoryInput, "theory13TeVmm.root h1cs_theory");

  std::map<TVaried_t,TString> inpHistoNames;

  TCSType_t csType= _csPreFsrFullSp;
  TString inpFName="dyee_test_dressed_El2skim3.root";
  inpHistoNames[_varYield] = "h1recoSel_MWPU";
  //inpHistoNames[_varBkg] = tag_UseZero;
  inpHistoNames[_varDetRes]= "rooUnf_detResRespPU";
  inpHistoNames[_varRho] = tag_UseUnity;
  inpHistoNames[_varEffAcc] = "h1EffPUAcc";
  inpHistoNames[_varFSRRes] = "rooUnf_fsrResp";
  inpHistoNames[_varLast] = "h1_preFsr_Mweighted";

  if (1) inpHistoNames[_varEffAcc] = "h1EffPUEleMatchAcc";

  if (dataUnfold) {
    inpHistoNames[_varYield] = "data";
    //inpHistoNames[_varRho] = "h1rho";
    //inpHistoNames[_varRho] = "h1rho_inPostFsrAcc";
  }

  TH1D *h1Zero= new TH1D("h1Zero","h1Zero",DYtools::nMassBins,DYtools::massBinEdges);
  h1Zero->SetDirectory(0);
  TH1D *h1Unity= cloneHisto(h1Zero,"h1Unity","h1Unity");
  for (int ibin=1; ibin<=h1Unity->GetNbinsX(); ibin++) {
    h1Unity->SetBinContent(ibin,1.);
  }

  CrossSection_t cs("cs",inpVerTag,csType,inpVer);


  for (std::map<TVaried_t,TString>::const_iterator it=inpHistoNames.begin();
       it!=inpHistoNames.end(); it++) {
    TH1D *h1tmp=NULL, *h1tmp2=NULL;
    std::cout << "load " << variedVarName(it->first) << " named " << it->second << "\n";
    switch(it->first) {
    case _varYield:
      if (it->second == tag_UseData) {
	h1tmp= loadHisto(dataInput[0],dataInput[1], "h1Yield",1,h1dummy);
	h1tmp2=loadHisto(dataInput[0],dataInput[2], "h1Bkg",1,h1dummy);
	double lumi=2316.97;
	if (h1tmp && h1tmp2) {
	  h1tmp->Scale(1/lumi);
	  h1tmp2->Scale(1/lumi);
	}
      }
      else {
	h1tmp= loadHisto(inpFName,it->second, "h1Yield",1,h1dummy);
	h1tmp2=h1Zero;
      }
      if (!h1tmp || !h1tmp2) return;
      cs.h1Yield(h1tmp);
      cs.h1Bkg(h1tmp2);
      break;
    case _varDetRes: {
      RooUnfoldResponse *r= loadRooUnfoldResponse(inpFName,it->second,"detRes");
      if (!r) return;
      cs.detRes(*r);
    }
      break;
    case _varRho:
      if (it->second == tag_UseUnity) {
	h1tmp= h1Unity;
      }
      else {
	h1tmp= loadHisto(inpFName,it->second, "h1Rho",1,h1dummy);
      }
      if (!h1tmp) return;
      cs.h1Rho(h1tmp);
      break;
    case _varEff:
      h1tmp= loadHisto(inpFName,it->second, "h1Eff",1,h1dummy);
      if (!h1tmp) return;
      cs.h1Eff(h1tmp);
      break;
    case _varEffAcc:
      h1tmp= loadHisto(inpFName,it->second, "h1EffAcc",1,h1dummy);
      if (!h1tmp) return;
      cs.h1EffAcc(h1tmp);
      break;
    case _varAcc:
      h1tmp= loadHisto(inpFName,it->second, "h1Acc",1,h1dummy);
      if (!h1tmp) return;
      cs.h1Acc(h1tmp);
      break;
    case _varFSRRes: {
      RooUnfoldResponse *r= loadRooUnfoldResponse(inpFName,it->second,"fsrRes");
      if (!r) return;
      cs.fsrRes(*r);
    }
      break;
    case _varLast:
      h1tmp2= loadHisto(inpFName,it->second, "h1theory",1,h1dummy);
      if (!h1tmp2) return;
      h1tmp= perMassBinWidth(h1tmp2);
      cs.h1Theory(h1tmp);
      break;
    default:
      std::cout << "unprepared case " << variedVarName(it->first) << "\n";
      return;
    }
  }


  TH1D* h1cs= cs.calcCrossSection();
  printRatio(h1cs,cs.h1Theory());
  cs.plotCrossSection();

  if (1) {
    // unfold in 1 step
    TH1D* h1reco=cloneHisto(cs.h1Yield(),"h1reco","h1reco");
    RooUnfoldResponse *rFull= loadRooUnfoldResponse(inpFName,"rooUnf_detResEffAccFsrResp","fullResponse");
    if (!h1reco || !rFull) return;
    RooUnfoldBayes bayesFull( rFull, h1reco, 4);
    TH1D *h1UnfFull=perMassBinWidth((TH1D*)bayesFull.Hreco());
    histoStyle(h1UnfFull,kOrange,5,1);
    plotHistoSame(h1UnfFull,"cs","LPE1","fullUnf");
    printRatio(h1UnfFull,cs.h1Theory());
    std::cout << "\n\ndifferent unfoldings\n";
    printRatio(h1UnfFull,cs.h1PreFsrCS());
  }

  //cs.save(csOutFName);
}
