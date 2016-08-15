#include "inputs.h"
#include "crossSection.h"
#include "../RooUnfold/src/RooUnfoldResponse.h"
#include <map>

typedef enum { _fn_main=0, _fn_cmp1, _fn_cmp2, _fn_last } TSource_t;

typedef enum { _iEff, _iRho, _iAcc, _iEffAcc,
	       _iSignal, _iUnf, _iUnfRho, _iUnfRhoEffAcc,
	       _iPostFsr, _iPreFsr
	       //_iRRmeas, _iRRtrue
} THisto_t;

// ---------------------------------------------------------

TString THistoName(THisto_t ih);
TH1D* loadHisto(const std::map<THisto_t,TString> &histoNames,
		THisto_t idxHisto,
		TString fileName, TString newHistoName);

// ---------------------------------------------------------

void compareVersions(int theCase=1)
{
  TString mainFName;
  TString cmpFName1, cmpFName2;
  std::map<THisto_t,TString> histoNames[_fn_last];

  if (theCase==1) {
    std::cout << "compare DYee corrections\n";
    mainFName="cs_DYee_13TeV_El2.root";
    histoNames[_fn_main][_iEff]="h1Eff";
    histoNames[_fn_main][_iRho]="h1Rho";
    histoNames[_fn_main][_iAcc]="h1Acc";
    histoNames[_fn_main][_iEffAcc]="h1EffAcc";
    histoNames[_fn_main][_iPostFsr]="RRmeas+FSRRes";
    histoNames[_fn_main][_iPreFsr]="RRtrue+FSRRes";
    cmpFName1="dyee_test_El2skim_excludeGap.root";
    histoNames[_fn_cmp1][_iEff]="h1Eff";
    histoNames[_fn_cmp1][_iEff]="h1EffPU";
    histoNames[_fn_cmp1][_iAcc]="h1Acc";
    histoNames[_fn_cmp1][_iEffAcc]="h1EffPUAcc";
    TString lumiScale=":scale2316.97";
    histoNames[_fn_cmp1][_iPostFsr]="h1_postFsr_Mweighted"+lumiScale;
    //histoNames[_fn_cmp1][_iPostFsr]="RRmeas+rooUnf_fsrResp"+lumiScale;
    histoNames[_fn_cmp1][_iPreFsr]="RRtrue+rooUnf_fsrResp"+lumiScale;
    cmpFName2="dyee_test_RECO_El2skim.root";
    histoNames[_fn_cmp2][_iRho]="h1rho";
  }
  else if ((theCase==2) || (theCase==3)) {
    std::cout << "compare DYee CS closure test\n";
    mainFName="csClosure_DYee_13TeV_El2.root";
    if (theCase==3) {
      std::cout << " -- A x eff replaced by my distribution\n";
      mainFName="csClosure_DYee_13TeV_El2_AJeffAcc.root";
    }
    histoNames[_fn_main][_iSignal]="h1Signal";
    histoNames[_fn_main][_iUnf]="h1Unf";
    histoNames[_fn_main][_iUnfRho]="h1UnfRho";
    histoNames[_fn_main][_iUnfRhoEffAcc]="h1UnfRhoEffAcc";
    cmpFName1="cs_DYee_13TeV_El2.root";
    histoNames[_fn_cmp1][_iSignal]="RRmeas+detRes";
    histoNames[_fn_cmp1][_iUnf]="RRtrue+detRes";
    histoNames[_fn_cmp1][_iUnfRho]="RRtrue+detRes";
    histoNames[_fn_cmp1][_iUnfRhoEffAcc]="RRmeas+FSRRes";
  }
  //else if (theCase==4)
  else {
    std::cout << "not ready for theCase=" << theCase << "\n";
  }


  for (std::map<THisto_t,TString>::const_iterator it= histoNames[_fn_main].begin();
       it!=histoNames[_fn_main].end(); it++) {
    std::cout << THistoName(it->first) << " " << it->second << "\n";

    TString h1nameBase= "h1" + THistoName(it->first);
    TH1D *h1main= loadHisto(histoNames[_fn_main],it->first,mainFName,h1nameBase+"_main");
    TH1D *h1cmp1= loadHisto(histoNames[_fn_cmp1],it->first,cmpFName1,h1nameBase+"_cmp1");
    TH1D* h1cmp2= loadHisto(histoNames[_fn_cmp2],it->first,cmpFName2,h1nameBase+"_cmp2");

    histoStyle(h1main,kBlue,24,1);
    if (h1cmp1) histoStyle(h1cmp1,kGreen+1,5,1);
    if (h1cmp2) histoStyle(h1cmp2,kRed,7,1);

    TString cName="c" + h1nameBase;
    plotHisto(h1main,cName,1,0,"LPE1","main");
    if (h1cmp1) plotHistoSame(h1cmp1,cName,"LPE1","cmp1");
    if (h1cmp2) plotHistoSame(h1cmp2,cName,"LPE1","cmp2");
    if (1) {
      if (h1cmp1) printRatio(h1main,h1cmp1);
      if (h1cmp2) printRatio(h1main,h1cmp2);
    }
  }
}

// ---------------------------------------------------------
// ---------------------------------------------------------

TString THistoName(THisto_t ih)
{
  TString name;
  switch(ih) {
    // corrections
  case _iEff: name="Eff"; break;
  case _iRho: name="Rho"; break;
  case _iAcc: name="Acc"; break;
  case _iEffAcc: name="EffAcc"; break;
    // distributions
  case _iSignal: name="Signal"; break;
  case _iUnf: name="Unf"; break;
  case _iUnfRho: name="UnfRho"; break;
  case _iUnfRhoEffAcc: name="UnfRhoEffAcc"; break;
  case _iPostFsr: name="postFSR"; break;
  case _iPreFsr: name="preFSR"; break;
  default:
    name="UNKNOWN";
  }
  return name;
}

// ---------------------------------------------------------

TH1D* loadHisto(const std::map<THisto_t,TString> &histoNames,
		THisto_t idxHisto,
		TString fileName, TString newHistoName)
{
  TH1D *h1=NULL;
  if (fileName.Length()==0) return h1;

  for (std::map<THisto_t,TString>::const_iterator it=histoNames.begin();
       it!=histoNames.end(); it++) {
    if (it->first == idxHisto) {
      TString hNameDisk=it->second;
      double scale=1.;
      if (hNameDisk.Index(":scale")!=-1) {
	TString valueStr=hNameDisk(hNameDisk.Index(":scale")+6,hNameDisk.Length());
	//std::cout << "valueStr=" << valueStr << "\n";
	if (valueStr.Length()) scale=atof(valueStr.Data());
	hNameDisk.ReplaceAll(":scale" + valueStr,"");
      }
      if (hNameDisk.Index("+")==-1) {
	h1= loadHisto(fileName,hNameDisk, newHistoName,h1dummy);
      }
      else if ((hNameDisk.Index("RRmeas+")!=-1) ||
	       (hNameDisk.Index("RRtrue+")!=-1)) {
	TString respName=hNameDisk;
	respName.Remove(0,respName.Index("+")+1);
	std::cout << "respName=" << respName << "\n";
	RooUnfoldResponse *r= loadRooUnfoldResponse(fileName,respName,"resp"+newHistoName);
	if (!r) {
	  std::cout << "failed to get response " << respName << " from "
		    << fileName << "\n";
	}
	else {
	  if (hNameDisk.Index("RRmeas+")!=-1)
	    h1=cloneHisto((TH1D*)r->Hmeasured(),newHistoName,newHistoName);
	  else if (hNameDisk.Index("RRtrue+")!=-1)
	    h1=cloneHisto((TH1D*)r->Htruth(),newHistoName,newHistoName);
	  delete r;
	}
      }
      if (!h1) {
	std::cout << "failed to get " << it->second << "\n";
      }
      else if (scale!=1.) h1->Scale(scale);
      break;
    }
  }
  return h1;
}

// ---------------------------------------------------------
