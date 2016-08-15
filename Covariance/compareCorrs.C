// Enhanced version of the code is compareVersions.C -- 2016.08.14

#include "inputs.h"
#include "crossSection.h"
#include <map>

typedef enum { _fn_main=0, _fn_cmp1, _fn_cmp2, _fn_last } TSource_t;

void compareCorrs(int theCase=1)
{
  TString mainFName;
  TString cmpFName1, cmpFName2;
  std::map<TVaried_t,TString> histoNames[_fn_last];

  if (theCase==1) {
    std::cout << "compare DYee corrections\n";
    mainFName="cs_DYee_13TeV_El2.root";
    histoNames[_fn_main][_varEff]="h1Eff";
    histoNames[_fn_main][_varRho]="h1Rho";
    histoNames[_fn_main][_varAcc]="h1Acc";
    histoNames[_fn_main][_varEffAcc]="h1EffAcc";
    cmpFName1="dyee_test_El2skim_excludeGap.root";
    histoNames[_fn_cmp1][_varEff]="h1Eff";
    histoNames[_fn_cmp1][_varEff]="h1EffPU";
    histoNames[_fn_cmp1][_varAcc]="h1Acc";
    histoNames[_fn_cmp1][_varEffAcc]="h1EffPUAcc";
    cmpFName2="dyee_test_RECO_El2skim.root";
    histoNames[_fn_cmp2][_varRho]="h1rho";
  }
  else {
    std::cout << "not ready for theCase=" << theCase << "\n";
  }

  for (std::map<TVaried_t,TString>::const_iterator it= histoNames[_fn_main].begin();
       it!=histoNames[_fn_main].end(); it++) {
    std::cout << variedVarName(it->first) << " " << it->second << "\n";

    TString h1nameBase= "h1" + variedVarName(it->first);
    TH1D *h1main= loadHisto(mainFName,it->second, h1nameBase+"_main",h1dummy);
    if (!h1main) {
      std::cout << "failed to get h1main " << it->second << "\n";
      return;
    }

    TH1D *h1cmp1= NULL, *h1cmp2= NULL;

    for (std::map<TVaried_t,TString>::const_iterator it1=histoNames[_fn_cmp1].begin();
	 it1!=histoNames[_fn_cmp1].end(); it1++) {
      if (it1->first == it->first) {
	h1cmp1= loadHisto(cmpFName1,it1->second, h1nameBase+"_cmp1",h1dummy);
	if (!h1cmp1) {
	  std::cout << "failed to get h1cmp1 " << it1->second << "\n";
	}
	break;
      }
    }

    for (std::map<TVaried_t,TString>::const_iterator it2=histoNames[_fn_cmp2].begin();
	 it2!=histoNames[_fn_cmp2].end(); it2++) {
      if (it2->first == it->first) {
	h1cmp2= loadHisto(cmpFName2,it2->second, h1nameBase+"_cmp2",h1dummy);
	if (!h1cmp2) {
	  std::cout << "failed to get h1cmp2 " << it2->second << "\n";
	}
	break;
      }
    }

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
