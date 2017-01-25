#include "inputs.h"

void plotYield(int ele=1)
{
  if (ele==1) {
    TString fname="cs_DYee_13TeV_El3.root";
    TH1D *h1Yield=loadHisto(fname,"h1Yield","h1Yield",h1dummy);
    TH1D *h1Signal=loadHisto(fname,"h1Signal","h1Signal",h1dummy);
    hsBlue.SetStyle(h1Yield);
    h1Yield->SetMarkerStyle(5);
    hsGreen.SetStyle(h1Signal);
    h1Yield->GetYaxis()->SetRangeUser(1e-2,1e6);
    h1Yield->GetYaxis()->SetTitleOffset(1.5);
    logAxis(h1Yield,1,niceMassAxisLabel(0,"",0),"count","");
    plotHisto(h1Yield,"cYield_EE",1,1,"LPE","yield");
    plotHistoSame(h1Signal,"cYield_EE","LPE1","signal");
    //printHisto(h1Yield);
    //printHisto(h1Signal);
  }
  else if ((ele==0) || (ele==-1)) {
    TString tag=(ele==0) ? "csA" : "csB";
    TString fname="cs_DYmm_13TeVMuApproved_" + tag + ".root";
    TH1D *h1Yield=loadHisto(fname,"h1Yield","h1Yield",h1dummy);
    TH1D *h1Signal=loadHisto(fname,"h1Signal","h1Signal",h1dummy);
    hsBlue.SetStyle(h1Yield);
    h1Yield->SetMarkerStyle(5);
    hsGreen.SetStyle(h1Signal);
    h1Yield->GetYaxis()->SetRangeUser(1e-2,1e6);
    h1Yield->GetYaxis()->SetTitleOffset(1.5);
    logAxis(h1Yield,1,niceMassAxisLabel(1,"",0),"count","");
    plotHisto(h1Yield,"cYield_MM_"+tag,1,1,"LPE","yield " + tag);
    plotHistoSame(h1Signal,"cYield_MM_"+tag,"LPE1","signal " + tag);
  }
}
