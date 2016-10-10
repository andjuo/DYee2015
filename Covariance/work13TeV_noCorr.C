#include "inputs.h"
#include "crossSection.h"
#include "Blue.h"

void work13TeV_noCorr()
{
  TString eeCSFName="cs_DYee_13TeV_El3.root";
  TString mmCSFName="cs_DYmm_13TeVMuApproved_cs.root";

  TH1D* h1csEE=loadHisto(eeCSFName, "h1PreFSRCS", "h1csEE", 1,h1dummy);
  TH1D* h1csMM=loadHisto(mmCSFName, "h1CS", "h1csMM", 1,h1dummy);
  TH1D* h1csTheory=perMassBinWidth(loadHisto("theory13TeVmm.root",
					     "h1cs_theory", "h1cs_theory",
					     1,h1dummy));

  TString massStr= niceMassAxisLabel(2,"");
  TString eemassStr= niceMassAxisLabel(0,"");
  TString mmmassStr= niceMassAxisLabel(1,"");
  TString sigmaStr= niceMassAxisLabel(2,"",1);
  TString axisStr= ";" + massStr + ";" + sigmaStr;

  histoStyle(h1csEE, kGreen+1, 24, 1, 0.8);
  histoStyle(h1csMM, 46, 20, 1, 0.8);
  h1csTheory->SetStats(0);

  if (1) {
    TH1D *h1frame= cloneHisto(h1csEE,"h1frame_inpCS","frame");
    h1frame->Reset();
    logAxis(h1frame,1+0*2,massStr,sigmaStr);
    h1frame->SetTitle("input cross section");
    h1frame->GetYaxis()->SetTitleOffset(1.5);
    h1frame->GetYaxis()->SetRangeUser(1e-8,1e3);
    TCanvas *cx=plotHisto(h1frame,"cCScheck",1,1, "");
    plotHistoSame(h1csEE,"cCScheck","LPE1", "DY#rightarrowee");
    setLeftMargin(cx,0.15); moveLegend(cx,0.05,0.);
    plotHistoSame(h1csMM,"cCScheck", "LPE1", "DY#rightarrow#mu#mu");
    plotHistoSame(h1csTheory,"cCScheck","LP", "theory");
    //return;
  }

  if (1) {
    TGraphErrors *grCSEE= createGraph(h1csEE,"grEE", -1);
    TGraphErrors *grCSMM= createGraph(h1csMM,"grMM",  1);
    graphStyle(grCSEE, kGreen+1, 24, 1, 0.8);
    graphStyle(grCSMM, 46, 20, 1, 0.8);
    for (int i=0; i<5; i++) {
      TString canvName="";
      TCanvas *cx= createMassFrame(i,"cCSCheckRange_",
				   "input cross section" + axisStr,
				   &canvName);
      //plotHistoSame(h1csEE,canvName,"LPE1", "DY#rightarrowee");
      //plotHistoSame(h1csMM,canvName,"LPE1", "DY#rightarrow#mu#mu");
      plotGraphSame(grCSEE,canvName,"PE1", "DY#rightarrowee");
      plotGraphSame(grCSMM,canvName,"PE1", "DY#rightarrow#mu#mu");
      plotHistoSame(h1csTheory,canvName,"LP", "theory");
      if (i==1) moveLegend(cx,0.,0.55);
      if (i==4) moveLegend(cx,0.05,0.);
    }
  }

}
