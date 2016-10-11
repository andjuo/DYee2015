#include "inputs.h"
#include "Blue.h"

void plotDYSignal(int print=0)
{
  TH1D *h1signalMM4p2= loadHisto("cs_DYmm_13TeVMuApproved_csA.root","h1Signal",
				 "h1signalMM4p2",1,h1dummy);
  TH1D *h1signalMM4p3= loadHisto("cs_DYmm_13TeVMuApproved_csB.root","h1Signal",
				 "h1signalMM4p2",1,h1dummy);
  TH1D *h1signalEE= loadHisto("cs_DYee_13TeV_El3.root","h1Signal",
			      "h1signalEE",1,h1dummy);

  double lumiMM4p2= 865.919;
  double lumiMM4p3= 2832.673 - lumiMM4p2;
  double lumiEE= 2316.969;

  h1signalMM4p2->Scale(1/lumiMM4p2);
  h1signalMM4p3->Scale(1/lumiMM4p3);
  h1signalEE->Scale(1/lumiEE);

  TString massStr= niceMassAxisLabel(2,"");
  //TString eemassStr= niceMassAxisLabel(0,"");
  //TString mmmassStr= niceMassAxisLabel(1,"");
  TString axisStr= ";" + massStr + ";signal yield/Lumi";

  TGraphErrors *grEE= createGraph(h1signalEE,"grEE", 0);
  TGraphErrors *grMM4p2= createGraph(h1signalMM4p2, "grMM4p2", -1);
  TGraphErrors *grMM4p3= createGraph(h1signalMM4p3, "grMM4p3",  1);

  printHisto(h1signalEE);
  printHisto(h1signalMM4p2);

  graphStyle(grEE, kGreen+1, 7, 1, 0.8);
  graphStyle(grMM4p2, 46, 24, 1, 0.8);
  graphStyle(grMM4p3, kBlack, 20, 1, 0.8);

  std::vector<TCanvas*> canvasV;
  for (int i=0; i<5; i++) {
    TString canvName="";
    TH2D *h2frame=NULL;
    TCanvas *cx= createMassFrame(i,"cSignal_",
				 "signal" + axisStr, 2,
				 &canvName,&h2frame);
    canvasV.push_back(cx);
    //plotHistoSame(h1csEE,canvName,"LPE1", "DY#rightarrowee");
    //plotHistoSame(h1csMM,canvName,"LPE1", "DY#rightarrow#mu#mu");
    plotGraphSame(grEE,canvName,"PE1", "DY#rightarrowee");
    plotGraphSame(grMM4p2,canvName,"PE1", "DY#rightarrow#mu#mu (HLT4.2)");
    plotGraphSame(grMM4p3,canvName,"PE1", "DY#rightarrow#mu#mu (HLT4.3)");
    if (i==1) {
      h2frame->GetYaxis()->SetTitleOffset(1.5);
      setLeftMargin(cx,0.11);
      moveLegend(cx,0.45,0.);
    }
    if (i==4) {
      h2frame->GetYaxis()->SetNoExponent(false);
      h2frame->GetYaxis()->SetTitleOffset(1.5);
      setLeftMargin(cx,0.15);
      moveLegend(cx,0.45,0.55);
    }
    if (i==5) moveLegend(cx,0.45,0.55);
  }


  if (print==1) {
    TFile fout("foutCanvas_DYSignal.root","RECREATE");
    SaveCanvases(canvasV,"dir-plot-DYSignal",&fout);
    writeTimeTag();
    fout.Close();
    std::cout << "file <" << fout.GetName() << "> created\n";
  }
}
