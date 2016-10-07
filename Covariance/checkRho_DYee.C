#include "inputs.h"
#include "DYmm13TeV_eff.h"

void checkRho_DYee()
{
  TString inpAvgRhoFName500="dir-RhoEl3/dyee_rhoRndVec_El3_500.root";
  TString inpAvgRhoFName1000="dir-RhoEl3/dyee_rhoRndVec_El3_1000.root";
  TString rhoInpName="h1rho_avg";
  TH1D *h1rho500= NULL;
  if (0) {
    h1rho500=loadHisto(inpAvgRhoFName500,rhoInpName,"h1rho_avg500",
		       h1dummy);
    if (!h1rho500) return;
  }
  TH1D *h1rho1000= NULL;
  if (1) {
    h1rho1000=loadHisto(inpAvgRhoFName1000,rhoInpName,"h1rho_avg1000",
			h1dummy);
    if (!h1rho1000) return;
  }

  TString inpFName="dyee_test_dressed_El3_chkRho.root";
  TString h1nomName="h1_postFsrInAccSel_MWPURho";
  TString h1denomName="h1_postFsrInAccSel_MWPU";
  TString h1rhoName="h1rho";
  //TString 
  DYTnPEff_t tnpEff(2);
  EventSpace_t esPostFsr;

  TFile fin(inpFName);
  if (!fin.IsOpen()) {
    std::cout << "failed to open the file <" << fin.GetName() << ">\n";
    return;
  }
  TH1D *h1rhoNom= loadHisto(fin,h1nomName,h1nomName,1,h1dummy);
  TH1D *h1rhoDenom= loadHisto(fin,h1denomName,h1denomName,1,h1dummy);
  TH1D *h1rho= loadHisto(fin,h1rhoName,h1rhoName,1,h1dummy);
  if (!h1rhoNom || !h1rhoDenom || !h1rho) {
    std::cout << "failed to load histos\n";
    return;
  }
  if (!esPostFsr.load(fin, "mainES") ||
      !tnpEff.load(fin,"DYTnPEff")) {
    std::cout << "failed to load esPostFsr or tnpEff\n";
    return;
  }
  fin.Close();

  TH1D *h1rhoHere= esPostFsr.calculateScaleFactor(tnpEff,0,"h1rhoHere",
						  "scale factor (here);M_{ee} [GeV];#rho");
  printHisto(h1rhoHere);

  TH1D *h1rhoFromDiv= cloneHisto(h1rhoNom,"h1rhoFromDiv","h1rhoFromDiv");
  h1rhoFromDiv->Divide(h1rhoDenom);
  h1rhoFromDiv->SetStats(0);

  if (1 && h1rho1000) {
    for (int ibin=1; ibin<=h1rho1000->GetNbinsX(); ibin++) {
      h1rho1000->SetBinContent(ibin,h1rho->GetBinContent(ibin));
    }
  }


  histoStyle(h1rho,46,7,1);
  histoStyle(h1rhoHere,1,20,1);
  histoStyle(h1rhoFromDiv,kBlue,24,1);
  if (h1rho500) histoStyle(h1rho500,kGreen+1,24,1);
  if (h1rho1000) histoStyle(h1rho1000,kRed,7,1);
  h1rhoFromDiv->SetTitle("Efficiency scale factor");
  h1rhoFromDiv->GetXaxis()->SetNoExponent(true);
  h1rhoFromDiv->GetXaxis()->SetMoreLogLabels(true);
  h1rhoFromDiv->GetYaxis()->SetTitle("#rho");
  h1rhoFromDiv->GetYaxis()->SetDecimals(true);

  TCanvas *c=plotHisto(h1rhoFromDiv,"cRho",1,0,"LPE1","rho from TH1D division");
  c->SetCanvasSize(800,800);
  if (h1rho500) plotHistoSame(h1rho500,"cRho","LPE1","rho avg 500");
  if (h1rho1000) plotHistoSame(h1rho1000,"cRho","LPE1","rho avg 1000");
  plotHistoSame(h1rhoHere,"cRho","LPE1","rho with err (from ES)");
  //plotHistoSame(h1rho,"cRho","LPE1","rho from macro w/o err");

  h1rhoHere->SetStats(0);
  plotHisto(h1rhoHere,"cRhoMCStat",1,0,"LPE1","rho with err (from ES)");
}
