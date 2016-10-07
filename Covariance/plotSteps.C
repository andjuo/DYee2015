#include "inputs.h"
//#include "latex_output.h"

void checkInputs();
void addRatio(const TH1D *h1a, const TH1D *h1b, TString label,
	      std::vector<TH1D*> &h1RatioV, std::vector<TString> &labelV);
void plotRatios(std::vector<TH1D*> &h1RatioV, const std::vector<TString> &labelV);

// --------------------------------------------

void plotSteps()
{

  //checkInputs(); return;

  TFile fin("/mnt/sdc/andriusj/DY13TeV/DYanalysis-20160817/ElectronNtupler/test/Analysis_Codes/AccEff/dyee_preFSR_forAccEff_v1steps.root");
  TH1D *h1effPass_new= loadHisto(fin,"h1_eff_sumPass","h1effPass_new",1,h1dummy);
  TH1D *h1effPass_noPU= loadHisto(fin,"h1_eff_sumPass_noPU","h1effPass_noPU",1,h1dummy);
  TH1D *h1effPass_noTrigObjMatching= loadHisto(fin,"h1_eff_sumPass_noTrigObjMatching","h1effPass_noTrigObjMatching",1,h1dummy);
  TH1D *h1effPass_matchDR3= loadHisto(fin,"h1_eff_sumPass_matchDR3","h1effPass_matchDR3",1,h1dummy);
  TH1D *h1effPass_noTrig= loadHisto(fin,"h1_eff_sumPass_noTrig","h1effPass_noTrig",1,h1dummy);
  TH1D *h1effPass_noTrigNoPU= loadHisto(fin,"h1_eff_sumPass_noTrigNoPU","h1effPass_noTrigNoPU",1,h1dummy);
  //TH1D *h1effTot_new= loadHisto(fin, "h1_eff_sumTot", "h1effTot_new",1,h1dummy);
  fin.Close();
  if (!h1effPass_new  || !h1effPass_noPU || !h1effPass_noTrigObjMatching ||
      !h1effPass_matchDR3 || !h1effPass_noTrig || !h1effPass_noTrigNoPU) {
    std::cout << "null ptr\n";
    return;
  }

  std::vector<TH1D*> h1InfoV;
  std::vector<TString> infoLabelV;

  histoStyle(h1effPass_noPU,kRed,5,1);
  histoStyle(h1effPass_noTrig,kBlue,7,2);
  histoStyle(h1effPass_matchDR3,kGreen+1,5,1);
  histoStyle(h1effPass_noTrigObjMatching,kOrange,20,2);

  std::cout << "\nwPU effect on new effPass distribution\n";
  plotHisto(h1effPass_new,"cPUnew",1,0,"LPE1","new wPU");
  plotHistoSame(h1effPass_noPU,"cPUnew","LPE1","new noPU");
  printRatio(h1effPass_new,h1effPass_noPU);
  addRatio(h1effPass_noPU,h1effPass_new,"no PU rew.",h1InfoV,infoLabelV);

  std::cout << "\nwPU effect on old effPass distribution\n";
  plotHisto(h1effPass_noTrig,"cPUold",1,0,"LPE1","old wPU");
  plotHistoSame(h1effPass_noTrigNoPU,"cPUold","LPE1","old noPU");
  printRatio(h1effPass_noTrig,h1effPass_noTrigNoPU);

  std::cout << "\ntriggering\n";
  plotHisto(h1effPass_new,"cTrig",1,0,"LPE1","new wPU");
  plotHistoSame(h1effPass_noTrig,"cTrig","LPE1","no event trigger");
  printRatio(h1effPass_new,h1effPass_noTrig);
  addRatio(h1effPass_noTrig,h1effPass_new,"no trigger",h1InfoV,infoLabelV);

  std::cout << "\ntrigger object matching\n";
  plotHisto(h1effPass_new,"cTrigObjMatch",1,0,"LPE1","new wPU");
  plotHistoSame(h1effPass_noTrigObjMatching,"cTrigObjMatch","LPE1","new, noTrigObjMatch");
  printRatio(h1effPass_new,h1effPass_noTrigObjMatching);
  addRatio(h1effPass_noTrigObjMatching,h1effPass_new,"no trig.obj.match.",h1InfoV,infoLabelV);

  std::cout << "\nDR effect on new effPass distribution\n";
  plotHisto(h1effPass_new,"cDR",1,0,"LPE1","new wPU, match DR1");
  plotHistoSame(h1effPass_matchDR3,"cDR","LPE1","new wPU, match DR3");
  printRatio(h1effPass_new,h1effPass_matchDR3);
  addRatio(h1effPass_matchDR3,h1effPass_new,"#DeltaR<0.3",h1InfoV,infoLabelV);

  //printLatex("effPassSteps.tex",h1InfoV,infoLabelV,1);
  plotRatios(h1InfoV,infoLabelV);
}

// --------------------------------------------

void checkInputs()
{
  TFile fin("/mnt/sdc/andriusj/DY13TeV/DYanalysis-20160817/ElectronNtupler/test/Analysis_Codes/AccEff/dyee_preFSR_forAccEff_v1steps.root");
  TH1D *h1effPass_new= loadHisto(fin,"h1_eff_sumPass","h1effPass_new",1,h1dummy);
  TH1D *h1effPass_chk= loadHisto(fin,"h1_eff_sumPass_chk","h1effPass_chk",1,h1dummy);
  TH1D *h1effPass_noTrigNoPU= loadHisto(fin,"h1_eff_sumPass_noTrigNoPU","h1effPass_noTrigNoPU",1,h1dummy);
  TH1D *h1effTot_new= loadHisto(fin, "h1_eff_sumTot", "h1effTot_new",1,h1dummy);
  fin.Close();
  if (!h1effPass_new || !h1effPass_chk || !h1effPass_noTrigNoPU || !h1effTot_new) return;

  TFile fin2("/mnt/sdc/andriusj/DY13TeV/DYanalysis-20160817/ElectronNtupler/test/Analysis_Codes/AccEff/dyee_preFSR_forAccEff_v1.root");
  TH1D *h1effPass_old= loadHisto(fin2,"h1_eff_sumPass","h1effPass_old",1,h1dummy);
  TH1D *h1effTot_old= loadHisto(fin2,"h1_eff_sumTot","h1effTot_old",1,h1dummy);
  fin2.Close();
  if (!h1effPass_old || !h1effTot_old) return;

  histoStyle(h1effPass_chk,kRed,5,1);
  histoStyle(h1effPass_old,kBlue,7,1);

  plotHisto(h1effPass_new,"cCheck",1,0,"LPE1","new effPass");
  plotHistoSame(h1effPass_chk,"cCheck","LPE1","new effPass (control)");
  printRatio(h1effPass_chk,h1effPass_new);

  //plotHisto(h1effPass_noTrigNoPU,"cCmpOld",1,0,"LPE1","noTrigNoPU");
  //plotHistoSame(h1effPass_old,"cCmpOld","LPE1","old effPass");
  //printRatio(h1effPass_old,h1effPass_noTrigNoPU);

  plotHisto(h1effTot_new,"cCmpTot",1,0,"LPE1","new effTot");
  plotHistoSame(h1effTot_old,"cCmpTot","LPE1","old effTot");
  printRatio(h1effTot_old,h1effTot_new);
}

// --------------------------------------------

void addRatio(const TH1D *h1a, const TH1D *h1b, TString label,
	      std::vector<TH1D*> &hRatioV, std::vector<TString> &labelV)
{
  TH1D *h1r= cloneHisto(h1a,TString(h1a->GetName()) + "clone",h1a->GetTitle());
  h1r->Divide(h1b);
  hRatioV.push_back(h1r);
  labelV.push_back(label);
}

// --------------------------------------------

void plotRatios(std::vector<TH1D*> &h1RatioV, const std::vector<TString> &labelV)
{
  for (unsigned int ih=0; ih<h1RatioV.size(); ih++) {
    removeError(h1RatioV[ih]);
    h1RatioV[ih]->SetStats(0);
  }

  TH1D *h1=h1RatioV[0];
  h1->SetTitle("Ratio of the modified selection to the new selection");
  h1->SetTitleSize(24);
  h1->GetXaxis()->SetNoExponent();
  h1->GetXaxis()->SetMoreLogLabels();
  h1->GetXaxis()->SetTitleSize(0.05);
  h1->GetYaxis()->SetRangeUser(0.95,1.15);
  h1->GetYaxis()->SetTitle("ratio");
  h1->GetYaxis()->SetTitleSize(0.05);
  h1->GetYaxis()->SetTitleOffset(1.38);

  TString cName="cRatio";
  TCanvas *c=plotHisto(h1RatioV[0],cName,1,0,"LP",labelV[0]);
  gPad->SetGridx();
  gPad->SetGridy();
  for (unsigned int ih=1; ih<h1RatioV.size(); ih++) {
    plotHistoSame(h1RatioV[ih],cName,"LP",labelV[ih]);
  }
}

// --------------------------------------------
