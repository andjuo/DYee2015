#include "inputs.h"
#include <TLine.h>
#include <fstream>

// ----------------------------------------------------------------

const int nPt=5;
const double Ele_pt[nPt+1] = { 10.,20.,30.,40.,50.,200. };
const int nPtHLT=12;
const double Ele_ptHLT[nPtHLT+1] = {10,22,23,24,25,26,28,30,35,40,50,60,100};
const int nAbsEta=5;
const double SC_abseta[nAbsEta+1] = { 0,0.8,1.4442,1.566,2.0,2.5 };
const int nEta=10;
const double SC_eta[nEta+1] = { -2.5,-2.0,-1.566,-1.4442,-0.8,0,0.8,1.4442,1.566,2.0,2.5 };

void plotEffs(std::vector<TH2D*> &hV,
	      const std::vector<HistoStyle_t> &hsV,
	      const std::vector<TString> &labels,
	      TString plotNameBase, TString effName,
	      int allInOne=0);

void plotEffs_separate(std::vector<TH2D*> &hV,
		       const std::vector<HistoStyle_t> &hsV,
		       const std::vector<TString> &labels,
		       TString plotNameBase, TString effName);

// ----------------------------------------------------------------

void plotDYeeESF(TString theSet="ID", int save=0)
{
  TH2D *h2bin= NULL;
  if (theSet=="ID")
    h2bin=new TH2D("h2bin","h2bin;|#eta|;p_{T} [GeV]",
		   nAbsEta,SC_abseta,nPt,Ele_pt);
  else if (theSet=="RECO")
    h2bin=new TH2D("h2bin","h2bin;#eta;p_{T} [GeV]",
		   nEta,SC_eta,nPt,Ele_pt);
  else if (theSet=="HLT_finePt")
    h2bin=new TH2D("h2bin","h2bin;#eta;p_{T} [GeV]",
		   nEta,SC_eta,nPtHLT,Ele_ptHLT);
  else if (theSet=="HLT")
    h2bin=new TH2D("h2bin","h2bin;#eta;p_{T} [GeV]",
		   nEta,SC_eta,nPt,Ele_pt);
  else {
    std::cout << "h2bin is not ready for theSet=" << theSet << "\n";
    return;
  }

  TH2D* h2effData= cloneHisto(h2bin,"h2effData","h2effData");
  TH2D* h2effMC= cloneHisto(h2bin,"h2effMC","h2effMC");
  TH2D* h2effData_bkgPdf= cloneHisto(h2bin,"h2effData_bkgPdf","h2effData_bkgPdf");
  TH2D* h2effData_sigPdf= cloneHisto(h2bin,"h2effData_sigPdf","h2effData_sigPdf");
  TH2D* h2effMC_NLOvsLO= cloneHisto(h2bin,"h2effMC_NLOvsLO","h2effMC_NLOvsLO");
  TH2D* h2effData_tag= cloneHisto(h2bin,"h2effData_tag","h2effData_tag");
  TH2D* h2effData_Zrange= cloneHisto(h2bin,"h2effData_Zrange","h2effData_Zrange");

  std::vector<TH2D*> allHistosV;
  allHistosV.push_back(h2effData);
  allHistosV.push_back(h2effMC);
  allHistosV.push_back(h2effData_bkgPdf);
  allHistosV.push_back(h2effData_sigPdf);
  allHistosV.push_back(h2effMC_NLOvsLO);
  allHistosV.push_back(h2effData_tag);
  allHistosV.push_back(h2effData_Zrange);

  TString fname_RECO="dyee_eleRECO_76X.txt";
  TString fname_ID="dyee_CutBasedID_MediumWP_76X_18Feb.txt";
  TString fname_HLT_finePt="dyee_triggerEff_76X_finePt.txt";
  TString fname_HLT="dyee_triggerEff_76X.txt";
  TString effName="RECO";
  TString fname;
  if (theSet=="ID") { fname=fname_ID; effName="ID"; }
  else if (theSet=="RECO") { fname=fname_RECO; effName="RECO"; }
  else if (theSet=="HLT") { fname=fname_HLT; effName="HLT"; }
  else if (theSet=="HLT_finePt") { fname=fname_HLT_finePt; effName="HLT_finePt"; }
  else {
    std::cout << "code is not ready for theSet=" << theSet << "\n";
    return;
  }

  std::ifstream fin(fname.Data());
  if (!fin.is_open()) {
    std::cout << "failed to open the file <" << fname << ">\n";
    return;
  }
  {
    std::string line;
    for (int i=0; i<2; i++) {
      getline(fin,line);
      if (line[0]!='#') {
	fin.close();
	std::cout << "skipped uncommented line. stopping\n";
	return;
      }
    }
  }

  double eta1,eta2,pt1,pt2, eff,deff;
  int iEta,iPt;
  while (!fin.eof()) {
    fin >> eta1 >> eta2 >> pt1 >> pt2;
    if (fin.rdstate()!=0) continue;
    iEta= h2bin->GetXaxis()->FindBin(0.5*(eta1+eta2));
    iPt = h2bin->GetYaxis()->FindBin(0.5*(pt1+pt2));
    std::cout << "eta1=" << eta1 << ", eta2=" << eta2 << ", iEta=" << iEta << "\n";
    std::cout << "pt1=" << pt1 << ", pt2=" << pt2 << ", iPt=" << iPt << "\n";

    fin >> eff >> deff;
    h2effData->SetBinContent(iEta,iPt, eff);
    h2effData->SetBinError(iEta,iPt, deff);

    fin >> eff >> deff;
    h2effMC->SetBinContent(iEta,iPt, eff);
    h2effMC->SetBinError(iEta,iPt, deff);

    std::cout << "data: " << h2effData->GetBinContent(iEta,iPt) << " +- "
	      << h2effData->GetBinError(iEta,iPt) << "\n";
    std::cout << "MC: " << h2effMC->GetBinContent(iEta,iPt) << " +- "
	      << h2effMC->GetBinError(iEta,iPt) << "\n";

    deff=0;

    if ((effName!="HLT") && (effName!="HLT_finePt")) {
      fin >> eff;
      h2effData_bkgPdf->SetBinContent(iEta,iPt, eff);
      h2effData_bkgPdf->SetBinError(iEta,iPt, deff);

      fin >> eff;
      h2effData_sigPdf->SetBinContent(iEta,iPt, eff);
      h2effData_sigPdf->SetBinError(iEta,iPt, deff);

      fin >> eff;
      h2effMC_NLOvsLO->SetBinContent(iEta,iPt, eff);
      h2effMC_NLOvsLO->SetBinError(iEta,iPt, deff);
    }

    fin >> eff;
    h2effData_tag->SetBinContent(iEta,iPt, eff);
    h2effData_tag->SetBinError(iEta,iPt, deff);

    std::cout << "dataTag: " << h2effData_tag->GetBinContent(iEta,iPt) << " +- "
	      << h2effData_tag->GetBinError(iEta,iPt) << "\n";

    if (effName=="ID") {
      fin >> eff >> deff; // load -1,-1
      deff=0;
      fin >> eff;
      h2effData_Zrange->SetBinContent(iEta,iPt, eff);
      h2effData_Zrange->SetBinError(iEta,iPt, deff);
    }
  }
  fin.close();
  
  HistoStyle_t hsData(kRed,  5,1,1.0,2.); hsData.SetStyle(h2effData);
  HistoStyle_t hsDataBkgPdf(kBlue,24,1,0.8,2.); hsDataBkgPdf.SetStyle(h2effData_bkgPdf);
  HistoStyle_t hsDataSigPdf(kGreen+1,20,1,0.8,2.); hsDataSigPdf.SetStyle(h2effData_sigPdf);
  HistoStyle_t hsMC(kBlack,3,1,1.0,2.); hsMC.SetStyle(h2effMC);
  HistoStyle_t hsMCNLO(kBlue,24,1,0.8,2.); hsMCNLO.SetStyle(h2effMC_NLOvsLO);
  HistoStyle_t hsDataTag(6 ,   26,1,0.8,2.); hsDataTag.SetStyle(h2effData_tag);
  HistoStyle_t hsDataZRange(46,   23,1,0.8,2.); hsDataZRange.SetStyle(h2effData_Zrange);

  if (1) {
    std::vector<TH2D*> h2_data_vs_mc;
    std::vector<HistoStyle_t> hs_data_vs_mc;
    std::vector<TString> label_data_vs_mc;
    h2_data_vs_mc.push_back(h2effMC); label_data_vs_mc.push_back("MC");
    hs_data_vs_mc.push_back(hsMC);
    h2_data_vs_mc.push_back(h2effData); label_data_vs_mc.push_back("data");
    hs_data_vs_mc.push_back(hsData);
    plotEffs(h2_data_vs_mc,hs_data_vs_mc,label_data_vs_mc,"Data_vs_MC",effName,0);
  }

  if (1) {
    std::vector<TH2D*> h2_data_var;
    std::vector<HistoStyle_t> hs_data_var;
    std::vector<TString> label_data_var;
    h2_data_var.push_back(h2effData); label_data_var.push_back("data");
    hs_data_var.push_back(hsData);
    h2_data_var.push_back(h2effData_bkgPdf); label_data_var.push_back("data: bkgPdf");
    hs_data_var.push_back(hsDataBkgPdf);
    h2_data_var.push_back(h2effData_sigPdf); label_data_var.push_back("data: sigPdf");
    hs_data_var.push_back(hsDataSigPdf);
    h2_data_var.push_back(h2effData_tag); label_data_var.push_back("data: tag def");
    hs_data_var.push_back(hsDataTag);
    if (effName=="ID") {
      h2_data_var.push_back(h2effData_Zrange); label_data_var.push_back("data: Zrange");
      hs_data_var.push_back(hsDataZRange);
    }
    plotEffs(h2_data_var,hs_data_var,label_data_var, "Data",effName,0);
  }

  if (1) {
    std::vector<TH2D*> h2_mc_var;
    std::vector<HistoStyle_t> hs_mc_var;
    std::vector<TString> label_mc_var;
    h2_mc_var.push_back(h2effMC); label_mc_var.push_back("MC");
    hs_mc_var.push_back(hsMC);
    h2_mc_var.push_back(h2effMC_NLOvsLO); label_mc_var.push_back("MC: NLO vs LO");
    hs_mc_var.push_back(hsMCNLO);
    plotEffs(h2_mc_var,hs_mc_var,label_mc_var, "MC",effName,0);
  }


  if (save) {
    if (theSet=="ID") {
      TH2D *h2signed= new TH2D("h2signed","h2signed",nEta,SC_eta,nPt,Ele_pt);
      h2signed->SetStats(0);
      unsigned int nHistos=allHistosV.size();
      for (unsigned int i=0; i<nHistos; i++) {
	const TH2D *h2eff= allHistosV[i];
	if (h2eff->Integral()==0) continue;
	std::cout << " working with " << h2eff->GetName() << "\n";
	TH2D *h2=cloneHisto(h2eff,h2eff->GetName() + TString("_inEta"),
			    h2eff->GetTitle());
	allHistosV.push_back(h2);
	for (int ibin=1; ibin<=h2signed->GetNbinsX(); ibin++) {
	  double eta= h2signed->GetXaxis()->GetBinLowEdge(ibin)
	    + 0.5*h2signed->GetXaxis()->GetBinWidth(ibin);
	  int iEta= h2->GetXaxis()->FindBin(abs(eta));
	  for (int jbin=1; jbin<=h2signed->GetNbinsY(); jbin++) {
	    h2signed->SetBinContent(ibin, jbin, h2->GetBinContent(iEta,jbin));
	    h2signed->SetBinError  (ibin, jbin, h2->GetBinError(iEta,jbin));
	    //std::cout << " set eta=" << eta << ", " << h2signed->GetYaxis()->GetBinLowEdge(jbin) << " .. " << (h2signed->GetYaxis()->GetBinLowEdge(jbin)+h2signed->GetYaxis()->GetBinWidth(jbin)) << " to  iEta=" << iEta << " eff=" << h2->GetBinContent(iEta,jbin) << " +- " << h2->GetBinError(iEta,jbin) << "\n";
	  }
	}
      }
    }

    TString foutName="dyee_effSyst_" + theSet + ".root";
    TFile fout(foutName,"RECREATE");
    if (!fout.IsOpen()) {
      std::cout << "failed to create a file <" << fout.GetName() << ">\n";
      return;
    }
    for (unsigned int ih=0; ih<allHistosV.size(); ih++) {
      TH2D *h2=allHistosV[ih];
      if (h2->Integral()==0) {
	std::cout << "histogram " << h2->GetName() << " is empty\n";
      }
      else {
	h2->Write();
	h2->SetDirectory(0);
	std::cout << "histogram " << h2->GetName() << " saved\n";
      }
    }
    fout.Close();
  }

}

// --------------------------------------------------------------

void plotEffs(std::vector<TH2D*> &h2V,
	      const std::vector<HistoStyle_t> &hsV,
	      const std::vector<TString> &labels,
	      TString plotNameBase, TString effName,
	      int allInOne)
{
  std::cout << "allInOne=" << allInOne << "\n";
  if (!allInOne) {
    plotEffs_separate(h2V,hsV,labels,plotNameBase,effName);
    return;
  }

  for (int iEta=0; iEta<h2V[0]->GetNbinsX(); iEta++) {
    TString tag=Form("_%s_iEta%d",effName.Data(),iEta);
    std::vector<TH1D*> h1V;
    for (unsigned int ih=0; ih<h2V.size(); ih++) {
      TH2D *h2= h2V[ih];
      if (h2V[ih]->Integral()==0) continue;
      TH1D *h1eff_iEta=	h2->ProjectionY(h2->GetName() + tag,iEta+1,iEta+1);
      h1eff_iEta->SetStats(0);
      hsV[ih].SetStyle(h1eff_iEta);
      h1V.push_back(h1eff_iEta);
    }

    TString cName=plotNameBase + tag;
    TH1D *h1= h1V[0];
    logAxis(h1,0,"","#varepsilon");
    h1->SetTitle(plotNameBase);
    h1->GetYaxis()->SetTitleSize(0.05);
    h1->GetYaxis()->SetTitleOffset(1.2);
    plotHisto(h1,cName,0,0,"LPE1",labels[0]);
    TCanvas *cx=findCanvas(cName);
    setLeftMargin(cx,0.15);
    moveLegend(cx,0.45,0.);
    for (unsigned int ih=1; ih<h1V.size(); ih++) {
      plotHistoSame(h1V[ih],cName,"LPE1",labels[ih]);
    }
  }
}

// --------------------------------------------------------------

void plotEffs_separate(std::vector<TH2D*> &h2V,
		       const std::vector<HistoStyle_t> &hsV,
		       const std::vector<TString> &labels,
		       TString plotNameBase, TString effName)
{
  std::cout << "plotEffs_separate: labels = ";
  for (unsigned int i=0; i<labels.size(); i++) {
    std::cout << " <" << labels[i] << ">";
  }
  std::cout << "\n";

  std::vector<TH1D*> h1mainV;
  std::vector<TH1D*> h1cmpV;
  std::vector<TH1D*> h1mecV; // for auto-range
  std::vector<TH1D*> h1ratioV;

  TString plotOutDir= "dir-plot-DYeeEff-" + effName;

  for (unsigned int ih=0; ih<h2V.size(); ih++) {

    TH2D *h2= h2V[ih];

    if (1 && (h2->Integral()!=0)) {
      TString cName=TString("canv") + h2->GetName();
      TH2D *h2Clone= cloneHisto(h2,
				h2->GetName() + TString("_log"),h2->GetTitle());
      h2Clone->SetStats(0);
      logAxis(h2Clone,4);
      TCanvas *cx= plotHisto(h2Clone,cName,0,1);
      cx->SetCanvasSize(800,600);
      cx->Update();
      if (gROOT->IsBatch()==true) {
	std::cout << "saving H2 canvas\n";
	SaveCanvas(cx,cx->GetName(),plotOutDir,1);
      }
      else {
	std::cout << "H2 canvas not saved\n";
      }
    }

    if ((ih>0) && (h2->Integral()!=0)) {
      TString h2name=effName + "__" + h2->GetName() + TString("_div_")
	+ h2V[0]->GetName();
      TString h2title=labels[ih];
      if (plotNameBase=="Data_vs_MC") h2title="Data vs MC";
      TH2D *h2ratio= cloneHisto(h2,h2name,h2title);
      h2ratio->Divide(h2V[0]);
      h2ratio->SetStats(0);
      logAxis(h2ratio,4);
      TCanvas *cx=plotHisto(h2ratio,"canv" + h2name,0,1);
      cx->SetCanvasSize(800,600);
      if (gROOT->IsBatch()==true) {
	std::cout << "saving H2 canvas\n";
	SaveCanvas(cx,cx->GetName(),plotOutDir);
      }
      else {
	std::cout << "H2 canvas not saved\n";
      }
    }


    for (int iEta=0; iEta<h2->GetNbinsX(); iEta++) {
      TString tag=Form("_%s_iEta%d",effName.Data(),iEta);
      std::vector<TH1D*> *h1V= (ih==0) ? &h1mainV : &h1cmpV;

      if ((ih>1) && (iEta==0)) {
	h1cmpV.clear(); h1ratioV.clear();
	h1mecV.clear();
	for (unsigned int ii=0; ii<h1mainV.size(); ii++)
	  h1mecV.push_back(h1mainV[ii]);
      }

      TH1D *h1eff_iEta=	h2->ProjectionY(h2->GetName() + tag,iEta+1,iEta+1);
      std::cout << "ih=" << ih << ", iEta=" << iEta << ", lineColor=" << h2->GetLineColor() << "\n";
      //copyStyle(h1eff_iEta, h2);
      hsV[ih].SetStyle(h1eff_iEta);
      h1eff_iEta->SetStats(0);
      h1V->push_back(h1eff_iEta);
      h1mecV.push_back(h1eff_iEta);
      //printHisto(h1eff_iEta);

      if (ih==0) continue;
      TString div= "_div_";
      TH1D *h1cmp= h1eff_iEta;
      TH1D *h1main= h1mainV[iEta];
      TH1D *h1ratio=cloneHisto(h1cmp,
			       h1cmp->GetName()+div+h1main->GetName(),
			       TString(";")+h1main->GetXaxis()->GetTitle()
			       +";"+labels[ih]+"/"+labels[0]);
      h1ratio->Divide(h1main);
      std::cout << "ih=" << ih << ", iEta=" << iEta << ", ratio lineColor=" << h1ratio->GetLineColor() << "\n";
      h1ratioV.push_back(h1ratio);
    }

    if (ih==0) continue;

    if (h2->Integral()==0) {
      std::cout << "histogram " << h2->GetName() << " is empty. Skipping\n";
      continue;
    }

    if (0) {
      std::cout << "work with label=" << labels[ih] << "\n";
      for (unsigned int i=0; i<h1mecV.size(); i++) {
	std::cout << " i = " << i << "  " << h1mecV[i]->GetName() << "\n";
      }
      for (unsigned int i=0; i<h1ratioV.size(); i++) {
	std::cout << " i = " << i << "  " << h1ratioV[i]->GetName() << "\n";
      }
      return;
    }

    double effRange_min=0.;
    double effRange_max=1.;
    double yrange_min=0.8;
    double yrange_max=1.2;

    if (1) {
      // auto range on ratio
      typedef std::pair<double,double> ddpair_t;
      std::vector<std::pair<double,double> > yRange, yRangeRatio;
      yRange.push_back(ddpair_t(0.90, 1.01));
      yRange.push_back(ddpair_t(0.80, 1.01));
      yRange.push_back(ddpair_t(0.70, 1.01));
      yRange.push_back(ddpair_t(0.60, 1.01));
      yRange.push_back(ddpair_t(0.50, 1.01));
      yRange.push_back(ddpair_t(0.40, 1.01));
      yRange.push_back(ddpair_t(0.30, 1.01));
      yRange.push_back(ddpair_t(0.20, 1.01));
      yRange.push_back(ddpair_t(0.10, 1.01));
      yRange.push_back(ddpair_t(0.00, 1.01));
      yRangeRatio.push_back(ddpair_t(0.95,1.05));
      yRangeRatio.push_back(ddpair_t(0.90,1.10));
      yRangeRatio.push_back(ddpair_t(0.85,1.15));
      yRangeRatio.push_back(ddpair_t(0.80,1.20));
      yRangeRatio.push_back(ddpair_t(0.75,1.25));
      yRangeRatio.push_back(ddpair_t(0.70,1.30));
      yRangeRatio.push_back(ddpair_t(0.65,1.35));
      yRangeRatio.push_back(ddpair_t(0.60,1.40));
      yRangeRatio.push_back(ddpair_t(0.50,1.50));

      if (!checkRange(h1mecV,effRange_min,effRange_max,yRange,1)) {
	std::cout << "auto range on efficiency failed\n";
      }
      if (!checkRange(h1ratioV,yrange_min,yrange_max,yRangeRatio,1)) {
	std::cout << "auto range on ratio failed\n";
      }
    }


    TString cmpTitle= labels[0] + " vs " + labels[ih];

    for (unsigned int iEta=0; iEta<h1mainV.size(); iEta++) {
      TString tag=Form("_%s_iEta%d",effName.Data(),iEta);
      TString cName=cmpTitle + tag;
      TCanvas *cx= new TCanvas(cName,cName,600,800);
      cx->Divide(1,2);
      TPad *pad1 = (TPad*)cx->GetPad(1);
      TPad *pad2 = (TPad*)cx->GetPad(2);

      const double padYDivPoint=0.29;
      const double dyPads=0.005;
      double ydiv= pad1->GetYlowNDC() - padYDivPoint * pad1->GetHNDC();
      pad1->SetPad(pad1->GetXlowNDC(), ydiv+dyPads,
		   pad1->GetXlowNDC() + pad1->GetWNDC(),
		   pad1->GetYlowNDC() + pad1->GetHNDC());
      pad2->SetPad(pad2->GetXlowNDC(), pad2->GetYlowNDC(),
		   pad2->GetXlowNDC() + pad2->GetWNDC(),
		   ydiv-dyPads);
      pad1->SetLogx();
      pad2->SetLogx();
      pad1->SetGrid(1,1);
      pad2->SetGrid(1,1);
      pad1->SetBottomMargin(0.015);
      pad2->SetTopMargin(0.1);
      pad2->SetBottomMargin(0.18);
      cx->SetBottomMargin(0.015);

      TString formatNoGap="%3.1lf #leq |#eta| #leq %3.1lf";
      TString formatPreGap="%3.1lf #leq |#eta| #leq %5.3lf";
      TString formatGap="%5.3lf #leq |#eta| #leq %5.3lf";
      TString formatPostGap="%5.3lf #leq |#eta| #leq %3.1lf";
      TString etaFormat=formatNoGap;
      double absEta1=h2->GetXaxis()->GetBinLowEdge(iEta+1);
      double absEta2=fabs(absEta1 + h2->GetXaxis()->GetBinWidth(iEta+1));
      if (absEta1<0) absEta1=-absEta1;
      std::cout << "absEta1=" << absEta1 << ", absEta2=" << absEta2 << "  ";
      if (((fabs(absEta1-1.444)<1e-3) && (fabs(absEta2-1.566)<1e-3)) ||
	  ((fabs(absEta1-1.566)<1e-3) && (fabs(absEta2-1.444)<1e-3)))
	{
	  std::cout << "gap\n";
	  etaFormat= formatGap;
	}
      else if ((fabs(absEta2-1.444)<1e-3) ||
	       (fabs(absEta2-1.566)<1e-3)){
	std::cout << "preGap\n";
	etaFormat= formatPreGap;
      }
      else if ((fabs(absEta1-1.444)<1e-3) ||
	       (fabs(absEta1-1.566)<1e-3)) {
	std::cout << "postGap\n";
	etaFormat= formatPostGap;
      }
      if (effName!="ID") etaFormat.ReplaceAll("|#eta|","#eta");
      TString etaRange=Form(etaFormat,
			    h2->GetXaxis()->GetBinLowEdge(iEta+1),
			    h2->GetXaxis()->GetBinLowEdge(iEta+1)
			    + h2->GetXaxis()->GetBinWidth(iEta+1));

      TH1D *h1main= h1mainV[iEta];
      logAxis(h1main,0,"","#varepsilon_{" + effName + "}");
      h1main->SetTitle(plotNameBase + "  " + etaRange);
      h1main->GetXaxis()->SetLabelOffset(0.15);
      h1main->GetYaxis()->SetTitleSize(0.06);
      h1main->GetYaxis()->SetTitleOffset(0.8);
      h1main->GetYaxis()->SetLabelSize(0.05);
      h1main->GetYaxis()->SetLabelOffset(0.01);
      h1main->GetYaxis()->SetRangeUser(effRange_min,effRange_max);

      TH1D *h1cmp= h1cmpV[iEta];
      copyStyle(h1cmp,h2); // if canvas does not fit, auto-resize changes color!
      //hsV[ih].SetStyle(h1cmp);
      std::cout << "h1cmp color=" << h1cmp->GetMarkerColor() << "\n";

      pad1->cd();
      h1main->Draw("LPE1");
      h1cmp->Draw("LPE1 same");
      checkRange(h1main, effRange_min,effRange_max);
      checkRange(h1cmp, effRange_min,effRange_max);

      double dxleg=0.45;
      TLegend *leg= new TLegend(0.15+dxleg,0.15,0.40+dxleg,0.21+0.06);
      leg->SetName("myLegend");
      leg->AddEntry(h1main,labels[0]);
      leg->AddEntry(h1cmp,labels[ih]);
      leg->Draw();

      pad2->cd();
      TH1D *h1ratio= h1ratioV[iEta];
      copyStyle(h1ratio,h2);
      logAxis(h1ratio,1);
      h1ratio->GetXaxis()->SetTitleSize(0.08);
      h1ratio->GetXaxis()->SetTitleOffset(1.);
      h1ratio->GetXaxis()->SetLabelSize(0.08);
      h1ratio->GetXaxis()->SetLabelOffset(0.02);
      h1ratio->GetYaxis()->SetTitleSize(0.08);
      h1ratio->GetYaxis()->SetLabelSize(0.08);
      h1ratio->GetYaxis()->SetTitleOffset(0.65);
      h1ratio->GetYaxis()->SetRangeUser(yrange_min,yrange_max);
      h1ratio->Draw("LPE1");
      std::cout << "h1ratio color=" << h1ratio->GetMarkerColor() << "\n";
      checkRange(h1ratio, yrange_min,yrange_max);


      std::cout << "names: "<< h1main->GetName() << ", " << h1cmp->GetName() << ", " << h1ratio->GetName() << "\n";

      int nBins=h1ratio->GetNbinsX();
      TLine *line= new TLine(h1ratio->GetBinLowEdge(1), 1.,
	     h1ratio->GetBinLowEdge(nBins)+h1ratio->GetBinWidth(nBins),1.);
      line->SetLineStyle(2);
      line->SetLineColor(kBlue);
      line->Draw();
      h1ratio->Draw("LPE1 same");

      cx->Update();
      //setLeftMargin(cx,0.15);
      //moveLegend(cx,0.45,0.);

      if (gROOT->IsBatch()==true) {
	std::cout << "saving canvas\n";
	TString outName=cx->GetName();
	outName.ReplaceAll(" ","_");
	outName.ReplaceAll(":","_");
	SaveCanvas(cx,outName,plotOutDir);
      }
      else {
	std::cout << "canvas not saved\n";
      }
    }
  }
}

// --------------------------------------------------------------
// --------------------------------------------------------------