#include "inputs.h"
#include <TLine.h>
#include <fstream>

// ----------------------------------------------------------------

const int nPt=5;
const double Ele_pt[nPt+1] = { 10.,20.,30.,40.,50.,200. };
const int nEta=5;
const double SC_abseta[nEta+1] = { 0,0.8,1.4442,1.566,2.0,2.5 };

void plotEffs(std::vector<TH2D*> &hV,
	      const std::vector<HistoStyle_t> &hsV,
	      const std::vector<TString> &labels,
	      TString plotNameBase,
	      int allInOne=0);

void plotEffs_separate(std::vector<TH2D*> &hV,
		       const std::vector<HistoStyle_t> &hsV,
		       const std::vector<TString> &labels,
		       TString plotNameBase);

void plotEffs(TH2D *h2main, TString labelMain,
	      TH2D *h2cmp, TString labelCmp,
	      TString plotNameBase);

// ----------------------------------------------------------------

void plotDYeeESF()
{
  TH2D *h2bin= new TH2D("h2bin","h2bin;|#eta|;p_{T} [GeV]",
			nEta,SC_abseta,nPt,Ele_pt);

  TH2D* h2effData= cloneHisto(h2bin,"h2effData","h2effData");
  TH2D* h2effMC= cloneHisto(h2bin,"h2effMC","h2effMC");
  TH2D* h2effData_bkgPdf= cloneHisto(h2bin,"h2effData_bkgPdf","h2effData_bkgPdf");
  TH2D* h2effData_sigPdf= cloneHisto(h2bin,"h2effData_sigPdf","h2effData_sigPdf");
  TH2D* h2effMC_NLOvsLO= cloneHisto(h2bin,"h2effMC_NLOvsLO","h2effMC_NLOvsLO");
  TH2D* h2effData_tag= cloneHisto(h2bin,"h2effData_tag","h2effData_tag");
  TH2D* h2effData_Zrange= cloneHisto(h2bin,"h2effData_Zrange","h2effData_Zrange");

  const char *fname="dyee_CutBasedID_MediumWP_76X_18Feb.txt";
  std::ifstream fin(fname);
  if (!fin.is_open()) {
    std::cout << "failed to open the file <" << fname << ">\n";
    return;
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

    deff=0;
    fin >> eff;
    h2effData_bkgPdf->SetBinContent(iEta,iPt, eff);
    h2effData_bkgPdf->SetBinError(iEta,iPt, deff);

    fin >> eff;
    h2effData_sigPdf->SetBinContent(iEta,iPt, eff);
    h2effData_sigPdf->SetBinError(iEta,iPt, deff);

    fin >> eff;
    h2effMC_NLOvsLO->SetBinContent(iEta,iPt, eff);
    h2effMC_NLOvsLO->SetBinError(iEta,iPt, deff);

    fin >> eff;
    h2effData_tag->SetBinContent(iEta,iPt, eff);
    h2effData_tag->SetBinError(iEta,iPt, deff);

    fin >> eff >> deff; // load -1,-1
    deff=0;
    fin >> eff;
    h2effData_Zrange->SetBinContent(iEta,iPt, eff);
    h2effData_Zrange->SetBinError(iEta,iPt, deff);
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
    plotEffs(h2_data_vs_mc,hs_data_vs_mc,label_data_vs_mc,"Data_vs_MC",0);
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
    h2_data_var.push_back(h2effData_Zrange); label_data_var.push_back("data: Zrange");
    hs_data_var.push_back(hsDataZRange);
    plotEffs(h2_data_var,hs_data_var,label_data_var, "Data",0);
  }

  if (1) {
    std::vector<TH2D*> h2_mc_var;
    std::vector<HistoStyle_t> hs_mc_var;
    std::vector<TString> label_mc_var;
    h2_mc_var.push_back(h2effMC); label_mc_var.push_back("MC");
    hs_mc_var.push_back(hsMC);
    h2_mc_var.push_back(h2effMC_NLOvsLO); label_mc_var.push_back("MC: NLO vs LO");
    hs_mc_var.push_back(hsMCNLO);
    plotEffs(h2_mc_var,hs_mc_var,label_mc_var, "MC",0);
  }
}

// --------------------------------------------------------------

void plotEffs(std::vector<TH2D*> &h2V,
	      const std::vector<HistoStyle_t> &hsV,
	      const std::vector<TString> &labels,
	      TString plotNameBase,
	      int allInOne)
{
  std::cout << "allInOne=" << allInOne << "\n";
  if (!allInOne) {
    plotEffs_separate(h2V,hsV,labels,plotNameBase);
    return;
  }

  for (int iEta=0; iEta<nEta; iEta++) {
    TString tag=Form("_iEta%d",iEta);
    std::vector<TH1D*> h1V;
    for (unsigned int ih=0; ih<h2V.size(); ih++) {
      TH2D *h2= h2V[ih];
      TH1D *h1eff_iEta=	h2->ProjectionY(h2->GetName() + tag,iEta+1,iEta+2);
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
		       TString plotNameBase)
{

  std::vector<TH1D*> h1mainV;
  std::vector<TH1D*> h1cmpV;
  std::vector<TH1D*> h1ratioV;

  for (unsigned int ih=0; ih<h2V.size(); ih++) {

    TH2D *h2= h2V[ih];
    for (int iEta=0; iEta<nEta; iEta++) {
      TString tag=Form("_iEta%d",iEta);
      std::vector<TH1D*> *h1V= (ih==0) ? &h1mainV : &h1cmpV;
      if (ih>1) { h1cmpV.clear(); h1ratioV.clear(); }

      TH1D *h1eff_iEta=	h2->ProjectionY(h2->GetName() + tag,iEta+1,iEta+1);
      std::cout << "ih=" << ih << ", iEta=" << iEta << ", lineColor=" << h2->GetLineColor() << "\n";
      //copyStyle(h1eff_iEta, h2);
      hsV[ih].SetStyle(h1eff_iEta);
      h1eff_iEta->SetStats(0);
      h1V->push_back(h1eff_iEta);
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

    double effRange_min=0.;
    double effRange_max=1.;
    double yrange_min=0.8;
    double yrange_max=1.2;

    if (1) {
      // auto range on ratio
      const int nRatioRanges=8;
      const double dRatio[nRatioRanges]= { 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5 };
      int irr=0, ok=1;
      for (irr=0; irr<nRatioRanges; irr++) {
	for (unsigned int iih=0; ok && (iih<h1ratioV.size()); iih++) {
	  ok= checkRange(h1ratioV[iih], 1-dRatio[irr],1+dRatio[irr],1);
	}
	if (ok) break;
      }
      if (ok) {
	if (irr+1<nRatioRanges) irr++;
	yrange_min= 1-dRatio[irr];
	yrange_max= 1+dRatio[irr];
      }
      else {
	std::cout << "auto range on ratio failed\n";
      }
    }


    TString cmpTitle= labels[0] + " vs " + labels[ih];

    for (int iEta=0; iEta<nEta; iEta++) {
      TString tag=Form("_iEta%d",iEta);
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
      if (iEta==1) etaFormat= formatPreGap;
      else if (iEta==2) etaFormat= formatGap;
      else if (iEta==3) etaFormat= formatPostGap;
      TString etaRange=Form(etaFormat,
			    h2->GetXaxis()->GetBinLowEdge(iEta+1),
			    h2->GetXaxis()->GetBinLowEdge(iEta+1)
			    + h2->GetXaxis()->GetBinWidth(iEta+1));

      TH1D *h1main= h1mainV[iEta];
      logAxis(h1main,0,"","#varepsilon");
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
	SaveCanvas(cx,outName,"dir-plot-DYeeEff/");
      }
      else {
	std::cout << "canvas not saved\n";
      }
    }
  }
}

// --------------------------------------------------------------

void plotEffs(TH2D *h2main, TString labelMain,
	      TH2D *h2cmp, TString labelCmp,
	      TString plotNameBase) {}

// --------------------------------------------------------------
