#include "inputs.h"
#include <fstream>

// ----------------------------------------------------------------

const int nPt=5;
const double Ele_pt[nPt+1] = { 10.,20.,30.,40.,50.,200. };
const int nEta=5;
const double SC_abseta[nEta+1] = { 0,0.8,1.4442,1.566,2.0,2.5 };

void plotEffs(std::vector<TH2D*> &hV,
	      const std::vector<TString> &labels,
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

  const char *fname="CutBasedID_MediumWP_76X_18Feb.txt";
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
  
  HistoStyle_t(kRed,  5,1,1.0,2.).SetStyle(h2effData);
  HistoStyle_t(kBlue,24,1,0.8,2.).SetStyle(h2effData_bkgPdf);
  HistoStyle_t(kGreen+1,20,1,0.8,2.).SetStyle(h2effData_sigPdf);
  HistoStyle_t(kBlack,3,1,1.0,2.).SetStyle(h2effMC);
  HistoStyle_t(kBlue,24,1,0.8,2.).SetStyle(h2effMC_NLOvsLO);
  HistoStyle_t(6 ,   26,1,0.8,2.).SetStyle(h2effData_tag);
  HistoStyle_t(46,   23,1,0.8,2.).SetStyle(h2effData_Zrange);

  if (0) {
    std::vector<TH2D*> h2_data_vs_mc;
    std::vector<TString> label_data_vs_mc;
    h2_data_vs_mc.push_back(h2effMC); label_data_vs_mc.push_back("MC");
    h2_data_vs_mc.push_back(h2effData); label_data_vs_mc.push_back("data");
    plotEffs(h2_data_vs_mc,label_data_vs_mc,"Data_vs_MC");
  }

  if (1) {
    std::vector<TH2D*> h2_data_var;
    std::vector<TString> label_data_var;
    h2_data_var.push_back(h2effData); label_data_var.push_back("data");
    h2_data_var.push_back(h2effData_bkgPdf); label_data_var.push_back("data: bkgPdf");
    h2_data_var.push_back(h2effData_sigPdf); label_data_var.push_back("data: sigPdf");
    h2_data_var.push_back(h2effData_tag); label_data_var.push_back("data: tag def");
    h2_data_var.push_back(h2effData_Zrange); label_data_var.push_back("data: Zrange");
    plotEffs(h2_data_var,label_data_var, "Data");
  }

  if (0) {
    std::vector<TH2D*> h2_mc_var;
    std::vector<TString> label_mc_var;
    h2_mc_var.push_back(h2effMC); label_mc_var.push_back("MC");
    h2_mc_var.push_back(h2effMC_NLOvsLO); label_mc_var.push_back("MC: NLO vs LO");
    plotEffs(h2_mc_var,label_mc_var, "MC");
  }
}

// --------------------------------------------------------------

void plotEffs(std::vector<TH2D*> &h2V,
	      const std::vector<TString> &labels,
	      TString plotNameBase)
{

  for (int iEta=0; iEta<nEta; iEta++) {
    TString tag=Form("_iEta%d",iEta);
    std::vector<TH1D*> h1V;
    for (unsigned int ih=0; ih<h2V.size(); ih++) {
      TH2D *h2= h2V[ih];
      TH1D *h1eff_iEta=	h2->ProjectionY(h2->GetName() + tag,iEta+1,iEta+2);
      h1eff_iEta->SetStats(0);
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
// --------------------------------------------------------------
