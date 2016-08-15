#include "inputs.h"
#include "DYbinning.h"
#include "DYmm13TeV_eff.h"

// ------------------------------------------------------

void studyAvgValues(const EventSpace_t &es);

// ------------------------------------------------------

void calcRhoRndVec_ee(int nSamples=100)
{

  //TH2D* h2PostrFsrEventSpace=NULL;
  DYTnPEff_t tnpEff(1); // 1 - electron channel
  TH1D* h1postFsrInAccSel_MWPU=NULL; // will provide mass binning

  EventSpace_t esPostFsr;

  TString fnameChk="/mnt/sdb/andriusj/v2_08082016_CovarianceMatrixInputs/";
  TVersion_t inpVersion=_verEl2skim;

  if (DYtools::nMassBins!=DYtools::nMassBins43) {
    std::cout << "a potential DYbinning.h problem\n";
    return;
  }
  fnameChk="/mnt/sdb/andriusj/v2_08082016_CovarianceMatrixInputs/ROOTFile_Input6_CrossCheck.root";

  // Load data
  TString fname="dyee_test_RECO_" + versionName(inpVersion) + TString(".root");
  TFile fin(fname);
  if (!fin.IsOpen()) {
    std::cout << "failed to open the file <" << fin.GetName() << ">\n";
    return;
  }
  h1postFsrInAccSel_MWPU=loadHisto(fin, "h1_postFsrInAccSel_MWPU","h1postFsrInAccSel_MWPU",1, h1dummy);
  if (!esPostFsr.load(fin, "mainES")) return;
  if (!tnpEff.load(fin,"DYTnPEff")) {
    std::cout << "failed to load efficiencies\n";
    return;
  }
  fin.Close();

  if (0) {
    if (!tnpEff.checkBinning()) std::cout << " binning check failed\n";
    else std::cout << " binning check ok\n";
    return;
  }
  if (0) {
    printHisto(tnpEff.h2Eff_RecoID_Data);
    tnpEff.printNumbers();
    return;
  }

  // check the event space
  if (0) {
    TH1D* h1postFsrCheck= cloneHisto(h1postFsrInAccSel_MWPU,"h1postFsrCheck","h1postFsrCheck");
    h1postFsrCheck->Reset();

    for (int im=0; im<DYtools::nMassBins; im++) {
      const TH2D* h2= esPostFsr.h2ES(im);
      double sum=0, sumErr=0;
      for (int ibin=1; ibin<=h2->GetNbinsX(); ibin++) {
	for (int jbin=1; jbin<=h2->GetNbinsY(); jbin++) {
	  sum+= h2->GetBinContent(ibin,jbin);
	  sumErr+= pow(h2->GetBinError(ibin,jbin),2);
	}
      }
      sumErr=sqrt(sumErr);
      h1postFsrCheck->SetBinContent(im+1, sum);
      h1postFsrCheck->SetBinError  (im+1, sumErr);
    }

    plotHisto(h1postFsrInAccSel_MWPU,"cSel",1,1,"LPE");
    histoStyle(h1postFsrCheck,kBlue,24);
    plotHistoSame(h1postFsrCheck,"cSel","LPE");
    printRatio(h1postFsrInAccSel_MWPU,h1postFsrCheck);
    return;
  }

  // calculate event scale factors
  std::cout << "calculate h1rho\n";
  TH1D* h1rho= esPostFsr.calculateScaleFactor(tnpEff,0,"h1rho","scale factor;M_{ee} [GeV];#rho");

  if (1) {
    h1rho->GetYaxis()->SetRangeUser(0.9,1.0);
  }

  h1rho->SetStats(0);
  histoStyle(h1rho,kBlack,5);
  plotHisto(h1rho,"cRho",1,0,"LPE1","my code");

  if (1) {
    TH1D *h1rhoChk= loadHisto(fnameChk,"h_EfficiencySF", "h1rhoChk",1,h1dummy);
    if (!h1rhoChk) return;
    histoStyle(h1rhoChk,kBlue,24);
    plotHistoSame(h1rhoChk,"cRho","LPE1","RC");
    if (0) {
      EventSpace_t esPostFsr_noPU;
      TFile fin2(fname);
      if (!esPostFsr_noPU.load(fin2,"mainES_noPU")) return;
      fin2.Close();
      TH1D* h1rho_noPU= esPostFsr_noPU.calculateScaleFactor(tnpEff,0,"h1rho_noPU","scale factor (no PU);M_{ee} [GeV];#rho (noPU)");
      histoStyle(h1rho_noPU,kRed,25);
      plotHistoSame(h1rho_noPU,"cRho","LPE1","my code (noPU)");
      //printRatio(h1rho,h1rho_noPU);
    }
    //return;
  }

  if (1 && (nSamples==1)) {
    studyAvgValues(esPostFsr);
  }

  TString outDir="dir-Rho" + versionName(inpVersion) + "/";
  gSystem->mkdir(outDir);
  TString foutName=outDir + "dyee_rhoRndVec_" + versionName(inpVersion) +
    Form("_%d.root",nSamples);
  TFile fout(foutName,"recreate");
  if (!fout.IsOpen()) {
    std::cout << "failed to create the file <" << fout.GetName() << ">\n";
    return;
  }

  // randomize
  DYTnPEff_t rndEff;
  TH1D* h1rho_avg= cloneHisto(h1rho,"h1rho_avg","h1rho_avg");
  h1rho_avg->Reset();
  TH1D* h1rho_sqrAvg= cloneHisto(h1rho,"h1rho_sqrAvg","h1rho_sqrAvg");
  h1rho_sqrAvg->Reset();

  int displayRnd=0;
  for (int i=0; i<nSamples; i++) {
    TString tag=Form("_var%d",i);
    rndEff.randomize(tnpEff,tag);
    if (nSamples<51) rndEff.save(fout);
    TH1D *h1rho_rnd= esPostFsr.calculateScaleFactor(rndEff,0,"h1rho"+tag, "scale factor" + tag);
    h1rho_rnd->Write();
    h1rho_avg->Add(h1rho_rnd,1/double(nSamples));
    if (displayRnd) {
      histoStyle(h1rho_rnd,kGreen,1,2);
      plotHistoSame(h1rho_rnd,"cRho","LPE");
    }

    TH1D* h1rho_rndSqr= (TH1D*)h1rho_rnd->Clone("h1rho_rndSqr");
    h1rho_rndSqr->Multiply(h1rho_rnd);
    h1rho_sqrAvg->Add(h1rho_rndSqr, 1/double(nSamples));
    delete h1rho_rndSqr;
  }
  h1rho_avg->SetDirectory(0);
  h1rho_sqrAvg->SetDirectory(0);
  fout.Close();

  for (int ibin=1; ibin<=h1rho_avg->GetNbinsX(); ibin++) {
    for (int jbin=1; jbin<=h1rho_avg->GetNbinsY(); jbin++) {
      double err=
	h1rho_sqrAvg->GetBinContent(ibin,jbin) -
	pow( h1rho_avg->GetBinContent(ibin,jbin), 2 );
      h1rho_avg->SetBinError(ibin,jbin, sqrt(err));
    }
  }

  histoStyle(h1rho_avg,kRed,24);
  plotHistoSame(h1rho_avg,"cRho","LPE1","avg");

  return;
}

// ------------------------------------------------------

void studyAvgValues(const EventSpace_t &es) {

  std::cout << "\nStudyAvgValues\n";

  plotHisto(es.h2EffBinDef(),"cEffBinDef",0,0);

  std::vector<TH2D*> h2ptSpace, h2etaSpace;
  std::vector<TH1D*> avgH1V=es.avgAxisValues(&h2ptSpace,&h2etaSpace);

  for (unsigned int im=0; im<h2ptSpace.size(); im++) {
    if (h2ptSpace[im]->Integral()!=0) {
      plotHisto(h2ptSpace[im],"cPtSpace"+DYtools::massStr(im),0,0);
    }
  }
  for (unsigned int im=0; im<h2etaSpace.size(); im++) {
    if (h2etaSpace[im]->Integral()!=0) {
      plotHisto(h2etaSpace[im],"cEtaSpace"+DYtools::massStr(im),0,0);
    }
  }

  histoStyle(avgH1V[0],kGreen+1,5);
  histoStyle(avgH1V[1],kGreen+1,5);
  histoStyle(avgH1V[2],kBlue,24);
  histoStyle(avgH1V[3],kBlue,24);
  plotHisto(avgH1V[1],"cAvgPt",1,0,"LPE","avgPt1");
  plotHistoSame(avgH1V[3],"cAvgPt","LPE","avgPt2");
  plotHisto(avgH1V[0],"cAvgEta",1,0,"LPE","avg|Eta1|");
  plotHistoSame(avgH1V[2],"cAvgEta","LPE","avg|Eta2|");

  if (1) {
    for (unsigned int ih=0; ih<avgH1V.size(); ih++) {
      printHisto(avgH1V[ih]);
    }
  }
}

// ------------------------------------------------------
