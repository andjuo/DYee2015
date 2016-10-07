#include "inputs.h"
#include "DYbinning.h"
#include "DYmm13TeV_eff.h"

// ------------------------------------------------------

void studyAvgValues(const EventSpace_t &es);

// ------------------------------------------------------

void calcRhoRndVec(int nSamples=100)
{

  //TH2D* h2PostrFsrEventSpace=NULL;
  DYTnPEff_t tnpEff;
  TH1D* h1postFsrInAccSel_MWPU=NULL; // will provide mass binning

  EventSpace_t esPostFsr;

  TString fnameChk="../../../v20160214_1st_CovarianceMatrixInputs/Input6/ROOTFile_Input6_CrossCheck.root";
  TVersion_t inpVersion=_verMu1;
  inpVersion=_verMu76X;
  inpVersion=_verMuApproved;

  if ((inpVersion!=_verMu1) && (DYtools::nMassBins!=DYtools::nMassBins43)) {
    std::cout << "a potential DYbinning.h problem\n";
    return;
  }

  if (inpVersion==_verMu76X) {
    fnameChk="/mnt/sdb/andriusj/v20160527_1st_CovarianceMatrixInputs_76X/Input6/ROOTFile_Input6_CrossCheck.root";
  }
  else if (inpVersion==_verMuApproved) {
    fnameChk="/mnt/sdb/andriusj/v20160915_CovInput_ApprovedResults/ROOTFile_Input6_CrossCheck.root";
  }

  // Load data
  TString fname="dymm_test_RECO_" + versionName(inpVersion) + TString(".root");
  TFile fin(fname);
  if (!fin.IsOpen()) {
    std::cout << "failed to open the file <" << fin.GetName() << ">\n";
    return;
  }
  h1postFsrInAccSel_MWPU=loadHisto(fin, "h1_postFsrInAccSel_MWPU","h1postFsrInAccSel_MW",1, h1dummy);
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
  std::cout << "calculate h1rho_4p2\n";
  TH1D* h1rho_4p2= esPostFsr.calculateScaleFactor(tnpEff,0,"h1rho_4p2","scale factor 4p2");
  std::cout << "calculate h1rho_4p3\n";
  TH1D* h1rho_4p3= esPostFsr.calculateScaleFactor(tnpEff,1,"h1rho_4p3","scale factor 4p3");
  //printHisto(h1rho_4p2); return;

  if (1) {
    h1rho_4p2->GetYaxis()->SetRangeUser(0.9,1.0);
    h1rho_4p3->GetYaxis()->SetRangeUser(0.9,1.0);
  }

  h1rho_4p2->SetStats(0);
  h1rho_4p3->SetStats(0);
  histoStyle(h1rho_4p2,kBlack,5);
  histoStyle(h1rho_4p3,kBlack,5);
  plotHisto(h1rho_4p2,"c4p2",1,0,"LPE1","my code");
  plotHisto(h1rho_4p3,"c4p3",1,0,"LPE1","my code");

  if (1) {
    TFile finChk(fnameChk);
    if (!finChk.IsOpen()) {
      std::cout << "file <" << finChk.GetName() << "> could not be found\n";
      return;
    }
    TGraphAsymmErrors *gr4p2= (TGraphAsymmErrors*)finChk.Get("g_EffSF_HLTv4p2");
    TGraphAsymmErrors *gr4p3= (TGraphAsymmErrors*)finChk.Get("g_EffSF_HLTv4p3");
    finChk.Close();
    if (!gr4p2 || !gr4p3) {
      std::cout << "failed to load graphs\n";
      return;
    }
    TH1D *h1KP4p2= convert(gr4p2,"h1KP4p2","HLTv4p2 sf",0);
    TH1D *h1KP4p3= convert(gr4p3,"h1KP4p3","HLTv4p3 sf",0);
    histoStyle(h1KP4p2,kBlue,24);
    histoStyle(h1KP4p3,kBlue,24);
    plotHistoSame(h1KP4p2,"c4p2","LPE1","KL");
    plotHistoSame(h1KP4p3,"c4p3","LPE1","KL");
    //return;
  }

  if (1 && (nSamples==1)) {
    studyAvgValues(esPostFsr);
  }

  TString outDir="dir-Rho" + versionName(inpVersion) + "/";
  gSystem->mkdir(outDir);
  TString foutName=outDir + "dymm_rhoRndVec_" + versionName(inpVersion) +
    Form("_%d.root",nSamples);
  TFile fout(foutName,"recreate");
  if (!fout.IsOpen()) {
    std::cout << "failed to create the file <" << fout.GetName() << ">\n";
    return;
  }

  // randomize
  DYTnPEff_t rndEff;
  TH1D* h1rho_4p2_avg= cloneHisto(h1rho_4p2,"h1rho_4p2_avg","h1rho_4p2_avg");
  TH1D* h1rho_4p3_avg= cloneHisto(h1rho_4p3,"h1rho_4p3_avg","h1rho_4p3_avg");
  h1rho_4p2_avg->Reset();
  h1rho_4p3_avg->Reset();
  TH1D* h1rho_4p2_sqrAvg= cloneHisto(h1rho_4p2,"h1rho_4p2_sqrAvg","h1rho_4p2_sqrAvg");
  TH1D* h1rho_4p3_sqrAvg= cloneHisto(h1rho_4p3,"h1rho_4p3_sqrAvg","h1rho_4p3_sqrAvg");
  h1rho_4p2_sqrAvg->Reset();
  h1rho_4p3_sqrAvg->Reset();

  const int displayRnd4p3=0;
  for (int i=0; i<nSamples; i++) {
    TString tag=Form("_var%d",i);
    rndEff.randomize(tnpEff,tag);
    if (nSamples<51) rndEff.save(fout);
    TH1D *h1rho_4p2_rnd= esPostFsr.calculateScaleFactor(rndEff,0,"h1rho_4p2"+tag, "scale factor 4p2" + tag);
    TH1D *h1rho_4p3_rnd= esPostFsr.calculateScaleFactor(rndEff,1,"h1rho_4p3"+tag, "scale factor 4p3" + tag);
    h1rho_4p2_rnd->Write();
    h1rho_4p3_rnd->Write();
    h1rho_4p2_avg->Add(h1rho_4p2_rnd,1/double(nSamples));
    h1rho_4p3_avg->Add(h1rho_4p3_rnd,1/double(nSamples));
    if (displayRnd4p3) {
      histoStyle(h1rho_4p3_rnd,kGreen,24,2);
      plotHistoSame(h1rho_4p3_rnd,"c4p3","LPE");
    }

    TH1D* h1rho_4p2_rndSqr= (TH1D*)h1rho_4p2_rnd->Clone("h1rho_4p2_rndSqr");
    h1rho_4p2_rndSqr->Multiply(h1rho_4p2_rnd);
    TH1D* h1rho_4p3_rndSqr= (TH1D*)h1rho_4p3_rnd->Clone("h1rho_4p3_rndSqr");
    h1rho_4p3_rndSqr->Multiply(h1rho_4p3_rnd);
    h1rho_4p2_sqrAvg->Add(h1rho_4p2_rndSqr, 1/double(nSamples));
    h1rho_4p3_sqrAvg->Add(h1rho_4p3_rndSqr, 1/double(nSamples));
    delete h1rho_4p2_rndSqr;
    delete h1rho_4p3_rndSqr;
  }
  h1rho_4p2_avg->SetDirectory(0);
  h1rho_4p3_avg->SetDirectory(0);
  h1rho_4p2_sqrAvg->SetDirectory(0);
  h1rho_4p3_sqrAvg->SetDirectory(0);
  fout.Close();

  for (int ibin=1; ibin<=h1rho_4p2_avg->GetNbinsX(); ibin++) {
    for (int jbin=1; jbin<=h1rho_4p2_avg->GetNbinsY(); jbin++) {
      double err4p2=
	h1rho_4p2_sqrAvg->GetBinContent(ibin,jbin) -
	pow( h1rho_4p2_avg->GetBinContent(ibin,jbin), 2 );
      h1rho_4p2_avg->SetBinError(ibin,jbin, sqrt(err4p2));
      double err4p3=
	h1rho_4p3_sqrAvg->GetBinContent(ibin,jbin) -
	pow( h1rho_4p3_avg->GetBinContent(ibin,jbin), 2 );
      h1rho_4p3_avg->SetBinError(ibin,jbin, sqrt(err4p3));
    }
  }

  histoStyle(h1rho_4p2_avg,kRed,24);
  histoStyle(h1rho_4p3_avg,kRed,24);
  plotHistoSame(h1rho_4p2_avg,"c4p2","LPE1","avg");
  plotHistoSame(h1rho_4p3_avg,"c4p3","LPE1","avg");
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
