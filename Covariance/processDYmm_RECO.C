#include "DYmm13TeV.h"
#include "DYbinning.h"
#include "inputs.h"
#include "../RooUnfold/src/RooUnfoldResponse.h"
#include "../RooUnfold/src/RooUnfoldBayes.h"

// --------------------------------------------------------------

inline
void prepareHisto(TH1D *h1) {
  h1->SetDirectory(0);
  //if (!h1->GetSumw2()) 
  h1->Sumw2();
}

// --------------------------------------------------------------

TH2D* loadEffHisto(TString hName) {
  TString srcPath="/media/ssd/v20160214_1st_CovarianceMatrixInputs/";
  srcPath="/mnt/sdb/andriusj/v20160214_1st_CovarianceMatrixInputs/";

  TString fname= srcPath + "Input5/ROOTFile_TagProbeEfficiency.root";
  TFile fin(fname);
  if (!fin.IsOpen()) {
    std::cout << "failed to open the file <" << fname << ">\n";
    return NULL;
  }
  TH2D *h2= loadHisto(fin, "h_2D_Eff_Iso_MC", hName,1,h2dummy);
  fin.Close();
  return h2;
}

// --------------------------------------------------------------

void processDYmm_RECO(Int_t maxEntries=-1)
{
  std::cout << "processDYmm_RECO\n";
  std::cout << "DY range=" << DYtools::minMass << " .. " << DYtools::maxMass << "\n";

  TString srcPath="/media/ssd/v20160214_1st_CovarianceMatrixInputs/";
  srcPath="/mnt/sdb/andriusj/v20160214_1st_CovarianceMatrixInputs/";

  TH2D* h2EffBinDef= loadEffHisto("h2EffBinDef");
  if (!h2EffBinDef) return;

  DYmm13TeV_t data(srcPath + "Input3/ROOTFile_ntuple_CovarianceMatrixInput.root");
  data.DeactivateBranches();
  data.ActivateBranches("Momentum_Reco_Lead Momentum_Reco_Sub Momentum_postFSR_Lead  Momentum_postFSR_Sub  Weight_Norm  Weight_PU  Weight_Gen");
  data.ActivateBranches("Flag_EventSelection");

  TH1D *h1recoSel_M= new TH1D("h1recoSel_M", "RECO selected;M_{reco} [GeV];count",DYtools::nMassBins,DYtools::massBinEdges);
  TH1D *h1recoSel_MW= new TH1D("h1recoSel_MW", "RECO selected (weighted);M_{reco} [GeV];weighted count",DYtools::nMassBins,DYtools::massBinEdges);

  TH1D *h1postFsrInAccSel_M= new TH1D("h1_postFsrInAccSel_M", "postFSR selected events in acceptance;M_{gen,postFSR} [GeV];count",DYtools::nMassBins,DYtools::massBinEdges);
  TH1D *h1postFsrInAccSel_MW= new TH1D("h1_postFsrInAccSel_MW", "postFSR selected events in acceptance (weighted);M_{gen,postFSR} [GeV];weighted count",DYtools::nMassBins,DYtools::massBinEdges);

  TH2D* h2PostFsrEventSpace= new TH2D("h2PostFsrEventSpace","post-FSR event space;(eta,pt) flat index 1;(eta,pt) flat index 2", DYtools::EtaPtFIMax+1, -0.5, DYtools::EtaPtFIMax+0.5, DYtools::EtaPtFIMax+1, -0.5, DYtools::EtaPtFIMax+0.5);
  h2PostFsrEventSpace->SetDirectory(0);
  h2PostFsrEventSpace->Sumw2();

  prepareHisto(h1recoSel_M);
  prepareHisto(h1recoSel_MW);
  prepareHisto(h1postFsrInAccSel_M);
  prepareHisto(h1postFsrInAccSel_MW);

  TH2D* h2DetResMig= new TH2D("h2DetResMig","Det.resolution migration", DYtools::nMassBins, DYtools::massBinEdges, DYtools::nMassBins, DYtools::massBinEdges);
  h2DetResMig->Sumw2();
  h2DetResMig->SetDirectory(0);
  RooUnfoldResponse detResResp(h1recoSel_MW,h1postFsrInAccSel_MW,h2DetResMig,"rooUnf_detResResp","detResResp;M_{#mu#mu} [GeV];det.res. unfolded yield");
  detResResp.UseOverflow();
  RooUnfoldResponse detResRespPU(h1recoSel_MW,h1postFsrInAccSel_MW,h2DetResMig,"rooUnf_detResRespPU","detResRespPU;M_{#mu#mu} [GeV];det.res.(PU) unfolded yield");
  detResRespPU.UseOverflow();

  UInt_t nEvents= data.GetEntries();
  UInt_t nProcessed=0, nSelected=0;

  std::cout << "process data\n";
  for (UInt_t iEntry=0; iEntry<nEvents; iEntry++) {
    if (data.GetEntry(iEntry)<0) break;
    if (iEntry%100000==0) std::cout << "iEntry=" << iEntry << Form("(%4.2lf%%)\n",iEntry*100/double(nEvents));
    if ((maxEntries>0) && (iEntry>UInt_t(maxEntries))) {
      std::cout << "debug run\n";
      break;
    }
    nProcessed++;

    if (!data.Flag_EventSelection) continue;
    nSelected++;

    double w= data.Weight_Norm * data.Weight_Gen;
    double wPU= w* data.Weight_PU;
    if (iEntry<1000) {
      std::cout << "w=" << w << ": weight_norm=" << data.Weight_Norm << ", gen=" << data.Weight_Gen << "\n";
      std::cout << "   wPU=" << wPU << ": weight_PU=" << data.Weight_PU << "\n";
    }
    double mReco= (*data.Momentum_Reco_Lead + *data.Momentum_Reco_Sub).M();
    double mPostFsr= (*data.Momentum_postFSR_Lead + *data.Momentum_postFSR_Sub).M();

    if (DYtools::InsideMassRange(mReco)) {
      h1recoSel_M->Fill( mReco, 1. );
      h1recoSel_MW->Fill( mReco, wPU );
    }

    if ( ! DYtools::InsideMassRange(mPostFsr) &&
	 data.Flag_EventSelection ) {
      if (mReco> DYtools::minMass) {
	std::cout << "fake " << mPostFsr << " (mReco=" << mReco << ")\n";
	detResResp.Fake(mReco, wPU);
	detResRespPU.Fake(mReco, wPU);
      }
    }
    else if ( DYtools::InsideMassRange(mPostFsr) &&
	      ! data.Flag_EventSelection ) {
      std::cout << "miss\n";
      detResResp.Miss(mPostFsr, w);
      detResRespPU.Miss(mPostFsr, wPU);
    }
    else {
      detResResp.Fill( mReco, mPostFsr, w );
      detResRespPU.Fill( mReco, mPostFsr, wPU );
      h2DetResMig->Fill( mReco, mPostFsr, w );
    }

    if (DYtools::InAcceptance_mm(data.Momentum_postFSR_Lead,data.Momentum_postFSR_Sub)) {
      h1postFsrInAccSel_M->Fill( mPostFsr, 1. );
      h1postFsrInAccSel_MW->Fill( mPostFsr, w );

      double fi1= DYtools::FlatIndex(h2EffBinDef, data.Momentum_postFSR_Lead);
      double fi2= DYtools::FlatIndex(h2EffBinDef, data.Momentum_postFSR_Sub);
      if ((fi1>=0) && (fi2>=0)) {
	h2PostFsrEventSpace->Fill(fi1,fi2, w);
      }
      if ((fi1>DYtools::EtaPtFIMax) || (fi2>DYtools::EtaPtFIMax)) {
	std::cout << "EtaPtFIMax=" << DYtools::EtaPtFIMax
		  << ", fi1=" << fi1 << ", fi2=" << fi2 << "\n";
	return;
      }
    }
  }
  std::cout << "\nProcessed " << nProcessed << " events, selected " << nSelected << "\n";

  RooUnfoldBayes bayesDetRes( &detResResp, h1recoSel_MW, 4 );
  TH1D *h1Unf= (TH1D*) bayesDetRes.Hreco()->Clone("h1Unf");
  h1Unf->SetTitle("unfolded reco->postFSR");
  h1Unf->SetDirectory(0);
  h1Unf->SetLineColor(kRed);
  h1Unf->SetMarkerStyle(24);
  h1Unf->SetMarkerColor(kRed);
  h1Unf->GetXaxis()->SetMoreLogLabels();
  h1Unf->GetXaxis()->SetNoExponent();
  TCanvas *cDetResTest=plotHisto(h1Unf,"cDetResTest",1,1,"LPE1");
  plotHistoSame(h1postFsrInAccSel_MW,"cDetResTest","LPE");

  RooUnfoldBayes bayesDetResPU( &detResRespPU, h1recoSel_MW, 4 );
  TH1D *h1UnfPU= (TH1D*) bayesDetResPU.Hreco()->Clone("h1UnfPU");
  h1UnfPU->SetTitle("unfolded (PU) reco->postFSR");
  h1UnfPU->SetDirectory(0);
  h1UnfPU->SetLineColor(kBlue);
  h1UnfPU->SetMarkerStyle(5);
  h1UnfPU->SetMarkerColor(kBlue);
  TCanvas *cDetResTestPU=plotHisto(h1UnfPU,"cDetResTestPU",1,1,"LPE1");
  plotHistoSame(h1postFsrInAccSel_MW,"cDetResTestPU","LPE");

  TH1D* h1UnfBinByBinCorr=(TH1D*)h1recoSel_MW->Clone("h1UnfBinByBinCorr");
  h1UnfBinByBinCorr->SetDirectory(0);
  if (!h1UnfBinByBinCorr->GetSumw2()) h1UnfBinByBinCorr->Sumw2();
  h1UnfBinByBinCorr->SetTitle("Unf corr bin-by-bin;M_{#mu#mu} [GeV];postFSR_inAcc_Sel/recoSel");
  h1UnfBinByBinCorr->Divide(h1postFsrInAccSel_MW,h1recoSel_MW,1,1,"B");

  TString fname="dymm_test_RECO.root";
  if (maxEntries>0) fname.ReplaceAll(".root","_debug.root");
  TFile fout(fname,"RECREATE");
  h1UnfBinByBinCorr->Write();
  detResResp.Write();
  detResRespPU.Write();
  h1recoSel_M->Write();
  h1recoSel_MW->Write();
  h1postFsrInAccSel_M->Write();
  h1postFsrInAccSel_MW->Write();
  h2PostFsrEventSpace->Write();
  h1Unf->Write();
  h1UnfPU->Write();
  h2DetResMig->Write();
  cDetResTest->Write();
  cDetResTestPU->Write();
  TObjString timeTag(DayAndTimeTag(0));
  timeTag.Write("timeTag");

  fout.Close();
  std::cout << "file <" << fout.GetName() << "> created\n";
}

