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
// --------------------------------------------------------------


void processDYee(Int_t maxEntries=-1, TString includeOnlyRange="")
{
  std::cout << "DY range=" << DYtools::minMass << " .. " << DYtools::maxMass << "\n";

  TVersion_t inpVersion=_verEl2skim;
  const int testSC=1;

  if (DYtools::nMassBins!=DYtools::nMassBins43) {
    std::cout << "a potential DYbinning.h problem\n";
    return;
  }

  TString srcPath="/mnt/sdb/andriusj/v2ee-skim1-20160812/";
  inpVersion=_verEl2skim2; // preFSR is dressedElectrons(0.1)
  srcPath="/media/sf_Share2/DYee_76X_Calibrated/mySkim-20160906/";

  TString dataFName;
  TString fileFormat="DYEE_M%dto%d_v1.root ";
  if (inpVersion==_verEl2skim2) fileFormat="DY_%dto%d_v2_orig.root ";

  if (includeOnlyRange.Length()>0) {
    std::stringstream ss(includeOnlyRange.Data());
    TString rStr;
    TString strFileFormat=fileFormat;
    strFileFormat.ReplaceAll("%dto%d","%s");
    while (!ss.eof()) {
      ss >> rStr;
      if (rStr.Length()) {
	std::cout << "rStr=" << rStr << "\n";
	dataFName.Append(srcPath + Form(strFileFormat,rStr.Data()));
      }
    }
  }
  else {
    for (int im=0; im<DYtools::nEEDatasetMassBins; im++) {
      dataFName.Append(srcPath + Form(fileFormat,
		      DYtools::eeDatasetMass[im],DYtools::eeDatasetMass[im+1]));
    }
  }
  //std::cout << "dataFName=<" << dataFName << ">\n";

  DYmm13TeV_t data(dataFName);
  //data.DeactivateBranches(); // does not work
  data.ActivateBranches("Momentum_postFSR_Lead  Momentum_postFSR_Sub  Momentum_preFSR_Lead  Momentum_preFSR_Sub  Weight_Norm  Weight_Gen Weight_PU");
  data.ActivateBranches("Flag_EventSelection");

  TH1D *h1postFsrInAccSel_M= new TH1D("h1_postFsrInAccSel_M", "postFSR selected events in acceptance;M_{gen,postFSR} [GeV];count",DYtools::nMassBins,DYtools::massBinEdges);
  TH1D *h1postFsrInAccSel_MW= new TH1D("h1_postFsrInAccSel_MW", "postFSR selected events in acceptance (weighted);M_{gen,postFSR} [GeV];weighted count",DYtools::nMassBins,DYtools::massBinEdges);
  TH1D *h1postFsrInAccSel_MWPU= new TH1D("h1_postFsrInAccSel_MWPU", "postFSR selected events in acceptance (weighted wPU);M_{gen,postFSR} [GeV];weighted (wPU) count",DYtools::nMassBins,DYtools::massBinEdges);
  TH1D *h1postFsrInAcc_M= new TH1D("h1_postFsrInAcc_M", "postFSR in acceptace;M_{gen,postFSR} [GeV];count",DYtools::nMassBins,
			    DYtools::massBinEdges);
  TH1D *h1postFsrInAcc_MW= new TH1D("h1_postFsrInAcc_Mweighted", "postFSR in acceptance;M_{gen,postFSR} [GeV];weighted count", DYtools::nMassBins,
			    DYtools::massBinEdges);
  TH1D *h1postFsrInAcc_MWPU= new TH1D("h1_postFsrInAcc_MweightedPU", "postFSR in acceptance (weighted wPU);M_{gen,postFSR} [GeV];weighted (wPU) count", DYtools::nMassBins,
			    DYtools::massBinEdges);
  TH1D *h1postFsr_M= new TH1D("h1_postFsr_M", "postFSR;M_{gen,postFSR} [GeV];count",DYtools::nMassBins,
			    DYtools::massBinEdges);
  TH1D *h1postFsr_MW= new TH1D("h1_postFsr_Mweighted", "postFSR;M_{gen,postFSR} [GeV];weighted count", DYtools::nMassBins,
			    DYtools::massBinEdges);
  TH1D *h1preFsr_M= new TH1D("h1_preFsr_M", "preFSR;M_{gen,preFSR} [GeV];count",DYtools::nMassBins,
			   DYtools::massBinEdges);
  TH1D *h1preFsr_MW= new TH1D("h1_preFsr_Mweighted", "preFSR;M_{gen,preFSR} [GeV];weighted count", DYtools::nMassBins,
			   DYtools::massBinEdges);
  TH1D *h1preFsrTest_MW= new TH1D("h1_preFsrTest_Mweighted", "preFSR test version (different from preFSR);M_{gen,preFSR} [GeV];weighted count", DYtools::nMassBins,
			   DYtools::massBinEdges);
  prepareHisto(h1postFsrInAccSel_M);
  prepareHisto(h1postFsrInAccSel_MW);
  prepareHisto(h1postFsrInAccSel_MWPU);
  prepareHisto(h1postFsrInAcc_M);
  prepareHisto(h1postFsrInAcc_MW);
  prepareHisto(h1postFsrInAcc_MWPU);
  prepareHisto(h1postFsr_M);
  prepareHisto(h1postFsr_MW);
  prepareHisto(h1preFsr_M);
  prepareHisto(h1preFsr_MW);
  prepareHisto(h1preFsrTest_MW);

  TH2D* h2FSRmig= new TH2D("h2FSRmig","FSR migration", DYtools::nMassBins, DYtools::massBinEdges, DYtools::nMassBins, DYtools::massBinEdges);
  h2FSRmig->Sumw2();
  RooUnfoldResponse fsrResp(h1postFsr_MW,h1preFsr_MW,h2FSRmig,"rooUnf_fsrResp","fsrResp;M_{ee} [GeV];FSR unfolded yield");
  fsrResp.UseOverflow(false);

  TH2D* h2effAccFSRmig= new TH2D("h2effAccFSRmig","eff x Acc x FSR migration", DYtools::nMassBins, DYtools::massBinEdges, DYtools::nMassBins, DYtools::massBinEdges);
  h2effAccFSRmig->Sumw2();
  RooUnfoldResponse effAccFsrResp(h1postFsrInAccSel_MW,h1preFsr_MW,h2effAccFSRmig,"rooUnf_effAccFsrResp","effAccFsrResp;M_{ee} [GeV];eff Acc FSR unfolded yield");
  effAccFsrResp.UseOverflow(false);


  UInt_t nEvents= data.GetEntries();

  std::cout << "process data\n";
  for (UInt_t iEntry=0; iEntry<nEvents; iEntry++) {
    if (data.GetEntry(iEntry)<0) break;
    if (iEntry%100000==0) std::cout << "iEntry=" << iEntry << Form("(%4.2lf%%)\n",iEntry*100./double(nEvents));
    if ((maxEntries>0) && (iEntry>UInt_t(maxEntries))) {
      std::cout << "debug run\n";
      break;
    }

    double w= data.Weight_Norm * data.Weight_Gen;
    if (iEntry<100) std::cout << "w=" << w << ": weight_norm=" << data.Weight_Norm << ", gen=" << data.Weight_Gen << "\n";
    double mPostFsr= (*data.Momentum_postFSR_Lead + *data.Momentum_postFSR_Sub).M();
    double mPreFsr= (*data.Momentum_preFSR_Lead + *data.Momentum_preFSR_Sub).M();

    //std::cout << "momentum postFSR: " << data.Momentum_postFSR_Lead << ", " << data.Momentum_postFSR_Sub << "\n";
    //std::cout << "mPostFsr=" << mPostFsr << ", mPreFsr=" << mPreFsr << "\n";

    if (DYtools::InAcceptance_ee(data.Momentum_postFSR_Lead,data.Momentum_postFSR_Sub,testSC)) {
      if (data.Flag_EventSelection) {
	h1postFsrInAccSel_M->Fill( mPostFsr, 1. );
	h1postFsrInAccSel_MW->Fill( mPostFsr, w );
	h1postFsrInAccSel_MWPU->Fill( mPostFsr, w * data.Weight_PU );
      }
      h1postFsrInAcc_M->Fill( mPostFsr, 1. );
      h1postFsrInAcc_MW->Fill( mPostFsr, w );
      h1postFsrInAcc_MWPU->Fill( mPostFsr, w * data.Weight_PU );
    }
    h1postFsr_M->Fill( mPostFsr, 1. );
    h1postFsr_MW->Fill( mPostFsr, w );
    h1preFsr_M->Fill( mPreFsr, 1. );
    h1preFsr_MW->Fill( mPreFsr, w );
    if (data.Flag_EventSelection) h1preFsrTest_MW->Fill( mPreFsr, w );
    h2FSRmig->Fill( mPostFsr, mPreFsr, w );

    if (0) {
      if ( ! DYtools::InsideMassRange(mPreFsr) )
	std::cout << "preFSR mass " << mPreFsr << ", postFSR mass " << mPostFsr << "\n";
    }

    if ( ! DYtools::InsideMassRange(mPreFsr) &&
	 DYtools::InsideMassRange(mPostFsr) ) {
      std::cout << "fake detected\n";
      fsrResp.Fake(mPostFsr,w);
    }
    else if ( DYtools::InsideMassRange(mPreFsr) &&
	      ! DYtools::InsideMassRange(mPostFsr) ) {
      fsrResp.Miss(mPreFsr,w);
    }
    else {
      fsrResp.Fill(mPostFsr,mPreFsr,w);
    }

    int sel=(DYtools::InAcceptance_ee(data.Momentum_postFSR_Lead,data.Momentum_postFSR_Sub,testSC) && data.Flag_EventSelection && DYtools::InsideMassRange(mPostFsr)) ? 1:0;
    if (sel) h2effAccFSRmig->Fill(mPostFsr,mPreFsr,w);
    if ( ! DYtools::InsideMassRange(mPreFsr) && sel ) {
      effAccFsrResp.Fake(mPostFsr,w);
    }
    else if ( DYtools::InsideMassRange(mPreFsr) && ! sel ) {
      effAccFsrResp.Miss(mPreFsr,w);
    }
    else {
      effAccFsrResp.Fill(mPostFsr,mPreFsr,w);
    }
  }

  // plots

  TCanvas *cPreFSR_WG=NULL, *cPreFSR_WGN=NULL;
  if (data.fH1PreFSR_WG.size()) {
    for (int iloop=1; iloop<2; iloop++) {
      std::vector<TH1D*> *h1V=
	(iloop==0) ? &data.fH1PreFSR_WG : &data.fH1PreFSR_WGN;
      TString cNameLoop= (iloop==0) ? "cPreFSR_WG" : "cPreFSR_WGN";
      for (unsigned int i=0; i<h1V->size(); i+=2) {
	h1V->at(i)->SetLineColor(kRed);
	h1V->at(i)->SetMarkerColor(kRed);
      }
      TCanvas **c= (iloop==0) ? &cPreFSR_WG : &cPreFSR_WGN;
      (*c)=plotHisto(h1V->at(0),cNameLoop,1,1,"hist", Form("M%dto%d",
		     DYtools::eeDatasetMass[0],DYtools::eeDatasetMass[1]));
      for (unsigned int i=1; i<h1V->size(); i++) {
	plotHistoSame(h1V->at(i),cNameLoop,"hist", Form("M%dto%d",
		      DYtools::eeDatasetMass[i],DYtools::eeDatasetMass[i+1]));
      }
    }
  }

  if (0) {
    histoStyle(h1postFsr_MW,kGreen+1,24);
    histoStyle(h1preFsrTest_MW,kRed,5);
    std::cout << "preFsrTest_MW->Integral()=" << h1preFsrTest_MW->Integral() << "\n";
    plotHisto(h1postFsr_MW,"cPreFsrTEST",1,1,"LPE","postFsr");
    plotHistoSame(h1preFsr_MW,"cPreFsrTEST","LPE","preFsr");
    plotHistoSame(h1preFsrTest_MW,"cPreFsrTEST","LPE","preFsrTEST");
  }

  TCanvas *cPreFSRcheck=NULL;
  if (data.fH1PreFSR_binned) {
    histoStyle(data.fH1PreFSR_binned,kOrange,5);
    cPreFSRcheck=plotHisto(data.fH1PreFSR_binned,"cPreFsrEmbeded",1,1,"LPE","preFSR embeded");
    plotHistoSame(h1preFsr_MW,"cPreFsrEmbeded","LPE","preFsr_MW");
    printRatio(data.fH1PreFSR_binned,h1preFsr_MW);
  }

  RooUnfoldBayes bayesFSR( &fsrResp, h1postFsr_MW, 4 );
  TH1D *h1preFsrUnf= (TH1D*) bayesFSR.Hreco()->Clone("h1preFsrUnf");
  h1preFsrUnf->SetTitle("unfolded postFSR->preFSR");
  h1preFsrUnf->SetDirectory(0);
  histoStyle(h1preFsrUnf,kRed,24);
  h1preFsrUnf->GetXaxis()->SetMoreLogLabels();
  h1preFsrUnf->GetXaxis()->SetNoExponent();
  TCanvas *cFSRTest=plotHisto(h1preFsrUnf,"cFSRTest",1,1,"LPE1","BayesUnf4");
  plotHistoSame(h1preFsr_MW,"cFSRTest","LPE","preFsr_MW");

  RooUnfoldBayes bayesEffAccFsr( &effAccFsrResp, h1postFsrInAccSel_MW, 4 );
  TH1D *h1EffAccFsrUnf= (TH1D*) bayesEffAccFsr.Hreco()->Clone("h1EffAccFsrUnf");
  h1EffAccFsrUnf->SetTitle("unfolded postFsrInAccSel -> preFSR");
  h1EffAccFsrUnf->SetDirectory(0);
  histoStyle(h1EffAccFsrUnf,kRed,24);
  h1EffAccFsrUnf->GetXaxis()->SetMoreLogLabels();
  h1EffAccFsrUnf->GetXaxis()->SetNoExponent();
  TCanvas *cScaleToPreFsrTest= plotHisto(h1EffAccFsrUnf,"cScaleToPreFsrTest",1,1,"LPE1","BayesUnf4");
  histoStyle(h1preFsr_MW,kBlue,5);
  plotHistoSame(h1preFsr_MW,"cScaleToPreFsrTest","LPE","preFsr_MW");

  TH1D* h1Eff=(TH1D*)h1postFsrInAccSel_MW->Clone("h1Eff");
  h1Eff->SetDirectory(0);
  if (!h1Eff->GetSumw2()) h1Eff->Sumw2();
  h1Eff->SetTitle("Efficiency;M_{ee,GEN postFSR} [GeV];postFSR_inAcc_Sel/postFSR_inAcc");
  h1Eff->Divide(h1postFsrInAccSel_MW,h1postFsrInAcc_MW,1,1,"B");

  TH1D* h1EffPU=(TH1D*)h1postFsrInAccSel_MWPU->Clone("h1EffPU");
  h1EffPU->SetDirectory(0);
  if (!h1EffPU->GetSumw2()) h1EffPU->Sumw2();
  h1EffPU->SetTitle("Efficiency (wPU/noPU);M_{ee,GEN postFSR} [GeV];postFSR_inAcc_Sel(wPU)/postFSR_inAcc(noPU)");
  h1EffPU->Divide(h1postFsrInAccSel_MWPU,h1postFsrInAcc_MW,1,1,"B");

  histoStyle(h1EffPU,kBlue,24);
  TCanvas *cEffCmp=plotHisto(h1Eff,"cEff_noPU_vs_wPU",1,0,"LPE1","Eff (no.PUw)");
  plotHistoSame(h1EffPU,"cEff_noPU_vs_wPU","LPE","Eff (w.PUw)");

  TH1D* h1Acc=(TH1D*)h1postFsrInAcc_MW->Clone("h1Acc");
  h1Acc->SetDirectory(0);
  if (!h1Acc->GetSumw2()) h1Acc->Sumw2();
  h1Acc->SetTitle("Acceptance;M_{ee,GEN postFSR} [GeV];postFSR_inAcc/postFSR");
  h1Acc->Divide(h1postFsrInAcc_MW,h1postFsr_MW,1,1,"B");

  TH1D* h1EffPUAcc=(TH1D*)h1postFsrInAccSel_MWPU->Clone("h1EffPUAcc");
  h1EffPUAcc->SetDirectory(0);
  if (!h1EffPUAcc->GetSumw2()) h1EffPUAcc->Sumw2();
  h1EffPUAcc->SetTitle("Efficiency(wPU) x Acc;M_{ee,GEN postFSR} [GeV];postFSR_inAccSel(wPU)/h1postFsr_MW(noPU)");
  h1EffPUAcc->Divide(h1postFsrInAccSel_MWPU,h1postFsr_MW,1,1,"B");

  TH1D *h1FSRCorr_binByBin=(TH1D*)h1preFsr_MW->Clone("h1FSRCorr_binByBin");
  h1FSRCorr_binByBin->SetDirectory(0);
  if (!h1FSRCorr_binByBin->GetSumw2()) h1FSRCorr_binByBin->Sumw2();
  h1FSRCorr_binByBin->SetTitle("FSR correction bin-by-bin;M_{ee} [GeV];preFSR/postFSR");
  h1FSRCorr_binByBin->Divide(h1preFsr_MW,h1postFsr_MW,1.,1.,"B");

  TString fname="dyee_test_" + versionName(inpVersion) + TString(".root");
  if ((maxEntries>0) || includeOnlyRange.Length()) {
    fname.ReplaceAll(".root","_debug.root");
  }
  if (testSC && (fname.Index("_debug")==-1))
    fname.ReplaceAll(".root","_excludeGap.root");

  TFile fout(fname,"RECREATE");
  h1Eff->Write();
  h1EffPU->Write();
  h1Acc->Write();
  h1EffPUAcc->Write();
  h1FSRCorr_binByBin->Write();
  fsrResp.Write();
  effAccFsrResp.Write();
  h1postFsrInAccSel_M->Write();
  h1postFsrInAccSel_MW->Write();
  h1postFsrInAccSel_MWPU->Write();
  h1postFsrInAcc_M->Write();
  h1postFsrInAcc_MW->Write();
  h1postFsr_M->Write();
  h1postFsr_MW->Write();
  h1preFsr_M->Write();
  h1preFsr_MW->Write();
  h1preFsrUnf->Write();
  if (data.fH1PreFSR_binned) data.fH1PreFSR_binned->Write();
  h2FSRmig->Write();
  h1EffAccFsrUnf->Write();
  h2effAccFSRmig->Write();
  if (cPreFSR_WG) cPreFSR_WG->Write();
  if (cPreFSR_WGN) cPreFSR_WGN->Write();
  if (cPreFSRcheck) cPreFSRcheck->Write();
  cEffCmp->Write();
  cFSRTest->Write();
  cScaleToPreFsrTest->Write();
  TObjString timeTag(DayAndTimeTag(0));
  timeTag.Write("timeTag");
  fout.Close();
  std::cout << "file <" << fout.GetName() << "> created\n";
}

