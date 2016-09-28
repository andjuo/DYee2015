#include "DYmm13TeV.h"
#include "DYmm13TeV_eff.h"
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

void processDYmm_RECO(Int_t maxEntries=-1)
{
  std::cout << "processDYmm_RECO\n";
  std::cout << "DY range=" << DYtools::minMass << " .. " << DYtools::maxMass << "\n";

  TVersion_t inpVersion=_verMu1;
  inpVersion=_verMu76X;
  inpVersion=_verMuApproved;

  if (inpVersion!=_verMu1) {
    if (DYtools::nMassBins!=DYtools::nMassBins43) {
      std::cout << "a potential DYbinning.h problem\n";
      return;
    }
  }
  else {
    if (DYtools::nMassBins!=DYtools::nMassBins45) {
      std::cout << "a potential DYbinning.h problem\n";
      return;
    }
  }


  TString srcPath="/media/ssd/v20160214_1st_CovarianceMatrixInputs/";
  srcPath="/mnt/sdb/andriusj/v20160214_1st_CovarianceMatrixInputs/";
  TString fnameEff= srcPath + "Input5/ROOTFile_TagProbeEfficiency.root";
  TString dataFName=srcPath + "Input3/ROOTFile_ntuple_CovarianceMatrixInput.root";

  if (inpVersion==_verMu76X) {
    srcPath="/mnt/sdb/andriusj/v20160527_1st_CovarianceMatrixInputs_76X/";
    fnameEff= srcPath + "Input5/ROOTFile_TagProbeEfficiency_76X_v20160502.root";
    dataFName=srcPath + "Input3/ROOTFile_Input_CovarianceMatrix.root";
  }
  else if (inpVersion==_verMuApproved) {
    TString srcPath76X="/mnt/sdb/andriusj/v20160527_1st_CovarianceMatrixInputs_76X/";
    srcPath="/mnt/sdb/andriusj/v20160915_CovInput_ApprovedResults/";
    fnameEff= srcPath76X + "Input5/ROOTFile_TagProbeEfficiency_76X_v20160502.root";
    dataFName=srcPath + "ROOTFile_Input_CovarianceMatrix.root";
  }
  std::cout << "dataFName=" << dataFName << ", fnameEff=" << fnameEff << "\n";


  DYTnPEff_t tnpEff;
  TFile finEff(fnameEff);
  if (!finEff.IsOpen()) {
    std::cout << "failed to open the file <" << fnameEff << ">\n";
    return;
  }
  tnpEff.h2Eff_RecoID_Data= loadHisto(finEff,"h_2D_Eff_RecoID_Data","h2Eff_RecoID_Data",1,h2dummy);
  tnpEff.h2Eff_Iso_Data= loadHisto(finEff,"h_2D_Eff_Iso_Data","h2Eff_Iso_Data",1,h2dummy);
  tnpEff.h2Eff_HLT4p2_Data= loadHisto(finEff,"h_2D_Eff_HLTv4p2_Data","h2Eff_HLT4p2_Data",1,h2dummy);
  tnpEff.h2Eff_HLT4p3_Data= loadHisto(finEff,"h_2D_Eff_HLTv4p3_Data","h2Eff_HLT4p3_Data",1,h2dummy);
  tnpEff.h2Eff_RecoID_MC= loadHisto(finEff,"h_2D_Eff_RecoID_MC","h2Eff_RecoID_MC",1,h2dummy);
  tnpEff.h2Eff_Iso_MC= loadHisto(finEff,"h_2D_Eff_Iso_MC","h2Eff_Iso_MC",1,h2dummy);
  tnpEff.h2Eff_HLT4p2_MC= loadHisto(finEff,"h_2D_Eff_HLTv4p2_MC","h2Eff_HLT4p2_MC",1,h2dummy);
  tnpEff.h2Eff_HLT4p3_MC= loadHisto(finEff,"h_2D_Eff_HLTv4p3_MC","h2Eff_HLT4p3_MC",1,h2dummy);
  finEff.Close();
  if (!tnpEff.updateVectors()) {
    std::cout << "failed to load TnP efficiencies from <" << finEff.GetName() << ">\n";
    return;
  }
  TH2D* h2EffBinDef= cloneHisto(tnpEff.h2Eff_RecoID_Data,"h2EffBinDef","h2EffBinDef");
  if (!h2EffBinDef) return;

  DYmm13TeV_t data(dataFName);
  if (!data.fChain || (data.fChain->GetEntries()==0)) {
    std::cout << "error opening the file\n";
    return;
  }
  data.DeactivateBranches();
  data.ActivateBranches("Momentum_Reco_Lead Momentum_Reco_Sub Momentum_postFSR_Lead  Momentum_postFSR_Sub  Weight_Norm  Weight_PU  Weight_Gen");
  data.ActivateBranches("Flag_EventSelection");

  TH1D *h1recoSel_M= new TH1D("h1recoSel_M", "RECO selected;M_{reco} [GeV];count",DYtools::nMassBins,DYtools::massBinEdges);
  TH1D *h1recoSel_MW= new TH1D("h1recoSel_MW", "RECO selected (weighted);M_{reco} [GeV];weighted count",DYtools::nMassBins,DYtools::massBinEdges);
  TH1D *h1recoSel_MWPU= new TH1D("h1recoSel_MWPU", "RECO selected (weighted,wPU);M_{reco} [GeV];weighted (wPU) count",DYtools::nMassBins,DYtools::massBinEdges);

  TH1D *h1postFsrInAccSel_M= new TH1D("h1_postFsrInAccSel_M", "postFSR selected events in acceptance;M_{gen,postFSR} [GeV];count",DYtools::nMassBins,DYtools::massBinEdges);
  TH1D *h1postFsrInAccSel_MW= new TH1D("h1_postFsrInAccSel_MW", "postFSR selected events in acceptance (weighted);M_{gen,postFSR} [GeV];weighted count",DYtools::nMassBins,DYtools::massBinEdges);
  TH1D *h1postFsrInAccSel_MWPU= new TH1D("h1_postFsrInAccSel_MWPU", "postFSR selected events in acceptance (weighted,wPU);M_{gen,postFSR} [GeV];weighted (wPU) count",DYtools::nMassBins,DYtools::massBinEdges);

  prepareHisto(h1recoSel_M);
  prepareHisto(h1recoSel_MW);
  prepareHisto(h1recoSel_MWPU);
  prepareHisto(h1postFsrInAccSel_M);
  prepareHisto(h1postFsrInAccSel_MW);
  prepareHisto(h1postFsrInAccSel_MWPU);

  TH2D* h2DetResMig= new TH2D("h2DetResMig","Det.resolution migration", DYtools::nMassBins, DYtools::massBinEdges, DYtools::nMassBins, DYtools::massBinEdges);
  h2DetResMig->Sumw2();
  h2DetResMig->SetDirectory(0);
  RooUnfoldResponse detResResp(h1recoSel_MW,h1postFsrInAccSel_MW,h2DetResMig,"rooUnf_detResResp","detResResp;M_{#mu#mu} [GeV];det.res. unfolded yield");
  detResResp.UseOverflow(false);
  RooUnfoldResponse detResRespPU(h1recoSel_MW,h1postFsrInAccSel_MW,h2DetResMig,"rooUnf_detResRespPU","detResRespPU;M_{#mu#mu} [GeV];det.res.(PU) unfolded yield");
  detResRespPU.UseOverflow(false);

  EventSpace_t postFsrES("mainES",tnpEff.h2Eff_RecoID_Data);
  EventSpace_t recoES("recoES",tnpEff.h2Eff_RecoID_Data);
  // expanded event space data
  std::vector<TH2D*> h2PostFsrEventSpaceV;
  h2PostFsrEventSpaceV.reserve(DYtools::nMassBins);

  for (int im=0; im<DYtools::nMassBins; im++) {
    TString mStr= DYtools::massStr(im);
    TString hTitle=TString("post-FSR event space for ") + mStr + TString(" GeV;(eta,pt) flat index 1;(eta,pt) flat index 2");
    mStr.ReplaceAll("-","_");
    TH2D* h2PostFsrEventSpace= new TH2D("h2PostFsrEventSpace"+mStr,hTitle,
	   DYtools::EtaPtFIMax+1, -0.5, DYtools::EtaPtFIMax+0.5,
	   DYtools::EtaPtFIMax+1, -0.5, DYtools::EtaPtFIMax+0.5);
    h2PostFsrEventSpace->SetDirectory(0);
    h2PostFsrEventSpace->Sumw2();
    h2PostFsrEventSpaceV.push_back(h2PostFsrEventSpace);
  }

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

    if (!DYtools::InAcceptance_mm(data.Momentum_Reco_Lead,data.Momentum_Reco_Sub)) {
      const TLorentzVector *v1= data.Momentum_Reco_Lead;
      const TLorentzVector *v2= data.Momentum_Reco_Sub;
      std::cout << " Not in acceptance reco (eta,pt)= ("
		<< v1->Eta() << "," << v1->Pt() << ") and ("
		<< v2->Eta() << "," << v2->Pt() << ")\n";
    }


    double w= data.Weight_Norm * data.Weight_Gen;
    double wPU= w* data.Weight_PU;
    if (iEntry<1000) {
      std::cout << "w=" << w << ": weight_norm=" << data.Weight_Norm << ", gen=" << data.Weight_Gen << "\n";
      std::cout << "   wPU=" << wPU << ": weight_PU=" << data.Weight_PU << "\n";
    }
    double mReco= (*data.Momentum_Reco_Lead + *data.Momentum_Reco_Sub).M();
    double mPostFsr= (*data.Momentum_postFSR_Lead + *data.Momentum_postFSR_Sub).M();

    if (!recoES.fill(data.Momentum_Reco_Lead,data.Momentum_Reco_Sub,wPU)) {
      //std::cout << "recoES: not added\n";
    }

    if (DYtools::InsideMassRange(mReco)) {
      h1recoSel_M->Fill( mReco, 1. );
      h1recoSel_MW->Fill( mReco, w );
      h1recoSel_MWPU->Fill( mReco, wPU );
    }

    // Fill detResResponse
    if ( ! DYtools::InsideMassRange(mPostFsr) &&
	 DYtools::InsideMassRange(mReco) &&
	 data.Flag_EventSelection ) {
      if (mReco> DYtools::minMass) {
	std::cout << "fake " << mPostFsr << " (mReco=" << mReco << ")\n";
	detResResp.Fake(mReco, w);
	detResRespPU.Fake(mReco, wPU);
      }
    }
    else if ( DYtools::InsideMassRange(mPostFsr) &&
	      ! DYtools::InsideMassRange(mReco)
	      ) {
      std::cout << "miss " << mPostFsr << " (mReco=" << mReco << "), is selected=" << data.Flag_EventSelection << "\n";
      detResResp.Miss(mPostFsr, w);
      detResRespPU.Miss(mPostFsr, wPU);
    }
    else {
      detResResp.Fill( mReco, mPostFsr, w );
      detResRespPU.Fill( mReco, mPostFsr, wPU );
      h2DetResMig->Fill( mReco, mPostFsr, w );
    }


    if (DYtools::InAcceptance_mm(data.Momentum_postFSR_Lead,data.Momentum_postFSR_Sub)) {
      if (DYtools::InsideMassRange(mPostFsr)) {
	h1postFsrInAccSel_M->Fill( mPostFsr, 1. );
	h1postFsrInAccSel_MW->Fill( mPostFsr, w );
	h1postFsrInAccSel_MWPU->Fill( mPostFsr, wPU );
	//std::cout << "\nadding mPostFsr=" << mPostFsr << " with wPU=" << wPU << "\n";
      }

      //std::cout << "mPostFsr=" << mPostFsr << "\n";
      if (!postFsrES.fill(data.Momentum_postFSR_Lead,data.Momentum_postFSR_Sub,wPU)) {
	//std::cout << "postFsrES: not added\n";
      }

      int fi1= DYtools::FlatIndex(h2EffBinDef, data.Momentum_postFSR_Lead,1);
      int fi2= DYtools::FlatIndex(h2EffBinDef, data.Momentum_postFSR_Sub,1);
      if ((fi1>=0) && (fi2>=0)) {
	int idx= DYtools::massIdx(mPostFsr);
	if (idx!=-1)
	  h2PostFsrEventSpaceV[idx]->Fill(fi1,fi2, wPU);
      }
      else {
	const TLorentzVector *v1= data.Momentum_postFSR_Lead;
	const TLorentzVector *v2= data.Momentum_postFSR_Sub;
	std::cout << " iEntry=" << iEntry << "\n";
	std::cout << " failed postFSR fidx : " << v1->Eta() << "," << v1->Pt() << "; "
		  << v2->Eta() << "," << v2->Pt() << "\n";
	const TLorentzVector *v1r= data.Momentum_Reco_Lead;
	const TLorentzVector *v2r= data.Momentum_Reco_Sub;
	std::cout << "        Reco      : " << v1r->Eta() << "," << v1r->Pt() << "; "
		  << v2r->Eta() << "," << v2r->Pt() << "\n";
	return;
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
  printRatio(h1postFsrInAccSel_MW, h1Unf);

  RooUnfoldBayes bayesDetResPU( &detResRespPU, h1recoSel_MWPU, 4 );
  TH1D *h1UnfPU= (TH1D*) bayesDetResPU.Hreco()->Clone("h1UnfPU");
  h1UnfPU->SetTitle("unfolded (PU) reco->postFSR");
  h1UnfPU->SetDirectory(0);
  h1UnfPU->SetLineColor(kBlue);
  h1UnfPU->SetMarkerStyle(5);
  h1UnfPU->SetMarkerColor(kBlue);
  TCanvas *cDetResTestPU=plotHisto(h1UnfPU,"cDetResTestPU",1,1,"LPE1");
  plotHistoSame(h1postFsrInAccSel_MWPU,"cDetResTestPU","LPE");
  printRatio(h1postFsrInAccSel_MWPU,h1UnfPU);

  TH1D* h1UnfBinByBinCorr=(TH1D*)h1recoSel_MW->Clone("h1UnfBinByBinCorr");
  h1UnfBinByBinCorr->SetDirectory(0);
  if (!h1UnfBinByBinCorr->GetSumw2()) h1UnfBinByBinCorr->Sumw2();
  h1UnfBinByBinCorr->SetTitle("Unf corr bin-by-bin;M_{#mu#mu} [GeV];postFSR_inAcc_Sel/recoSel");
  h1UnfBinByBinCorr->Divide(h1postFsrInAccSel_MW,h1recoSel_MW,1,1,"B");

  // check postFSR event space
  if (1) {
    TH1D* h1esChk= cloneHisto(h1recoSel_MWPU,"h1esChk","h1esChk");
    h1esChk->Reset();
    for (int iM=0; iM<DYtools::nMassBins; iM++) {
      const TH2D *h2es= postFsrES.h2ES(iM);
      double sum=0;
      for (int ibin=1; ibin<=h2es->GetNbinsX(); ibin++) {
	for (int jbin=1; jbin<=h2es->GetNbinsY(); jbin++) {
	  sum+= h2es->GetBinContent(ibin,jbin);
	}
      }
      h1esChk->SetBinContent(iM+1, sum);
    }
    histoStyle(h1esChk,kRed,24);
    plotHisto(h1postFsrInAccSel_MWPU,"cChkES",1,1,"LPE1","postFSRInAccSel_MWPU");
    plotHistoSame(h1esChk,"cChkES","LPE","h1esChk");
    printRatio(h1postFsrInAccSel_MWPU, h1esChk);
  }




  TString fname="dymm_test_RECO_" + versionName(inpVersion) + TString(".root");
  //fname.ReplaceAll(".root","_withOverflow.root");
  if (maxEntries>0) fname.ReplaceAll(".root","_debug.root");
  TFile fout(fname,"RECREATE");
  tnpEff.save(fout);
  h1UnfBinByBinCorr->Write();
  detResResp.Write();
  detResRespPU.Write();
  h1recoSel_M->Write();
  h1recoSel_MW->Write();
  h1recoSel_MWPU->Write();
  h1postFsrInAccSel_M->Write();
  h1postFsrInAccSel_MW->Write();
  h1postFsrInAccSel_MWPU->Write();
  h1Unf->Write();
  h1UnfPU->Write();
  h2DetResMig->Write();
  cDetResTest->Write();
  cDetResTestPU->Write();
  recoES.save(fout);
  postFsrES.save(fout);
  fout.mkdir("eventSpace");
  fout.cd("eventSpace");
  for (unsigned int i=0; i<h2PostFsrEventSpaceV.size(); i++)
    h2PostFsrEventSpaceV[i]->Write();
  fout.cd();
  TObjString timeTag(DayAndTimeTag(0));
  timeTag.Write("timeTag");

  fout.Close();
  std::cout << "file <" << fout.GetName() << "> created\n";
}

