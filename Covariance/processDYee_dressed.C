//#include "DYmm13TeV.h"
#include "DYee13TeV.h"
#include "DYmm13TeV_eff.h"
#include "DYbinning.h"
#include "inputs.h"
#include "../RooUnfold/src/RooUnfoldResponse.h"
#include "../RooUnfold/src/RooUnfoldBayes.h"

// --------------------------------------------------------------

template<class histo_t>
inline
void prepareHisto(histo_t *h1) {
  h1->SetDirectory(0);
  //if (!h1->GetSumw2())
  h1->Sumw2();
}

// --------------------------------------------------------------
// --------------------------------------------------------------

void processDYee_dressed(Int_t maxEntries=100, TString includeOnlyRange="")
{
  std::cout << "processDYee_dressed\n";
  std::cout << "DY range=" << DYtools::minMass << " .. " << DYtools::maxMass << "\n";

  TVersion_t inpVersion=_verEl2skim;
  const int testSC_postFSR=1;

  if (DYtools::nMassBins!=DYtools::nMassBins43) {
    std::cout << "a potential DYbinning.h problem\n";
    return;
  }

  TString srcPath="/mnt/sdb/andriusj/v2ee-skim1-20160812/";
  TString fnameEff="/mnt/sdb/andriusj/v2_08082016_CovarianceMatrixInputs/ROOTFile_Input5_TagProbeEfficiency.root";

  inpVersion=_verEl2skim3;
  srcPath="/media/sf_Share2/DYee_76X_Calibrated/mySkim/";

  inpVersion=_verEl3;
  fnameEff="/mnt/sdb/andriusj/v3_09092016_CovarianceMatrixInputs/ROOTFile_Input5_TagProbeEfficiency.root";

  TString dataFName;
  TString fileFormat="DYEE_M%dto%d_v1.root ";
  if ((inpVersion==_verEl2skim2) ||
      (inpVersion==_verEl2skim3) ||
      (inpVersion==_verEl3)) fileFormat="DY_%dto%d_v2_orig.root ";

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

  DYTnPEff_t tnpEff(2); // 0 - muon channel, 1,2 - electron channel
  TFile finEff(fnameEff);
  if (!finEff.IsOpen()) {
    std::cout << "failed to open the file <" << fnameEff << ">\n";
    return;
  }
  if (!tnpEff.load(finEff,"","")) {
    std::cout << "failed to load tnpEff from file <" << fnameEff << ">\n";
    std::cout << "isElChannel= " << tnpEff.isElChannel() << "\n";
    return;
  }
  finEff.Close();

  // h2EffBinDef needed by flat idx
  TH2D* h2EffBinDef= cloneHisto(tnpEff.h2Eff_RecoID_Data,"h2EffBinDef","h2EffBinDef");
  if (!h2EffBinDef) return;

  DYee13TeV_t data(dataFName);
  if (!data.fChain || (data.fChain->GetEntries()==0)) {
    std::cout << "error opening the file\n";
    return;
  }
  //data.DeactivateBranches(); // does not work
  //data.ActivateBranches("Momentum_Reco_Lead Momentum_Reco_Sub Momentum_postFSR_Lead  Momentum_postFSR_Sub  Momentum_preFSR_Lead Momentum_preFSR_Sub  Weight_Norm  Weight_PU  Weight_Gen");
  //data.ActivateBranches("Flag_EventSelection Flag_EventSelectionExceptSCGap");
  //data.ActivateBranches("SCEta_Lead SCEta_Sub");

  TH1D *h1recoSel_M= new TH1D("h1recoSel_M", "RECO selected;M_{reco} [GeV];count",DYtools::nMassBins,DYtools::massBinEdges);
  TH1D *h1recoSel_MW= new TH1D("h1recoSel_MW", "RECO selected (weighted);M_{reco} [GeV];weighted count",DYtools::nMassBins,DYtools::massBinEdges);
  TH1D *h1recoSel_MWPU= new TH1D("h1recoSel_MWPU", "RECO selected (wPU);M_{reco} [GeV];weighted (wPU) count",DYtools::nMassBins,DYtools::massBinEdges);
  TH1D *h1recoSelInPostFsrAcc_MWPU= new TH1D("h1recoSelInPostFsrAcc_MWPU","RECO selected, in postFSR acc (wPU);M_{reco} [GeV];weighted (wPU) count",DYtools::nMassBins,DYtools::massBinEdges);

  TH1D *h1postFsrInAccSel_M= new TH1D("h1_postFsrInAccSel_M", "postFSR selected events in acceptance;M_{gen,postFSR} [GeV];count",DYtools::nMassBins,DYtools::massBinEdges);
  TH1D *h1postFsrInAccSel_MW= new TH1D("h1_postFsrInAccSel_MW", "postFSR selected events in acceptance (weighted);M_{gen,postFSR} [GeV];weighted count",DYtools::nMassBins,DYtools::massBinEdges);
  TH1D *h1postFsrInAccSel_MWPU= new TH1D("h1_postFsrInAccSel_MWPU", "postFSR selected events in acceptance (weighted,wPU);M_{gen,postFSR} [GeV];weighted (wPU) count",DYtools::nMassBins,DYtools::massBinEdges);
  TH1D *h1postFsrInAccSelNoRecoGap_MWPU= new TH1D("h1_postFsrInAccSelNoRecoGap_MWPU", "postFSR selected events in acceptance (ignoring recoSCGap) (weighted,wPU);M_{gen,postFSR} [GeV];weighted (wPU) count",DYtools::nMassBins,DYtools::massBinEdges);
  TH1D *h1postFsrInAccSel_MissAndMig_MWPU= new TH1D("h1_postFsrInAccSel_MissAndMig_MWPU", "postFSR selected events in acceptance to check detResRespPU (weighted,wPU);M_{gen,postFSR} [GeV];weighted (wPU) count",DYtools::nMassBins,DYtools::massBinEdges);
  TH1D *h1postFsrInAccSel_Mig_MWPU= new TH1D("h1_postFsrInAccSel_Mig_MWPU", "postFSR selected events in acceptance to check detResRespPU (weighted,wPU);M_{gen,postFSR} [GeV];weighted (wPU) count",DYtools::nMassBins,DYtools::massBinEdges);
  TH1D *h1postFsrInAccSel_DieleOk_MWPU= new TH1D("h1_postFsrInAccSel_DieleOk_MWPU","postFSR selected events in acceptance when dielectron exists (weighted,wPU);M_{gen,postFSR} [GeV];weighted (wPU) count",DYtools::nMassBins,DYtools::massBinEdges);

  prepareHisto(h1recoSel_M);
  prepareHisto(h1recoSel_MW);
  prepareHisto(h1recoSel_MWPU);
  prepareHisto(h1recoSelInPostFsrAcc_MWPU);
  prepareHisto(h1postFsrInAccSel_M);
  prepareHisto(h1postFsrInAccSel_MW);
  prepareHisto(h1postFsrInAccSel_MWPU);
  prepareHisto(h1postFsrInAccSelNoRecoGap_MWPU);
  prepareHisto(h1postFsrInAccSel_MissAndMig_MWPU);
  prepareHisto(h1postFsrInAccSel_Mig_MWPU);
  prepareHisto(h1postFsrInAccSel_DieleOk_MWPU);

  TH1D *h1detRes_Miss= new TH1D("h1detRes_Miss","h1detRes_Miss;M_{postFSR};count",DYtools::nMassBins,DYtools::massBinEdges);
  TH1D *h1detRes_Fake= new TH1D("h1detRes_Fake","h1detRes_Fake;M_{RECO};count",DYtools::nMassBins,DYtools::massBinEdges);
  TH1D *h1detRes_Meas= new TH1D("h1detRes_Meas","h1detRes_Meas;M_{RECO};count",DYtools::nMassBins,DYtools::massBinEdges);
  TH1D *h1detRes_True= new TH1D("h1detRes_True","h1detRes_True;M_{postFSR};count",DYtools::nMassBins,DYtools::massBinEdges);
  TH2D *h2detRes_Mig= new TH2D("h2detRes_Mig","h2detRes_Mig;M_{RECO};M_{postFSR}",DYtools::nMassBins,DYtools::massBinEdges, DYtools::nMassBins,DYtools::massBinEdges);

  prepareHisto(h1detRes_Miss);
  prepareHisto(h1detRes_Fake);
  prepareHisto(h1detRes_Meas);
  prepareHisto(h1detRes_True);
  prepareHisto(h2detRes_Mig);

  TH2D* h2DetResMig= new TH2D("h2DetResMig","Det.resolution migration", DYtools::nMassBins, DYtools::massBinEdges, DYtools::nMassBins, DYtools::massBinEdges);
  h2DetResMig->Sumw2();
  h2DetResMig->SetDirectory(0);
  RooUnfoldResponse detResResp(h1recoSel_MW,h1postFsrInAccSel_MW,h2DetResMig,"rooUnf_detResResp","detResResp;M_{ee} [GeV];det.res. unfolded yield");
  detResResp.UseOverflow(false);
  RooUnfoldResponse detResRespPU(h1recoSel_MWPU,h1postFsrInAccSel_MWPU,h2DetResMig,"rooUnf_detResRespPU","detResRespPU;M_{ee} [GeV];det.res.(PU) unfolded yield");
  detResRespPU.UseOverflow(false);

  TH1D *h1postFsrInAcc_M= new TH1D("h1_postFsrInAcc_M", "postFSR events in acceptance;M_{gen,postFSR} [GeV];count",DYtools::nMassBins,DYtools::massBinEdges);
  TH1D *h1postFsrInAcc_MW= new TH1D("h1_postFsrInAcc_MW", "postFSR events in acceptance (weighted);M_{gen,postFSR} [GeV];weighted count",DYtools::nMassBins,DYtools::massBinEdges);
  TH1D *h1postFsrInAcc_MWPU= new TH1D("h1_postFsrInAcc_MWPU", "postFSR events in acceptance (weighted,wPU);M_{gen,postFSR} [GeV];weighted (wPU) count",DYtools::nMassBins,DYtools::massBinEdges);

  prepareHisto(h1postFsrInAcc_M);
  prepareHisto(h1postFsrInAcc_MW);
  prepareHisto(h1postFsrInAcc_MWPU);

  TH1D *h1postFsr_M= new TH1D("h1_postFsr_M", "postFSR;M_{gen,postFSR} [GeV];count",DYtools::nMassBins, DYtools::massBinEdges);
  TH1D *h1postFsr_MW= new TH1D("h1_postFsr_Mweighted", "postFSR;M_{gen,postFSR} [GeV];weighted count", DYtools::nMassBins, DYtools::massBinEdges);
  TH1D *h1preFsr_M= new TH1D("h1_preFsr_M", "preFSR;M_{gen,preFSR} [GeV];count",DYtools::nMassBins, DYtools::massBinEdges);
  TH1D *h1preFsr_MW= new TH1D("h1_preFsr_Mweighted", "preFSR;M_{gen,preFSR} [GeV];weighted count", DYtools::nMassBins, DYtools::massBinEdges);

  prepareHisto(h1postFsr_M);
  prepareHisto(h1postFsr_MW);
  prepareHisto(h1preFsr_M);
  prepareHisto(h1preFsr_MW);

  TH2D* h2FSRmig= new TH2D("h2FSRmig","FSR migration", DYtools::nMassBins, DYtools::massBinEdges, DYtools::nMassBins, DYtools::massBinEdges);
  h2FSRmig->Sumw2();
  RooUnfoldResponse fsrResp(h1postFsr_MW,h1preFsr_MW,h2FSRmig,"rooUnf_fsrResp","fsrResp;M_{ee} [GeV];FSR unfolded yield");
  fsrResp.UseOverflow(false);

  TH2D* h2detResEffAccFSRmig= new TH2D("h2detResEffAccFSRmig","detRes x eff x Acc x FSR migration", DYtools::nMassBins, DYtools::massBinEdges, DYtools::nMassBins, DYtools::massBinEdges);
  h2detResEffAccFSRmig->Sumw2();
  RooUnfoldResponse detResEffAccFsrResp(h1recoSel_MW,h1preFsr_MW,h2detResEffAccFSRmig,"rooUnf_detResEffAccFsrResp","detResEffAccFsrResp;M_{ee} [GeV];detRes eff Acc FSR unfolded yield");
  detResEffAccFsrResp.UseOverflow(false);

  EventSpace_t postFsrES("mainES",tnpEff.h2Eff_RecoID_Data);
  EventSpace_t postFsrES_trigOnly("mainES_trigOnly",tnpEff.h2Eff_RecoID_Data);
  EventSpace_t postFsrES_trigOnly_noPU("mainES_trigOnly_noPU",tnpEff.h2Eff_RecoID_Data);
  EventSpace_t recoES("recoES",tnpEff.h2Eff_RecoID_Data);
  EventSpace_t recoESinPostFsrAcc("recoESinPostFsrAcc",tnpEff.h2Eff_RecoID_Data);

  UInt_t nEvents= data.GetEntries();
  UInt_t nProcessed=0, nSelected=0;

  std::cout << "process data\n";
  for (UInt_t iEntry=0; iEntry<nEvents; iEntry++) {
    if (data.GetEntry(iEntry)<0) break;
    if (iEntry%100000==0) std::cout << "iEntry=" << iEntry << Form("(%4.2lf%%)\n",iEntry*100./double(nEvents));
    if ((maxEntries>0) && (iEntry>UInt_t(maxEntries))) {
      std::cout << "debug run\n";
      break;
    }
    nProcessed++;


    double w= data.Weight_Norm * data.Weight_Gen;
    double wPU= w* data.Weight_PU;

    // Flag_RecoEleSelection = Flag_EventSelection + dielectron in acceptance
    int inAccReco=data.Flag_RecoEleSelection;

    //if (data.Flag_EventSelection && !data.Flag_RecoEleSelection) {
    //  std::cout << "Flag_EventSelection on, RecoEleSelection off\n";
    //}
    //if (!data.Flag_EventSelection && data.Flag_RecoEleSelection) {
    //  std::cout << "Flag_EventSelection OFF, RecoEleSelection on\n";
    //}

    int inAccPostFsr=DYtools::InAcceptance_ee(data.Momentum_postFSR_Lead,data.Momentum_postFSR_Sub, testSC_postFSR);
    double mReco= (*data.Momentum_Reco_Lead + *data.Momentum_Reco_Sub).M();
    double mPostFsr= (*data.Momentum_postFSR_Lead + *data.Momentum_postFSR_Sub).M();
    double mPreFsr= (*data.Momentum_preFSR_Lead + *data.Momentum_preFSR_Sub).M();

    if (0) if (iEntry<1000) {
      std::cout << "w=" << w << ": weight_norm=" << data.Weight_Norm << ", gen=" << data.Weight_Gen << "\n";
      std::cout << "   wPU=" << wPU << ": weight_PU=" << data.Weight_PU << "\n";
    }

    if (inAccPostFsr) {
      if (data.Flag_EventSelectionExceptSCGap)
	h1postFsrInAccSelNoRecoGap_MWPU->Fill( mPostFsr, wPU );
      if (data.Flag_EventSelection) {
	h1postFsrInAccSel_M->Fill( mPostFsr, 1. );
	h1postFsrInAccSel_MW->Fill( mPostFsr, w );
	h1postFsrInAccSel_MWPU->Fill( mPostFsr, wPU );
	if (data.Flag_RecoEleSelection)
	  h1postFsrInAccSel_DieleOk_MWPU->Fill( mPostFsr, wPU );
      }
      h1postFsrInAcc_M->Fill( mPostFsr, 1. );
      h1postFsrInAcc_MW->Fill( mPostFsr, w );
      h1postFsrInAcc_MWPU->Fill( mPostFsr, wPU );
    }
    h1postFsr_M->Fill( mPostFsr, 1. );
    h1postFsr_MW->Fill( mPostFsr, w );
    h1preFsr_M->Fill( mPreFsr, 1. );
    h1preFsr_MW->Fill( mPreFsr, w );
    h2FSRmig->Fill( mPostFsr, mPreFsr, w );
    if ( ! DYtools::InsideMassRange(mPreFsr) &&
	 DYtools::InsideMassRange(mPostFsr) ) {
      //std::cout << "fake detected\n";
      fsrResp.Fake(mPostFsr,w);
    }
    else if ( DYtools::InsideMassRange(mPreFsr) &&
	      ! DYtools::InsideMassRange(mPostFsr) ) {
      fsrResp.Miss(mPreFsr,w);
    }
    else {
      fsrResp.Fill(mPostFsr,mPreFsr,w);
    }

    int recoRequirements= (data.Flag_EventSelection &&
			   data.Flag_RecoEleSelection &&
			   DYtools::InsideMassRange(mReco) && inAccReco);
    if (! DYtools::InsideMassRange(mPreFsr) &&
	recoRequirements ) {
      detResEffAccFsrResp.Fake(mReco, w);
    }
    else if ( DYtools::InsideMassRange(mPreFsr) &&
	      ! recoRequirements ) {
      detResEffAccFsrResp.Miss(mPreFsr, w);
    }
    else if ( DYtools::InsideMassRange(mPreFsr) &&
	      recoRequirements ) {
      detResEffAccFsrResp.Fill( mReco, mPreFsr, w );
    }

    // Fill detResResponse
    if (data.Flag_EventSelection) {
      if ( ! (DYtools::InsideMassRange(mPostFsr) && inAccPostFsr) &&
	   DYtools::InsideMassRange(mReco) && inAccReco) {
	if (mReco> DYtools::minMass) {
	  //std::cout << "fake " << mPostFsr << " (mReco=" << mReco << ")\n";
	  detResResp.Fake(mReco, w);
	  detResRespPU.Fake(mReco, wPU);
	  h1detRes_Fake->Fill(mReco, wPU);
	}
      }
      else if ( DYtools::InsideMassRange(mPostFsr) && inAccPostFsr &&
		! ( DYtools::InsideMassRange(mReco)  && inAccReco )
		) {
	std::cout << "miss " << mPostFsr << " (mReco=" << mReco << "), is selected=" << data.Flag_EventSelection << "\n";
	detResResp.Miss(mPostFsr, w);
	detResRespPU.Miss(mPostFsr, wPU);
	h1postFsrInAccSel_MissAndMig_MWPU->Fill(mPostFsr, wPU);
	h1detRes_Miss->Fill(mPostFsr, wPU);
      }
      else {
	if (inAccReco && inAccPostFsr) {
	  if (( DYtools::InsideMassRange(mPostFsr) && !DYtools::InsideMassRange(mReco)) ||
	      (!DYtools::InsideMassRange(mPostFsr) &&  DYtools::InsideMassRange(mReco)) ){
	    std::cout << "in acceptances mReco=" << mReco << ", mPostFsr=" << mPostFsr << "\n";
	  }
	}
	if (inAccReco && inAccPostFsr && DYtools::InsideMassRange(mPostFsr) && DYtools::InsideMassRange(mReco)) {
	  detResResp.Fill( mReco, mPostFsr, w );
	  detResRespPU.Fill( mReco, mPostFsr, wPU );
	  h2DetResMig->Fill( mReco, mPostFsr, w );
	  h1postFsrInAccSel_MissAndMig_MWPU->Fill(mPostFsr, wPU);
	  h1postFsrInAccSel_Mig_MWPU->Fill(mPostFsr, wPU);
	  h1detRes_Meas->Fill( mReco, wPU );
	  h1detRes_True->Fill( mPostFsr, wPU );
	  h2detRes_Mig->Fill( mReco, mPostFsr, wPU);
	}
      }
    }


    if (!data.Flag_EventSelection) continue;
    nSelected++;

    // Introduction of Flag_RecoEleSelection took out this check from
    // Flag_EventSelection
    //// event selection flag guarantees that testSC was enabled
    //if (!DYtools::InAcceptance_ee(data.Momentum_Reco_Lead,data.Momentum_Reco_Sub, 0)) {
    // std::cout << "event selected but not in acceptance reco\n";
    //  return;
    //}

    //
    // probably not needed
    //
    if (1) {
      const TLorentzVector *v1= data.Momentum_Reco_Lead;
      const TLorentzVector *v2= data.Momentum_Reco_Sub;
      int ptOk=
	( (((v1->Pt()>30) && (v2->Pt()>10)) ||
	   ((v1->Pt()>10) && (v2->Pt()>30))) ) ? 1 :0;

      if (!DYtools::InAcceptance_ee(data.Momentum_Reco_Lead,data.Momentum_Reco_Sub, 0)) {
	if (ptOk) {
	  std::cout << " Not in acceptance reco (eta,pt)= ("
		    << v1->Eta() << "," << v1->Pt() << ") and ("
		    << v2->Eta() << "," << v2->Pt() << ")\n";
	}
      }
    }

    if (data.Flag_RecoEleSelection) {
      h1recoSel_M->Fill( mReco, 1. );
      h1recoSel_MW->Fill( mReco, w );
      h1recoSel_MWPU->Fill( mReco, wPU );
      recoES.fill(data.Momentum_Reco_Lead,data.Momentum_Reco_Sub,wPU);
      if (inAccPostFsr) {
	h1recoSelInPostFsrAcc_MWPU->Fill( mReco, wPU );
      }
    }

    if (inAccPostFsr) {
      recoESinPostFsrAcc.fill(data.Momentum_Reco_Lead,data.Momentum_Reco_Sub,wPU);
    }

    //std::cout << "mPostFsr=" << mPostFsr << "\n";
    postFsrES_trigOnly.fill(data.Momentum_postFSR_Lead,data.Momentum_postFSR_Sub,wPU);
    postFsrES_trigOnly_noPU.fill(data.Momentum_postFSR_Lead,data.Momentum_postFSR_Sub,w);

    //
    // It is needed
    //
    if (inAccPostFsr) {
      postFsrES.fill(data.Momentum_postFSR_Lead,data.Momentum_postFSR_Sub,wPU);
      int fi1= DYtools::FlatIndex(h2EffBinDef, data.Momentum_postFSR_Lead,1);
      int fi2= DYtools::FlatIndex(h2EffBinDef, data.Momentum_postFSR_Sub,1);
      if ((fi1<0) || (fi2<0)) {
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

  TCanvas *cPreFSR_WG=NULL, *cPreFSR_WGN=NULL;
  if (data.fH1PreFSR_WGN.size()) {
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


  TH1D* h1rho= postFsrES.calculateScaleFactor(tnpEff,0,"h1rho","scale factor;M_{ee] [GeV];#rho");
  TH1D* h1rhoTrigOnly= postFsrES_trigOnly.calculateScaleFactor(tnpEff,0,"h1rho_trigOnly","scale factor (trigOnly)");
  TH1D* h1rhoTrigOnlyNoPU= postFsrES_trigOnly_noPU.calculateScaleFactor(tnpEff,0,"h1rho_trigOnly_noPU","scale factor (no PU)");
  TH1D* h1rho_reco= recoES.calculateScaleFactor(tnpEff,0,"h1rho_reco","scale factor (reco space)");
  TH1D* h1rho_recoInPostFsrAcc= recoESinPostFsrAcc.calculateScaleFactor(tnpEff,0,"h1rho_recoInPostFsrAcc","scale factor (reco space, in postFSR acc)");
  histoStyle(h1rho,kBlack,24,1);
  histoStyle(h1rhoTrigOnly,kGreen+1,5,3);
  histoStyle(h1rhoTrigOnlyNoPU,kBlue,7,1);
  histoStyle(h1rho_reco,kOrange,5,1);
  histoStyle(h1rho_recoInPostFsrAcc,kRed,7,1);
  plotHisto(h1rho,"cRho",1,0,"LPE1","rho (postFSRES)");
  plotHistoSame(h1rhoTrigOnly,"cRho","LPE1","rho (trigOnly)");
  plotHistoSame(h1rhoTrigOnlyNoPU,"cRho","LPE1","rho (trigOnly_noPU)");
  plotHistoSame(h1rho_reco,"cRho","LPE1","rho (reco)");
  plotHistoSame(h1rho_recoInPostFsrAcc,"cRho","LPE1","rho (reco inPostFsrAcc)");

  RooUnfoldBayes bayesDetRes( &detResResp, h1recoSel_MW, 4 );
  TH1D *h1Unf= (TH1D*) bayesDetRes.Hreco()->Clone("h1Unf");
  h1Unf->SetTitle("unfolded reco->postFSR");
  h1Unf->SetDirectory(0);
  h1Unf->SetLineColor(kRed);
  h1Unf->SetMarkerStyle(24);
  h1Unf->SetMarkerColor(kRed);
  h1Unf->GetXaxis()->SetMoreLogLabels();
  h1Unf->GetXaxis()->SetNoExponent();
  TCanvas *cDetResTest=plotHisto(h1Unf,"cDetResTest",1,1,"LPE1","BayesUnf4");
  plotHistoSame(h1postFsrInAccSel_MW,"cDetResTest","LPE","postFsrInAccSel_MW");
  printRatio(h1postFsrInAccSel_MW, h1Unf);

  RooUnfoldBayes bayesDetResPU( &detResRespPU, h1recoSel_MWPU, 4 );
  TH1D *h1UnfPU= (TH1D*) bayesDetResPU.Hreco()->Clone("h1UnfPU");
  h1UnfPU->SetTitle("unfolded (PU) reco->postFSR");
  h1UnfPU->SetDirectory(0);
  h1UnfPU->SetLineColor(kBlue);
  h1UnfPU->SetMarkerStyle(5);
  h1UnfPU->SetMarkerColor(kBlue);
  TCanvas *cDetResTestPU=plotHisto(h1UnfPU,"cDetResTestPU",1,1,"LPE1","BayesUnf4PU");
  plotHistoSame(h1postFsrInAccSel_MWPU,"cDetResTestPU","LPE","postFsrInAccSel_MWPU");
  printRatio(h1postFsrInAccSel_MWPU,h1UnfPU);

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

  RooUnfoldBayes bayesFSR( &fsrResp, h1postFsr_MW, 4 );
  TH1D *h1preFsrUnf= (TH1D*) bayesFSR.Hreco()->Clone("h1preFsrUnf");
  h1preFsrUnf->SetTitle("unfolded postFSR->preFSR");
  h1preFsrUnf->SetDirectory(0);
  histoStyle(h1preFsrUnf,kRed,24);
  h1preFsrUnf->GetXaxis()->SetMoreLogLabels();
  h1preFsrUnf->GetXaxis()->SetNoExponent();
  TCanvas *cFSRTest=plotHisto(h1preFsrUnf,"cFSRTest",1,1,"LPE1","BayesUnf4");
  plotHistoSame(h1preFsr_MW,"cFSRTest","LPE","preFsr_MW");

  RooUnfoldBayes bayesDetResEffAccFsr( &detResEffAccFsrResp, h1recoSel_MWPU, 4 );
  TH1D *h1DetResEffAccFsrUnf= (TH1D*) bayesDetResEffAccFsr.Hreco()->Clone("h1DetResEffAccFsrUnf");
  h1DetResEffAccFsrUnf->SetTitle("unfolded recoSel -> preFSR");
  h1DetResEffAccFsrUnf->SetDirectory(0);
  histoStyle(h1DetResEffAccFsrUnf,kRed,24);
  h1DetResEffAccFsrUnf->GetXaxis()->SetMoreLogLabels();
  h1DetResEffAccFsrUnf->GetXaxis()->SetNoExponent();
  TCanvas *cScaleRecoSelToPreFsrTest= plotHisto(h1DetResEffAccFsrUnf,"cScaleRecoSelToPreFsrTest",1,1,"LPE1","RecoSel->preFSR(BayesUnf4)");
  histoStyle(h1preFsr_MW,kBlue,5);
  plotHistoSame(h1preFsr_MW,"cScaleRecoSelToPreFsrTest","LPE","preFsr_MW");

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

  TH1D* h1EffPUNoRecoGap=(TH1D*)h1postFsrInAccSelNoRecoGap_MWPU->Clone("h1EffPUNoRecoGap");
  h1EffPUNoRecoGap->SetDirectory(0);
  if (!h1EffPUNoRecoGap->GetSumw2()) h1EffPUNoRecoGap->Sumw2();
  h1EffPUNoRecoGap->SetTitle("Efficiency (wPU noRecoGap/noPU);M_{ee,GEN postFSR} [GeV];postFSR_inAcc_SelNoRecoGap(wPU)/postFSR_inAcc(noPU)");
  h1EffPUNoRecoGap->Divide(h1postFsrInAccSelNoRecoGap_MWPU,h1postFsrInAcc_MW,1,1,"B");

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

  TH1D *h1EffPUEleMatchAcc=(TH1D*)h1postFsrInAccSel_DieleOk_MWPU->Clone("h1EffPUEleMatchAcc");
  h1EffPUEleMatchAcc->SetDirectory(0);
  if (!h1EffPUEleMatchAcc->GetSumw2()) h1EffPUEleMatchAcc->Sumw2();
  h1EffPUEleMatchAcc->SetTitle("Efficiency(wPU,eleMatch) #times Acc;M_{ee,GEN postFSR} [GeV];postFSR_inAccSel_DieleOk(wPU)/h1postFsr_MW(noPU)");
  h1EffPUEleMatchAcc->Divide(h1EffPUEleMatchAcc,h1postFsr_MW,1,1,"B");


  TString fname="dyee_test_dressed_" + versionName(inpVersion)+TString(".root");
  //fname.ReplaceAll(".root","_withOverflow.root");
  if ((maxEntries>0) || includeOnlyRange.Length()) {
    fname.ReplaceAll(".root","_debug.root");
  }
  TFile foutH(fname,"RECREATE");
  tnpEff.save(foutH);
  h1rho->Write();
  h1rhoTrigOnly->Write();
  h1rhoTrigOnlyNoPU->Write();
  h1rho_reco->Write();
  h1rho_recoInPostFsrAcc->Write();
  detResResp.Write();
  detResRespPU.Write();
  h1detRes_Miss->Write();
  h1detRes_Fake->Write();
  h1detRes_Meas->Write();
  h1detRes_True->Write();
  h2detRes_Mig->Write();
  h1recoSel_M->Write();
  h1recoSel_MW->Write();
  h1recoSel_MWPU->Write();
  h1recoSelInPostFsrAcc_MWPU->Write();
  h1postFsrInAccSel_M->Write();
  h1postFsrInAccSel_MW->Write();
  h1postFsrInAccSel_MWPU->Write();
  h1postFsrInAccSelNoRecoGap_MWPU->Write();
  h1postFsrInAccSel_MissAndMig_MWPU->Write();
  h1postFsrInAccSel_Mig_MWPU->Write();
  h1postFsrInAccSel_DieleOk_MWPU->Write();
  h1postFsrInAcc_M->Write();
  h1postFsrInAcc_MW->Write();
  h1postFsr_M->Write();
  h1postFsr_MW->Write();
  h1preFsr_M->Write();
  h1preFsr_MW->Write();
  h1Unf->Write();
  h1UnfPU->Write();
  h2DetResMig->Write();
  h1Eff->Write();
  h1EffPU->Write();
  h1EffPUNoRecoGap->Write();
  h1Acc->Write();
  h1EffPUAcc->Write();
  h1EffPUEleMatchAcc->Write();
  fsrResp.Write();
  detResEffAccFsrResp.Write();
  h1preFsrUnf->Write();
  h1DetResEffAccFsrUnf->Write();
  if (data.fH1PreFSR_binned) data.fH1PreFSR_binned->Write();
  h2FSRmig->Write();
  h1DetResEffAccFsrUnf->Write();
  //h2detResEffAccFsrMig->Write();
  //if (cPreFSR_WG) cPreFSR_WG->Write();
  if (cPreFSR_WGN) cPreFSR_WGN->Write();
  cDetResTest->Write();
  cDetResTestPU->Write();
  cEffCmp->Write();
  cFSRTest->Write();
  cScaleRecoSelToPreFsrTest->Write();
  recoES.save(foutH);
  recoESinPostFsrAcc.save(foutH);
  postFsrES.save(foutH);
  postFsrES_trigOnly.save(foutH);
  postFsrES_trigOnly_noPU.save(foutH);
  TObjString timeTag(DayAndTimeTag(0));
  timeTag.Write("timeTag");

  foutH.Close();
  std::cout << "file <" << foutH.GetName() << "> created\n";
}


// ------------------------------------------------------------

#ifdef __CXX__
int main(int argc, char **argv)
{
  int maxEntries=100;
  TString includeOnlyRange;
  if (argc>1) {
    for (int i=1; i<argc; i++) {
      if (PosOk(argv[i],"-max-events")!=-1) {
	maxEntries=atoi(argv[i+1]);
	std::cout << " ** maxEntries=" << maxEntries << "\n";
      }
      else if (PosOk(argv[i],"-include-only")!=-1) {
	includeOnlyRange=argv[i+1];
	std::cout << " ** includeOnlyRange=" << includeOnlyRange << "\n";
      }
    }
  }
  processDYee_dressed(maxEntries,includeOnlyRange);
  return 0;
}
#endif

// ------------------------------------------------------------
