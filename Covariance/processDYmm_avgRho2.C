#include "DYmm13TeV.h"
#include "DYmm13TeV_eff.h"
#include "DYbinning.h"
#include "inputs.h"

// --------------------------------------------------------------
// --------------------------------------------------------------

void processDYmm_avgRho2(Int_t maxEntries=100)
{
  std::cout << "processDYmm_avgRho\n";
  std::cout << "DY range=" << DYtools::minMass << " .. " << DYtools::maxMass << "\n";

  //TVersion_t inpVersion=_verMu76X;

  if (DYtools::nMassBins!=DYtools::nMassBins43) {
    std::cout << "a potential DYbinning.h problem\n";
    return;
  }

  TString srcPath="/mnt/sdb/andriusj/v20160527_1st_CovarianceMatrixInputs_76X/";
  TString fnameEff= srcPath + "Input5/ROOTFile_TagProbeEfficiency_76X_v20160502.root";
  TString dataFName=srcPath + "Input3/ROOTFile_Input_CovarianceMatrix.root";
  TString avgRhoFName=srcPath + "Input6/ROOTFile_Input6_CrossCheck.root";

  TFile finEff(fnameEff);
  if (!finEff.IsOpen()) {
    std::cout << "failed to open the file <" << fnameEff << ">\n";
    return;
  }
  TH2D *h2Eff_RecoID_Data= loadHisto(finEff,"h_2D_Eff_RecoID_Data","h2Eff_RecoID_Data",1,h2dummy);
  h2Eff_RecoID_Data->SetDirectory(0);
  TH2D* h2Eff_Iso_Data= loadHisto(finEff,"h_2D_Eff_Iso_Data","h2Eff_Iso_Data",1,h2dummy);
  h2Eff_Iso_Data->SetDirectory(0);
  TH2D* h2Eff_HLT4p2_Data= loadHisto(finEff,"h_2D_Eff_HLTv4p2_Data","h2Eff_HLT4p2_Data",1,h2dummy);
  h2Eff_HLT4p2_Data->SetDirectory(0);
  TH2D* h2Eff_HLT4p3_Data= loadHisto(finEff,"h_2D_Eff_HLTv4p3_Data","h2Eff_HLT4p3_Data",1,h2dummy);
  h2Eff_HLT4p3_Data->SetDirectory(0);
  TH2D* h2Eff_RecoID_MC= loadHisto(finEff,"h_2D_Eff_RecoID_MC","h2Eff_RecoID_MC",1,h2dummy);
  h2Eff_RecoID_MC->SetDirectory(0);
  TH2D* h2Eff_Iso_MC= loadHisto(finEff,"h_2D_Eff_Iso_MC","h2Eff_Iso_MC",1,h2dummy);
  h2Eff_Iso_MC->SetDirectory(0);
  TH2D* h2Eff_HLT4p2_MC= loadHisto(finEff,"h_2D_Eff_HLTv4p2_MC","h2Eff_HLT4p2_MC",1,h2dummy);
  h2Eff_HLT4p2_MC->SetDirectory(0);
  TH2D* h2Eff_HLT4p3_MC= loadHisto(finEff,"h_2D_Eff_HLTv4p3_MC","h2Eff_HLT4p3_MC",1,h2dummy);
  h2Eff_HLT4p3_MC->SetDirectory(0);
  finEff.Close();

  if (!DYtools::compareBinning(8, h2Eff_RecoID_Data, h2Eff_Iso_Data,
			       h2Eff_HLT4p2_Data, h2Eff_HLT4p3_Data,
			       h2Eff_RecoID_MC, h2Eff_Iso_MC,
			       h2Eff_HLT4p2_MC, h2Eff_HLT4p3_MC)) {
    std::cout << "unforeseen binning problem\n";
    return;
  }

  DYTnPEff_t tnpEff;
  tnpEff.h2Eff_RecoID_Data= h2Eff_RecoID_Data;
  tnpEff.h2Eff_Iso_Data= h2Eff_Iso_Data;
  tnpEff.h2Eff_HLT4p2_Data= h2Eff_HLT4p2_Data;
  tnpEff.h2Eff_HLT4p3_Data= h2Eff_HLT4p3_Data;
  tnpEff.h2Eff_RecoID_MC= h2Eff_RecoID_MC;
  tnpEff.h2Eff_Iso_MC= h2Eff_Iso_MC;
  tnpEff.h2Eff_HLT4p2_MC= h2Eff_HLT4p2_MC;
  tnpEff.h2Eff_HLT4p3_MC= h2Eff_HLT4p3_MC;
  if (!tnpEff.updateVectors()) {
    std::cout << "tnpEff.updateVectors failed\n";
    return;
  }
  EventSpace_t postFsrES("mainES",tnpEff.h2Eff_RecoID_Data);

  TH2D *h2SF_RecoID = (TH2D*)h2Eff_RecoID_Data->Clone("h2SF_RecoID");
  h2SF_RecoID->Divide(h2Eff_RecoID_Data,h2Eff_RecoID_MC);
  h2SF_RecoID->SetTitle("h2SF_RecoID");
  TH2D *h2SF_Iso = (TH2D*)h2Eff_Iso_Data->Clone("h2SF_Iso");
  h2SF_Iso->Divide(h2Eff_Iso_Data,h2Eff_Iso_MC);
  h2SF_Iso->SetTitle("h2SF_Iso");


  TH2D* h2EffBinDef= cloneHisto(h2Eff_RecoID_Data,"h2EffBinDef","h2EffBinDef");
  if (!h2EffBinDef) return;
  const int iPtMax= h2EffBinDef->GetNbinsY();
  const double maxPt= (h2EffBinDef->GetYaxis()->GetBinLowEdge(iPtMax) +
		       h2EffBinDef->GetYaxis()->GetBinWidth(iPtMax));
  std::cout << "maxPt=" << maxPt << "\n";


  DYmm13TeV_t data(dataFName);
  if (!data.fChain || (data.fChain->GetEntries()==0)) {
    std::cout << "error opening the file\n";
    return;
  }
  data.DeactivateBranches();
  data.ActivateBranches("Momentum_Reco_Lead Momentum_Reco_Sub Momentum_postFSR_Lead  Momentum_postFSR_Sub  Weight_Norm  Weight_PU  Weight_Gen");
  data.ActivateBranches("Flag_EventSelection");

  TH1D *h1recoSel4p2_MWPU= new TH1D("h1recoSel4p2_MWPU", "RECO 4p2 selected (weighted,wPU);M_{reco} [GeV];weighted (wPU) 4.2 count",DYtools::nMassBins,DYtools::massBinEdges);
  TH1D *h1recoSel4p3_MWPU= new TH1D("h1recoSel_MWPU", "RECO 4p3 selected (weighted,wPU);M_{reco} [GeV];weighted (wPU) 4.3 count",DYtools::nMassBins,DYtools::massBinEdges);
  TH1D *h1recoSel4p2_Rho_MWPU= new TH1D("h1recoSel4p2_Rho_MWPU", "Rho_4p2 in RECO selected (weighted,wPU);M_{reco} [GeV];#rho_{v4.2} weighted (wPU) count",DYtools::nMassBins,DYtools::massBinEdges);
  TH1D *h1recoSel4p3_Rho_MWPU= new TH1D("h1recoSel4p3_Rho_MWPU", "Rho_4p3 in RECO selected (weighted,wPU);M_{reco} [GeV];#rho_{v4.3} weighted (wPU) count",DYtools::nMassBins,DYtools::massBinEdges);

  TH1D *h1postFsrInAccSel4p2_MWPU= new TH1D("h1_postFsrInAccSel4p2_MWPU", "postFSR 4p2 selected events in acceptance (weighted,wPU);M_{gen,postFSR} [GeV];weighted (wPU) 4.2 count",DYtools::nMassBins,DYtools::massBinEdges);
  TH1D *h1postFsrInAccSel4p3_MWPU= new TH1D("h1_postFsrInAccSel4p3_MWPU", "postFSR 4p3 selected events in acceptance (weighted,wPU);M_{gen,postFSR} [GeV];weighted (wPU) 4.3 count",DYtools::nMassBins,DYtools::massBinEdges);
  TH1D *h1postFsrInAccSel4p2_Rho_MWPU= new TH1D("h1_postFsrInAccSel4p2_Rho_MWPU", "Rho_4p2 in postFSR selected events in acceptance (weighted,wPU);M_{gen,postFSR} [GeV];#rho_{4p2} weighted (wPU) count",DYtools::nMassBins,DYtools::massBinEdges);
  TH1D *h1postFsrInAccSel4p3_Rho_MWPU= new TH1D("h1_postFsrInAccSel4p3_Rho_MWPU", "Rho_4p3 in postFSR selected events in acceptance (weighted,wPU);M_{gen,postFSR} [GeV];#rho_{4p3} weighted (wPU) count",DYtools::nMassBins,DYtools::massBinEdges);

  h1recoSel4p2_MWPU->SetDirectory(0); h1recoSel4p2_MWPU->Sumw2();
  h1postFsrInAccSel4p2_MWPU->SetDirectory(0); h1postFsrInAccSel4p2_MWPU->Sumw2();
  h1recoSel4p3_MWPU->SetDirectory(0); h1recoSel4p3_MWPU->Sumw2();
  h1postFsrInAccSel4p3_MWPU->SetDirectory(0); h1postFsrInAccSel4p3_MWPU->Sumw2();
  h1recoSel4p2_Rho_MWPU->SetDirectory(0); h1recoSel4p2_Rho_MWPU->Sumw2();
  h1postFsrInAccSel4p2_Rho_MWPU->SetDirectory(0); h1postFsrInAccSel4p2_Rho_MWPU->Sumw2();
  h1recoSel4p3_Rho_MWPU->SetDirectory(0); h1recoSel4p3_Rho_MWPU->Sumw2();
  h1postFsrInAccSel4p3_Rho_MWPU->SetDirectory(0); h1postFsrInAccSel4p3_Rho_MWPU->Sumw2();

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

    const int debug_print=0;

    // RECO mass binning
    if (DYtools::InAcceptance_mm(data.Momentum_Reco_Lead,data.Momentum_Reco_Sub)) {
      //std::cout << "mReco=" << mReco << "\n";

      const TLorentzVector *p1= data.Momentum_Reco_Lead;
      const TLorentzVector *p2= data.Momentum_Reco_Sub;
      double pt1tmp= (p1->Pt() >= maxPt) ? maxPt-1 : p1->Pt();
      double pt2tmp= (p2->Pt() >= maxPt) ? maxPt-1 : p2->Pt();
      int iEta1=h2EffBinDef->GetXaxis()->FindBin(p1->Eta());
      int iPt1 =h2EffBinDef->GetYaxis()->FindBin(pt1tmp);
      int iEta2=h2EffBinDef->GetXaxis()->FindBin(p2->Eta());
      int iPt2 =h2EffBinDef->GetYaxis()->FindBin(pt2tmp);
      double sf_recoID_1= h2SF_RecoID->GetBinContent(iEta1,iPt1);
      double sf_iso_1= h2SF_Iso->GetBinContent(iEta1,iPt1);
      double sf_recoID_2= h2SF_RecoID->GetBinContent(iEta2,iPt2);
      double sf_iso_2= h2SF_Iso->GetBinContent(iEta2,iPt2);

      //debug_print=(p1->Pt()>250) ? 1:0;

      double eff_HLT4p2_data1= h2Eff_HLT4p2_Data->GetBinContent(iEta1,iPt1);
      double eff_HLT4p3_data1= h2Eff_HLT4p3_Data->GetBinContent(iEta1,iPt1);
      double eff_HLT4p2_data2= h2Eff_HLT4p2_Data->GetBinContent(iEta2,iPt2);
      double eff_HLT4p3_data2= h2Eff_HLT4p3_Data->GetBinContent(iEta2,iPt2);

      double eff_HLT4p2_mc1= h2Eff_HLT4p2_MC->GetBinContent(iEta1,iPt1);
      double eff_HLT4p3_mc1= h2Eff_HLT4p3_MC->GetBinContent(iEta1,iPt1);
      double eff_HLT4p2_mc2= h2Eff_HLT4p2_MC->GetBinContent(iEta2,iPt2);
      double eff_HLT4p3_mc2= h2Eff_HLT4p3_MC->GetBinContent(iEta2,iPt2);

      double evEff_HLT4p2_data= 1-(1 - eff_HLT4p2_data1)*(1 - eff_HLT4p2_data2);
      double evEff_HLT4p2_mc  = 1-(1 - eff_HLT4p2_mc1  )*(1 - eff_HLT4p2_mc2);
      double evEff_HLT4p3_data= 1-(1 - eff_HLT4p3_data1)*(1 - eff_HLT4p3_data2);
      double evEff_HLT4p3_mc  = 1-(1 - eff_HLT4p3_mc1  )*(1 - eff_HLT4p3_mc2);
      double sf_HLT4p2= evEff_HLT4p2_data/evEff_HLT4p2_mc;
      double sf_HLT4p3= evEff_HLT4p3_data/evEff_HLT4p3_mc;

      double rhoReco4p2= sf_recoID_1*sf_recoID_2*sf_iso_1*sf_iso_2*sf_HLT4p2;
      double rhoReco4p3= sf_recoID_1*sf_recoID_2*sf_iso_1*sf_iso_2*sf_HLT4p3;

      if (debug_print) {
	std::cout << Form("pt1=%6.2lf (%d), eta1=%5.2lf (%d) ",p1->Pt(),iPt1,p1->Eta(),iEta1);
	std::cout << Form("pt2=%6.2lf (%d), eta2=%5.2lf (%d) ",p2->Pt(),iPt2,p2->Eta(),iEta2);
	std::cout << "\n";
	std::cout << " " << sf_recoID_1 << " " << sf_iso_1 << " " << eff_HLT4p2_data1/eff_HLT4p2_mc1 << "," << eff_HLT4p3_data1/eff_HLT4p2_mc1 << "\n";
	std::cout << " " << sf_recoID_2 << " " << sf_iso_2 << " " << eff_HLT4p2_data2/eff_HLT4p2_mc2 << "," << eff_HLT4p3_data2/eff_HLT4p2_mc2 << "\n";

	std::cout << "rhoReco4p2=" << rhoReco4p2 << ", 4p3=" << rhoReco4p3 << "\n";
      }

      if (rhoReco4p2 == rhoReco4p2) {
	h1recoSel4p2_MWPU->Fill( mReco, wPU );
	h1recoSel4p2_Rho_MWPU->Fill( mReco, rhoReco4p2 * wPU);
      }
      if (rhoReco4p3 == rhoReco4p3) {
	h1recoSel4p3_MWPU->Fill( mReco, wPU );
	h1recoSel4p3_Rho_MWPU->Fill( mReco, rhoReco4p3 * wPU);
      }
    }

    // post-FSR mass binning
    if (DYtools::InAcceptance_mm(data.Momentum_postFSR_Lead,data.Momentum_postFSR_Sub)) {
      //std::cout << "mPostFsr=" << mPostFsr << "\n";

      const TLorentzVector *p1= data.Momentum_postFSR_Lead;
      const TLorentzVector *p2= data.Momentum_postFSR_Sub;
      double pt1tmp= (p1->Pt() >= maxPt) ? maxPt-1 : p1->Pt();
      double pt2tmp= (p2->Pt() >= maxPt) ? maxPt-1 : p2->Pt();
      int iEta1=h2EffBinDef->GetXaxis()->FindBin(p1->Eta());
      int iPt1 =h2EffBinDef->GetYaxis()->FindBin(pt1tmp);
      int iEta2=h2EffBinDef->GetXaxis()->FindBin(p2->Eta());
      int iPt2 =h2EffBinDef->GetYaxis()->FindBin(pt2tmp);
      double sf_recoID_1= h2SF_RecoID->GetBinContent(iEta1,iPt1);
      double sf_iso_1= h2SF_Iso->GetBinContent(iEta1,iPt1);
      double sf_recoID_2= h2SF_RecoID->GetBinContent(iEta2,iPt2);
      double sf_iso_2= h2SF_Iso->GetBinContent(iEta2,iPt2);

      double eff_HLT4p2_data1= h2Eff_HLT4p2_Data->GetBinContent(iEta1,iPt1);
      double eff_HLT4p3_data1= h2Eff_HLT4p3_Data->GetBinContent(iEta1,iPt1);
      double eff_HLT4p2_data2= h2Eff_HLT4p2_Data->GetBinContent(iEta2,iPt2);
      double eff_HLT4p3_data2= h2Eff_HLT4p3_Data->GetBinContent(iEta2,iPt2);

      double eff_HLT4p2_mc1= h2Eff_HLT4p2_MC->GetBinContent(iEta1,iPt1);
      double eff_HLT4p3_mc1= h2Eff_HLT4p3_MC->GetBinContent(iEta1,iPt1);
      double eff_HLT4p2_mc2= h2Eff_HLT4p2_MC->GetBinContent(iEta2,iPt2);
      double eff_HLT4p3_mc2= h2Eff_HLT4p3_MC->GetBinContent(iEta2,iPt2);

      double evEff_HLT4p2_data= 1-(1 - eff_HLT4p2_data1)*(1 - eff_HLT4p2_data2);
      double evEff_HLT4p2_mc  = 1-(1 - eff_HLT4p2_mc1  )*(1 - eff_HLT4p2_mc2);
      double evEff_HLT4p3_data= 1-(1 - eff_HLT4p3_data1)*(1 - eff_HLT4p3_data2);
      double evEff_HLT4p3_mc  = 1-(1 - eff_HLT4p3_mc1  )*(1 - eff_HLT4p3_mc2);
      double sf_HLT4p2= evEff_HLT4p2_data/evEff_HLT4p2_mc;
      double sf_HLT4p3= evEff_HLT4p3_data/evEff_HLT4p3_mc;

      double rhoPostFsr4p2= sf_recoID_1*sf_recoID_2*sf_iso_1*sf_iso_2*sf_HLT4p2;
      double rhoPostFsr4p3= sf_recoID_1*sf_recoID_2*sf_iso_1*sf_iso_2*sf_HLT4p3;

      HERE("debug vvvvvvvv");

      if (debug_print) {
	std::cout << Form("pt1=%6.2lf (%d), eta1=%5.2lf (%d) ",p1->Pt(),iPt1,p1->Eta(),iEta1);
	std::cout << Form("pt2=%6.2lf (%d), eta2=%5.2lf (%d) ",p2->Pt(),iPt2,p2->Eta(),iEta2);
	std::cout << "\n";
	std::cout << " " << sf_recoID_1 << " " << sf_iso_1 << " " << eff_HLT4p2_data1/eff_HLT4p2_mc1 << "," << eff_HLT4p3_data1/eff_HLT4p2_mc1 << "\n";
	std::cout << " " << sf_recoID_2 << " " << sf_iso_2 << " " << eff_HLT4p2_data2/eff_HLT4p2_mc2 << "," << eff_HLT4p3_data2/eff_HLT4p2_mc2 << "\n";

	std::cout << "rhoPostFsr4p2=" << rhoPostFsr4p2 << ", 4p3=" << rhoPostFsr4p3 << "\n";
      }

      int fi1_flatIdx= DYtools::FlatIndex(h2EffBinDef, data.Momentum_postFSR_Lead,1);
      int fi2_flatIdx= DYtools::FlatIndex(h2EffBinDef, data.Momentum_postFSR_Sub,1);
      if ((fi1_flatIdx<0) || (fi2_flatIdx<0)) {
	std::cout << " bad fi_flatIdx values: " << fi1_flatIdx << "," << fi2_flatIdx << "\n";
	return;
      }
      int iMass= DYtools::massIdx(mPostFsr);
      std::cout << "iMass=" << iMass << std::endl;
      const TH2D* h2fromES= (iMass>=0) ? postFsrES.h2ES(iMass) : NULL;
      int hasERROR=0;
      if (h2fromES) {
      int fi1= postFsrES.h2ES(iMass)->GetXaxis()->FindBin(fi1_flatIdx);
      int fi2= postFsrES.h2ES(iMass)->GetYaxis()->FindBin(fi2_flatIdx);
      std::cout << "iEntry=" << iEntry << ", fi1=" << fi1 << ", fi2=" << fi2 << std::endl;
      double iniVal1= postFsrES.h2ES(iMass)->GetBinContent(fi1,fi2);
      double iniVal2= postFsrES.h2ES(iMass)->GetBinContent(fi2,fi1);
      postFsrES.fill(data.Momentum_postFSR_Lead,data.Momentum_postFSR_Sub, wPU);
      double finVal1= postFsrES.h2ES(iMass)->GetBinContent(fi1,fi2);
      double finVal2= postFsrES.h2ES(iMass)->GetBinContent(fi2,fi1);
      if ((fabs(finVal1-iniVal1-wPU)>1e-5) && (fabs(finVal2-iniVal2-wPU)>1e-5)) {
	TH1D* h1rho= NULL;// postFsrES.calculateScaleFactor(tnpEff,0,"h1rhoTest","h1rhoTest");
	std::cout << "postFsrES.fill problem\n";
	std::cout << "fi_flat = "<< fi1_flatIdx << "," << fi2_flatIdx << "\n";
	std::cout << "fi1="<<fi1<<",fi2="<<fi2<<", iniVal="<<iniVal1
		  << "/"<<iniVal2<<", finVal="<<finVal1<<"/"<<finVal2
		  << " : " << (finVal1-iniVal1) << "/" << (finVal2-iniVal2)
		  << ", wPU=" << wPU <<"\n";
	//printHisto(postFsrES.h2ES(iMass));
	const TH2D* h2tmp= postFsrES.h2ES(iMass);
	for (int ibin=1; ibin<=h2tmp->GetNbinsX(); ibin++) {
	  for (int jbin=1; jbin<=h2tmp->GetNbinsY(); jbin++) {
	    if (h2tmp->GetBinContent(ibin,jbin)) {
	      std::cout << "ibin=" << ibin << ",jbin=" << jbin
		<< "  " << h2tmp->GetXaxis()->GetBinLowEdge(ibin)
		<< "-" << (h2tmp->GetXaxis()->GetBinLowEdge(ibin)+
			   h2tmp->GetXaxis()->GetBinWidth(ibin))
		<< "  " << h2tmp->GetYaxis()->GetBinLowEdge(jbin)
		<< "-" << (h2tmp->GetYaxis()->GetBinLowEdge(jbin)+
			   h2tmp->GetYaxis()->GetBinWidth(jbin))
		<< "  " << h2tmp->GetBinContent(ibin,jbin) << " +- "
		<< h2tmp->GetBinError(ibin,jbin) << "\n";
	    }
	  }
	}
	if (h1rho) printHisto(h1rho);
	//DYtools::compareBinning(2,postFsrES.h2EffBinDef(),h2EffBinDef);
	//return;
	hasERROR=1;
      }
      }
      else std::cout << "h2fromES is null for mass " << mPostFsr << "\n";


      double rho4p2Chk= tnpEff.scaleFactor(data.Momentum_postFSR_Lead->Eta(),
					   data.Momentum_postFSR_Lead->Pt(),
					   data.Momentum_postFSR_Sub->Eta(),
					   data.Momentum_postFSR_Sub->Pt(),0);
      if (fabs(rho4p2Chk-rhoPostFsr4p2)>1e-6) {
	std::cout << "rho4p2Chk!=rhoPostFsr4p2: " << rho4p2Chk << ","
		  << rhoPostFsr4p2 << "\n";
	hasERROR=1;
	//return;
      }
      if (hasERROR) return;

      HERE("debug ^^^^^^^^^");

      if (rhoPostFsr4p2 == rhoPostFsr4p2) {
	h1postFsrInAccSel4p2_MWPU->Fill( mPostFsr, wPU );
	h1postFsrInAccSel4p2_Rho_MWPU->Fill( mPostFsr, rhoPostFsr4p2 * wPU);
      }
      if (rhoPostFsr4p3 == rhoPostFsr4p3) {
	h1postFsrInAccSel4p3_MWPU->Fill( mPostFsr, wPU );
	h1postFsrInAccSel4p3_Rho_MWPU->Fill( mPostFsr, rhoPostFsr4p3 * wPU);
      }
    }

  }

  std::cout << "\nProcessed " << nProcessed << " events, selected " << nSelected << "\n";

  TH1D *h1avgRho4p2_reco=(TH1D*) h1recoSel4p2_Rho_MWPU->Clone("h1avgRho4p2_reco");
  h1avgRho4p2_reco->SetTitle("avg rho 4p2;M_{reco} [GeV];<#rho_{4p2}>");
  h1avgRho4p2_reco->Divide(h1recoSel4p2_MWPU);
  TH1D *h1avgRho4p3_reco=(TH1D*) h1recoSel4p3_Rho_MWPU->Clone("h1avgRho4p3_reco");
  h1avgRho4p3_reco->SetTitle("avg rho 4p3;M_{reco} [GeV];<#rho_{4p3}>");
  h1avgRho4p3_reco->Divide(h1recoSel4p3_MWPU);

  TH1D *h1avgRho4p2_postFSR=(TH1D*) h1postFsrInAccSel4p2_Rho_MWPU->Clone("h1avgRho4p2_postFSR");
  h1avgRho4p2_postFSR->SetTitle("avg rho 4p2;M_{postFSR} [GeV];<#rho_{4p2}>");
  h1avgRho4p2_postFSR->Divide(h1postFsrInAccSel4p2_MWPU);
  TH1D *h1avgRho4p3_postFSR=(TH1D*) h1postFsrInAccSel4p3_Rho_MWPU->Clone("h1avgRho4p3_postFSR");
  h1avgRho4p3_postFSR->SetTitle("avg rho 4p3;M_{postFSR} [GeV];<#rho_{4p3}>");
  h1avgRho4p3_postFSR->Divide(h1postFsrInAccSel4p3_MWPU);


  TH1D *h1avgRho4p2_KL=NULL, *h1avgRho4p3_KL=NULL;

  if (avgRhoFName.Length()) {
    TFile fileAvgRho(avgRhoFName);
    if (fileAvgRho.IsOpen()) {
      TGraphAsymmErrors* g_EffSF_HLTv4p2= (TGraphAsymmErrors*)fileAvgRho.Get("g_EffSF_HLTv4p2");
      TGraphAsymmErrors* g_EffSF_HLTv4p3= (TGraphAsymmErrors*)fileAvgRho.Get("g_EffSF_HLTv4p3");
      fileAvgRho.Close();
      h1avgRho4p2_KL=convert(g_EffSF_HLTv4p2,"h1avgRho4p2_KL","KL g_EffSF_HLTv4p2",0);
      h1avgRho4p3_KL=convert(g_EffSF_HLTv4p3,"h1avgRho4p3_KL","KL g_EffSF_HLTv4p3",0);
    }
    else {
      std::cout << "failed to load from <" << avgRhoFName << ">\n";
    }
  }

  if (1) {
    h1avgRho4p2_postFSR->SetLineColor(kGreen);
    h1avgRho4p2_postFSR->SetMarkerColor(kGreen);
    h1avgRho4p3_postFSR->SetLineColor(kGreen);
    h1avgRho4p3_postFSR->SetMarkerColor(kGreen);
    plotHisto(h1avgRho4p2_reco,"cAvgRho4p2",1,0,"LPE1","AJ reco");
    plotHistoSame(h1avgRho4p2_postFSR,"cAvgRho4p2","LPE1","AJ postFSR");
    plotHisto(h1avgRho4p3_reco,"cAvgRho4p3",1,0,"LPE1","AJ reco");
    plotHistoSame(h1avgRho4p3_postFSR,"cAvgRho4p3","LPE1","AJ postFSR");

    if (h1avgRho4p2_KL) {
      h1avgRho4p2_KL->SetLineColor(kBlue);
      h1avgRho4p2_KL->SetMarkerColor(kBlue);
      h1avgRho4p3_KL->SetLineColor(kBlue);
      h1avgRho4p3_KL->SetMarkerColor(kBlue);
      plotHistoSame(h1avgRho4p2_KL,"cAvgRho4p2","LPE1","KL");
      plotHistoSame(h1avgRho4p3_KL,"cAvgRho4p3","LPE1","KL");
    }
  }

  /*
  TString fname="dymm_test_avgRho_" + versionName(inpVersion) + TString(".root");
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
  */
}

