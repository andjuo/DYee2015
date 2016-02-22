#include </Users/KyeongPil_Lee/Codes/DYAnalysisCodes/MyCanvas.C>
#include </Users/KyeongPil_Lee/Codes/DYAnalysisCodes/UncertantyCalcTool.h>
// #define Lumi 2595.974 // -- Up to Run260627 (Full 2015 Data), MuonPhys JSON. unit: /pb -- //
// #define Lumi 2656.777 // -- Up to Run260627 (Full 2015 Data), MuonPhys JSON. unit: /pb, Updated at 2015.12.07-- //
#define Lumi 2765.229 // -- Up to Run260627 (Full 2015 Data), MuonPhys_v2 JSON. unit: /pb, Updated at 2015.12.19 -- //
#define Lumi_HLTv4p2 848.104 // -- integrated luminosity before Run 257933 -- //
// gSystem->Load("/Users/KyeongPil_Lee/ROOT5/Unfolding/RooUnfold/libRooUnfold.so");

void ObtainYieldHistogram(TString Type, TFile *f_data, TFile *f_MC, TH1D* h_yield);

void DiffXsec_dataDrivenBkg()
{
	TFile *f_MC = (TFile*)gROOT->GetListOfFiles()->At(0);
	TFile *f_data = (TFile*)gROOT->GetListOfFiles()->At(1);
	TFile *f_bkg_dataDriven = (TFile*)gROOT->GetListOfFiles()->At(2);
	TFile *f_AccEff = (TFile*)gROOT->GetListOfFiles()->At(3);
	TFile *f_Unfold = (TFile*)gROOT->GetListOfFiles()->At(4);
	TFile *f_FSR = (TFile*)gROOT->GetListOfFiles()->At(5);
	TFile *f_theory = (TFile*)gROOT->GetListOfFiles()->At(6);

	TFile *f_output = new TFile("ROOTFile_Results_DYAnalysis.root", "RECREATE");

	vector< TString > ntupleDirectory; vector< TString > Tag; vector< Double_t > Xsec; vector< Double_t > nEvents; vector< Int_t > color;
	SetupMCsamples(&ntupleDirectory, &Tag, &Xsec, &nEvents);	

	const Int_t nMassBin = 45;
	Double_t MassBinEdges[nMassBin+1] = {15, 20, 25, 30, 35, 40, 45, 50, 55, 60,
										 64, 68, 72, 76, 81, 86, 91, 96, 101, 106,
										 110, 115, 120, 126, 133, 141, 150, 160, 171, 185,
										 200, 220, 243, 273, 320, 380, 440, 510, 600, 700,
										 830, 1000, 1200, 1500, 2000, 3000};

	///////////////////////////////////////////
 	// -- Make reco-level total signal MC -- //
	///////////////////////////////////////////
 	TH1D *h_totSignalMC = NULL;
 	f_MC->cd();
 	Int_t nTag = Tag.size();
 	for(Int_t i_tag=0; i_tag < nTag; i_tag++)
 	{
 		TH1D *h_temp = (TH1D*)f_MC->Get( "h_mass_OS_"+Tag[i_tag] )->Clone();

 		// -- Rebin using DY analysis binning -- //
 		h_temp = (TH1D*)h_temp->Rebin(nMassBin, h_temp->GetName(), MassBinEdges);

 		
 		// -- Normalized to the integrated luminosity of the data -- //
 		Double_t Norm = (Lumi * Xsec[i_tag]) / nEvents[i_tag];
 		h_temp->Scale( Norm );

 		if( Tag[i_tag].Contains("DYMuMu") ) // -- if this is signal MC -- //
 		{
 			// if( Tag[i_tag] != "DYMuMu_M50to200" ) continue;
 			// cout << "Tag = " << Tag[i_tag] << endl;

 			// -- Sum of all signal MC histogram -- //
 			if( h_totSignalMC == NULL )
 				h_totSignalMC = (TH1D*)h_temp->Clone();
 			else
 				h_totSignalMC->Add( h_temp );
 			
 		}
 	}
 	TH1D* h_totSignalMC_HLTv4p2 = (TH1D*)h_totSignalMC->Clone();
 	h_totSignalMC_HLTv4p2->Scale( Lumi_HLTv4p2 / Lumi );

 	TH1D* h_totSignalMC_HLTv4p3 = (TH1D*)h_totSignalMC->Clone();
 	h_totSignalMC_HLTv4p3->Scale( (Lumi - Lumi_HLTv4p2) / Lumi );

 	////////////////////////////////////////////////////////
 	// -- Get the gen-level MC distribution (post-FSR) -- //
 	////////////////////////////////////////////////////////
 	TH1D *h_totSignalMC_GenLevel = NULL;
 	Int_t nTag = Tag.size();
 	for(Int_t i_tag=0; i_tag < nTag; i_tag++)
 	{
 		if( Tag[i_tag].Contains("DYMuMu") )
 		{
 			TH1D *h_temp = (TH1D*)f_MC->Get( "h_GenMass_"+Tag[i_tag] )->Clone();

 			// -- Rebin using DY analysis binning -- //
 			h_temp = (TH1D*)h_temp->Rebin(nMassBin, h_temp->GetName(), MassBinEdges);
 			
 			// -- Normalized to the integrated luminosity of the data -- //
 			Double_t Norm = (Lumi * Xsec[i_tag]) / nEvents[i_tag];
 			h_temp->Scale( Norm );

 			// -- Sum of all background MC histogram -- //
 			if( h_totSignalMC_GenLevel == NULL )
 				h_totSignalMC_GenLevel = (TH1D*)h_temp->Clone();
 			else
 				h_totSignalMC_GenLevel->Add( h_temp );
 		}
 	}
 	TH1D* h_totSignalMC_GenLevel_HLTv4p2 = (TH1D*)h_totSignalMC_GenLevel->Clone();
 	h_totSignalMC_GenLevel_HLTv4p2->Scale( Lumi_HLTv4p2 / Lumi );

 	TH1D* h_totSignalMC_GenLevel_HLTv4p3 = (TH1D*)h_totSignalMC_GenLevel->Clone();
 	h_totSignalMC_GenLevel_HLTv4p3->Scale( (Lumi - Lumi_HLTv4p2) / Lumi );

	
	////////////////////////////////////////////////////////////////
	// -- Background subtraction for each HLTv4.2 and 4.3 data -- //
	////////////////////////////////////////////////////////////////
	TH1D *h_yield_HLTv4p2 = (TH1D*)f_data->Get("h_mass_OS_HLTv4p2_Data")->Clone();
	h_yield_HLTv4p2 = (TH1D*)h_yield_HLTv4p2->Rebin(nMassBin, h_yield_HLTv4p2->GetName(), MassBinEdges);
	ObtainYieldHistogram("HLTv4p2", f_data, f_MC, f_bkg_dataDriven, h_yield_HLTv4p2);

	TH1D *h_yield_HLTv4p3 = (TH1D*)f_data->Get("h_mass_OS_HLTv4p3_Data")->Clone();
	h_yield_HLTv4p3 = (TH1D*)h_yield_HLTv4p3->Rebin(nMassBin, h_yield_HLTv4p3->GetName(), MassBinEdges);
	ObtainYieldHistogram("HLTv4p3", f_data, f_MC, f_bkg_dataDriven, h_yield_HLTv4p3);

	SaveCanvas_RecoLevel_Data_vs_recoMC_HLTv4p2( h_yield_HLTv4p2, h_totSignalMC_HLTv4p2 );
	SaveCanvas_RecoLevel_Data_vs_recoMC_HLTv4p3( h_yield_HLTv4p3, h_totSignalMC_HLTv4p3 );


	TH1D* h_yield_Raw = (TH1D*)h_yield_HLTv4p2->Clone();
	h_yield_Raw->Add( h_yield_HLTv4p3 );
	SaveCanvas_Yield_BeforeEffCorr_vs_RecoMC( h_yield_Raw, h_totSignalMC );


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// -- Unfolding correction for detector resolution: apply on the data before applying Acc*Eff correction -- //
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	f_Unfold->cd();
	TH1D *h_Truth_WithinAcc = (TH1D*)f_Unfold->Get("h_Truth_RooUnfold")->Clone();
	// -- HLTv4.2 -- //
	RooUnfoldResponse *UnfoldRes_HLTv4p2 = f_Unfold->Get("h_RecoMass_h_GenMass")->Clone();
	RooUnfoldBayes *UnfoldBayes_HLTv4p2 = new RooUnfoldBayes(UnfoldRes_HLTv4p2, h_yield_HLTv4p2, 4);
	TH1D *h_yield_HLTv4p2_Unfolded = (TH1D*)UnfoldBayes_HLTv4p2->Hreco();
	h_yield_HLTv4p2_Unfolded->SetName("h_yield_HLTv4p2_Unfolded");

	SaveCanvas_preUnfold_vs_Unfold( "HLTv4p2", h_yield_HLTv4p2, h_yield_HLTv4p2_Unfolded );
	
	TH1D *h_Truth_WithinAcc_HLTv4p2 = h_Truth_WithinAcc->Clone();
	h_Truth_WithinAcc_HLTv4p2->Scale( Lumi_HLTv4p2 / Lumi );

	SaveCanvas_preUnfold_vs_genMC( "HLTv4p2", h_yield_HLTv4p2, h_Truth_WithinAcc_HLTv4p2 );
	SaveCanvas_Unfold_vs_genMC( "HLTv4p2", h_yield_HLTv4p2_Unfolded, h_Truth_WithinAcc_HLTv4p2 );

	// -- HLTv4.3 -- //
	RooUnfoldResponse *UnfoldRes_HLTv4p3 = f_Unfold->Get("h_RecoMass_h_GenMass")->Clone();
	RooUnfoldBayes *UnfoldBayes_HLTv4p3 = new RooUnfoldBayes(UnfoldRes_HLTv4p3, h_yield_HLTv4p3, 4);
	TH1D *h_yield_HLTv4p3_Unfolded = (TH1D*)UnfoldBayes_HLTv4p3->Hreco();
	h_yield_HLTv4p3_Unfolded->SetName("h_yield_HLTv4p3_Unfolded");

	SaveCanvas_preUnfold_vs_Unfold( "HLTv4p3", h_yield_HLTv4p3, h_yield_HLTv4p3_Unfolded );

	TH1D *h_Truth_WithinAcc_HLTv4p3 = h_Truth_WithinAcc->Clone();
	h_Truth_WithinAcc_HLTv4p3->Scale( (Lumi - Lumi_HLTv4p2) / Lumi );

	SaveCanvas_preUnfold_vs_genMC( "HLTv4p3", h_yield_HLTv4p3, h_Truth_WithinAcc_HLTv4p3 );
	SaveCanvas_Unfold_vs_genMC( "HLTv4p3", h_yield_HLTv4p3_Unfolded, h_Truth_WithinAcc_HLTv4p3 );

	h_yield_HLTv4p2_Unfolded->Sumw2(); h_yield_HLTv4p3_Unfolded->Sumw2();
	TH1D* h_yield_Unfolded = (TH1D*)h_yield_HLTv4p2_Unfolded->Clone();
	h_yield_Unfolded->Add( h_yield_HLTv4p3_Unfolded );


	// // -- Inversion Method -- //
	// RooUnfoldInvert *UnfoldInvert = new RooUnfoldInvert(UnfoldRes, h_yield);
	// TH1D *h_yield_Unfold_Invert = (TH1D*)UnfoldInvert->Hreco();
	// h_yield_Unfold_Invert->SetName("h_yield_Unfold_Invert");

	// SaveCanvas_UnfoldInvert_vs_genMC( h_yield_Unfold_Invert, h_Truth_WithinAcc );
	// SaveCanvas_Unfold_Bayes_vs_Invert( h_yield_Unfolded, h_yield_Unfold_Invert );


	//////////////////////////////
	// -- Acc*Eff Correction -- //
	//////////////////////////////
	f_AccEff->cd();
	TGraphAsymmErrors *g_AccEff = (TGraphAsymmErrors*)f_AccEff->Get("g_AccEff");

	// -- HLTv4p2 -- //
	TH1D* h_yield_HLTv4p2_Unfolded_AccEff = h_yield_HLTv4p2_Unfolded->Clone();
	h_yield_HLTv4p2_Unfolded_AccEff->SetName("h_yield_HLTv4p2_Unfolded_AccEff");
	Correction_AccEff(h_yield_HLTv4p2_Unfolded_AccEff, h_yield_HLTv4p2_Unfolded, g_AccEff);
	SaveCanvas_YieldAccEff_Unfold_vs_genMC( "HLTv4p2", h_yield_HLTv4p2_Unfolded_AccEff, h_totSignalMC_GenLevel_HLTv4p2 );

	// -- Acc*Eff correction on pre-unfolded data -- //
	TH1D* h_yield_HLTv4p2_AccEff = h_yield_HLTv4p2->Clone();
	h_yield_HLTv4p2_AccEff->SetName("h_yield_HLTv4p2_AccEff");
	Correction_AccEff(h_yield_HLTv4p2_AccEff, h_yield_HLTv4p2, g_AccEff);
	SaveCanvas_YieldAccEff_preUnfold_vs_genMC( "HLTv4p2", h_yield_HLTv4p2_AccEff, h_totSignalMC_GenLevel_HLTv4p2 );

	// -- HLTv4p3 -- //
	TH1D* h_yield_HLTv4p3_Unfolded_AccEff = h_yield_HLTv4p3_Unfolded->Clone();
	h_yield_HLTv4p3_Unfolded_AccEff->SetName("h_yield_HLTv4p3_Unfolded_AccEff");
	Correction_AccEff(h_yield_HLTv4p3_Unfolded_AccEff, h_yield_HLTv4p3_Unfolded, g_AccEff);
	SaveCanvas_YieldAccEff_Unfold_vs_genMC( "HLTv4p3", h_yield_HLTv4p3_Unfolded_AccEff, h_totSignalMC_GenLevel_HLTv4p3 );

	// -- Acc*Eff correction on pre-unfolded data -- //
	TH1D* h_yield_HLTv4p3_AccEff = h_yield_HLTv4p3->Clone();
	h_yield_HLTv4p3_AccEff->SetName("h_yield_HLTv4p3_AccEff");
	Correction_AccEff(h_yield_HLTv4p3_AccEff, h_yield_HLTv4p3, g_AccEff);
	SaveCanvas_YieldAccEff_preUnfold_vs_genMC( "HLTv4p3", h_yield_HLTv4p3_AccEff, h_totSignalMC_GenLevel_HLTv4p3 );

	h_yield_HLTv4p2_Unfolded_AccEff->Sumw2(); h_yield_HLTv4p3_Unfolded_AccEff->Sumw2();
	TH1D* h_yield_Unfolded_AccEff = (TH1D*)h_yield_HLTv4p2_Unfolded_AccEff->Clone();
	h_yield_Unfolded_AccEff->Add( h_yield_HLTv4p3_Unfolded_AccEff );
	
	///////////////////////////////////////////////////////
	// -- Calculation of efficiency correction factor -- //
	///////////////////////////////////////////////////////
	TGraphAsymmErrors *g_Eff = (TGraphAsymmErrors*)f_AccEff->Get("g_Eff");
	TGraphAsymmErrors *g_Eff_Corr_HLTv4p2 = (TGraphAsymmErrors*)f_AccEff->Get("g_Eff_Corr_HLTv4p2");
	TGraphAsymmErrors *g_Eff_Corr_HLTv4p3 = (TGraphAsymmErrors*)f_AccEff->Get("g_Eff_Corr_HLTv4p3");

	// -- Calculate efficiency scale factor for each mass bin: SF = Corrected Eff / Un-corrected Eff -- //
	TGraphAsymmErrors *g_EffCorr_HLTv4p2 = (TGraphAsymmErrors*)g_Eff->Clone();
	MakeRatioGraph(g_EffCorr_HLTv4p2, g_Eff_Corr_HLTv4p2, g_Eff);

	TGraphAsymmErrors *g_EffCorr_HLTv4p3 = (TGraphAsymmErrors*)g_Eff->Clone();
	MakeRatioGraph(g_EffCorr_HLTv4p3, g_Eff_Corr_HLTv4p3, g_Eff);

	SaveCanvas_EffCorr_HLTv4p2_vs_HLTv4p3( g_EffCorr_HLTv4p2, g_EffCorr_HLTv4p3 );


	/////////////////////////////////////////////////////////////
	// -- Apply Efficiency correction factors to each Yield -- //
	/////////////////////////////////////////////////////////////
	TH1D *h_yield_HLTv4p2_EffCorr = (TH1D*)h_yield_HLTv4p2_Unfolded_AccEff->Clone();
	ApplyEffCorr_Yield(h_yield_HLTv4p2_EffCorr, h_yield_HLTv4p2_Unfolded_AccEff, g_EffCorr_HLTv4p2);

	TH1D *h_yield_HLTv4p3_EffCorr = (TH1D*)h_yield_HLTv4p3_Unfolded_AccEff->Clone();
	ApplyEffCorr_Yield(h_yield_HLTv4p3_EffCorr, h_yield_HLTv4p3_Unfolded_AccEff, g_EffCorr_HLTv4p3);

	SaveCanvas_Yield_HLTv4p2_UnCorr_vs_EffCorr( h_yield_HLTv4p2_EffCorr, h_yield_HLTv4p2_Unfolded_AccEff );
	SaveCanvas_Yield_HLTv4p3_UnCorr_vs_EffCorr( h_yield_HLTv4p3_EffCorr, h_yield_HLTv4p3_Unfolded_AccEff );

	TH1D *h_totSingalMC_HLTv4p2 = (TH1D*)h_totSignalMC_GenLevel->Clone();
	h_totSingalMC_HLTv4p2->Scale( Lumi_HLTv4p2 / Lumi );
	SaveCanvas_Yield_HLTv4p2_vs_genMC( h_yield_HLTv4p2_EffCorr, h_totSingalMC_HLTv4p2 );

	TH1D *h_totSingalMC_HLTv4p3 = (TH1D*)h_totSignalMC_GenLevel->Clone();
	h_totSingalMC_HLTv4p3->Scale( (Lumi - Lumi_HLTv4p2) / Lumi );
	SaveCanvas_Yield_HLTv4p3_vs_genMC( h_yield_HLTv4p3_EffCorr, h_totSingalMC_HLTv4p3 );


	/////////////////////////////
	// -- Combine the yield -- //
	/////////////////////////////
	h_yield_HLTv4p2_EffCorr->Sumw2(); h_yield_HLTv4p3_EffCorr->Sumw2();
	TH1D *h_yield_EffCorr = (TH1D*)h_yield_HLTv4p2_EffCorr->Clone();
	h_yield_EffCorr->Add( h_yield_HLTv4p3_EffCorr );

	//////////////////////////
	// -- FSR Correction -- //
	//////////////////////////
	f_FSR->cd();
	TH1D *h_Truth_preFSR = (TH1D*)f_FSR->Get("h_mass_preFSR");
	RooUnfoldResponse *UnfoldRes_FSR = (RooUnfoldResponse*)f_FSR->Get("UnfoldRes");
	RooUnfoldBayes *UnfoldBayes_FSR = new RooUnfoldBayes(UnfoldRes_FSR, h_yield_EffCorr, 4);
	TH1D *h_FSRCorrected = UnfoldBayes_FSR->Hreco();

	SaveCanvas_beforeFSRCorr_vs_genMC_preFSR( h_yield_EffCorr, h_Truth_preFSR );
	SaveCanvas_afterFSRCorr_vs_genMC_preFSR( h_FSRCorrected, h_Truth_preFSR );
	SaveCanvas_beforeFSRCorr_vs_AfterFSRCorr( h_yield_EffCorr, h_FSRCorrected );

	////////////////////////////////////
	// -- Make x-section histogram -- //
	////////////////////////////////////
	TH1D* h_xSec_Raw = (TH1D*)h_yield_Raw->Clone();
	h_xSec_Raw->Sumw2();
	h_xSec_Raw->Scale( 1 / Lumi );

	TH1D* h_xSec_Unfolded = (TH1D*)h_yield_Unfolded->Clone();
	h_xSec_Unfolded->Sumw2();
	h_xSec_Unfolded->Scale( 1 / Lumi );

	TH1D* h_xSec_Unfolded_AccEff = (TH1D*)h_yield_Unfolded_AccEff->Clone();
	h_xSec_Unfolded_AccEff->Sumw2();
	h_xSec_Unfolded_AccEff->Scale( 1 / Lumi );

	TH1D* h_xSec_FSRCorr = (TH1D*)h_FSRCorrected->Clone();
	h_xSec_FSRCorr->Sumw2();
	h_xSec_FSRCorr->Scale( 1 / Lumi );

	TH1D* h_xSec_MC = (TH1D*)h_totSignalMC_GenLevel->Clone();
	h_xSec_MC->Sumw2();
	h_xSec_MC->Scale( 1 / Lumi );

	TH1D* h_xSec_MC_preFSR = (TH1D*)h_Truth_preFSR->Clone();
	h_xSec_MC_preFSR->Sumw2();
	h_xSec_MC_preFSR->Scale( 1 / Lumi );

	SaveCanvas_xSec_Data_vs_aMCNLO( h_xSec_FSRCorr, h_xSec_MC_preFSR );


	////////////////////////////////
	// -- X-sec/dM Calculation -- //
	////////////////////////////////
	TH1D* h_xSec_dM_Raw = (TH1D*)h_xSec_Raw->Clone();
	h_xSec_dM_Raw->Sumw2();
	Obtain_dSigma_dM(h_xSec_dM_Raw);

	TH1D* h_xSec_dM_Unfolded = (TH1D*)h_xSec_Unfolded->Clone();
	h_xSec_dM_Unfolded->Sumw2();
	Obtain_dSigma_dM(h_xSec_dM_Unfolded);

	TH1D* h_xSec_dM_Unfolded_AccEff = (TH1D*)h_xSec_Unfolded_AccEff->Clone();
	h_xSec_dM_Unfolded_AccEff->Sumw2();
	Obtain_dSigma_dM(h_xSec_dM_Unfolded_AccEff);

	TH1D* h_xSec_dM_FSRCorr = (TH1D*)h_xSec_FSRCorr->Clone();
	h_xSec_dM_FSRCorr->Sumw2();
	Obtain_dSigma_dM(h_xSec_dM_FSRCorr);

	TH1D* h_xSec_MC_dM_preFSR = (TH1D*)h_xSec_MC_preFSR->Clone();
	h_xSec_MC_dM_preFSR->Sumw2();
	Obtain_dSigma_dM(h_xSec_MC_dM_preFSR);

	f_theory->cd();
	TH1D *h_xSec_FEWZ = (TH1D*)f_theory->Get("h_NNPDF30_nlo_as_0118");
	TH1D *h_xSec_dM_FEWZ = (TH1D*)h_xSec_FEWZ->Clone();
	Obtain_dSigma_dM( h_xSec_dM_FEWZ );

	SaveCanvas_DiffXsec_Data_AllCorrStep( h_xSec_dM_FSRCorr, h_xSec_dM_Unfolded_AccEff, h_xSec_dM_Unfolded, h_xSec_dM_Raw );

	SaveCanvas_DiffXsec_Data_vs_FEWZ( h_xSec_dM_FSRCorr, h_xSec_dM_FEWZ );
	SaveCanvas_DiffXsec_Data_vs_aMCNLO( h_xSec_dM_FSRCorr, h_xSec_MC_dM_preFSR );


	//////////////////////////
	// -- Save histogram -- //
	//////////////////////////
	f_output->cd();
	// h_FSRCorrected->SetName("h_FSRCorrected");
	// h_FSRCorrected->Write();

	h_xSec_FSRCorr->SetName("h_xSec_FSRCorr");
	h_xSec_FSRCorr->Write();

	h_xSec_Unfolded->SetName("h_xSec_postFSR");
	h_xSec_Unfolded->Write();
	
	h_xSec_MC_preFSR->SetName("h_xSec_MC_preFSR");
	h_xSec_MC_preFSR->Write();



	f_output->cd();
	h_xSec_dM_FSRCorr->SetName("h_xSec_dM_FSRCorr");
	h_xSec_dM_FSRCorr->Write();

	h_xSec_dM_FEWZ->SetName("h_xSec_dM_FEWZ");
	h_xSec_dM_FEWZ->Write();

	h_xSec_MC_dM_preFSR->SetName("h_xSec_dM_aMCNLO");
	h_xSec_MC_dM_preFSR->Write();

	h_yield_HLTv4p2->SetName("h_yield_HLTv4p2");
	h_yield_HLTv4p2->Write();

	h_yield_HLTv4p3->SetName("h_yield_HLTv4p3");
	h_yield_HLTv4p3->Write();

	// h_yield->SetName("h_yield");
	// h_yield->Write();

	h_yield_Unfolded->Write();

	h_yield_Unfolded_AccEff->SetName("h_yield_Unfolded_AccEff");
	h_yield_Unfolded_AccEff->Write();

	h_FSRCorrected->SetName("h_FSRCorrected");
	h_FSRCorrected->Write();

	f_data->cd();
	TH1D *h_mass_OS_HLTv4p2_Data = (TH1D*)f_data->Get("h_mass_OS_HLTv4p2_Data")->Clone();
	h_mass_OS_HLTv4p2_Data = (TH1D*)h_mass_OS_HLTv4p2_Data->Rebin(nMassBin, h_mass_OS_HLTv4p2_Data->GetName(), MassBinEdges);
	
	TH1D *h_mass_OS_HLTv4p3_Data = (TH1D*)f_data->Get("h_mass_OS_HLTv4p3_Data")->Clone();
	h_mass_OS_HLTv4p3_Data = (TH1D*)h_mass_OS_HLTv4p3_Data->Rebin(nMassBin, h_mass_OS_HLTv4p3_Data->GetName(), MassBinEdges);

	f_output->cd();
	h_mass_OS_HLTv4p2_Data->Write();
	h_mass_OS_HLTv4p3_Data->Write();

	//////////////////////////////////
	// -- Uncertainty Estimation -- //
	//////////////////////////////////
	UncertantyCalcTool *UncTool = new UncertantyCalcTool();

	///////////////////////////////////
	// -- Statistical Uncertainty -- //
	///////////////////////////////////
	f_data->cd();
	TH1D *h_data = (TH1D*)f_data->Get("h_mass_OS_Data")->Clone();
	h_data = (TH1D*)h_data->Rebin(nMassBin, h_data->GetName(), MassBinEdges);
	UncTool->StatUnc( h_yield_Unfolded, h_data );

	/////////////////////////////////////////////////
	// -- Uncertainty from background estimation-- //
	/////////////////////////////////////////////////
	f_bkg_dataDriven->cd();
	TH1D* h_ttbar = (TH1D*)f_bkg_dataDriven->Get("ttbar")->Clone();
	TH1D* h_DYTauTau = (TH1D*)f_bkg_dataDriven->Get("DYtautau")->Clone();
	TH1D* h_tW = (TH1D*)f_bkg_dataDriven->Get("tW")->Clone();
	TH1D* h_WJets = (TH1D*)f_bkg_dataDriven->Get("wjets")->Clone();
	TH1D* h_QCD = (TH1D*)f_bkg_dataDriven->Get("dijet")->Clone();

	f_MC->cd();
	TH1D *h_ZZ = (TH1D*)f_MC->Get("h_mass_OS_ZZ")->Clone();
	h_ZZ = (TH1D*)h_ZZ->Rebin(nMassBin, h_ZZ->GetName(), MassBinEdges);
	TH1D *h_WZ = (TH1D*)f_MC->Get("h_mass_OS_WZ")->Clone();
	h_WZ = (TH1D*)h_WZ->Rebin(nMassBin, h_WZ->GetName(), MassBinEdges);
	TH1D *h_WW = (TH1D*)f_MC->Get("h_mass_OS_WW")->Clone();
	h_WW = (TH1D*)h_WW->Rebin(nMassBin, h_WW->GetName(), MassBinEdges);

	// Int_t nTag = (Int_t)Tag.size();
	for(Int_t i_tag=0; i_tag<nTag; i_tag++)
	{
		Double_t norm = ( Lumi * Xsec[i_tag] ) / nEvents[i_tag];
		if( Tag[i_tag] == "ZZ" )
			h_ZZ->Scale( norm );
		else if( Tag[i_tag] == "WZ" )
			h_WZ->Scale( norm );
		else if( Tag[i_tag] == "WW" )
			h_WW->Scale( norm );
	}

	UncTool->BkgUnc( h_yield_Unfolded, h_ttbar, h_DYTauTau, h_tW, h_WJets, h_QCD, h_ZZ, h_WZ, h_WW);

	vector< TH1D* > Histos_Unc; vector< TString > Names_Unc;
	TH1D *h_Unc_Stat = UncTool->h_RelUnc_Stat->Clone(); Histos_Unc.push_back( h_Unc_Stat ); Names_Unc.push_back( "Stat." );
	TH1D *h_Unc_Bkg = UncTool->h_RelUnc_Bkg->Clone(); Histos_Unc.push_back( h_Unc_Bkg ); Names_Unc.push_back( "Bkg.Est." );

	h_Unc_Stat->Scale( 100 );
	h_Unc_Bkg->Scale( 100 );

	MyCanvas *myc_Unc = new MyCanvas("c_Unc", "Dimuon Mass [GeV]", "Uncertainty (%)");
	myc_Unc->isLogX = kTRUE;
	myc_Unc->Legend_x1 = 0.20; myc_Unc->Legend_x2 = 0.50;
	myc_Unc->LowerEdge_X = 15; myc_Unc->UpperEdge_X = 200;
	// myc_Unc->LowerEdge_Y = 0.001; myc_Unc->UpperEdge_Y = 25;

	myc_Unc->CanvasWithMultipleHistograms( Histos_Unc, Names_Unc, "LP" );
	myc_Unc->PrintCanvas();

	//////////////////////////////////////////
	// -- Store Uncertainty distribution -- //
	//////////////////////////////////////////
	f_output->cd();
	h_Unc_Stat->Write();
	h_Unc_Bkg->Write();

	// TFile *f_DiffXsec_Theory = new TFile("ROOTFile_DiffXsec_Theory.root", "RECREATE");
	// f_DiffXsec_Theory->cd();
	// h_xSec_dM_FEWZ->SetName("h_DiffXsec_FEWZ");
	// h_xSec_dM_FEWZ->Write();
	// h_xSec_MC_dM_preFSR->SetName("h_DiffXsec_aMCNLO");
	// h_xSec_MC_dM_preFSR->Write();

}

void CompareHistogram(TH1D *h1, TH1D *h2)
{
	Int_t nBins1 = h1->GetNbinsX();
	Int_t nBins2 = h2->GetNbinsX();

	if( nBins1 != nBins2 )
	{
		printf("Different number of bins! ... (nBins1, nBins2) = (%d, %d)\n", nBins1, nBins2);
		return;
	}

	for(Int_t i=0; i<nBins1; i++)
	{
		Int_t i_bin = i+1;

		Double_t value1 = h1->GetBinContent(i_bin);
		Double_t error1 = h1->GetBinError(i_bin);
		Double_t RelError1 = error1 / value1;

		Double_t value2 = h2->GetBinContent(i_bin);
		Double_t error2 = h2->GetBinError(i_bin);
		Double_t RelError2 = error2 / value2;

		printf("[%2d bin] ", i_bin);
		printf("(value1, value2) = (%12.3lf, %12.3lf), (error1, error2) = (%12.5lf, %12.5lf), (RelError1, RelError2) = (%12.5lf, %12.5lf)\n",
			value1, value2, error1, error2, RelError1, RelError2);
	}

}

void PrintHistogram( TH1D *h )
{
	cout << "========================================================================================================================" << endl;
	TString HistName = h->GetName();
	cout << "HistName = " << HistName << endl;

	Int_t nBins = h->GetNbinsX();
	for(Int_t i=0; i<nBins; i++)
	{
		Int_t i_bin = i+1;

		Double_t value = h->GetBinContent(i_bin);
		Double_t error = h->GetBinError(i_bin);

		Double_t RelError = 0;
		if( value == 0 )
			RelError = 0;
		else
			RelError = error / value;

		printf("[%.2d bin] ", i_bin);
		printf("(value, error, RelError) = (%12.3lf, %12.5lf, %12.5lf)\n", value, error, RelError);
	}
	cout << "========================================================================================================================" << endl;
}

void Obtain_dSigma_dM(TH1D *h)
{
	Int_t nBins = h->GetNbinsX();
	for(Int_t i=0; i<nBins; i++)
	{
		Int_t i_bin = i+1;
		Double_t BinWidth = h->GetBinWidth(i_bin);

		Double_t low = h->GetBinLowEdge(i_bin);
		Double_t high = h->GetBinLowEdge(i_bin + 1);

		Double_t xSec = h->GetBinContent(i_bin);
		Double_t xSec_dM = xSec / BinWidth;

		Double_t error_before = h->GetBinError(i_bin);
		Double_t error_after = error_before / BinWidth;

		h->SetBinContent(i_bin, xSec_dM);
		h->SetBinError(i_bin, error_after);

		// printf("%2dth bin [%5.lf, %5.lf] (xSec, BinWidth, dSigma/dM) = (%15.9lf, %6.1lf, %15.9lf), (error_before, error_after) = (%8.5lf, %8.5lf)\n", 
			// i_bin, low, high, xSec, BinWidth, xSec_dM, error_before, error_after );
	}
}

void Correction_AccEff(TH1D *h_yield_AccEff, TH1D *h_yield, TGraphAsymmErrors *g_AccEff)
{
	Int_t nBins = h_yield->GetNbinsX();
	for(Int_t i=0; i<nBins; i++)
	{
		Int_t i_bin = i+1;

		Double_t x_AccEff, y_AccEff;
		g_AccEff->GetPoint(i, x_AccEff, y_AccEff);

		Double_t Yield = h_yield->GetBinContent(i_bin);
		Double_t Yield_AfterAccEff = Yield / y_AccEff;

		// -- Set the central value -- //
		h_yield_AccEff->SetBinContent( i_bin, Yield_AfterAccEff );

		// -- Calculate the error -- //
		Double_t Error = CalcError_Yield_AccEff(Yield_AfterAccEff, 
												Yield, sqrt(Yield), 
												y_AccEff, g_AccEff->GetErrorYhigh(i) );
		if( Error != 0 )
			h_yield_AccEff->SetBinError(i_bin, Error);
	}
}

Double_t CalcError_Yield_AccEff(Double_t Yield_AfterAccEff, Double_t Yield, Double_t sigma_Yield, Double_t AccEff, Double_t sigma_AccEff)
{
	if( Yield <= 0 || AccEff == 0 )
		return 0;
	Double_t Partial_Yield = ( sigma_Yield / Yield ) * ( sigma_Yield / Yield );
	Double_t Partial_AccEff = ( sigma_AccEff / AccEff ) * ( sigma_AccEff / AccEff );

	Double_t error = fabs(Yield_AfterAccEff) * sqrt( Partial_Yield + Partial_AccEff );
	return error;
}

void ObtainYieldHistogram(TString Type, TFile *f_data, TFile *f_MC, TFile *f_bkg_dataDriven, TH1D *h_yield)
{
	Double_t lumi = 0;
	if( Type == "HLTv4p2")
		lumi = Lumi_HLTv4p2;
	else if( Type == "HLTv4p3" )
		lumi = Lumi - Lumi_HLTv4p2;
	else
	{
		cout << "ERROR! Wrong Type" << endl;
		return;
	}

	vector< TString> ntupleDirectory; vector< TString > Tag; vector< Double_t > Xsec; vector< Double_t > nEvents; vector< Int_t > color;
	SetupMCsamples(&ntupleDirectory, &Tag, &Xsec, &nEvents);

	const Int_t nMassBin = 45;
	Double_t MassBinEdges[nMassBin+1] = {15, 20, 25, 30, 35, 40, 45, 50, 55, 60,
										 64, 68, 72, 76, 81, 86, 91, 96, 101, 106,
										 110, 115, 120, 126, 133, 141, 150, 160, 171, 185,
										 200, 220, 243, 273, 320, 380, 440, 510, 600, 700,
										 830, 1000, 1200, 1500, 2000, 3000};

	// -- Get Data histogram -- //									 
	f_data->cd();
	TH1D *h_data = f_data->Get("h_mass_OS_" + Type + "_Data");
	h_data = (TH1D*)h_data->Rebin(nMassBin, h_data->GetName(), MassBinEdges);

	// -- Get MC histograms -- //
	f_MC->cd();
	vector< TH1D* > h_MC;
	vector< TH1D* > h_bkgMC;
	TH1D *h_totBkgMC;
	TH1D *h_totSignalMC;
	Int_t nTag = Tag.size();
	for(Int_t i_tag=0; i_tag < nTag; i_tag++)
	{
		if( Tag[i_tag] == "ttbar" || Tag[i_tag].Contains("DYTauTau") || Tag[i_tag] == "WJets") // -- estimated by data-driven method -- //
			continue;

		TH1D *h_temp = (TH1D*)f_MC->Get( "h_mass_OS_"+Tag[i_tag] )->Clone();

		// -- Rebin using DY analysis binning -- //
		h_temp = (TH1D*)h_temp->Rebin(nMassBin, h_temp->GetName(), MassBinEdges);

		
		// -- Normalized to the integrated luminosity of the data -- //
		Double_t Norm = (lumi * Xsec[i_tag]) / nEvents[i_tag];
		h_temp->Scale( Norm );

		h_MC.push_back( h_temp );

		if( !Tag[i_tag].Contains("DYMuMu") ) // -- if this is not signal MC -- //
		{
			h_bkgMC.push_back( h_temp );


			// -- Sum of all background MC histogram -- //
			if( h_totBkgMC == NULL )
				h_totBkgMC = (TH1D*)h_temp->Clone();
			else
				h_totBkgMC->Add( h_temp );
		}
		else
		{
			// -- Sum of all background MC histogram -- //
			if( h_totSignalMC == NULL )
				h_totSignalMC = (TH1D*)h_temp->Clone();
			else
				h_totSignalMC->Add( h_temp );
		}
	}

	// -- background events from data-driven method -- //
	TH1D* h_ttbar_emu = (TH1D*)f_bkg_dataDriven->Get("ttbar")->Clone();
	h_ttbar_emu->Scale( lumi / Lumi );

	TH1D* h_tW_emu = (TH1D*)f_bkg_dataDriven->Get("tW")->Clone();
	h_tW_emu->Scale( lumi / Lumi );

	TH1D* h_DYTauTau_emu = (TH1D*)f_bkg_dataDriven->Get("DYtautau")->Clone();
	h_DYTauTau_emu->Scale( lumi / Lumi );

	TH1D* h_WJets_emu = (TH1D*)f_bkg_dataDriven->Get("wjets")->Clone();
	h_WJets_emu->Scale( lumi / Lumi );

	TH1D* h_QCD_FR = (TH1D*)f_bkg_dataDriven->Get("dijet")->Clone();
	h_QCD_FR->Scale( lumi / Lumi );
	
	h_totBkgMC->Add( h_ttbar_emu );
	h_totBkgMC->Add( h_tW_emu );
	h_totBkgMC->Add( h_DYTauTau_emu );
	h_totBkgMC->Add( h_WJets_emu );
	h_totBkgMC->Add( h_QCD_FR );

	h_totBkgMC->SetName("h_totBkgMC_"+Type);
	PrintHistogram( h_yield );
	PrintHistogram( h_totBkgMC );
	
	h_yield->Sumw2();
	h_totBkgMC->Sumw2();

	// -- Remove background MC from the data: Yield for DY events -- //
	h_yield->Add( h_totBkgMC, -1);

	for(Int_t i=0; i<h_yield->GetNbinsX(); i++)
	{
		Double_t content = h_yield->GetBinContent(i+1);

		if( content < 0 )
			h_yield->SetBinContent(i+1, 0);
	}

	PrintHistogram( h_yield );


	// // -- Comparison: Data(Before Background Subtraction) vs. Data(After Background Subtraction) -- // 
	// MyCanvas *myc1 = new MyCanvas("c_All_vs_Yield_"+Type, "Dimuon Mass [GeV]", "Number of Events");
	// myc1->isLogX = kTRUE;
	// myc1->isLogY = kTRUE;
	// myc1->isSetNoExpo_MoreLogLabels_X = kTRUE;
	// myc1->isSetNoExpo_MoreLogLabels_Y = kFALSE;
	// myc1->LowerEdge_Y = 2e-2; myc1->UpperEdge_Y = 5e6; 

	// myc1->CanvasWithHistogramsRatioPlot( (TH1D*)h_data->Clone(), (TH1D*)h_yield->Clone(),
	// 										"Data (before BkgSubtraction)", "Data (After BkgSubtraction)", "Before/After",
	// 										kBlack, kRed,
	// 										kFALSE, kFALSE,
	// 										"EP", "EPSAME" );
	// myc1->c->SaveAs("c_All_vs_Yield_"+Type+".pdf");


	// // -- Comparison: Data(After background Subtraction) vs. MC(Reco-level, signal MC)
	// MyCanvas *myc2 = new MyCanvas("c_Data_vs_MC_reco_"+Type, "Reconstructed Dimuon Mass [GeV]", "Number of events");
	// myc2->isLogX = kTRUE;
	// myc2->isLogY = kTRUE;
	// myc2->isSetNoExpo_MoreLogLabels_X = kTRUE;
	// myc2->isSetNoExpo_MoreLogLabels_Y = kFALSE;
	// myc2->LowerEdge_Y = 2e-2; myc2->UpperEdge_Y = 5e6;

	// myc2->CanvasWithHistogramsRatioPlot( (TH1D*)h_yield->Clone(), (TH1D*)h_totSignalMC->Clone(),
	// 										"Yield (Data)", "MC(DY->#mu#mu, reco-level)", "Yield/MC",
	// 										kBlack, kOrange,
	// 										kFALSE, kTRUE,
	// 										"EP", "HISTSAME" );
	// myc2->c->SaveAs("c_Data_vs_MC_reco_"+Type+".pdf");

}


void ApplyEffCorr_Yield(TH1D *h_yield_EffCorr, TH1D *h_yield, TGraphAsymmErrors* g_EffCorr)
{
	Int_t NPoints = g_EffCorr->GetN();
	Int_t nBins = h_yield->GetNbinsX();

	if( NPoints != nBins )
	{
		cout << "# Points in EffCorr != # Bins of yield histogram!" << endl;
		return;
	}

	for(Int_t i_p=0; i_p<NPoints; i_p++)
	{
		// -- Get g_EffCorr point -- //
		Double_t mass, EffCorr;
		g_EffCorr->GetPoint(i_p, mass, EffCorr);
		Double_t error_EffCorr = ReturnLargerValue( g_EffCorr->GetErrorYhigh(i_p), g_EffCorr->GetErrorYlow(i_p) );

		Int_t i_bin = i_p + 1;
		Double_t yield = h_yield->GetBinContent(i_bin);
		Double_t error_yield = h_yield->GetBinError(i_bin);

		Double_t yield_EffCorr = yield / EffCorr;
		Double_t error_yield_EffCorr = CalcError_Yield_AccEff(yield_EffCorr, yield, error_yield, EffCorr, error_EffCorr);

		// printf("[%d bin] (yield, EffCorr, yield_EffCorr) = (%.3lf, %.3lf, %.3lf)\n", i_bin, yield, EffCorr, yield_EffCorr);

		h_yield_EffCorr->SetBinContent(i_bin, yield_EffCorr);
		h_yield_EffCorr->SetBinError(i_bin, error_yield_EffCorr);
	}

}

Double_t ReturnLargerValue(Double_t a, Double_t b)
{
	if( a > b )
		return a;
	else
		return b;
}

void SetupMCsamples( vector<TString> *ntupleDirectory, vector<TString> *Tag, vector<Double_t> *Xsec, vector<Double_t> *nEvents )
{
	// -- Background Samples -- //
	ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_ZZ_25ns" ); Tag->push_back( "ZZ" ); Xsec->push_back( 15.4 ); nEvents->push_back( 996944 );
	ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_WZ_25ns" ); Tag->push_back( "WZ" ); Xsec->push_back( 66.1 ); nEvents->push_back( 978512 );
	ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_WW_25ns" ); Tag->push_back( "WW" ); Xsec->push_back( 118.7 ); nEvents->push_back( 993640 );
	ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_WJets_25ns" ); Tag->push_back( "WJets" ); Xsec->push_back( 6.15e4 ); nEvents->push_back( 16541203.0 ); //nEvents: sum of weights
	ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M10to50_25ns" ); Tag->push_back( "DYTauTau_M10to50" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 7255646.0 ); //nEvents: sum of DYTauTau weights
	ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M50toInf_25ns" ); Tag->push_back( "DYTauTau" ); Xsec->push_back( 6104/3.0 ); nEvents->push_back( 6419292.0 ); //nEvents: sum of DYTauTau weights
	ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_ttbar_25ns" ); Tag->push_back( "ttbar" ); Xsec->push_back( 831.76 ); nEvents->push_back( 19757182 );
	
	// -- Signal binned samples -- //
	ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M10to50_25ns" ); Tag->push_back( "DYMuMu_M10to50" ); Xsec->push_back( 18610.0/3.0 ); nEvents->push_back( 7293818.0 ); //nEvents: sum of weights within 10<M<50
	ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M50toInf_25ns" ); Tag->push_back( "DYMuMu_M50to200" ); Xsec->push_back( 6095.58346/3.0 ); nEvents->push_back( 6413327.0 ); //nEvents: sum of DYMuMu weights
	ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M200to400_25ns" ); Tag->push_back( "DYMuMu_M200to400" ); Xsec->push_back( 7.67/3.0 ); nEvents->push_back( 18497.0 ); //nEvents: sum of DYMuMu weights 
	ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M400to500_25ns" ); Tag->push_back( "DYMuMu_M400to500" ); Xsec->push_back( 0.423/3.0 ); nEvents->push_back( 17143.0 ); //nEvents: sum of DYMuMu weights 
	ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M500to700_25ns" ); Tag->push_back( "DYMuMu_M500to700" ); Xsec->push_back( 0.24/3.0 ); nEvents->push_back( 17397.0 ); //nEvents: sum of DYMuMu weights 
	ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M700to800_25ns" ); Tag->push_back( "DYMuMu_M700to800" ); Xsec->push_back( 0.035/3.0 ); nEvents->push_back( 15827.0 ); //nEvents: sum of DYMuMu weights 
	ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M800to1000_25ns" ); Tag->push_back( "DYMuMu_M800to1000" ); Xsec->push_back( 0.03/3.0 ); nEvents->push_back( 14742.0 ); //nEvents: sum of DYMuMu weights 
	ntupleDirectory->push_back( "Spring15DR/25ns/v20160130_MINIAODv2_DYLL_M1000to1500_25ns_Resubmit" ); Tag->push_back( "DYMuMu_M1000to1500" ); Xsec->push_back( 0.016/3.0 ); nEvents->push_back( 14381.0 ); //nEvents: sum of DYMuMu weights 
	ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M1500to2000_25ns" ); Tag->push_back( "DYMuMu_M1500to2000" ); Xsec->push_back( 0.002/3.0 ); nEvents->push_back( 13855.0 ); //nEvents: sum of DYMuMu weights 
	ntupleDirectory->push_back( "Spring15DR/25ns/v20160123_MINIAODv2_DYLL_M2000to3000_25ns" ); Tag->push_back( "DYMuMu_M2000to3000" ); Xsec->push_back( 0.00054/3.0 ); nEvents->push_back( 12376.0 ); //nEvents: sum of DYMuMu weights 
}

void MakeRatioGraph(TGraphAsymmErrors *g_ratio, TGraphAsymmErrors *g1, TGraphAsymmErrors *g2)
{
	g_ratio->Set(0); // -- Remove all points (reset) -- //

	Int_t NPoints = g1->GetN();
	for(Int_t i_p=0; i_p<NPoints; i_p++)
	{
		// -- Get g1 point -- //
		Double_t x1, y1;
		g1->GetPoint(i_p, x1, y1);
		Double_t error1 = ReturnLargerValue( g1->GetErrorYhigh(i_p), g1->GetErrorYlow(i_p) );

		// -- Get g2 point -- //
		Double_t x2, y2;
		g2->GetPoint(i_p, x2, y2);
		Double_t error2 = ReturnLargerValue( g2->GetErrorYhigh(i_p), g2->GetErrorYlow(i_p) );

		Double_t ratio;
		Double_t ratio_error;
		if(y1 != 0 && error1 != 0 && y2 != 0 && error2 != 0)
		{
			// -- calculate ratio & error -- //
			ratio = y1 / y2;
			ratio_error = Error_PropagatedAoverB(y1, error1, y2, error2);
		}
		else if( y1 != 0 && y2 != 0 && (error1 == 0 || error2 == 0) )
		{
			ratio = y1 / y2;
			ratio_error = 0;
		}
		else
		{
			ratio = 0;
			ratio_error = 0;
		}

		// -- Set Central value -- //
		g_ratio->SetPoint(i_p, x1, ratio);

		// -- Set the error -- //
		Double_t error_XLow = g1->GetErrorXlow(i_p);
		Double_t error_Xhigh = g1->GetErrorXhigh(i_p);
		g_ratio->SetPointError(i_p, error_XLow, error_Xhigh, ratio_error, ratio_error);

	}
}


void SaveCanvas_EffCorr_HLTv4p2_vs_HLTv4p3( TGraphAsymmErrors* g_EffCorr_HLTv4p2, TGraphAsymmErrors *g_EffCorr_HLTv4p3 )
{
	MyCanvas *myc_EffCorr = new MyCanvas("c_EffCorr_HLTv4p2_vs_HLTv4p3", "Dimuon Mass [GeV]", "Correction");
	myc_EffCorr->isLogX = kTRUE;
	myc_EffCorr->LowerEdge_Ratio = 0.95; myc_EffCorr->UpperEdge_Ratio = 1.05;
	myc_EffCorr->LowerEdge_Y = 0.8; myc_EffCorr->UpperEdge_Y = 1.1;
	myc_EffCorr->Legend_x1 = 0.55;

	myc_EffCorr->CanvasWithGraphRatioPlot(g_EffCorr_HLTv4p2, g_EffCorr_HLTv4p3,
									"EffCorr (HLTv4.2)", "EffCorr (HLTv4.3)", "HLTv4.2/HLTv4.3",
									kOrange+1, kGreen+1 );

	myc_EffCorr->c->SaveAs("c_EffCorr_HLTv4p2_vs_HLTv4p3.pdf");
}

void SaveCanvas_Yield_HLTv4p2_UnCorr_vs_EffCorr( TH1D* h_yield_HLTv4p2_EffCorr, TH1D* h_yield_HLTv4p2 )
{
	MyCanvas *myc_yield_EffCorr_HLTv4p2 = new MyCanvas("c_Yield_UnCorr_vs_EffCorr_HLTv4p2", "Dimuon Mass [GeV]", "Events");
	myc_yield_EffCorr_HLTv4p2->isLogX = kTRUE;
	myc_yield_EffCorr_HLTv4p2->isLogY = kTRUE;
	myc_yield_EffCorr_HLTv4p2->isSetNoExpo_MoreLogLabels_X = kTRUE;
	myc_yield_EffCorr_HLTv4p2->isSetNoExpo_MoreLogLabels_Y = kFALSE;
	myc_yield_EffCorr_HLTv4p2->LowerEdge_Y = 2e-2; myc_yield_EffCorr_HLTv4p2->UpperEdge_Y = 5e6;
	myc_yield_EffCorr_HLTv4p2->LowerEdge_Ratio = 0.95; myc_yield_EffCorr_HLTv4p2->UpperEdge_Ratio = 1.1;

	myc_yield_EffCorr_HLTv4p2->CanvasWithHistogramsRatioPlot((TH1D*)h_yield_HLTv4p2_EffCorr->Clone(), (TH1D*)h_yield_HLTv4p2->Clone(),
											"Yield(HLTv4p2, EffCorr)", "Yield(HLTv4p2, UnCorr)", "Corr/UnCorr",
											kRed, kBlack,
											kFALSE, kFALSE, 
											"EP", "EPSAME");
	myc_yield_EffCorr_HLTv4p2->c->SaveAs("c_Yield_UnCorr_vs_EffCorr_HLTv4p2.pdf");
}

void SaveCanvas_Yield_HLTv4p3_UnCorr_vs_EffCorr( TH1D* h_yield_HLTv4p3_EffCorr, TH1D* h_yield_HLTv4p3)
{
	MyCanvas *myc_yield_EffCorr_HLTv4p3 = new MyCanvas("c_Yield_UnCorr_vs_EffCorr_HLTv4p3", "Dimuon Mass [GeV]", "Events");
	myc_yield_EffCorr_HLTv4p3->isLogX = kTRUE;
	myc_yield_EffCorr_HLTv4p3->isLogY = kTRUE;
	myc_yield_EffCorr_HLTv4p3->isSetNoExpo_MoreLogLabels_X = kTRUE;
	myc_yield_EffCorr_HLTv4p3->isSetNoExpo_MoreLogLabels_Y = kFALSE;
	myc_yield_EffCorr_HLTv4p3->LowerEdge_Y = 2e-2; myc_yield_EffCorr_HLTv4p3->UpperEdge_Y = 5e6;
	myc_yield_EffCorr_HLTv4p3->LowerEdge_Ratio = 0.95; myc_yield_EffCorr_HLTv4p3->UpperEdge_Ratio = 1.1;

	myc_yield_EffCorr_HLTv4p3->CanvasWithHistogramsRatioPlot((TH1D*)h_yield_HLTv4p3_EffCorr->Clone(), (TH1D*)h_yield_HLTv4p3->Clone(),
											"Yield(HLTv4p3, EffCorr)", "Yield(HLTv4p3, UnCorr)", "Corr/UnCorr",
											kRed, kBlack,
											kFALSE, kFALSE,
											"EP", "EPSAME" );

	myc_yield_EffCorr_HLTv4p3->c->SaveAs("c_Yield_UnCorr_vs_EffCorr_HLTv4p3.pdf");
}

void SaveCanvas_Yield_HLTv4p2_vs_genMC( TH1D *h_yield_HLTv4p2_EffCorr, TH1D *h_totSingalMC_HLTv4p2  )
{
	MyCanvas *myc_Yield_HLTv4p2_vs_recoMC = new MyCanvas("c_Yield_HLTv4p2_vs_genMC", "Dimuon Mass [GeV]", "Events");
	myc_Yield_HLTv4p2_vs_recoMC->isLogX = kTRUE;
	myc_Yield_HLTv4p2_vs_recoMC->isLogY = kTRUE;
	myc_Yield_HLTv4p2_vs_recoMC->isSetNoExpo_MoreLogLabels_X = kTRUE;
	myc_Yield_HLTv4p2_vs_recoMC->isSetNoExpo_MoreLogLabels_Y = kFALSE;
	myc_Yield_HLTv4p2_vs_recoMC->LowerEdge_Y = 2e-2; myc_Yield_HLTv4p2_vs_recoMC->UpperEdge_Y = 5e6;
	myc_Yield_HLTv4p2_vs_recoMC->LowerEdge_Ratio = 0.5; myc_Yield_HLTv4p2_vs_recoMC->UpperEdge_Ratio = 1.5;

	myc_Yield_HLTv4p2_vs_recoMC->CanvasWithHistogramsRatioPlot((TH1D*)h_yield_HLTv4p2_EffCorr->Clone(), (TH1D*)h_totSingalMC_HLTv4p2->Clone(),
											"Yield(HLTv4p2, EffCorr)", "Signal MC(gen-level)", "Data/MC",
											kBlack, kOrange,
											kFALSE, kTRUE,
											"EP", "HISTSAME" );
	myc_Yield_HLTv4p2_vs_recoMC->PrintCanvas();
}
void SaveCanvas_Yield_HLTv4p3_vs_genMC( TH1D* h_yield_HLTv4p3_EffCorr, TH1D* h_totSingalMC_HLTv4p3 )
{
	MyCanvas *myc_Yield_HLTv4p3_vs_recoMC = new MyCanvas("c_Yield_HLTv4p3_vs_genMC", "Dimuon Mass [GeV]", "Events");
	myc_Yield_HLTv4p3_vs_recoMC->isLogX = kTRUE;
	myc_Yield_HLTv4p3_vs_recoMC->isLogY = kTRUE;
	myc_Yield_HLTv4p3_vs_recoMC->isSetNoExpo_MoreLogLabels_X = kTRUE;
	myc_Yield_HLTv4p3_vs_recoMC->isSetNoExpo_MoreLogLabels_Y = kFALSE;
	myc_Yield_HLTv4p3_vs_recoMC->LowerEdge_Y = 2e-2; myc_Yield_HLTv4p3_vs_recoMC->UpperEdge_Y = 5e6;
	myc_Yield_HLTv4p3_vs_recoMC->LowerEdge_Ratio = 0.5; myc_Yield_HLTv4p3_vs_recoMC->UpperEdge_Ratio = 1.5;

	myc_Yield_HLTv4p3_vs_recoMC->CanvasWithHistogramsRatioPlot((TH1D*)h_yield_HLTv4p3_EffCorr->Clone(), (TH1D*)h_totSingalMC_HLTv4p3->Clone(),
											"Yield(HLTv4p3, EffCorr)", "Signal MC(gen-level)", "Data/MC",
											kBlack, kOrange,
											kFALSE, kTRUE,
											"EP", "HISTSAME" );
	myc_Yield_HLTv4p3_vs_recoMC->PrintCanvas();
}

void SaveCanvas_Yield_BeforeEffCorr_vs_RecoMC( TH1D* h_yield_beforeEffCorr, TH1D* h_totSignalMC )
{
	MyCanvas *myc_yield_beforeEffCorr_vs_RecoMC = new MyCanvas("c_yield_beforeEffCorr_vs_RecoMC", "Dimuon Mass [GeV]", "Events");
	myc_yield_beforeEffCorr_vs_RecoMC->isLogX = kTRUE;
	myc_yield_beforeEffCorr_vs_RecoMC->isLogY = kTRUE;
	myc_yield_beforeEffCorr_vs_RecoMC->isSetNoExpo_MoreLogLabels_X = kTRUE;
	myc_yield_beforeEffCorr_vs_RecoMC->isSetNoExpo_MoreLogLabels_Y = kFALSE;
	myc_yield_beforeEffCorr_vs_RecoMC->LowerEdge_Y = 2e-2; myc_yield_beforeEffCorr_vs_RecoMC->UpperEdge_Y = 5e6;
	myc_yield_beforeEffCorr_vs_RecoMC->LowerEdge_Ratio = 0.5; myc_yield_beforeEffCorr_vs_RecoMC->UpperEdge_Ratio = 1.5;

	myc_yield_beforeEffCorr_vs_RecoMC->CanvasWithHistogramsRatioPlot((TH1D*)h_yield_beforeEffCorr->Clone(), (TH1D*)h_totSignalMC->Clone(),
											"Total Yield", "Signal MC(reco-level)", "Data/MC",
											kBlack, kOrange,
											kFALSE, kTRUE,
											"EP", "HISTSAME" );
	myc_yield_beforeEffCorr_vs_RecoMC->PrintCanvas();
}

void SaveCanvas_yield_vs_RecoMC( TH1D* h_yield, TH1D* h_totSignalMC )
{

	MyCanvas *myc_yield_vs_RecoMC = new MyCanvas("c_yield_vs_RecoMC", "Dimuon Mass [GeV]", "Events");
	myc_yield_vs_RecoMC->isLogX = kTRUE;
	myc_yield_vs_RecoMC->isLogY = kTRUE;
	myc_yield_vs_RecoMC->isSetNoExpo_MoreLogLabels_X = kTRUE;
	myc_yield_vs_RecoMC->isSetNoExpo_MoreLogLabels_Y = kFALSE;
	myc_yield_vs_RecoMC->LowerEdge_Y = 2e-2; myc_yield_vs_RecoMC->UpperEdge_Y = 5e6;
	myc_yield_vs_RecoMC->LowerEdge_Ratio = 0.5; myc_yield_vs_RecoMC->UpperEdge_Ratio = 1.5;

	myc_yield_vs_RecoMC->CanvasWithHistogramsRatioPlot((TH1D*)h_yield->Clone(), (TH1D*)h_totSignalMC->Clone(),
											"Total Yield", "Signal MC(reco-level)", "Data/MC",
											kBlack, kOrange,
											kFALSE, kTRUE,
											"EP", "HISTSAME" );
	myc_yield_vs_RecoMC->PrintCanvas();
}

void SaveCanvas_YieldAccEff_preUnfold_vs_genMC( TString Type, TH1D* h_yield_AccEff, TH1D* h_totSignalMC_GenLevel)
{
	MyCanvas *myc3 = new MyCanvas("c_YieldAccEff_preUnfold_vs_genMC_"+Type, "Dimuon Mass [GeV]", "Number of events");
	myc3->isLogX = kTRUE;
	myc3->isLogY = kTRUE;
	myc3->isSetNoExpo_MoreLogLabels_X = kTRUE;
	myc3->isSetNoExpo_MoreLogLabels_Y = kFALSE;
	myc3->LowerEdge_Y = 2e-2; myc3->UpperEdge_Y = 5e6;

	myc3->CanvasWithHistogramsRatioPlot((TH1D*)h_yield_AccEff->Clone(), (TH1D*)h_totSignalMC_GenLevel->Clone(),
											"Yield("+Type+",AccEff)", "MC(gen-level)", "Data/MC",
											kBlack, kOrange,
											kFALSE, kTRUE,
											"EP", "HISTSAME" );
	myc3->PrintCanvas();
}

void SaveCanvas_preUnfold_vs_Unfold( TString Type, TH1D* h_yield, TH1D* h_yield_Unfolded )
{
	MyCanvas *myc4 = new MyCanvas("c_preUnfold_vs_Unfold_"+Type, "Dimuon Mass [GeV]", "Number of events");
	myc4->isLogX = kTRUE;
	myc4->isLogY = kTRUE;
	myc4->Legend_x1 = 0.55;
	myc4->isSetNoExpo_MoreLogLabels_X = kTRUE;
	myc4->isSetNoExpo_MoreLogLabels_Y = kFALSE;
	myc4->LowerEdge_Y = 2e-2; myc4->UpperEdge_Y = 5e6;
	myc4->LowerEdge_Ratio = 0.5; myc4->UpperEdge_Ratio = 1.5;

	myc4->CanvasWithHistogramsRatioPlot((TH1D*)h_yield->Clone(), (TH1D*)h_yield_Unfolded->Clone(),
											"Yield("+Type+", Pre-Unfolded)", "Yield("+Type+", Unfolded)", "PreUnfolded/Unfolded",
											kBlack, kRed,
											kFALSE, kFALSE,
											"HIST", "HISTSAME" );
	myc4->PrintCanvas();
}

void SaveCanvas_Unfold_Bayes_vs_Invert( TH1D* h_yield_Bayes, TH1D* h_yield_Invert )
{
	MyCanvas *myc4 = new MyCanvas("c_Unfold_Bayes_vs_Invert", "Dimuon Mass [GeV]", "Number of events");
	myc4->isLogX = kTRUE;
	myc4->isLogY = kTRUE;
	myc4->Legend_x1 = 0.55;
	myc4->isSetNoExpo_MoreLogLabels_X = kTRUE;
	myc4->isSetNoExpo_MoreLogLabels_Y = kFALSE;
	myc4->LowerEdge_Y = 2e-2; myc4->UpperEdge_Y = 5e6;
	myc4->LowerEdge_Ratio = 0; myc4->UpperEdge_Ratio = 2;

	myc4->CanvasWithHistogramsRatioPlot((TH1D*)h_yield_Invert->Clone(), (TH1D*)h_yield_Bayes->Clone(),
											"Yield(UnfoldInvert)", "Yield(UnfoldBayes)", "Invert/Bayes",
											kBlack, kRed,
											kFALSE, kFALSE,
											"HIST", "HISTSAME" );
	myc4->PrintCanvas();
}

void SaveCanvas_preUnfold_vs_genMC( TString Type, TH1D* h_yield, TH1D* h_Truth_WithinAcc )
{
	MyCanvas *myc_BeforeAccEff_preUnfold_vs_MC = new MyCanvas("c_preUnfold_vs_genMC_"+Type, "Dimuon Mass [GeV]", "Number of events");
	myc_BeforeAccEff_preUnfold_vs_MC->isLogX = kTRUE;
	myc_BeforeAccEff_preUnfold_vs_MC->isLogY = kTRUE;
	myc_BeforeAccEff_preUnfold_vs_MC->isSetNoExpo_MoreLogLabels_X = kTRUE;
	myc_BeforeAccEff_preUnfold_vs_MC->isSetNoExpo_MoreLogLabels_Y = kFALSE;
	myc_BeforeAccEff_preUnfold_vs_MC->LowerEdge_Y = 2e-2; myc_BeforeAccEff_preUnfold_vs_MC->UpperEdge_Y = 5e6;

	myc_BeforeAccEff_preUnfold_vs_MC->CanvasWithHistogramsRatioPlot((TH1D*)h_yield->Clone(), (TH1D*)h_Truth_WithinAcc->Clone(),
											"Yield(" + Type + ")", "MC(gen-level)", "Data/MC",
											kBlack, kOrange,
											kFALSE, kTRUE,
											"EP", "HISTSAME" );
	myc_BeforeAccEff_preUnfold_vs_MC->PrintCanvas();
}

void SaveCanvas_Unfold_vs_genMC( TString Type, TH1D* h_yield_Unfolded, TH1D* h_Truth_WithinAcc )
{
	MyCanvas *myc_BeforeAccEff_Unfold_vs_MC = new MyCanvas("c_Unfold_vs_genMC_"+Type, "Dimuon Mass [GeV]", "Number of events");
	myc_BeforeAccEff_Unfold_vs_MC->isLogX = kTRUE;
	myc_BeforeAccEff_Unfold_vs_MC->isLogY = kTRUE;
	myc_BeforeAccEff_Unfold_vs_MC->isSetNoExpo_MoreLogLabels_X = kTRUE;
	myc_BeforeAccEff_Unfold_vs_MC->isSetNoExpo_MoreLogLabels_Y = kFALSE;
	myc_BeforeAccEff_Unfold_vs_MC->LowerEdge_Y = 2e-2; myc_BeforeAccEff_Unfold_vs_MC->UpperEdge_Y = 5e6;

	myc_BeforeAccEff_Unfold_vs_MC->CanvasWithHistogramsRatioPlot((TH1D*)h_yield_Unfolded->Clone(), (TH1D*)h_Truth_WithinAcc->Clone(),
											"Yield("+Type+",Unfolded)", "MC(gen-level)", "Data/MC",
											kBlack, kOrange,
											kFALSE, kTRUE,
											"EP", "HISTSAME" );
	myc_BeforeAccEff_Unfold_vs_MC->PrintCanvas();
}

void SaveCanvas_UnfoldInvert_vs_genMC( TH1D* h_yield_Unfolded, TH1D* h_Truth_WithinAcc )
{
	MyCanvas *myc_BeforeAccEff_Unfold_vs_MC = new MyCanvas("c_UnfoldInvert_vs_genMC", "Dimuon Mass [GeV]", "Number of events");
	myc_BeforeAccEff_Unfold_vs_MC->isLogX = kTRUE;
	myc_BeforeAccEff_Unfold_vs_MC->isLogY = kTRUE;
	myc_BeforeAccEff_Unfold_vs_MC->isSetNoExpo_MoreLogLabels_X = kTRUE;
	myc_BeforeAccEff_Unfold_vs_MC->isSetNoExpo_MoreLogLabels_Y = kFALSE;
	myc_BeforeAccEff_Unfold_vs_MC->LowerEdge_Y = 2e-2; myc_BeforeAccEff_Unfold_vs_MC->UpperEdge_Y = 5e6;

	myc_BeforeAccEff_Unfold_vs_MC->CanvasWithHistogramsRatioPlot((TH1D*)h_yield_Unfolded->Clone(), (TH1D*)h_Truth_WithinAcc->Clone(),
											"Yield(Data,UnfoldInvert)", "MC(gen-level)", "Data/MC",
											kBlack, kOrange,
											kFALSE, kTRUE,
											"EP", "HISTSAME" );
	myc_BeforeAccEff_Unfold_vs_MC->PrintCanvas();
}

void SaveCanvas_YieldAccEff_Unfold_vs_genMC( TString Type, TH1D* h_yield_Unfolded_AccEff, TH1D* h_totSignalMC_GenLevel )
{
	MyCanvas *myc5 = new MyCanvas("c_YieldAccEff_Unfold_vs_genMC_"+Type, "Dimuon Mass [GeV]", "Number of events");
	myc5->isLogX = kTRUE;
	myc5->isLogY = kTRUE;
	myc5->isSetNoExpo_MoreLogLabels_X = kTRUE;
	myc5->isSetNoExpo_MoreLogLabels_Y = kFALSE;
	myc5->LowerEdge_Y = 2e-2; myc5->UpperEdge_Y = 5e6;

	myc5->CanvasWithHistogramsRatioPlot((TH1D*)h_yield_Unfolded_AccEff->Clone(), (TH1D*)h_totSignalMC_GenLevel->Clone(),
											"Yield("+Type+",Unfolded,AccEff)", "MC(gen-level)", "Data/MC",
											kBlack, kOrange,
											kFALSE, kTRUE,
											"EP", "HISTSAME" );
	myc5->PrintCanvas();
}

void SaveCanvas_beforeFSRCorr_vs_AfterFSRCorr( TH1D* h_yield_Unfolded_AccEff, TH1D* h_FSRCorrected )
{
	MyCanvas *myc6 = new MyCanvas("c_beforeFSRCorr_vs_AfterFSRCorr", "Dimuon Mass [GeV]", "Number of events");
	myc6->isLogX = kTRUE;
	myc6->isLogY = kTRUE;
	myc6->Legend_x1 = 0.55;
	myc6->isSetNoExpo_MoreLogLabels_X = kTRUE;
	myc6->isSetNoExpo_MoreLogLabels_Y = kFALSE;
	myc6->LowerEdge_Y = 2e-2; myc6->UpperEdge_Y = 5e6;

	myc6->CanvasWithHistogramsRatioPlot((TH1D*)h_yield_Unfolded_AccEff->Clone(), (TH1D*)h_FSRCorrected->Clone(),
											"Yield(Data,Unfolded,AccEff)", "Yield(Data,Unfolded,AccEff,FSRCorr)", "Before/After",
											kBlack, kRed,
											kFALSE, kFALSE,
											"HIST", "HISTSAME" );
	myc6->PrintCanvas();
}

void SaveCanvas_beforeFSRCorr_vs_genMC_preFSR( TH1D* h_yield_Unfolded_AccEff, TH1D* h_Truth_preFSR )
{
	MyCanvas *myc7 = new MyCanvas("c_beforeFSRCorr_vs_genMC_preFSR", "Dimuon Mass [GeV]", "Number of events");
	myc7->isLogX = kTRUE;
	myc7->isLogY = kTRUE;
	myc7->isSetNoExpo_MoreLogLabels_X = kTRUE;
	myc7->isSetNoExpo_MoreLogLabels_Y = kFALSE;
	myc7->LowerEdge_Y = 2e-2; myc7->UpperEdge_Y = 5e6;

	myc7->CanvasWithHistogramsRatioPlot((TH1D*)h_yield_Unfolded_AccEff->Clone(), (TH1D*)h_Truth_preFSR->Clone(),
											"Yield(Data,Unfolded,AccEff)", "MC(gen-level, preFSR)", "Data/MC",
											kBlack, kOrange,
											kFALSE, kTRUE,
											"EP", "HISTSAME" );

	myc7->PrintCanvas();
}

void SaveCanvas_afterFSRCorr_vs_genMC_preFSR( TH1D* h_FSRCorrected, TH1D* h_Truth_preFSR )
{
	MyCanvas *myc8 = new MyCanvas("c_afterFSRCorr_vs_genMC_preFSR", "Dimuon Mass [GeV]", "Number of events");
	myc8->isLogX = kTRUE;
	myc8->isSetNoExpo_MoreLogLabels_X = kTRUE;
	myc8->isLogY = kTRUE;
	myc8->isSetNoExpo_MoreLogLabels_Y = kFALSE;
	myc8->LowerEdge_Y = 2e-2; myc8->UpperEdge_Y = 5e6;
	
	myc8->CanvasWithHistogramsRatioPlot((TH1D*)h_FSRCorrected->Clone(), (TH1D*)h_Truth_preFSR->Clone(),
											"Yield(Data, Unfolded, AfterFSRCorr)", "MC(gen-level, preFSR)", "Data/MC",
											kBlack, kOrange,
											kFALSE, kTRUE,
											"EP", "HISTSAME" );

	myc8->PrintCanvas();
}

void SaveCanvas_xSec_Data_vs_aMCNLO( TH1D* h_xSec_FSRCorr, TH1D* h_xSec_MC_preFSR )
{
	MyCanvas *myc9 = new MyCanvas("c_xSec_Data_vs_aMCNLO", "Dimuon Mass [GeV]", "Cross Section [pb]");
	myc9->isLogX = kTRUE;
	myc9->isSetNoExpo_MoreLogLabels_X = kTRUE;
	myc9->isLogY = kTRUE;
	myc9->isSetNoExpo_MoreLogLabels_Y = kFALSE;


	// -- X-section comparison between data and MC -- //
	myc9->CanvasWithHistogramsRatioPlot( h_xSec_FSRCorr, (TH1D*)h_xSec_MC_preFSR->Clone(),
											"Data (Unfold+FSRCorr)", "aMC@NLO (pre-FSR)", "Data/aMC@NLO",
											kBlack, kRed,
											kFALSE, kFALSE,
											"EP", "EPSAME" );

	myc9->PrintCanvas();
}

void SaveCanvas_DiffXsec_Data_AllCorrStep( TH1D* h_xSec_dM_FSRCorr, TH1D* h_xSec_dM_Unfolded_AccEff, TH1D* h_xSec_dM_Unfolded, TH1D* h_xSec_dM_Raw )
{
	MyCanvas *myc10 = new MyCanvas("c_DiffXsec_Data_AllCorrStep", "Dimuon Mass [GeV]", "d#sigma/dM [pb/GeV]");
	myc10->isLogX = kTRUE;
	myc10->isSetNoExpo_MoreLogLabels_X = kTRUE;
	myc10->isLogY = kTRUE;
	myc10->isSetNoExpo_MoreLogLabels_Y = kFALSE;

	// vector< TH1D* > Histo_xSec_dM; vector< TString > Names_xSec_dM;
	// Histo_xSec_dM.push_back( h_xSec_dM_FSRCorr ); Names_xSec_dM.push_back( "Unfold + FSRCorr" );
	// Histo_xSec_dM.push_back( h_xSec_dM_Unfolded ); Names_xSec_dM.push_back( "Unfold" );
	// Histo_xSec_dM.push_back( h_xSec_dM_Raw ); Names_xSec_dM.push_back( "Raw" );

	// myc10->CanvasWithMultipleHistograms( Histo_xSec_dM, Names_xSec_dM );

	myc10->CanvasWithHistogramsRatioPlot( (TH1D*)h_xSec_dM_FSRCorr->Clone(), (TH1D*)h_xSec_dM_Raw->Clone(),
											"Data (FSR)", "Data (Raw)", "FSR/Raw",
											kBlack, kRed,
											kFALSE, kFALSE,
											"EP", "EPSAME" );


	h_xSec_dM_FSRCorr->SetMarkerSize(1);
	h_xSec_dM_FSRCorr->SetMarkerStyle(20);
	h_xSec_dM_FSRCorr->SetLineColor(kBlack);
	h_xSec_dM_FSRCorr->SetMarkerColor(kBlack);
	h_xSec_dM_FSRCorr->SetFillColorAlpha(kWhite, 0);

	h_xSec_dM_Unfolded_AccEff->SetMarkerSize(1);
	h_xSec_dM_Unfolded_AccEff->SetMarkerStyle(20);
	h_xSec_dM_Unfolded_AccEff->SetLineColor(kBlue);
	h_xSec_dM_Unfolded_AccEff->SetMarkerColor(kBlue);
	h_xSec_dM_Unfolded_AccEff->SetFillColorAlpha(kWhite, 0);

	h_xSec_dM_Unfolded->SetMarkerSize(1);
	h_xSec_dM_Unfolded->SetMarkerStyle(20);
	h_xSec_dM_Unfolded->SetLineColor(kGreen+1);
	h_xSec_dM_Unfolded->SetMarkerColor(kGreen+1);
	h_xSec_dM_Unfolded->SetFillColorAlpha(kWhite, 0);

	h_xSec_dM_Raw->SetMarkerSize(1);
	h_xSec_dM_Raw->SetMarkerStyle(20);
	h_xSec_dM_Raw->SetLineColor(kRed);
	h_xSec_dM_Raw->SetMarkerColor(kRed);
	h_xSec_dM_Raw->SetFillColorAlpha(kWhite, 0);


	myc10->legend->Clear();
	TLegend *legend2 = new TLegend(0.55, 0.75, 0.95, 0.95);
	legend2->SetBorderSize(0);
	legend2->SetFillStyle(0);
	legend2->AddEntry(h_xSec_dM_FSRCorr, "Data (FSR)");
	legend2->AddEntry(h_xSec_dM_Unfolded_AccEff, "Data (Acc*Eff)");
	legend2->AddEntry(h_xSec_dM_Unfolded, "Data (Unfolding)");
	legend2->AddEntry(h_xSec_dM_Raw, "Data (Raw)");


	myc10->TopPad->cd();
	h_xSec_dM_Unfolded->Draw("EPSAME");
	h_xSec_dM_Unfolded_AccEff->Draw("EPSAME");
	h_xSec_dM_FSRCorr->Draw("EPSAME");
	legend2->Draw();
	gPad->Update();

	myc10->PrintCanvas();
}

void SaveCanvas_DiffXsec_Data_vs_FEWZ( TH1D* h_xSec_dM_FSRCorr, TH1D* h_xSec_dM_FEWZ )
{
	MyCanvas *myc_xSec_dM_data_FEWZ = new MyCanvas("c_xSec_dM_data_FEWZ", "Dimuon Mass [GeV]", "d#sigma/dM [pb/GeV]");
	myc_xSec_dM_data_FEWZ->isLogX = kTRUE;
	myc_xSec_dM_data_FEWZ->isSetNoExpo_MoreLogLabels_X = kTRUE;
	myc_xSec_dM_data_FEWZ->isLogY = kTRUE;
	myc_xSec_dM_data_FEWZ->isSetNoExpo_MoreLogLabels_Y = kFALSE;
	myc_xSec_dM_data_FEWZ->Legend_x1 = 0.55;
	myc_xSec_dM_data_FEWZ->Legend_y1 = 0.75;
	myc_xSec_dM_data_FEWZ->Legend_x2 = 0.95;
	myc_xSec_dM_data_FEWZ->Legend_y2 = 0.95;

	myc_xSec_dM_data_FEWZ->CanvasWithHistogramsRatioPlot( (TH1D*)h_xSec_dM_FSRCorr->Clone(), (TH1D*)h_xSec_dM_FEWZ->Clone(),
											"Data", "FEWZ (NLO)", "Data/FEWZ",
											kBlack, kRed,
											kFALSE, kFALSE,
											"EP", "EPSAME" );


	myc_xSec_dM_data_FEWZ->PrintCanvas();
}

void SaveCanvas_DiffXsec_Data_vs_aMCNLO( TH1D* h_xSec_dM_FSRCorr, TH1D* h_xSec_MC_dM_preFSR )
{
	MyCanvas *myc_xSec_dM_data_aMCNLO = new MyCanvas("c_DiffXsec_Data_vs_aMCNLO", "Dimuon Mass [GeV]", "d#sigma/dM [pb/GeV]");
	myc_xSec_dM_data_aMCNLO->isLogX = kTRUE;
	myc_xSec_dM_data_aMCNLO->isSetNoExpo_MoreLogLabels_X = kTRUE;
	myc_xSec_dM_data_aMCNLO->isLogY = kTRUE;
	myc_xSec_dM_data_aMCNLO->isSetNoExpo_MoreLogLabels_Y = kFALSE;
	myc_xSec_dM_data_aMCNLO->Legend_x1 = 0.55;
	myc_xSec_dM_data_aMCNLO->Legend_y1 = 0.75;
	myc_xSec_dM_data_aMCNLO->Legend_x2 = 0.95;
	myc_xSec_dM_data_aMCNLO->Legend_y2 = 0.95;

	myc_xSec_dM_data_aMCNLO->CanvasWithHistogramsRatioPlot( (TH1D*)h_xSec_dM_FSRCorr->Clone(), (TH1D*)h_xSec_MC_dM_preFSR->Clone(),
											"Data", "aMC@NLO", "Data/aMC@NLO",
											kBlack, kRed,
											kFALSE, kFALSE,
											"EP", "EPSAME" );

	myc_xSec_dM_data_aMCNLO->PrintCanvas();
}

void SaveCanvas_RecoLevel_Data_vs_recoMC_HLTv4p2( TH1D* h_yield_HLTv4p2, TH1D* h_totSignalMC_GenLevel_HLTv4p2 )
{
	MyCanvas *myc = new MyCanvas("c_RecoLevel_Data_vs_recoMC_HLTv4p2", "Dimuon Mass [GeV]", "Events");
	myc->isLogX = kTRUE;
	myc->isLogY = kTRUE;
	myc->isSetNoExpo_MoreLogLabels_X = kTRUE;
	myc->isSetNoExpo_MoreLogLabels_Y = kFALSE;
	myc->LowerEdge_Y = 2e-2; myc->UpperEdge_Y = 5e6;
	myc->LowerEdge_Ratio = 0.5; myc->UpperEdge_Ratio = 1.5;

	myc->CanvasWithHistogramsRatioPlot((TH1D*)h_yield_HLTv4p2->Clone(), (TH1D*)h_totSignalMC_GenLevel_HLTv4p2->Clone(),
											"Data (HLTv4.2)", "Signal MC(reco-level)", "Data/MC",
											kBlack, kOrange,
											kFALSE, kTRUE,
											"EP", "HISTSAME" );
	myc->PrintCanvas();
}

void SaveCanvas_RecoLevel_Data_vs_recoMC_HLTv4p3( TH1D* h_yield_HLTv4p3, TH1D* h_totSignalMC_GenLevel_HLTv4p3 )
{
	MyCanvas *myc = new MyCanvas("c_RecoLevel_Data_vs_recoMC_HLTv4p3", "Dimuon Mass [GeV]", "Events");
	myc->isLogX = kTRUE;
	myc->isLogY = kTRUE;
	myc->isSetNoExpo_MoreLogLabels_X = kTRUE;
	myc->isSetNoExpo_MoreLogLabels_Y = kFALSE;
	myc->LowerEdge_Y = 2e-2; myc->UpperEdge_Y = 5e6;
	myc->LowerEdge_Ratio = 0.5; myc->UpperEdge_Ratio = 1.5;

	myc->CanvasWithHistogramsRatioPlot((TH1D*)h_yield_HLTv4p3->Clone(), (TH1D*)h_totSignalMC_GenLevel_HLTv4p3->Clone(),
											"Data (HLTv4.3)", "Signal MC(reco-level)", "Data/MC",
											kBlack, kOrange,
											kFALSE, kTRUE,
											"EP", "HISTSAME" );
	myc->PrintCanvas();
}

Double_t Error_PropagatedAoverB(Double_t A, Double_t sigma_A, Double_t B, Double_t sigma_B)
{
	Double_t ratio_A = (sigma_A) / A;
	Double_t ratio_B = (sigma_B) / B;

	Double_t errorSquare = ratio_A * ratio_A + ratio_B * ratio_B;

	return (A/B) * sqrt(errorSquare);
}

Double_t ReturnLargerValue(Double_t a, Double_t b)
{
	if( a > b )
		return a;
	else
		return b;
}

