#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLorentzVector.h>
#include <TStopwatch.h>
#include <TTimeStamp.h>
#include <TString.h>
#include <TLegend.h>
#include <THStack.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TAttMarker.h>
#include <TF1.h>
#include <TStyle.h>
#include <TEfficiency.h>
#include <stdio.h>

#include <vector>

static inline void loadBar(int x, int n, int r, int w);
void TestInputTree_signalMC()
{
	TTimeStamp ts_start;
	cout << "[Start Time(local time): " << ts_start.AsString("l") << "]" << endl;

	TStopwatch totaltime;
	totaltime.Start();

	TLorentzVector* Momentum_Reco_Lead_BeforeMomCorr = new TLorentzVector();
	TLorentzVector* Momentum_Reco_Sub_BeforeMomCorr = new TLorentzVector();
	TLorentzVector* Momentum_Reco_Lead = new TLorentzVector();
	TLorentzVector* Momentum_Reco_Sub = new TLorentzVector();


	Int_t Charge_Reco_Lead;
	Int_t Charge_Reco_Sub;
	Int_t TrackerLayers_Reco_Lead;
	Int_t TrackerLayers_Reco_Sub;
	TLorentzVector* Momentum_postFSR_Lead = new TLorentzVector();
	TLorentzVector* Momentum_postFSR_Sub = new TLorentzVector();
	TLorentzVector* Momentum_preFSR_Lead = new TLorentzVector();
	TLorentzVector* Momentum_preFSR_Sub = new TLorentzVector();
	Bool_t Flag_EventSelection;
	Double_t Weight_Norm;
	Double_t Weight_PU;
	Double_t Weight_Gen;

	// -- ntuple Setting -- //
	TChain *chain = new TChain("DYTree");
	chain->Add("/media/ssd/v20160214_1st_CovarianceMatrixInputs/Input3/ROOTFile_nutple_CovarianceMatrixInput.root");

	chain->SetBranchAddress("Momentum_Reco_Lead_BeforeMomCorr", &Momentum_Reco_Lead_BeforeMomCorr);
	chain->SetBranchAddress("Momentum_Reco_Sub_BeforeMomCorr", &Momentum_Reco_Sub_BeforeMomCorr);
	chain->SetBranchAddress("Momentum_Reco_Lead", &Momentum_Reco_Lead);
	chain->SetBranchAddress("Momentum_Reco_Sub", &Momentum_Reco_Sub);

	chain->SetBranchAddress("Charge_Reco_Lead", &Charge_Reco_Lead);
	chain->SetBranchAddress("Charge_Reco_Sub", &Charge_Reco_Sub);
	chain->SetBranchAddress("TrackerLayers_Reco_Lead", &TrackerLayers_Reco_Lead);
	chain->SetBranchAddress("TrackerLayers_Reco_Sub", &TrackerLayers_Reco_Sub);

	chain->SetBranchAddress("Momentum_postFSR_Lead", &Momentum_postFSR_Lead);
	chain->SetBranchAddress("Momentum_postFSR_Sub", &Momentum_postFSR_Sub);
	chain->SetBranchAddress("Momentum_preFSR_Lead", &Momentum_preFSR_Lead);
	chain->SetBranchAddress("Momentum_preFSR_Sub", &Momentum_preFSR_Sub);

	chain->SetBranchAddress("Flag_EventSelection", &Flag_EventSelection);
	chain->SetBranchAddress("Weight_Norm", &Weight_Norm);
	chain->SetBranchAddress("Weight_PU", &Weight_PU);
	chain->SetBranchAddress("Weight_Gen", &Weight_Gen);

	const Int_t nMassBin = 45;
	Double_t MassBinEdges[nMassBin+1] = {15, 20, 25, 30, 35, 40, 45, 50, 55, 60,
										 64, 68, 72, 76, 81, 86, 91, 96, 101, 106,
										 110, 115, 120, 126, 133, 141, 150, 160, 171, 185,
										 200, 220, 243, 273, 320, 380, 440, 510, 600, 700,
										 830, 1000, 1200, 1500, 2000, 3000};

	TFile *f_output = new TFile("ROOTFile_TestInputTree_signalMC.root", "RECREATE");
	TH1D *h_GenMass_preFSR = new TH1D("h_GenMass_preFSR", "", nMassBin, MassBinEdges);
	TH1D *h_GenMass_postFSR = new TH1D("h_GenMass_postFSR", "", nMassBin, MassBinEdges);
	TH1D *h_RecoMass = new TH1D("h_RecoMass", "", nMassBin, MassBinEdges);

	Int_t NEvents = chain->GetEntries();
	cout << "\t[Total Events: " << NEvents << "]" << endl;

	// NEvents = 1000;
	for(Int_t i=0; i<NEvents; i++)
	{
		loadBar(i+1, NEvents, 100, 100);

		// printf("[%d event]\n", i);

		chain->GetEntry(i);

		Double_t gen_M_preFSR = (*Momentum_preFSR_Lead + *Momentum_preFSR_Sub).M();
		Double_t gen_M_postFSR = (*Momentum_postFSR_Lead + *Momentum_postFSR_Sub).M();

		// printf("\t(gen_M_preFSR, gen_M_postFSR) = (%12.5lf, %12.5lf)\n", gen_M_preFSR, gen_M_postFSR);

		h_GenMass_preFSR->Fill( gen_M_preFSR, Weight_Gen * Weight_Norm );
		h_GenMass_postFSR->Fill( gen_M_postFSR, Weight_Gen * Weight_Norm );

		// printf("\t[Flag_EventSelection = %d]\n", (Int_t)Flag_EventSelection);
		if( Flag_EventSelection == kTRUE )
		{
			Double_t reco_M = (*Momentum_Reco_Lead + *Momentum_Reco_Sub).M();
			h_RecoMass->Fill( reco_M, Weight_Gen * Weight_Norm * Weight_PU );

			// printf("\t\treco_M = %12.5lf\n", reco_M);
		}
	}

	h_GenMass_preFSR->Write();
	h_GenMass_postFSR->Write();
	h_RecoMass->Write();

	Double_t TotalRunTime = totaltime.CpuTime();
	cout << "Total RunTime: " << TotalRunTime << " seconds" << endl;

	TTimeStamp ts_end;
	cout << "[End Time(local time): " << ts_end.AsString("l") << "]" << endl;
}

static inline void loadBar(int x, int n, int r, int w)
{
    // Only update r times.
    if( x == n )
    	cout << endl;

    if ( x % (n/r +1) != 0 ) return;

 
    // Calculuate the ratio of complete-to-incomplete.
    float ratio = x/(float)n;
    int   c     = ratio * w;
 
    // Show the percentage complete.
    printf("%3d%% [", (int)(ratio*100) );
 
    // Show the load bar.
    for (int x=0; x<c; x++) cout << "=";
 
    for (int x=c; x<w; x++) cout << " ";
 
    // ANSI Control codes to go back to the
    // previous line and clear it.
	cout << "]\r" << flush;

}


