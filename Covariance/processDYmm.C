#include "DYmm13TeV.h"
#include "DYbinning.h"
#include "inputs.h"

// --------------------------------------------------------------

inline
void prepareHisto(TH1D *h1) {
  h1->SetDirectory(0);
  if (!h1->GetSumw2()) h1->Sumw2();
}

// --------------------------------------------------------------

void processDYmm()
{
  DYmm13TeV_t data("/media/ssd/v20160214_1st_CovarianceMatrixInputs/Input3/ROOTFile_ntuple_CovarianceMatrixInput.root");
  data.DeactivateBranches();
  data.ActivateBranches("Momentum_postFSR_Lead  Momentum_postFSR_Sub  Momentum_preFSR_Lead  Momentum_preFSR_Sub  Weight_Norm  Weight_Gen");

  TH1D *h1postFsr_M= new TH1D("h1_postFsr_M", "postFSR;M_{gen,postFSR} [GeV];count",DYtools::nMassBins,
			    DYtools::massBinEdges);
  TH1D *h1preFsr_M= new TH1D("h1_preFsr_M", "preFSR;M_{gen,preFSR} [GeV];count",DYtools::nMassBins,
			   DYtools::massBinEdges);
  TH1D *h1postFsr_MW= new TH1D("h1_postFsr_Mweighted", "postFSR;M_{gen,postFSR} [GeV];weighted count", DYtools::nMassBins,
			    DYtools::massBinEdges);
  TH1D *h1preFsr_MW= new TH1D("h1_preFsr_Mweighted", "preFSR;M_{gen,preFSR} [GeV];weighted count", DYtools::nMassBins,
			   DYtools::massBinEdges);
  prepareHisto(h1postFsr_M);
  prepareHisto(h1preFsr_M);
  prepareHisto(h1postFsr_MW);
  prepareHisto(h1preFsr_MW);

  UInt_t nEvents= data.GetEntries();

  std::cout << "process data\n";
  for (UInt_t iEntry=0; iEntry<nEvents; iEntry++) {
    if (data.GetEntry(iEntry)<0) break;
    if (iEntry%100000==0) std::cout << "iEntry=" << iEntry << Form("(%4.2lf%%)\n",iEntry*100/double(nEvents));
    double w= data.Weight_Norm * data.Weight_Gen;
    double mPostFsr= (*data.Momentum_postFSR_Lead + *data.Momentum_postFSR_Sub).M();
    double mPreFsr= (*data.Momentum_preFSR_Lead + *data.Momentum_preFSR_Sub).M();
    h1postFsr_M->Fill( mPostFsr, 1. );
    h1postFsr_MW->Fill( mPostFsr, w );
    h1preFsr_M->Fill( mPreFsr, 1. );
    h1preFsr_MW->Fill( mPreFsr, w );
  }

  TFile fout("dymm_test.root","RECREATE");
  h1postFsr_M->Write();
  h1postFsr_MW->Write();
  h1preFsr_M->Write();
  h1preFsr_MW->Write();
  fout.Close();
  std::cout << "file <" << fout.GetName() << "> created\n";
}

