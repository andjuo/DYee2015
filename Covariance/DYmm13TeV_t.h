//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Feb 17 13:06:06 2016 by ROOT version 5.34/09
// from TTree DYTree/ntuple for signal MC
// found on file: ROOTFile_nutple_CovarianceMatrixInput.root
//////////////////////////////////////////////////////////

#ifndef DYmm13TeV_t_h
#define DYmm13TeV_t_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <TLorentzVector.h>


// Fixed size dimensions of array or collections stored in the TTree if any.

class DYmm13TeV_t {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   TLorentzVector  *Momentum_Reco_Lead_BeforeMomCorr;
   TLorentzVector  *Momentum_Reco_Sub_BeforeMomCorr;
   TLorentzVector  *Momentum_Reco_Lead;
   TLorentzVector  *Momentum_Reco_Sub;
   Int_t           Charge_Reco_Lead;
   Int_t           Charge_Reco_Sub;
   Int_t           TrackerLayers_Reco_Lead;
   Int_t           TrackerLayers_Reco_Sub;
   TLorentzVector  *Momentum_postFSR_Lead;
   TLorentzVector  *Momentum_postFSR_Sub;
   TLorentzVector  *Momentum_preFSR_Lead;
   TLorentzVector  *Momentum_preFSR_Sub;
   Bool_t          Flag_EventSelection;
   Double_t        Weight_Norm;
   Double_t        Weight_PU;
   Double_t        Weight_Gen;

   // List of branches
   TBranch        *b_Momentum_Reco_Lead_BeforeMomCorr;   //!
   TBranch        *b_Momentum_Reco_Sub_BeforeMomCorr;   //!
   TBranch        *b_Momentum_Reco_Lead;   //!
   TBranch        *b_Momentum_Reco_Sub;   //!
   TBranch        *b_Charge_Reco_Lead;   //!
   TBranch        *b_Charge_Reco_Sub;   //!
   TBranch        *b_TrackerLayers_Reco_Lead;   //!
   TBranch        *b_TrackerLayers_Reco_Sub;   //!
   TBranch        *b_Momentum_postFSR_Lead;   //!
   TBranch        *b_Momentum_postFSR_Sub;   //!
   TBranch        *b_Momentum_preFSR_Lead;   //!
   TBranch        *b_Momentum_preFSR_Sub;   //!
   TBranch        *b_Flag_EventSelection;   //!
   TBranch        *b_Weight_Norm;   //!
   TBranch        *b_Weight_PU;   //!
   TBranch        *b_Weight_Gen;   //!

   DYmm13TeV_t(TTree *tree=0);
   virtual ~DYmm13TeV_t();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef DYmm13TeV_t_cxx
DYmm13TeV_t::DYmm13TeV_t(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("ROOTFile_nutple_CovarianceMatrixInput.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("ROOTFile_nutple_CovarianceMatrixInput.root");
      }
      f->GetObject("DYTree",tree);

   }
   Init(tree);
}

DYmm13TeV_t::~DYmm13TeV_t()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t DYmm13TeV_t::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t DYmm13TeV_t::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void DYmm13TeV_t::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   Momentum_Reco_Lead_BeforeMomCorr = 0;
   Momentum_Reco_Sub_BeforeMomCorr = 0;
   Momentum_Reco_Lead = 0;
   Momentum_Reco_Sub = 0;
   Momentum_postFSR_Lead = 0;
   Momentum_postFSR_Sub = 0;
   Momentum_preFSR_Lead = 0;
   Momentum_preFSR_Sub = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Momentum_Reco_Lead_BeforeMomCorr", &Momentum_Reco_Lead_BeforeMomCorr, &b_Momentum_Reco_Lead_BeforeMomCorr);
   fChain->SetBranchAddress("Momentum_Reco_Sub_BeforeMomCorr", &Momentum_Reco_Sub_BeforeMomCorr, &b_Momentum_Reco_Sub_BeforeMomCorr);
   fChain->SetBranchAddress("Momentum_Reco_Lead", &Momentum_Reco_Lead, &b_Momentum_Reco_Lead);
   fChain->SetBranchAddress("Momentum_Reco_Sub", &Momentum_Reco_Sub, &b_Momentum_Reco_Sub);
   fChain->SetBranchAddress("Charge_Reco_Lead", &Charge_Reco_Lead, &b_Charge_Reco_Lead);
   fChain->SetBranchAddress("Charge_Reco_Sub", &Charge_Reco_Sub, &b_Charge_Reco_Sub);
   fChain->SetBranchAddress("TrackerLayers_Reco_Lead", &TrackerLayers_Reco_Lead, &b_TrackerLayers_Reco_Lead);
   fChain->SetBranchAddress("TrackerLayers_Reco_Sub", &TrackerLayers_Reco_Sub, &b_TrackerLayers_Reco_Sub);
   fChain->SetBranchAddress("Momentum_postFSR_Lead", &Momentum_postFSR_Lead, &b_Momentum_postFSR_Lead);
   fChain->SetBranchAddress("Momentum_postFSR_Sub", &Momentum_postFSR_Sub, &b_Momentum_postFSR_Sub);
   fChain->SetBranchAddress("Momentum_preFSR_Lead", &Momentum_preFSR_Lead, &b_Momentum_preFSR_Lead);
   fChain->SetBranchAddress("Momentum_preFSR_Sub", &Momentum_preFSR_Sub, &b_Momentum_preFSR_Sub);
   fChain->SetBranchAddress("Flag_EventSelection", &Flag_EventSelection, &b_Flag_EventSelection);
   fChain->SetBranchAddress("Weight_Norm", &Weight_Norm, &b_Weight_Norm);
   fChain->SetBranchAddress("Weight_PU", &Weight_PU, &b_Weight_PU);
   fChain->SetBranchAddress("Weight_Gen", &Weight_Gen, &b_Weight_Gen);
   Notify();
}

Bool_t DYmm13TeV_t::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void DYmm13TeV_t::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t DYmm13TeV_t::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef DYmm13TeV_t_cxx
