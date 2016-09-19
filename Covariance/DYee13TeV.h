//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Feb 17 13:06:06 2016 by ROOT version 5.34/09
// from TTree DYTree/ntuple for signal MC
// found on file: ROOTFile_nutple_CovarianceMatrixInput.root
//////////////////////////////////////////////////////////

#ifndef DYee13TeV_t_h
#define DYee13TeV_t_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TObjArray.h>

// Header file for the classes stored in the TTree if any.
#include <TLorentzVector.h>
#include <iostream>
#include <sstream>
#include <stdarg.h>
#include <TString.h>
#include <TH1D.h>
#include <vector>
#include "DYbinning.h"

// ------------------------------------------------------------------

// Fixed size dimensions of array or collections stored in the TTree if any.

class DYee13TeV_t {
public :
  //TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  TChain           *fChain;
  TTree            *fOutTree; // output tree
  TFile            *fOutFile; // output file
  Int_t           fCurrent; //!current Tree number in a TChain
  Int_t           fCurrentOld;
  TH1D            *fH1PreFSR_binned; // Weight = Gen x Norm
  std::vector<TH1D*> fH1PreFSR_WG, fH1PreFSR_WGN; // pre-FSR mass histograms

   // Declaration of leaf types
  Double_t RunNo;
  Double_t EvtNo;
  TLorentzVector  *Momentum_Reco_Lead_BeforeEnCorr;
  TLorentzVector  *Momentum_Reco_Sub_BeforeEnCorr;
  TLorentzVector  *Momentum_Reco_Lead;
  TLorentzVector  *Momentum_Reco_Sub;
  Double_t        SCEta_Lead;
  Double_t        SCEta_Sub;
  Int_t           Charge_Reco_Lead;
  Int_t           Charge_Reco_Sub;
  TLorentzVector  *Momentum_postFSR_Lead;
  TLorentzVector  *Momentum_postFSR_Sub;
  TLorentzVector  *Momentum_preFSR_Lead;
  TLorentzVector  *Momentum_preFSR_Sub;
  Bool_t          Flag_EventSelection;
  Bool_t          Flag_EventSelectionExceptSCGap;
  Bool_t          Flag_RecoEleSelection;
  Double_t        Weight_Norm;
  Double_t        Weight_PU;
  Double_t        Weight_Gen;

  // List of branches
  TBranch   *b_RunNo;
  TBranch   *b_EvtNo;
  TBranch        *b_Momentum_Reco_Lead_BeforeEnCorr;   //!
  TBranch        *b_Momentum_Reco_Sub_BeforeEnCorr;   //!
  TBranch        *b_Momentum_Reco_Lead;   //!
  TBranch        *b_Momentum_Reco_Sub;   //!
  TBranch        *b_SCEta_Lead;
  TBranch        *b_SCEta_Sub;
  TBranch        *b_Charge_Reco_Lead;   //!
  TBranch        *b_Charge_Reco_Sub;   //!
  TBranch        *b_Momentum_postFSR_Lead;   //!
  TBranch        *b_Momentum_postFSR_Sub;   //!
  TBranch        *b_Momentum_preFSR_Lead;   //!
  TBranch        *b_Momentum_preFSR_Sub;   //!
  TBranch        *b_Flag_EventSelection;   //!
  TBranch        *b_Flag_EventSelectionExceptSCGap;   //!
  TBranch        *b_Flag_RecoEleSelection;
  TBranch        *b_Weight_Norm;   //!
  TBranch        *b_Weight_PU;   //!
  TBranch        *b_Weight_Gen;   //!

  DYee13TeV_t(TString fname="", TString treeName="DYTree");
  virtual ~DYee13TeV_t();
  virtual Int_t    Cut(Long64_t entry);
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual int      Init(TString fname);
  virtual void     Loop();
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);

  // Methods for file creation
  void Zero();
  int CreateNew(TString fname, TString treeName);
  void Fill() { fOutTree->Fill(); }

  // Added methods
  UInt_t GetEntries() { return fChain->GetEntries(); }
  void DeactivateBranches();
  void ActivateBranches(TString brNames);
  int BranchExists(TString brName);

  friend
    std::ostream& operator<<(std::ostream &out, DYee13TeV_t &obj);
};

#ifndef Inputs_H
std::ostream& operator<<(std::ostream &out, const TLorentzVector &v);
std::ostream& operator<<(std::ostream &out, const TLorentzVector *v);
#endif

#endif

#ifdef DYee13TeV_t_cxx
DYee13TeV_t::DYee13TeV_t(TString fname, TString treeName) :
  fChain(new TChain(treeName)),
  fOutTree(NULL), fOutFile(NULL),
    fCurrent(-1), fCurrentOld(-1),
    fH1PreFSR_binned(NULL),
    fH1PreFSR_WG(), fH1PreFSR_WGN()
{
  if (fname.Length()==0) return;
  int ok=1;
  if (fname.Index("<new>")!=-1) {
    if (!this->CreateNew(fname,treeName)) {
      ok=0;
      std::cout << "CreateNew failed in constructor" << std::endl;
    }
  }
  else if (!this->Init(fname)) {
    ok=0;
    std::cout << "Initialization failed in constructor" << std::endl;
  }
  if (ok) {
    fH1PreFSR_binned=new TH1D("h1preFSR_DYee13TeV_binned",
	 "preFSR DYee13TeV from file;M_{preFSR} [GeV];sum(W_{norm}*W_{gen})",
			      DYtools::nMassBins,DYtools::massBinEdges);
    fH1PreFSR_binned->SetDirectory(0);
    fH1PreFSR_binned->Sumw2();
  }
}

DYee13TeV_t::~DYee13TeV_t()
{
  if (fOutTree) fOutTree->Write();
  if (fOutFile && fOutFile->IsOpen()) fOutFile->Close();
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t DYee13TeV_t::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   Int_t res= fChain->GetEntry(entry);
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   if (res && fH1PreFSR_WG.size()) {
     double m=( (*Momentum_preFSR_Lead) + (*Momentum_preFSR_Sub) ).M();
     fH1PreFSR_WG.back()->Fill(m, Weight_Gen);
     fH1PreFSR_WGN.back()->Fill(m, Weight_Gen * Weight_Norm);
     if (fH1PreFSR_binned) fH1PreFSR_binned->Fill(m, Weight_Gen * Weight_Norm);
   }
   return res;
}

Long64_t DYee13TeV_t::LoadTree(Long64_t entry)
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

int DYee13TeV_t::Init(TString fname)
{
  if (fname.Length()==0) {
    std::cout << "DYee13TeV_t::Init non-empty fname is expected\n";
    return 0;
  }
  if (fname.Index(" ")==-1) {
    // only one string defining files
    fChain->Add(fname);
  }
  else {
    std::stringstream ss(fname.Data());
    TString fn;
    while (!ss.eof()) {
      ss >> fn;
      if (fn.Length()) {
	std::cout << "adding file <" << fn << ">\n";
	if (fn=="<new>") {
	  std::cout << "keyword 'new' detected. Should call CreateNew(fname)\n";
	  return 0;
	}
	fChain->Add(fn);
      }
    }
  }

  // Clear fields
  RunNo=0;
  EvtNo=0;
  Charge_Reco_Lead=0;
  Charge_Reco_Sub=0;
  SCEta_Lead=0;
  SCEta_Sub=0;
  Flag_EventSelection=false;
  Flag_EventSelectionExceptSCGap=false;
  Flag_RecoEleSelection=false;
  Weight_Norm=0;
  Weight_PU=0;
  Weight_Gen=0;

   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   Momentum_Reco_Lead_BeforeEnCorr = 0;
   Momentum_Reco_Sub_BeforeEnCorr = 0;
   Momentum_Reco_Lead = 0;
   Momentum_Reco_Sub = 0;
   Momentum_postFSR_Lead = 0;
   Momentum_postFSR_Sub = 0;
   Momentum_preFSR_Lead = 0;
   Momentum_preFSR_Sub = 0;
   //// Set branch addresses and branch pointers
   //if (!tree) return;
   //fChain = tree;
   fCurrent = -1;
   fCurrentOld=-1;
   //fChain->SetMakeClass(1); // causes trouble with the skim from lxplus

   if ( BranchExists("RunNo") + BranchExists("EvtNo") == 2 ) {
     std::cout << "RunNo and EvtNo branches exist\n";
     fChain->SetBranchAddress("RunNo", &RunNo, &b_RunNo);
     fChain->SetBranchAddress("EvtNo", &EvtNo, &b_EvtNo);
   }
   else {
     b_RunNo=0;
     b_EvtNo=0;
   }

   fChain->SetBranchAddress("Momentum_Reco_Lead_BeforeEnCorr", &Momentum_Reco_Lead_BeforeEnCorr, &b_Momentum_Reco_Lead_BeforeEnCorr);
   fChain->SetBranchAddress("Momentum_Reco_Sub_BeforeEnCorr", &Momentum_Reco_Sub_BeforeEnCorr, &b_Momentum_Reco_Sub_BeforeEnCorr);
   fChain->SetBranchAddress("Momentum_Reco_Lead", &Momentum_Reco_Lead, &b_Momentum_Reco_Lead);
   fChain->SetBranchAddress("Momentum_Reco_Sub", &Momentum_Reco_Sub, &b_Momentum_Reco_Sub);
   fChain->SetBranchAddress("SCEta_Lead", &SCEta_Lead, &b_SCEta_Lead);
   fChain->SetBranchAddress("SCEta_Sub", &SCEta_Sub, &b_SCEta_Sub);
   fChain->SetBranchAddress("Charge_Reco_Lead", &Charge_Reco_Lead, &b_Charge_Reco_Lead);
   fChain->SetBranchAddress("Charge_Reco_Sub", &Charge_Reco_Sub, &b_Charge_Reco_Sub);
   fChain->SetBranchAddress("Momentum_postFSR_Lead", &Momentum_postFSR_Lead, &b_Momentum_postFSR_Lead);
   fChain->SetBranchAddress("Momentum_postFSR_Sub", &Momentum_postFSR_Sub, &b_Momentum_postFSR_Sub);
   fChain->SetBranchAddress("Momentum_preFSR_Lead", &Momentum_preFSR_Lead, &b_Momentum_preFSR_Lead);
   fChain->SetBranchAddress("Momentum_preFSR_Sub", &Momentum_preFSR_Sub, &b_Momentum_preFSR_Sub);
   fChain->SetBranchAddress("Flag_EventSelection", &Flag_EventSelection, &b_Flag_EventSelection);
   fChain->SetBranchAddress("Flag_EventSelectionExceptSCGap", &Flag_EventSelectionExceptSCGap, &b_Flag_EventSelectionExceptSCGap);
   fChain->SetBranchAddress("Flag_RecoEleSelection", &Flag_RecoEleSelection, &b_Flag_RecoEleSelection);
   fChain->SetBranchAddress("Weight_Norm", &Weight_Norm, &b_Weight_Norm);
   fChain->SetBranchAddress("Weight_PU", &Weight_PU, &b_Weight_PU);
   fChain->SetBranchAddress("Weight_Gen", &Weight_Gen, &b_Weight_Gen);
   Notify();
   return 1;
}

Bool_t DYee13TeV_t::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

  if (fCurrent!=fCurrentOld) {
    fCurrentOld=fCurrent;
    TString h1name=Form("h1preFSR_DYee13TeV_WG_%d",int(fH1PreFSR_WG.size()));
    TH1D *h1= new TH1D(h1name,h1name+TString(";M_{preFSR};weighted count(Gen)"),
		       3100,0,3100);
    h1->Sumw2();
    h1->SetDirectory(0);
    fH1PreFSR_WG.push_back(h1);
    h1name=Form("h1preFSR_DYee13TeV_WGN_%d",int(fH1PreFSR_WGN.size()));
    TH1D *h1gn= (TH1D*)h1->Clone("h1name");
    h1gn->SetTitle(h1name+TString(";M_{preFSR};weighted count (Gen#timesNorm)"));
    h1gn->SetDirectory(0);
    fH1PreFSR_WGN.push_back(h1gn);
  }
  return kTRUE;
}

void DYee13TeV_t::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
   std::cout << *this << "\n";
}

Int_t DYee13TeV_t::Cut(Long64_t )
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

// a few useful methods
void DYee13TeV_t::DeactivateBranches() {
  fChain->SetBranchStatus("*",0);
}

void DYee13TeV_t::ActivateBranches(TString brNames)
{
  if (brNames.Length()==0) {
    std::cout << "DYee13TeV_t::ActivateBranches non-empty string is expected\n";
    return;
  }
  std::stringstream ss(brNames.Data());
  TString fn;
  while (!ss.eof()) {
    ss >> fn;
    if (fn.Length()>1) {
      std::cout << "adding branch <" << fn << ">\n";
      fChain->SetBranchStatus(fn,1);
    }
  }
}

int DYee13TeV_t::BranchExists(TString brName)
{
  int yes=0;
  TObjArray *brObjArr= fChain->GetListOfBranches();
  for (int ibr=0; !yes && (ibr<brObjArr->GetEntries()); ibr++) {
    if (TString(brObjArr->At(ibr)->GetName()) == brName) yes=1;
  }
  return yes;
}

#undef DYee13TeV_t_cxx
#endif // #ifdef DYee13TeV_t_cxx
