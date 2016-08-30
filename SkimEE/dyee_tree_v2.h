//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Aug 10 20:26:10 2016 by ROOT version 5.34/36
// from TTree tree/ after preselections tree
// found on file: ../../DYEE_M400to500.root
//////////////////////////////////////////////////////////

#ifndef dyee_tree_v2_h
#define dyee_tree_v2_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <sstream>
#include <iostream>

// Fixed size dimensions of array or collections stored in the TTree if any.

class dyee_tree_v2 {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   vector<float>   *genPostFSR_Pt;
   vector<float>   *genPostFSR_Eta;
   vector<float>   *genPostFSR_Rap;
   vector<float>   *genPostFSR_Phi;
   vector<float>   *genPostFSR_En;
   vector<float>   *genPreFSR_Pt;
   vector<float>   *genPreFSR_Eta;
   vector<float>   *genPreFSR_Rap;
   vector<float>   *genPreFSR_Phi;
   vector<float>   *genPreFSR_En;
   vector<float>   *genPhoton_Pt;
   vector<float>   *genPhoton_Eta;
   vector<float>   *genPhoton_Rap;
   vector<float>   *genPhoton_Phi;
   vector<float>   *genPhoton_En;
   vector<float>   *ptElec;
   vector<float>   *etaElec;
   vector<float>   *rapElec;
   vector<float>   *phiElec;
   vector<float>   *energyElec;
   vector<float>   *chargeElec;
   vector<float>   *etaSC;
   vector<int>     *passMediumId;
   vector<double>  *pt_Ele23;
   vector<double>  *eta_Ele23;
   vector<double>  *phi_Ele23;
   Int_t           tauFlag;
   Int_t           nPV;
   Int_t           nPUTrue;
   Double_t        theWeight;
   Char_t          Ele23_WPLoose;
   vector<float>   *ptMuon;
   vector<float>   *etaMuon;
   vector<float>   *phiMuon;
   vector<float>   *energyMuon;
   vector<float>   *chargeMuon;
   vector<float>   *isoPFMuon;
   vector<bool>    *isTightMuon;
   Char_t          Mu8_Ele17;

   // List of branches
   TBranch        *b_genPostFSR_Pt;   //!
   TBranch        *b_genPostFSR_Eta;   //!
   TBranch        *b_genPostFSR_Rap;   //!
   TBranch        *b_genPostFSR_Phi;   //!
   TBranch        *b_genPostFSR_En;   //!
   TBranch        *b_genPreFSR_Pt;   //!
   TBranch        *b_genPreFSR_Eta;   //!
   TBranch        *b_genPreFSR_Rap;   //!
   TBranch        *b_genPreFSR_Phi;   //!
   TBranch        *b_genPreFSR_En;   //!
   TBranch        *b_genPhoton_Pt;   //!
   TBranch        *b_genPhoton_Eta;   //!
   TBranch        *b_genPhoton_Rap;   //!
   TBranch        *b_genPhoton_Phi;   //!
   TBranch        *b_genPhoton_En;   //!
   TBranch        *b_ptElec;   //!
   TBranch        *b_etaElec;   //!
   TBranch        *b_rapElec;   //!
   TBranch        *b_phiElec;   //!
   TBranch        *b_energyElec;   //!
   TBranch        *b_chargeElec;   //!
   TBranch        *b_etaSC;   //!
   TBranch        *b_passMediumId;   //!
   TBranch        *b_pt_Ele23;   //!
   TBranch        *b_eta_Ele23;   //!
   TBranch        *b_phi_Ele23;   //!
   TBranch        *b_tauFlag;   //!
   TBranch        *b_nPV;   //!
   TBranch        *b_nPUTrue;   //!
   TBranch        *b_theWeight;   //!
   TBranch        *b_Ele23_WPLoose;   //!
   TBranch        *b_ptMuon;   //!
   TBranch        *b_etaMuon;   //!
   TBranch        *b_phiMuon;   //!
   TBranch        *b_energyMuon;   //!
   TBranch        *b_chargeMuon;   //!
   TBranch        *b_isoPFMuon;   //!
   TBranch        *b_isTightMuon;   //!
   TBranch        *b_Mu8_Ele17;   //!

   dyee_tree_v2(TTree *tree=0);
   virtual ~dyee_tree_v2();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   // Added methods
   UInt_t GetEntries() { return fChain->GetEntries(); }
   void DeactivateAllBranches();
   void ActivateBranches(TString brNames);
};

#endif

#ifdef dyee_tree_v2_cxx
dyee_tree_v2::dyee_tree_v2(TTree *tree) : fChain(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../../DYEE_M400to500.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../../DYEE_M400to500.root");
      }
      f->GetObject("tree",tree);

   }
   Init(tree);
}

dyee_tree_v2::~dyee_tree_v2()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t dyee_tree_v2::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t dyee_tree_v2::LoadTree(Long64_t entry)
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

void dyee_tree_v2::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   genPostFSR_Pt = 0;
   genPostFSR_Eta = 0;
   genPostFSR_Rap = 0;
   genPostFSR_Phi = 0;
   genPostFSR_En = 0;
   genPreFSR_Pt = 0;
   genPreFSR_Eta = 0;
   genPreFSR_Rap = 0;
   genPreFSR_Phi = 0;
   genPreFSR_En = 0;
   genPhoton_Pt = 0;
   genPhoton_Eta = 0;
   genPhoton_Rap = 0;
   genPhoton_Phi = 0;
   genPhoton_En = 0;
   ptElec = 0;
   etaElec = 0;
   rapElec = 0;
   phiElec = 0;
   energyElec = 0;
   chargeElec = 0;
   etaSC = 0;
   passMediumId = 0;
   pt_Ele23 = 0;
   eta_Ele23 = 0;
   phi_Ele23 = 0;
   ptMuon = 0;
   etaMuon = 0;
   phiMuon = 0;
   energyMuon = 0;
   chargeMuon = 0;
   isoPFMuon = 0;
   isTightMuon = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("genPostFSR_Pt", &genPostFSR_Pt, &b_genPostFSR_Pt);
   fChain->SetBranchAddress("genPostFSR_Eta", &genPostFSR_Eta, &b_genPostFSR_Eta);
   fChain->SetBranchAddress("genPostFSR_Rap", &genPostFSR_Rap, &b_genPostFSR_Rap);
   fChain->SetBranchAddress("genPostFSR_Phi", &genPostFSR_Phi, &b_genPostFSR_Phi);
   fChain->SetBranchAddress("genPostFSR_En", &genPostFSR_En, &b_genPostFSR_En);
   fChain->SetBranchAddress("genPreFSR_Pt", &genPreFSR_Pt, &b_genPreFSR_Pt);
   fChain->SetBranchAddress("genPreFSR_Eta", &genPreFSR_Eta, &b_genPreFSR_Eta);
   fChain->SetBranchAddress("genPreFSR_Rap", &genPreFSR_Rap, &b_genPreFSR_Rap);
   fChain->SetBranchAddress("genPreFSR_Phi", &genPreFSR_Phi, &b_genPreFSR_Phi);
   fChain->SetBranchAddress("genPreFSR_En", &genPreFSR_En, &b_genPreFSR_En);
   fChain->SetBranchAddress("genPhoton_Pt", &genPhoton_Pt, &b_genPhoton_Pt);
   fChain->SetBranchAddress("genPhoton_Eta", &genPhoton_Eta, &b_genPhoton_Eta);
   fChain->SetBranchAddress("genPhoton_Rap", &genPhoton_Rap, &b_genPhoton_Rap);
   fChain->SetBranchAddress("genPhoton_Phi", &genPhoton_Phi, &b_genPhoton_Phi);
   fChain->SetBranchAddress("genPhoton_En", &genPhoton_En, &b_genPhoton_En);
   fChain->SetBranchAddress("ptElec", &ptElec, &b_ptElec);
   fChain->SetBranchAddress("etaElec", &etaElec, &b_etaElec);
   fChain->SetBranchAddress("rapElec", &rapElec, &b_rapElec);
   fChain->SetBranchAddress("phiElec", &phiElec, &b_phiElec);
   fChain->SetBranchAddress("energyElec", &energyElec, &b_energyElec);
   fChain->SetBranchAddress("chargeElec", &chargeElec, &b_chargeElec);
   fChain->SetBranchAddress("etaSC", &etaSC, &b_etaSC);
   fChain->SetBranchAddress("passMediumId", &passMediumId, &b_passMediumId);
   fChain->SetBranchAddress("pt_Ele23", &pt_Ele23, &b_pt_Ele23);
   fChain->SetBranchAddress("eta_Ele23", &eta_Ele23, &b_eta_Ele23);
   fChain->SetBranchAddress("phi_Ele23", &phi_Ele23, &b_phi_Ele23);
   fChain->SetBranchAddress("tauFlag", &tauFlag, &b_tauFlag);
   fChain->SetBranchAddress("nPV", &nPV, &b_nPV);
   fChain->SetBranchAddress("nPUTrue", &nPUTrue, &b_nPUTrue);
   fChain->SetBranchAddress("theWeight", &theWeight, &b_theWeight);
   fChain->SetBranchAddress("Ele23_WPLoose", &Ele23_WPLoose, &b_Ele23_WPLoose);
   fChain->SetBranchAddress("ptMuon", &ptMuon, &b_ptMuon);
   fChain->SetBranchAddress("etaMuon", &etaMuon, &b_etaMuon);
   fChain->SetBranchAddress("phiMuon", &phiMuon, &b_phiMuon);
   fChain->SetBranchAddress("energyMuon", &energyMuon, &b_energyMuon);
   fChain->SetBranchAddress("chargeMuon", &chargeMuon, &b_chargeMuon);
   fChain->SetBranchAddress("isoPFMuon", &isoPFMuon, &b_isoPFMuon);
   fChain->SetBranchAddress("isTightMuon", &isTightMuon, &b_isTightMuon);
   fChain->SetBranchAddress("Mu8_Ele17", &Mu8_Ele17, &b_Mu8_Ele17);
   Notify();
}

Bool_t dyee_tree_v2::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void dyee_tree_v2::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

Int_t dyee_tree_v2::Cut(Long64_t entry)
{
  if (0) entry=0;
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

// a few useful methods
void dyee_tree_v2::DeactivateAllBranches() {
  fChain->SetBranchStatus("*",0);
}

void dyee_tree_v2::ActivateBranches(TString brNames)
{
  if (brNames.Length()==0) {
    std::cout << "dyee_tree_v2::ActivateBranches non-empty string is expected\n";
    return;
  }
  std::stringstream ss(brNames.Data());
  TString fn;
  while (!ss.eof()) {
    ss >> fn;
    std::cout << "adding branch <" << fn << ">\n";
    fChain->SetBranchStatus(fn,1);
  }
}



#undef dyee_tree_v2_cxx
#endif // #ifdef dyee_tree_v2_cxx
