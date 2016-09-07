#define DYee13TeV_t_cxx
#include "DYee13TeV.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

// -------------------------------------------------------------

void DYee13TeV_t::Zero()
{
  TLorentzVector zero;
  zero.SetPxPyPzE(0,0,0,0);
  RunNo=0; EvtNo=0;
  *Momentum_Reco_Lead_BeforeEnCorr=zero;
  *Momentum_Reco_Sub_BeforeEnCorr=zero;
  *Momentum_Reco_Lead=zero;
  *Momentum_Reco_Sub=zero;
  SCEta_Lead=0;
  SCEta_Sub=0;
  Charge_Reco_Lead=0;
  Charge_Reco_Sub=0;
  *Momentum_postFSR_Lead=zero;
  *Momentum_postFSR_Sub=zero;
  *Momentum_preFSR_Lead=zero;
  *Momentum_preFSR_Sub=zero;
  Flag_EventSelection=false;
  Flag_EventSelectionExceptSCGap=false;
  Weight_Norm=1;
  Weight_PU=1;
  Weight_Gen=1;
}

// -------------------------------------------------------------

int DYee13TeV_t::CreateNew(TString fnameInp, TString treeName)
{
  std::stringstream ss(fnameInp.Data());
  TString fname;
  ss >> fname;
  delete fChain; fChain=NULL;

  // TChain cannot be used to fill a tree. Use TTree directly.
  fOutFile= new TFile(fname,"recreate");
  if (!fOutFile || !fOutFile->IsOpen()) {
    std::cout << "failed to open the file <" << fnameInp << ">\n";
    return 0;
  }
  fOutTree= new TTree(treeName,"semi-selected events");


  Momentum_Reco_Lead_BeforeEnCorr= new TLorentzVector();
  Momentum_Reco_Sub_BeforeEnCorr= new TLorentzVector();
  Momentum_Reco_Lead= new TLorentzVector();
  Momentum_Reco_Sub= new TLorentzVector();
  Momentum_postFSR_Lead= new TLorentzVector();
  Momentum_postFSR_Sub= new TLorentzVector();
  Momentum_preFSR_Lead= new TLorentzVector();
  Momentum_preFSR_Sub= new TLorentzVector();

  fOutTree->Branch("RunNo", &RunNo, "RunNo/D");
  fOutTree->Branch("EvtNo", &EvtNo, "EvtNo/D");
  fOutTree->Branch("Momentum_Reco_Lead_BeforeEnCorr", &Momentum_Reco_Lead_BeforeEnCorr);
  fOutTree->Branch("Momentum_Reco_Sub_BeforeEnCorr", &Momentum_Reco_Sub_BeforeEnCorr);
  fOutTree->Branch("Momentum_Reco_Lead", Momentum_Reco_Lead);
  fOutTree->Branch("Momentum_Reco_Sub", Momentum_Reco_Sub);
  fOutTree->Branch("SCEta_Lead", &SCEta_Lead, "SCEta_Lead/D");
  fOutTree->Branch("SCEta_Sub", &SCEta_Sub, "SCEta_Sub/D");
  fOutTree->Branch("Charge_Reco_Lead", &Charge_Reco_Lead);
  fOutTree->Branch("Charge_Reco_Sub", &Charge_Reco_Sub);
  fOutTree->Branch("Momentum_postFSR_Lead", &Momentum_postFSR_Lead);
  fOutTree->Branch("Momentum_postFSR_Sub", &Momentum_postFSR_Sub);
  fOutTree->Branch("Momentum_preFSR_Lead", &Momentum_preFSR_Lead);
  fOutTree->Branch("Momentum_preFSR_Sub", &Momentum_preFSR_Sub);
  fOutTree->Branch("Flag_EventSelection", &Flag_EventSelection, "Flag_EventSelection/O");
  fOutTree->Branch("Flag_EventSelectionExceptSCGap", &Flag_EventSelectionExceptSCGap, "Flag_EventSelectionExceptSCGap/O");
  fOutTree->Branch("Weight_Norm", &Weight_Norm, "Weight_Norm/D");
  fOutTree->Branch("Weight_PU", &Weight_PU, "Weight_PU/D");
  fOutTree->Branch("Weight_Gen", &Weight_Gen, "Weight_Gen/D");
  return 1;
}

// -------------------------------------------------------------

void DYee13TeV_t::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L DYee13TeV_t.C
//      Root > DYee13TeV_t t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
   }
}

// -------------------------------------------------------------

#ifndef Inputs_H
std::ostream& operator<<(std::ostream &out, const TLorentzVector &v)
{
  out << Form("(Pt,Eta,Phi,En)=(%lf,%lf,%lf,%lf)",v.Pt(),v.Eta(),v.Phi(),v.E());
  return out;
}
#endif

// -------------------------------------------------------------

#ifndef Inputs_H
std::ostream& operator<<(std::ostream &out, const TLorentzVector *v)
{
  if (!v) out << "(null ptr to TLorentzVector)";
  else out << (*v);
  return out;
}
#endif

// -------------------------------------------------------------

std::ostream& operator<<(std::ostream &out, DYee13TeV_t &obj)
{
  out << "DYee13TeV: \n";
  if ((obj.RunNo!=0) && (obj.EvtNo!=0))
    out << "Run,Evt No=" << obj.RunNo << "," << obj.EvtNo << "\n";
  out << "Mom.before corr. " << obj.Momentum_Reco_Lead_BeforeEnCorr
      << ", " << obj.Momentum_Reco_Sub_BeforeEnCorr << "\n";
  out << "Mom. " << obj.Momentum_Reco_Lead << " "
      << obj.Momentum_Reco_Sub << "\n";
  out << "ScEta " << obj.SCEta_Lead << " " << obj.SCEta_Sub << "\n";
  out << "Charge: " << obj.Charge_Reco_Lead << ", " << obj.Charge_Reco_Sub << "\n";
  out << "Mom. postFSR: " << (*obj.Momentum_postFSR_Lead)
      << ", " << (*obj.Momentum_postFSR_Sub) << "\n";
  out << "Mom. preFSR : " << (*obj.Momentum_preFSR_Lead)
      << ", " << (*obj.Momentum_preFSR_Sub) << "\n";
  out << "Event ";
  if (!obj.Flag_EventSelection) out << "NOT ";
  out << "selected";
  if (!obj.Flag_EventSelection && obj.Flag_EventSelectionExceptSCGap)
    out << " (because of SC gap)";
  out << ", Weights: lumi=" << obj.Weight_Norm
      << ", PU=" << obj.Weight_PU << ", Gen=" << obj.Weight_Gen << "\n";
  return out;
}

// -------------------------------------------------------------
