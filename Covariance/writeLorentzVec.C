#include <TROOT.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TBranch.h>
#include <TFile.h>
#include <iostream>

void writeLorentzVec(int read=0)
{
  TLorentzVector v;

  if (read==0) {
    TFile fout("test.root","RECREATE");
    TTree *tree= new TTree("tree","test");
    //tree->Branch("vec","TLorentzVector",&v);
    tree->Branch("vec",&v);
    v.SetPtEtaPhiE(20,1.5,2.1,100);
    tree->Fill();
    v.SetPtEtaPhiE(25,-1.2,1.3,40);
    tree->Fill();
    v.SetPtEtaPhiE(25,-1.2,1.3,140);
    tree->Fill();
    tree->Write();
    fout.Close();
  }

  v.SetPtEtaPhiE(0,0,0,0);

  if (read==1) {
    TLorentzVector *vin= new TLorentzVector();
    TFile fin("test.root");
    TTree *tin= (TTree*)fin.Get("tree");
    tin->SetBranchAddress("vec",&vin);
    int nentries= tin->GetEntries();
    for (int iEntry=0; iEntry<nentries; iEntry++) {
      tin->GetEntry(iEntry);
      std::cout << "got: " << vin->Pt() << ", " << vin->Eta() << ", " << vin->Phi() << ", " << vin->E() << "\n";
    }
    fin.Close();
  }

  if (read==2) {
    TLorentzVector *vin ;//= new TLorentzVector();
    TFile fin("../../../v2ee-skim1-20160812/DYEE_M400to500_v1.root");
    TTree *tin= (TTree*)fin.Get("DYTree");
    tin->SetBranchAddress("Momentum_postFSR_Lead",&vin);
    int nentries= tin->GetEntries();
    for (int iEntry=0; iEntry<nentries; iEntry++) {
      if (iEntry>20) break;
      tin->GetEntry(iEntry);
      std::cout << "got: " << vin->Pt() << ", " << vin->Eta() << ", " << vin->Phi() << ", " << vin->E() << "\n";
    }
    fin.Close();
  }
}
