// -*- C++ -*-
//
// Package:    DYee/DYeeNtuple
// Class:      DYeeNtuple
// 
/**\class DYeeNtuple DYeeNtuple.cc DYee/DYeeNtuple/plugins/DYeeNtuple.cc

 Description: N-tuplizer for DYee analysis at Run2

 Implementation:
     First version
*/
//
// Original Author:  Andrius Juodagalvis
//         Created:  Sat, 17 Oct 2015 08:43:46 GMT
//
//


// system include files
#include <memory>
#include <iostream>
#include <string>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
//#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
//#include "DataFormats/METReco/interface/METCollection.h"
//#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
//#include "DataFormats/TrackReco/interface/TrackFwd.h"
//#include "DataFormats/TrackReco/interface/Track.h"
//#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
//#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
//#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#ifndef DYee8TeV_reg
#define DYee8TeV_reg
#endif

//#include "DYee/Analysis/Include/TDielectron.hh"
#include "DYee/Analysis/Include/TElectron.hh"
#include "DYee/Analysis/Include/TEventInfo.hh"
//#include "DYee/DYeeNtuple/Include/TElectron.hh"
//#include "DYee/DYeeNtuple/Include/TEventInfo.hh"
//#include "DYee/Analysis/Include/TElectron.hh"
//#include "DYee/Analysis/Include/TEventInfo.hh"
//#include "DYee/Analysis/Include/TGenInfo.hh"
//#include "DYee/Analysis/Include/TGenPhoton.hh"
//#include "DYee/Analysis/Include/TJet.hh"
#include "DYee/Analysis/Include/TMuon.hh"
//#include "DYee/Analysis/Include/TPhoton.hh"
//#include "DYee/Analysis/Include/TVertex.hh"
//#include "DYee/Analysis/Include/ZeeData.hh"

#include <TClonesArray.h>
//#include <TFile.h>
#include <TTree.h>

//
// class declaration
//

class DYeeNtuple : public edm::EDAnalyzer {
  public:
  explicit DYeeNtuple(const edm::ParameterSet&);
  ~DYeeNtuple();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  // ----------member data ---------------------------
protected:
  std::string outFileName_;
  // Flags
  int doGenJets_; // process MC data
  int debug_; // whether to print debug info

  // input collections
  std::string electronCollName_;
  std::string photonCollName_;
  std::string jetCollName_;
  std::string muonCollName_;
  std::string vertexCollName_;
  std::string conversionCollName_;
  std::string offlineBSName_;
  edm::InputTag rhoCollTag_;
  edm::InputTag pfMETTag_;
  //std::string genParticleCollName_;
  //std::string genEventInfoName_;

  double electronPtMin_;
  std::vector<std::string> electronTrigNamesV_;
  std::vector<std::string> muonTrigNamesV_;

  // Tokens
  edm::EDGetTokenT<reco::GsfElectronCollection> tok_Elec;
  edm::EDGetTokenT<reco::PhotonCollection> tok_Photon;
  edm::EDGetTokenT<reco::PFJetCollection> tok_PFJet;
  edm::EDGetTokenT<reco::MuonCollection> tok_Muon;
  edm::EDGetTokenT<reco::VertexCollection> tok_Vertex;
  edm::EDGetTokenT<reco::ConversionCollection> tok_Conv;
  edm::EDGetTokenT<reco::BeamSpot> tok_BS;
  edm::EDGetTokenT<double> tok_Rho;
  edm::EDGetTokenT<reco::PFMETCollection> tok_PFMET;

  // Data to keep
  //TFile *outfile;
  TTree *outtree;
  mithep::TEventInfo *evInfo;
  TClonesArray *electronArr;

  // Misc info
  UInt_t nInputEvts, nPassEvts;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
DYeeNtuple::DYeeNtuple(const edm::ParameterSet& iConfig) :
  outtree(NULL),
  evInfo(NULL),
  electronArr(NULL),
  nInputEvts(0),
  nPassEvts(0)
{

  outFileName_ = iConfig.getUntrackedParameter<std::string>("outFileName");
  doGenJets_   = iConfig.getUntrackedParameter<int>("doGenJets");
  debug_       = iConfig.getUntrackedParameter<int>("debug");

  // input collections
  electronCollName_ = iConfig.getUntrackedParameter<std::string>("electronCollName");
  photonCollName_ = iConfig.getUntrackedParameter<std::string>("photonCollName");
  jetCollName_ = iConfig.getUntrackedParameter<std::string>("jetCollName");
  muonCollName_ = iConfig.getUntrackedParameter<std::string>("muonCollName");
  vertexCollName_ = iConfig.getUntrackedParameter<std::string>("vertexCollName");
  conversionCollName_ = iConfig.getUntrackedParameter<std::string>("conversionCollName");
  offlineBSName_ = iConfig.getUntrackedParameter<std::string>("offlineBSName");
  rhoCollTag_ = iConfig.getUntrackedParameter<edm::InputTag>("rhoInpTag");
  pfMETTag_ = iConfig.getUntrackedParameter<edm::InputTag>("pfMETTag");

  //genParticleCollName_ = iConfig.getUntrackedParameter<std::string>("genParticleCollName");
  //genEventInfoCollName_ = iConfig.getUntrackedParameter<std::string>("genEventInfoCollName");

  // Cut parameters
  electronPtMin_ = iConfig.getUntrackedParameter<double>("electronPtMin");
  electronTrigNamesV_ = iConfig.getUntrackedParameter<std::vector<std::string>>("electronTrigNames");
  muonTrigNamesV_ = iConfig.getUntrackedParameter<std::vector<std::string>>("muonTrigNames");

  // Register tokens
  tok_Elec = consumes<reco::GsfElectronCollection>(electronCollName_);
  tok_Photon= consumes<reco::PhotonCollection>(photonCollName_);
  //tok_PFJet = consumes<reco::PFJetCollection>(pfJetCollName_);
  tok_Muon = consumes<reco::MuonCollection>(muonCollName_);
  tok_Vertex = consumes<reco::VertexCollection>(vertexCollName_);
  //tok_Conv= consumes<reco::ConversionCollection>(convCollName_);
  tok_BS= consumes<reco::BeamSpot>(offlineBSName_);
  //tok_Rho= consumes<double>(rhoInpTag_);
  tok_PFMET= consumes<reco::PFMETCollection>(pfMETTag_);

  // Create structures
  evInfo = new mithep::TEventInfo();
  electronArr = new TClonesArray("mithep::TElectron");
  mithep::TElectron dummyE;
  mithep::TMuon dummyM;
}


DYeeNtuple::~DYeeNtuple()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
DYeeNtuple::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

}


// ------------ method called once each job just before starting event loop  ------------
void DYeeNtuple::beginJob()
{

  // Check structures
  if (!evInfo || !electronArr) {
    edm::LogError("DYee") << "structures contain null pointer\n";
    throw edm::Exception(edm::errors::NullPointerError);
  }

  edm::Service<TFileService> fs;
  outtree = fs->make<TTree> ("Events","Selected events");

  if (!outtree) {
    edm::LogError("DYee") << "failed to create the output tree";
    throw edm::Exception(edm::errors::FileOpenError);
  }
}

// ------------ method called once each job just after ending the event loop  ------------
void DYeeNtuple::endJob()
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
DYeeNtuple::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
DYeeNtuple::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
DYeeNtuple::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
DYeeNtuple::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DYeeNtuple::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DYeeNtuple);
