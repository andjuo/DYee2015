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

//#define OUTMSG edm::LogWarning("DYee")
#define OUTMSG std::cout


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
  //std::string outFileName_;
  // Flags
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
  edm::EDGetTokenT<edm::View<reco::GsfElectron> > tok_Elec;
  edm::EDGetTokenT<edm::View<reco::Photon> > tok_Photon;
  edm::EDGetTokenT<edm::View<reco::PFJet> > tok_PFJet;
  edm::EDGetTokenT<edm::View<reco::Muon> > tok_Muon;
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
  double rho;

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
  electronPtMin_(10),
  outtree(NULL),
  evInfo(NULL),
  electronArr(NULL),
  nInputEvts(0),
  nPassEvts(0)
{

  //outFileName_ = iConfig.getUntrackedParameter<std::string>("outFileName");
  debug_       = iConfig.getUntrackedParameter<int>("debug");
  electronPtMin_ = iConfig.getUntrackedParameter<double>("electronPtMin");

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
  tok_Rho= consumes<double>(rhoInpTag_);
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
  nInputEvts++;

  edm::Handle<edm::View<reco::GsfElectron> > elHandle;
  iEvent.getByToken(tok_Elec, elHandle);
  if (!elHandle.isValid()) {
    OUTMSG << "DYeeNtuple: electron collection is not valid";
    return;
  }

  edm::Handle<edm::View<reco::Photon> > phoHandle;
  iEvent.getByToken(tok_Photon, phoHandle);
  if (!phoHandle.isValid()) {
    OUTMSG << "DYeeNtuple: photon collection is not valid";
    return;
  }

  edm::Handle<edm::View<reco::Muon> > muonHandle;
  iEvent.getByToken(tok_Muon, muonHandle);
  if (!muonHandle.isValid()) {
    OUTMSG << "DYeeNtuple: muon collection is not valid";
    return;
  }

  edm::Handle<reco::VertexCollection> vertexHandle;
  iEvent.getByToken(tok_Vertex, vertexHandle);
  if (!vertexHandle.isValid()) {
    OUTMSG << "DYeeNtuple: vertex handle is not valid";
    return;
  }

  // Find the first vertex in the collection that passes
  // good quality criteria
  VertexCollection::const_iterator firstGoodVertex = vertices->end();
  int firstGoodVertexIdx = 0;
  int nGoodPVs=0;
  for (VertexCollection::const_iterator vtx = vertices->begin();
       vtx != vertices->end(); ++vtx, ++firstGoodVertexIdx) {
    // Replace isFake() for miniAOD because it requires tracks and miniAOD vertices don't have tracks:
    // Vertex.h: bool isFake() const {return (chi2_==0 && ndof_==0 && tracks_.empty());}
    bool isFake = vtx->isFake();
    //if( !isAOD )
    isFake = (vtx->chi2()==0 && vtx->ndof()==0);
    // Check the goodness
    if ( !isFake
	 && vtx->ndof()>=4. && vtx->position().Rho()<=2.0
	 && fabs(vtx->position().Z())<=24.0) {
      if (nGoodPVs==0) firstGoodVertex = vtx;
      nGoodPVs++;
      break;
    }
  }
  if ( firstGoodVertex==vertices->end() ) {
    OUTMSG << "event has no good PVs";
    return; // skip event if there are no good PVs
  }

  edm::Handle<double> rhoHandle;
  iEvent.getByToken(tok_Rho, rhoHandle);
  if (!rhoHandle.isValid()) {
    OUTMSG << "DYeeNtuple: rho handle is not valid";
    return;
  }

  // Good event. Collect info

  evInfo->Clear();
  electronArr->Clear();


  // fill event info
  evInfo->runNum= iEvent.id().run();
  evInfo->evtNum= iEvent.id().event();
  evInfo->lumiSec=iEvent.id().luminosityBlock();
  evInfo->pvx = firstGoodVertex->px();
  evInfo->pvy = firstGoodVertex->py();
  evInfo->pvz = firstGoodVertex->pz();
  evInfo->rho = (*rhoHandle.product());
  evInfo->hasGoodPV = kTRUE;

  // fill electron info
  unsigned int nElectrons=0;
  for (size_t i = 0; i < elHandle->size(); ++i){
    const auto cEl = elHandle->ptrAt(i);
    if( cEl->pt() < electronPtMin_ ) // keep only electrons above the threshold
      continue;

    new((*electronArr)[nElectrons]) mithep::TElectron();
    mithep::TElectron *el= (mithep::TElectron*)(*electronArr[nElectrons]);
    nElectrons++;

    el->pt = cEl->pt();
    el->eta= cEl->eta();
    el->phi= cEl->phi();
    el->scEt = cEl->superCluster()->scEt();
    el->scEta= cEl->superCluster()->eta();
    el->scPhi= cEl->superCluster()->phi();
    // ID and matching
    el->deltaEtaIn = cEl->deltaEtaSuperClusterTrackAtVtx();
    el->deltaPhiIn = cEl->deltaPhiSuperClusterTrackAtVtx();
    el->HoverE =     cEl->hcalOverEcal();
    el->sigiEtaiEta= cEl->full5x5_sigmaIetaIeta();
    /*
    // |1/E-1/p| = |1/E - EoverPinner/E| is computed below
    // The if protects against ecalEnergy == inf or zero
    // (always the case for miniAOD for electrons <5 GeV)
    if( cEl->ecalEnergy() == 0 ){
      printf("Electron energy is zero!\n");
      ooEmooP_.push_back( 1e30 );
      }else if( !std::isfinite(cEl->ecalEnergy())){
      printf("Electron energy is not finite!\n");
      ooEmooP_.push_back( 1e30 );
      }else{
      ooEmooP_.push_back( fabs(1.0/cEl->ecalEnergy() - cEl->eSuperClusterOverP()/cEl->ecalEnergy() ) );
      }
    */

    // Isolation
    GsfElectron::PflowIsolationVariables pfIso = cEl->pfIsolationVariables();
    // Compute individual PF isolations
    el->hadIso03 = pfIso.sumChargedHadronPt;
    el->neuHadIso_00_01 = pfIso.sumNeutralHadronEt;
    el->gammaIso_00_01  = pfIso.sumPhotonEt;
    //isoChargedFromPU_.push_back( pfIso.sumPUPt );
    // Compute combined relative PF isolation with the effective area correction for pile-up
    //float abseta = abs(el->superCluster()->eta());
    //float eA = effectiveAreas_.getEffectiveArea(abseta);
    //relCombIsoWithEA_.push_back( ( pfIso.sumChargedHadronPt
    //				   + std::max( 0.0f, pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - eA*rho_) )
    //				 / el->pt() );
    // Impact parameter
    reco::GsfTrackRef theTrack = cEl->gsfTrack();
    el->d0 = (-1) * theTrack->dxy(firstGoodVertex->position() );
    el->dz = theTrack->dz( firstGoodVertex->position() );
    // Conversion rejection
    el->nExpHitsInner = cEl->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);
    //el->isConv = ConversionTools::hasMatchedConversion(*el,conversions,
    //					       theBeamSpot->position());
    // Match to generator level truth
    //isTrue_.push_back( matchToTruth( el, genParticles) );
  }

  nPassEvts++;
  outtree->Fill();
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

  outtree->Branch("Info", &evInfo);
  outtree->Branch("Electron", &electronArr);
}

// ------------ method called once each job just after ending the event loop  ------------
void DYeeNtuple::endJob()
{
  edm::LogWarning("DYee") << "nInputEvts=" << nInputEvts
			  << ", nPassEvts=" << nPassEvts;
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
