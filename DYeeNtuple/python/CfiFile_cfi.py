import FWCore.ParameterSet.Config as cms

DYeeNtuplizer = cms.EDAnalyzer('DYeeNtuple',
    (
    #outFileName= cms.untracked.string("outFileName"),
    doGenJets= cms.untracked.int32(0),

    electronCollName= cms.untracked.string("gedGsfElectrons"),
    photonCollName= cms.untracked.string("gedPhotons"),
    jetCollName= cms.untracked.string("ak4PFJetsCHS"),
    muonCollName= cms.untracked.string("muons"),
    vertexCollName= cms.untracked.string("offlinePrimaryVertices"),
    conversionNollName= cms.untracked.string("allConversions"),
    offlineBSName= cms.untracked.string("offlineBeamSpot"),
    rhoInpTag= cms.untracked.InputTag("fixedGridRhoFastjetAll"),
    pfMETTag= cms.untracked.InputTag("pfMet"),

    electronPtMin= cms.untracked.double(10.),
    electronTrigNames= cms.untracked.vstring(""),
    muonTrigNames= cms.untracked.vstring("")
)
