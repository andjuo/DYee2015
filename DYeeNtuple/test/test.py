import FWCore.ParameterSet.Config as cms
process = cms.Process("Demo")

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/tmp/andriusj/FileA_468.root',
        'file:/tmp/andriusj/FileB_0D6.root'
        #'file:myfile.root'
    )
)

# Specify IdealMagneticField ESSource (needed for CMSSW 730)
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("MagneticField.Engine.autoMagneticFieldProducer_cfi")
from Configuration.AlCa.autoCond import autoCond
#process.GlobalTag.globaltag=autoCond['startup']
process.GlobalTag.globaltag = 'PHYS14_25_V1::All'


process.load("DYee.DYeeNtuple.CfiFile_cfi")

process.p = cms.Path(process.DYeeNtuplizer)

outputFile="test.root"
print "outputFile=%s" : outputFile

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string( outputFile )
)

# To check what collections are present, define the output module
process.Output = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('keep *'),
    fileName = cms.untracked.string('gjetOut.root')
)

