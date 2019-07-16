import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing ('analysis') # keep this name
options.outputFile = 'RecoSimDumper.root'
options.inputFiles = 'file:test/step3.root'
options.maxEvents = -1 # -1 means all events, maxEvents considers the total over files considered
options.parseArguments()

process = cms.Process("RecoSimAnalysis")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff") # gives deprecated message in 80X but still runs
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag,'106X_mc2017_realistic_v3','')

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1 )
                                                                       
process.source = cms.Source("PoolSource",
    skipEvents = cms.untracked.uint32(0),                       
    fileNames = cms.untracked.vstring (options.inputFiles),
    secondaryFileNames = cms.untracked.vstring()
    ) 

process.load('RecoSimStudies.Dumpers.Cfg_RecoSimDumper_cfi')

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(options.outputFile)
)

process.p = cms.Path(
    process.recosimdumper
)
