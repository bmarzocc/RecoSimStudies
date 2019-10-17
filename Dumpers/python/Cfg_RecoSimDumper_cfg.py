#example of utilisation: cmsRun Cfg_RecoSimDumper_cfg.py outputFile=../test/outputfiles/dumpedFiles/dumped_singlePhoton_5k_EB.root inputFiles=file:../test/outputfiles/singlePhoton_5k_EB/step3.root 

#cmsRun Cfg_RecoSimDumper_cfg.py outputFile=../test/outputfiles/dumpedFiles/dumped_singlePhoton_150k_EB.root inputFiles=file:root://t3dcachedb.psi.ch:1094//pnfs/psi.ch/cms/trivcat/store/user/anlyon/EcalProd/singlePhoton_closeECAL_0to100GeV_150k_EB/step3.root 

import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing ('analysis') # keep this name
options.outputFile = 'RecoSimDumper.root'
#options.inputFiles = 'file:test/step3.root'
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
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 100 )
                                                                       
process.source = cms.Source("PoolSource",
    skipEvents = cms.untracked.uint32(0),                       
    fileNames = cms.untracked.vstring (options.inputFiles),
    secondaryFileNames = cms.untracked.vstring(),
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck")
) 

process.load('RecoSimStudies.Dumpers.Cfg_RecoSimDumper_cfi')

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(options.outputFile)
)

#Setup FWK for multithreaded 
process.options = cms.untracked.PSet( 
 
 )
# uncomment these lines for multithreaded
#process.options.numberOfThreads=cms.untracked.uint32(8)
#process.options.numberOfStreams=cms.untracked.uint32(0)
#process.options.numberOfConcurrentLuminosityBlocks=cms.untracked.uint32(1)

process.p = cms.Path(
    process.recosimdumper
)

