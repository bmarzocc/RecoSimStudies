import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

process = cms.Process("RecoSimAnalysis")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff") # gives deprecated message in 80X but still runs
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag,'106X_mc2017_realistic_v3','')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( 1 ) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1 )
                                                                       
process.source = cms.Source("PoolSource",
    skipEvents = cms.untracked.uint32(0),                       
    fileNames = cms.untracked.vstring("file:/afs/cern.ch/work/b/bmarzocc/Clustering/CMSSW_10_6_0/src/SimGeneral/CaloAnalysis/10005.0_SingleGammaPt35+SingleGammaPt35_pythia8_2017_GenSimFull+DigiFull_2017+RecoFull_2017+HARVESTFull_2017+ALCAFull_2017/step3.root"),
    secondaryFileNames = cms.untracked.vstring()
    ) 

process.load('RecoSimStudies.Dumpers.RecoSimDumper_cfi')

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('RecoSimDumper.root')
)

process.p = cms.Path(
    process.recosimdumper
)
