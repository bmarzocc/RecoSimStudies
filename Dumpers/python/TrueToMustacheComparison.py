import FWCore.ParameterSet.Config as cms

process = cms.PSet()

process.ioFilesOpt = cms.PSet(

    ##input file
    inputFile = cms.string('/eos/cms/store/group/dpg_ecal/alca_ecalcalib/bmarzocc/Clustering/FourGammasGunPt1-100_pythia8_withPU_withTracker_106X_mcRun3_2021_realistic_v3_RAW_StdSeedingGathering_Mustache_Dumper_v5_Total_1.root'),

    ## output directory
    outDir = cms.string('~/www/Clustering_TrueToMustache/'),

    ## fitFunction
    fitFunction = cms.string('cruijff'),
    #fitFunction = cms.string('doubleCB'), 

    ## maxEvents
    maxEvents = cms.untracked.int32(100000),  

    useMustacheWindows = cms.bool(False),

    etBinning = cms.string('0 10 20 30 40 50 60 70 80 90 100'),
    etaBinning = cms.string('0.0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.479 1.75 2.0 2.25 3.0')
)
