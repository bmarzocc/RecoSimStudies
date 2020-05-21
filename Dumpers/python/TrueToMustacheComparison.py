import FWCore.ParameterSet.Config as cms

process = cms.PSet()

process.ioFilesOpt = cms.PSet(

    ##input file
    inputFile = cms.string('/eos/cms/store/group/dpg_ecal/alca_ecalcalib/bmarzocc/Clustering/FourElectronsGunPt1-100_pythia8_withPU_withTracker_106X_mcRun3_2021_realistic_v3_RAW_StdSeedingGathering_Mustache_Dumper_v4_Total_1.root'),
    
    ## output directory
    outDir = cms.string('~/www/Clustering_TrueToMustache_simFraction/'),

    ## fitFunction
    fitFunction = cms.string('cruijff'),
    #fitFunction = cms.string('doubleCB'), 

    ## maxEvents
    maxEvents = cms.untracked.int32(-1),  

    ## caloParticle matching
    matching = cms.string('sim_fraction'),
    #matching = cms.string('sim_fraction_old'),
    
    useMustacheWindows = cms.bool(False),

    etBinning = cms.string('0 10 20 30 40 50 60 70 80 90 100'),
    etaBinning = cms.string('0.0 0.2 0.4 0.6 0.8 1. 1.2 1.4 1.479 1.75 2.0 2.25 3.0'),
    scoreBinning = cms.string('0.0001 0.0002 0.0003 0.0004 0.0005 0.0006 0.0007 0.0008 0.0009 0.001 0.0015 0.002 0.0025 0.003 0.004 0.005 0.01 0.02 0.03 0.04 0.05' )
    #scoreBinning = cms.string('0.005 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5' )
)
