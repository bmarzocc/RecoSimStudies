import FWCore.ParameterSet.Config as cms

process = cms.PSet()

process.ioFilesOpt = cms.PSet(

    ##input file
    inputFile = cms.string('/eos/cms/store/group/dpg_ecal/alca_ecalcalib/bmarzocc/Clustering/FourElectronsGunPt1-100_pythia8_withPU_withTracker_106X_mcRun3_2021_realistic_v3_RAW_StdSeedingGathering_Mustache_Dumper_v4_Total_1.root'),
    
    ## maxEvents
    maxEvents = cms.untracked.int32(-1),
 
    ##fitFunction
    fitFunction = cms.string('cruijff'),
    #fitFunction = cms.string('doubleCB'),   

    matching = cms.string('sim_fraction'),
 
    etBinning = cms.string('0 10 20 30 40 50 60 70 80 90 100'),
    etaBinning = cms.string('-3.0 -2.5 -1.75 -1.479 -0.75 0.0 0.75 1.479 1.75 2.5 3.0'),
    scoreBinning = cms.string('0.0001 0.0002 0.0003 0.0004 0.0005 0.0006 0.0007 0.0008 0.0009 0.001 0.0015 0.002 0.0025 0.003 0.004 0.005 0.01 0.02 0.03 0.04 0.05' )

)
