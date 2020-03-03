import FWCore.ParameterSet.Config as cms

process = cms.PSet()

process.ioFilesOpt = cms.PSet(

    ##input file
    inputFiles = cms.string('/eos/cms/store/group/dpg_ecal/alca_ecalcalib/bmarzocc/Clustering/FourElectronsGunPt1-100_pythia8_withPU_withTracker_106X_mcRun3_2021_realistic_v3_RAW_StdSeedingGathering_Mustache_and_DeepSC_v2_showervars_Dumper_1.root'),
    
    ## base output directory: default output/
    outputDir = cms.string(''),

    ## maxEvents
    maxEvents = cms.untracked.int32(100000),

    ##fitFunction
    fitFunction = cms.string('cruijff'),
    #fitFunction = cms.string('doubleCB'),   
)
