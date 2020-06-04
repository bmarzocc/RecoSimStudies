import FWCore.ParameterSet.Config as cms

process = cms.PSet()

process.ioFilesOpt = cms.PSet(

    ##input file
    inputFiles = cms.string('/eos/cms/store/group/dpg_ecal/alca_ecalcalib/bmarzocc/Clustering/FourElectronsGunPt1-100_pythia8_withPU_withTracker_106X_mcRun3_2021_realistic_v3_RAW_StdSeedingGathering_Mustache_Dumper_v5_Total_1.root,/eos/cms/store/group/dpg_ecal/alca_ecalcalib/bmarzocc/Clustering/FourElectronsGunPt1-100_pythia8_withPU_withTracker_106X_mcRun3_2021_realistic_v3_RAW_StdSeedingGathering_Mustache_Dumper_v5_Total_2.root,/eos/cms/store/group/dpg_ecal/alca_ecalcalib/bmarzocc/Clustering/FourElectronsGunPt1-100_pythia8_withPU_withTracker_106X_mcRun3_2021_realistic_v3_RAW_StdSeedingGathering_Mustache_Dumper_v5_Total_3.root,/eos/cms/store/group/dpg_ecal/alca_ecalcalib/bmarzocc/Clustering/FourElectronsGunPt1-100_pythia8_withPU_withTracker_106X_mcRun3_2021_realistic_v3_RAW_StdSeedingGathering_Mustache_Dumper_v5_Total_4.root,/eos/cms/store/group/dpg_ecal/alca_ecalcalib/bmarzocc/Clustering/FourElectronsGunPt1-100_pythia8_withPU_withTracker_106X_mcRun3_2021_realistic_v3_RAW_StdSeedingGathering_Mustache_Dumper_v5_Total_5.root,/eos/cms/store/group/dpg_ecal/alca_ecalcalib/bmarzocc/Clustering/FourElectronsGunPt1-100_pythia8_withPU_withTracker_106X_mcRun3_2021_realistic_v3_RAW_StdSeedingGathering_Mustache_Dumper_v5_Total_6.root,/eos/cms/store/group/dpg_ecal/alca_ecalcalib/bmarzocc/Clustering/FourElectronsGunPt1-100_pythia8_withPU_withTracker_106X_mcRun3_2021_realistic_v3_RAW_StdSeedingGathering_Mustache_Dumper_v5_Total_7.root,/eos/cms/store/group/dpg_ecal/alca_ecalcalib/bmarzocc/Clustering/FourElectronsGunPt1-100_pythia8_withPU_withTracker_106X_mcRun3_2021_realistic_v3_RAW_StdSeedingGathering_Mustache_Dumper_v5_Total_8.root,/eos/cms/store/group/dpg_ecal/alca_ecalcalib/bmarzocc/Clustering/FourElectronsGunPt1-100_pythia8_withPU_withTracker_106X_mcRun3_2021_realistic_v3_RAW_StdSeedingGathering_Mustache_Dumper_v5_Total_9.root'),

    ## output directory
    outDir = cms.string('~/www/Clustering_simFraction_withFraction_ElectronsOnly/'),

    ## fitFunction
    fitFunction = cms.string('cruijff'),
    #fitFunction = cms.string('doubleCB'), 

    ## maxEvents
    maxEvents = cms.untracked.int32(-1),  

    ## caloParticle matching
    #matching = cms.string('sim_fraction'),
    matching = cms.string('sim_fraction_withFraction'),

    ## simEnergy cut
    simEnergyCut = cms.double(0.0),
    
    etBinning = cms.string('0 10 20 30 40 50 60 70 80 90 100'),
    etaBinning = cms.string('0.0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.479 1.75 2.0 2.25 3.0'),
    scoreBinning = cms.string('0.0001 0.00025 0.0005 0.00075 0.001 0.0015 0.002 0.0025 0.003 0.004 0.005 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 0.15 0.20 0.25 0.30 0.40 0.50 0.60 0.70 0.80' )
)
