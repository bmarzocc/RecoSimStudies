import FWCore.ParameterSet.Config as cms

process = cms.PSet()

process.ioFilesOpt = cms.PSet(

    ##input file
    inputFiles = cms.string('/eos/cms/store/group/dpg_ecal/alca_ecalcalib/bmarzocc/Clustering/FourElectronsGunPt1-100_pythia8_withPU_withTracker_106X_mcRun3_2021_realistic_v3_RAW_StdSeedingGathering_Mustache_optimizedDeepSC_v14_finalscore_v2_Total_1.root,/eos/cms/store/group/dpg_ecal/alca_ecalcalib/bmarzocc/Clustering/FourElectronsGunPt1-100_pythia8_withPU_withTracker_106X_mcRun3_2021_realistic_v3_RAW_StdSeedingGathering_Mustache_optimizedDeepSC_v14_finalscore_v2_Total_2.root,/eos/cms/store/group/dpg_ecal/alca_ecalcalib/bmarzocc/Clustering/FourElectronsGunPt1-100_pythia8_withPU_withTracker_106X_mcRun3_2021_realistic_v3_RAW_StdSeedingGathering_Mustache_optimizedDeepSC_v14_finalscore_v2_Total_3.root,/eos/cms/store/group/dpg_ecal/alca_ecalcalib/bmarzocc/Clustering/FourElectronsGunPt1-100_pythia8_withPU_withTracker_106X_mcRun3_2021_realistic_v3_RAW_StdSeedingGathering_Mustache_optimizedDeepSC_v14_finalscore_v2_Total_4.root,/eos/cms/store/group/dpg_ecal/alca_ecalcalib/bmarzocc/Clustering/FourElectronsGunPt1-100_pythia8_withPU_withTracker_106X_mcRun3_2021_realistic_v3_RAW_StdSeedingGathering_Mustache_optimizedDeepSC_v14_finalscore_v2_Total_5.root,/eos/cms/store/group/dpg_ecal/alca_ecalcalib/bmarzocc/Clustering/FourElectronsGunPt1-100_pythia8_withPU_withTracker_106X_mcRun3_2021_realistic_v3_RAW_StdSeedingGathering_Mustache_optimizedDeepSC_v14_finalscore_v2_Total_6.root,/eos/cms/store/group/dpg_ecal/alca_ecalcalib/bmarzocc/Clustering/FourElectronsGunPt1-100_pythia8_withPU_withTracker_106X_mcRun3_2021_realistic_v3_RAW_StdSeedingGathering_Mustache_optimizedDeepSC_v14_finalscore_v2_Total_7.root,/eos/cms/store/group/dpg_ecal/alca_ecalcalib/bmarzocc/Clustering/FourElectronsGunPt1-100_pythia8_withPU_withTracker_106X_mcRun3_2021_realistic_v3_RAW_StdSeedingGathering_Mustache_optimizedDeepSC_v14_finalscore_v2_Total_8.root,/eos/cms/store/group/dpg_ecal/alca_ecalcalib/bmarzocc/Clustering/FourElectronsGunPt1-100_pythia8_withPU_withTracker_106X_mcRun3_2021_realistic_v3_RAW_StdSeedingGathering_Mustache_optimizedDeepSC_v14_finalscore_v2_Total_9.root'),

    ## output directory
    outDir = cms.string('~/www/Clustering_MustacheToDeepSCComparison_Electrons/'),

    ## maxEvents
    maxEvents = cms.untracked.int32(-1),  

    etBinning = cms.string('0 10 20 30 40 50 60 70 80 90 100'),
    etaBinning = cms.string('0.0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.479 1.75 2.0 2.25 3.0')
)
