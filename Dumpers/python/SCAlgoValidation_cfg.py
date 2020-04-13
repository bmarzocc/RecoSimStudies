import FWCore.ParameterSet.Config as cms

process = cms.PSet()

process.ioFilesOpt = cms.PSet(

    ##input file
    inputFiles = cms.string('/eos/cms/store/group/dpg_ecal/alca_ecalcalib/bmarzocc/Clustering/FourElectronsGunPt1-100_pythia8_withPU_withTracker_106X_mcRun3_2021_realistic_v3_RAW_StdSeedingGathering_Mustache_RetunedMustache_DeepSC_0954_0898_0523_v3_showervars_Dumper_Test.root'),
    
    ## base output directory: default output/
    outputDir = cms.string(''),

    ## maxEvents
    maxEvents = cms.untracked.int32(250000),

    ##fitFunction
    fitFunction = cms.string('cruijff'),
    #fitFunction = cms.string('doubleCB'),   
 
    #superClusterVal = cms.string('superCluster'),
    #superClusterVal = cms.string('retunedSuperCluster'),
    #superClusterVal = cms.string('deepSuperCluster'),
    #superClusterVal = cms.string('deepSuperClusterLWP'),
    superClusterVal = cms.string('deepSuperClusterTWP'),

    superClusterRef = cms.string('superCluster'),
    #superClusterRef = cms.string('retunedSuperCluster'),
    #superClusterRef = cms.string('deepSuperCluster'),
    #superClusterRef = cms.string('deepSuperClusterLWP'),
    #superClusterRef = cms.string('deepSuperClusterTWP'),
)

process.histOpt = cms.PSet(
    nVtxBins_Barrel = cms.vdouble(60,0.,120), 
    nVtxBins_Endcap = cms.vdouble(60,0.,120),  
    RhoBins_Barrel = cms.vdouble(40,0.,80.), 
    RhoBins_Endcap = cms.vdouble(40,0.,80.),  
    nPFClustersBins = cms.vdouble(20,0.,20.), 
    nPFClustersBins_Barrel = cms.vdouble(20,0.,20.), 
    nPFClustersBins_Endcap = cms.vdouble(20,0.,20.),  
    EnergyBins_Barrel = cms.vdouble(150,0.,300.), 
    EnergyBins_Endcap = cms.vdouble(100,0.,1000.), 
    EoEtrueBins_Barrel = cms.vdouble(150,0.5,1.5), 
    EoEtrueBins_Endcap = cms.vdouble(150,0.5,1.5), 
    EoEgenBins_Barrel = cms.vdouble(150,0.5,1.5), 
    EoEgenBins_Endcap = cms.vdouble(150,0.5,1.5), 
    EtBins_Barrel = cms.vdouble(50,0.,100.), 
    EtBins_Endcap = cms.vdouble(50,0.,100.), 
    EtaBins = cms.vdouble(100,-3.2,3.2), 
    PhiBins_Barrel = cms.vdouble(100,-3.3,3.3), 
    PhiBins_Endcap = cms.vdouble(100,-3.3,3.3),  
    EtaWidthBins_Barrel = cms.vdouble(100,0.,0.1), 
    EtaWidthBins_Endcap = cms.vdouble(100,0.,0.1), 
    PhiWidthBins_Barrel = cms.vdouble(100,0.,0.5), 
    PhiWidthBins_Endcap = cms.vdouble(100,0.,0.5), 
    R9Bins_Barrel = cms.vdouble(200,0.,1.2),
    R9Bins_Endcap = cms.vdouble(200,0.,1.2), 
    full5x5_R9Bins_Barrel = cms.vdouble(200,0.,1.2),
    full5x5_R9Bins_Endcap = cms.vdouble(200,0.,1.2), 
    SigmaIetaIetaBins_Barrel = cms.vdouble(100,0.,0.02),
    SigmaIetaIetaBins_Endcap = cms.vdouble(400,0.,0.08), 
    full5x5_SigmaIetaIetaBins_Barrel = cms.vdouble(100,0.,0.02),
    full5x5_SigmaIetaIetaBins_Endcap = cms.vdouble(400,0.,0.08), 
    SigmaIetaIphiBins_Barrel = cms.vdouble(200,0.,0.01),
    SigmaIetaIphiBins_Endcap = cms.vdouble(400,0.,0.05), 
    full5x5_SigmaIetaIphiBins_Barrel = cms.vdouble(200,0.,0.01),
    full5x5_SigmaIetaIphiBins_Endcap = cms.vdouble(400,0.,0.05), 
    SigmaIphiIphiBins_Barrel = cms.vdouble(200,0.,0.04),
    SigmaIphiIphiBins_Endcap = cms.vdouble(400,0.,0.08), 
    full5x5_SigmaIphiIphiBins_Barrel = cms.vdouble(200,0.,0.04),
    full5x5_SigmaIphiIphiBins_Endcap = cms.vdouble(400,0.,0.08), 
)

