import FWCore.ParameterSet.Config as cms

import Geometry.CaloEventSetup.caloTowerConstituents_cfi 
CaloTowerConstituentsMapBuilder = cms.ESProducer("CaloTowerConstituentsMapBuilder",
   MapFile = cms.untracked.string('Geometry/CaloTopology/data/CaloTowerEEGeometric.map.gz')
   )

recosimdumper = cms.EDAnalyzer("RecoSimDumper",

    rhoCollection                   = cms.InputTag("fixedGridRhoAll"),
    pileupSummary                   = cms.InputTag("addPileupInfo"),
    vertexCollection                = cms.InputTag("offlinePrimaryVertices"),
    genParticleCollection           = cms.InputTag("genParticles",""),
    caloParticleCollection          = cms.InputTag("signalCaloParticles"),
    puCaloParticleCollection        = cms.InputTag("reducedCaloParticlesPU"),
    ootpuCaloParticleCollection     = cms.InputTag("reducedCaloParticlesOOTPU"),
    ebRechitCollection              = cms.InputTag("ecalRecHit","EcalRecHitsEB","RECO"),
    eeRechitCollection              = cms.InputTag("ecalRecHit","EcalRecHitsEE","RECO"),
    pfRechitCollection              = cms.InputTag("particleFlowRecHitECAL","","RECO"),
    pfClusterCollection             = cms.InputTag("particleFlowClusterECAL","","RECO"),
    ebSuperClusterCollection        = cms.InputTag("particleFlowSuperClusterECAL","particleFlowSuperClusterECALBarrel","RECO"), 
    eeSuperClusterCollection        = cms.InputTag("particleFlowSuperClusterECAL","particleFlowSuperClusterECALEndcapWithPreshower","RECO"), 
    useRetunedSC                    = cms.bool(False),  #run on new RetunedSCs
    useDeepSC                       = cms.bool(False),  #run on new DeepSCs
    ebRetunedSuperClusterCollection = cms.InputTag("particleFlowSuperClusterECALNewParams","particleFlowSuperClusterECALBarrelMustacheNewParams","RECO"), 
    eeRetunedSuperClusterCollection = cms.InputTag("particleFlowSuperClusterECALNewParams","particleFlowSuperClusterECALEndcapWithPreshowerMustacheNewParams","RECO"),
    ebDeepSuperClusterCollection    = cms.InputTag("particleFlowDeepSuperClusterECAL","particleFlowDeepSuperClusterECALBarrel","RECO"), 
    eeDeepSuperClusterCollection    = cms.InputTag("particleFlowDeepSuperClusterECAL","particleFlowDeepSuperClusterECALEndcapWithPreshower","RECO"),
    
    doCompression                   = cms.bool(True),  #do the compression of floats
    nBits                           = cms.int32(23),   #nbits for float compression (<=23)
    
    saveGenParticles                = cms.bool(True),  #save genParticles information   
    saveCaloParticles               = cms.bool(True),  #save caloParticles information
    saveCaloParticlesPU             = cms.bool(True),  #save PU caloParticles information
    saveCaloParticlesOOTPU          = cms.bool(False),  #save OOT PU caloParticles information
    saveSimhits                     = cms.bool(True), #save simHits information
    saveSimhitsPU                   = cms.bool(False), #save simHits of PU information
    saveRechits                     = cms.bool(True), #save recHits information
    savePFRechits                   = cms.bool(True), #save pfRecHits information
    savePFCluster                   = cms.bool(True),  #save pfClusters information
    savePFClusterhits               = cms.bool(True), #save pfClustershits information
    saveSuperCluster                = cms.bool(True),  #save superClusters information
    saveShowerShapes                = cms.bool(True),  #save showerShapes information
)
