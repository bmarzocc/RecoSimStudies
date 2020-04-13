import FWCore.ParameterSet.Config as cms

import Geometry.CaloEventSetup.caloTowerConstituents_cfi 
CaloTowerConstituentsMapBuilder = cms.ESProducer("CaloTowerConstituentsMapBuilder",
   MapFile = cms.untracked.string('Geometry/CaloTopology/data/CaloTowerEEGeometric.map.gz')
   )

recosimdumper = cms.EDAnalyzer("RecoSimDumper",

    rhoCollection                   = cms.InputTag("fixedGridRhoAll"),
    vertexCollection                = cms.InputTag("offlinePrimaryVertices"),
    genParticleCollection           = cms.InputTag("genParticles",""),
    caloParticleCollection          = cms.InputTag("mix","MergedCaloTruth"),
    ebRechitCollection              = cms.InputTag("ecalRecHit","EcalRecHitsEB","RECO"),
    eeRechitCollection              = cms.InputTag("ecalRecHit","EcalRecHitsEE","RECO"),
    pfRechitCollection              = cms.InputTag("particleFlowRecHitECAL","","RECO"),
    pfClusterCollection             = cms.InputTag("particleFlowClusterECAL","","RECO"),
    ebSuperClusterCollection        = cms.InputTag("particleFlowSuperClusterECAL","particleFlowSuperClusterECALBarrel","RECO"), 
    eeSuperClusterCollection        = cms.InputTag("particleFlowSuperClusterECAL","particleFlowSuperClusterECALEndcapWithPreshower","RECO"), 
    useHcalTowers                   = cms.bool(True),  #compute HoE
    hcalTowersCollection            = cms.InputTag("towerMaker"), #compute HoE
    useRetunedSC                    = cms.bool(True),  #run on new RetunedSCs
    useDeepSC                       = cms.bool(True),  #run on new DeepSCs
    ebRetunedSuperClusterCollection = cms.InputTag("particleFlowSuperClusterECALNewParams","particleFlowSuperClusterECALBarrelMustacheNewParams","RECO"), 
    eeRetunedSuperClusterCollection = cms.InputTag("particleFlowSuperClusterECALNewParams","particleFlowSuperClusterECALEndcapWithPreshowerMustacheNewParams","RECO"),
    ebDeepSuperClusterCollection    = cms.InputTag("particleFlowDeepSuperClusterECAL","particleFlowDeepSuperClusterECALBarrel","RECO"), 
    eeDeepSuperClusterCollection    = cms.InputTag("particleFlowDeepSuperClusterECAL","particleFlowDeepSuperClusterECALEndcapWithPreshower","RECO"),
    ebDeepSuperClusterLWPCollection = cms.InputTag("particleFlowDeepSuperClusterLWPECAL","particleFlowDeepSuperClusterLWPECALBarrel","RECO"), 
    eeDeepSuperClusterLWPCollection = cms.InputTag("particleFlowDeepSuperClusterLWPECAL","particleFlowDeepSuperClusterLWPECALEndcapWithPreshower","RECO"),
    ebDeepSuperClusterTWPCollection = cms.InputTag("particleFlowDeepSuperClusterTWPECAL","particleFlowDeepSuperClusterTWPECALBarrel","RECO"), 
    eeDeepSuperClusterTWPCollection = cms.InputTag("particleFlowDeepSuperClusterTWPECAL","particleFlowDeepSuperClusterTWPECALEndcapWithPreshower","RECO"),
    
    doCompression                   = cms.bool(True),  #do the compression of floats
    nBits                           = cms.int32(23),   #nbits for float compression (<=23)
    
    saveGenParticles                = cms.bool(True),  #save genParticles information   
    saveCaloParticles               = cms.bool(True),  #save caloParticles information
    saveSimhits                     = cms.bool(True), #save simHits information
    saveRechits                     = cms.bool(True), #save recHits information
    savePFRechits                   = cms.bool(True), #save pfRecHits information
    savePFCluster                   = cms.bool(True),  #save pfClusters information
    savePFClusterhits               = cms.bool(True), #save pfClustershits information
    saveSuperCluster                = cms.bool(True),  #save superClusters information
    saveShowerShapes                = cms.bool(True),  #save showerShapes information

    #scoreType                       = cms.string("simScore_final_combination"),  #score to be used for caloParticle matching
    scoreType                      = cms.string("sim_fraction"),  #score to be used for caloParticle matching
    #scoreType                      = cms.string("n_shared_xtals"),  #score to be used for caloParticle matching
    #scoreType                      = cms.string("sim_fraction_1MeVCut"),  #score to be used for caloParticle matching
    #scoreType                      = cms.string("sim_fraction_5MeVCut"),  #score to be used for caloParticle matching
    #scoreType                      = cms.string("sim_fraction_10MeVCut"),  #score to be used for caloParticle matching
    #scoreType                      = cms.string("sim_fraction_50MeVCut"),  #score to be used for caloParticle matching
    #scoreType                      = cms.string("sim_fraction_100MeVCut"),  #score to be used for caloParticle matching
    #scoreType                      = cms.string("sim_fraction_500MeVCut"),  #score to be used for caloParticle matching
    #scoreType                      = cms.string("sim_fraction_1GeVCut"),  #score to be used for caloParticle matching
    #scoreType                      = cms.string("sim_rechit_diff"),  #score to be used for caloParticle matching
    #scoreType                      = cms.string("sim_rechit_fraction"),  #score to be used for caloParticle matching
    #scoreType                      = cms.string("global_sim_rechit_fraction"),  #score to be used for caloParticle matching
    #scoreType                      = cms.string("hgcal_caloToCluster"),  #score to be used for caloParticle matching
    #scoreType                      = cms.string("hgcal_clusterToCalo"),  #score to be used for caloParticle matching
    #scoreType                      = cms.string("sim_rechit_combined_fraction"),  #score to be used for caloParticle matching
    #scoreType                      = cms.string("rechit_sim_combined_fraction"),  #score to be used for caloParticle matching
   
    saveScores                      = cms.bool(False),  #save scores information
    genID                           = cms.vint32(22,11, -11), #save only caloParticles with this pdgId 
    #genID                          = cms.vdouble(0),  #save only caloParticles with this pdgId 

)
