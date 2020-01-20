import FWCore.ParameterSet.Config as cms

recosimdumper = cms.EDAnalyzer("RecoSimDumper",

    rhoCollection                = cms.InputTag("fixedGridRhoAll"),
    vertexCollection             = cms.InputTag("offlinePrimaryVertices"),
    genParticleCollection        = cms.InputTag("genParticles",""),
    caloParticleCollection       = cms.InputTag("mix","MergedCaloTruth"),
    ebRechitCollection           = cms.InputTag("ecalRecHit","EcalRecHitsEB","RECO"),
    eeRechitCollection           = cms.InputTag("ecalRecHit","EcalRecHitsEE","RECO"),
    pfRechitCollection           = cms.InputTag("particleFlowRecHitECAL","","RECO"),
    pfClusterCollection          = cms.InputTag("particleFlowClusterECAL","","RECO"),
    ebSuperClusterCollection     = cms.InputTag("particleFlowSuperClusterECAL","particleFlowSuperClusterECALBarrel","RECO"), 
    eeSuperClusterCollection     = cms.InputTag("particleFlowSuperClusterECAL","particleFlowSuperClusterECALEndcapWithPreshower","RECO"), 

    useDeepSC                    = cms.bool(True),  #run on new DeepSCs
    ebDeepSuperClusterCollection = cms.InputTag("particleFlowSuperClusterECALFromDeepSC","particleFlowDeepSuperClusterECALBarrel","RECO"), 
    eeDeepSuperClusterCollection = cms.InputTag("particleFlowSuperClusterECALFromDeepSC","particleFlowDeepSuperClusterECALEndcapWithPreshower","RECO"),
    
    doCompression                = cms.bool(True),  #do the compression of floats
    nBits                        = cms.int32(23),   #nbits for float compression (<=23)
    
    saveGenParticles             = cms.bool(True),  #save genParticles information   
    saveCaloParticles            = cms.bool(True),  #save caloParticles information
    saveSimhits                  = cms.bool(False), #save simHits information
    saveRechits                  = cms.bool(False), #save recHits information
    savePFRechits                = cms.bool(False), #save pfRecHits information
    savePFCluster                = cms.bool(True),  #save pfClusters information
    savePFClusterhits            = cms.bool(False), #save pfClustershits information
    saveSuperCluster             = cms.bool(True),  #save superClusters information
    saveShowerShapes             = cms.bool(True),  #save showerShapes information

    scoreType                    = cms.string("sim_fraction_min1"),  #score to be used for caloParticle matching
    #scoreType                   = cms.string("n_shared_xtals"),  #score to be used for caloParticle matching
    #scoreType                   = cms.string("sim_fraction"),  #score to be used for caloParticle matching
    #scoreType                   = cms.string("sim_fraction_1MeVCut"),  #score to be used for caloParticle matching
    #scoreType                   = cms.string("sim_fraction_5MeVCut"),  #score to be used for caloParticle matching
    #scoreType                   = cms.string("sim_fraction_10MeVCut"),  #score to be used for caloParticle matching
    #scoreType                   = cms.string("sim_fraction_50MeVCut"),  #score to be used for caloParticle matching
    #scoreType                   = cms.string("sim_fraction_100MeVCut"),  #score to be used for caloParticle matching
    #scoreType                   = cms.string("sim_fraction_500MeVCut"),  #score to be used for caloParticle matching
    #scoreType                   = cms.string("sim_fraction_1GeVCut"),  #score to be used for caloParticle matching
    #scoreType                   = cms.string("sim_fraction_min3"),  #score to be used for caloParticle matching
    #scoreType                   = cms.string("sim_fraction_min3_1MeVCut"),  #score to be used for caloParticle matching
    #scoreType                   = cms.string("sim_fraction_min3_5MeVCut"),  #score to be used for caloParticle matching
    #scoreType                   = cms.string("sim_fraction_min3_10MeVCut"),  #score to be used for caloParticle matching
    #scoreType                   = cms.string("sim_fraction_min3_50MeVCut"),  #score to be used for caloParticle matching
    #scoreType                   = cms.string("sim_fraction_min3_100MeVCut"),  #score to be used for caloParticle matching
    #scoreType                   = cms.string("sim_rechit_diff"),  #score to be used for caloParticle matching
    #scoreType                   = cms.string("sim_rechit_fraction"),  #score to be used for caloParticle matching
    #scoreType                   = cms.string("global_sim_rechit_fraction"),  #score to be used for caloParticle matching
    #scoreType                   = cms.string("hgcal_caloToCluster"),  #score to be used for caloParticle matching
    #scoreType                   = cms.string("hgcal_clusterToCalo"),  #score to be used for caloParticle matching
   
    saveScores                   = cms.bool(True),  #save scores information
    genID                        = cms.vint32(22,11, -11), #save only caloParticles with this pdgId 
    #genID                       = cms.vdouble(0),  #save only caloParticles with this pdgId 

)
