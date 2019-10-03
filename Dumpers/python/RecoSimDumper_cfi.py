import FWCore.ParameterSet.Config as cms

recosimdumper = cms.EDAnalyzer("RecoSimDumper",

    genParticleCollection             = cms.InputTag("genParticles","","HLT"),
    caloParticleCollection            = cms.InputTag("mix","MergedCaloTruth","HLT"),
    PCaloHitEBCollection              = cms.InputTag("g4SimHits","EcalHitsEB","SIM"),
    PCaloHitEECollection              = cms.InputTag("g4SimHits","EcalHitsEE","SIM"),
    ebRechitCollection                = cms.InputTag("ecalRecHit","EcalRecHitsEB","RECO"),
    eeRechitCollection                = cms.InputTag("ecalRecHit","EcalRecHitsEE","RECO"),
    pfRechitCollection                = cms.InputTag("particleFlowRecHitECAL","","RECO"),
    pfClusterCollection               = cms.InputTag("particleFlowClusterECAL","","RECO"),
    ebSuperClusterCollection          = cms.InputTag("particleFlowSuperClusterECAL","particleFlowSuperClusterECALBarrel","RECO"), 
    eeSuperClusterCollection          = cms.InputTag("particleFlowSuperClusterECAL","particleFlowSuperClusterECALEndcapWithPreshower","RECO"), 
    
    doCompression                     = cms.bool(True),  #do the compression of floats
    nBits                             = cms.int32(12),   #nbits for float compression (<=23)

    saveCalohits                      = cms.bool(False), #save pCaloHits information
    saveSimhits                       = cms.bool(True),  #save simHits information
    saveRechits                       = cms.bool(True),  #save recHits information
    savePFRechits                     = cms.bool(True),  #save pfRecHits information
    savePFCluster                     = cms.bool(True),  #save pfClusters information
    saveSuperCluster                  = cms.bool(False),  #save superClusters information
    saveShowerShapes                  = cms.bool(True),  #save saveShowerShapes information
    useEnergyRegression               = cms.bool(False), #save corrected energy
    genID                             = cms.vint32(22,11)  #save only caloParticles with this pdgId 
    #genID                            = cms.vdouble(0),  #save only caloParticles with this pdgId 
)
