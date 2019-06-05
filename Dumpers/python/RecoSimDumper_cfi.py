import FWCore.ParameterSet.Config as cms

recosimdumper = cms.EDAnalyzer("RecoSimDumper",

    genParticleCollection             = cms.InputTag("genParticles","","HLT"),
    caloParticleCollection            = cms.InputTag("mix","MergedCaloTruth","HLT"),
    ebRechitCollection                = cms.InputTag("ecalRecHit","EcalRecHitsEB","RECO"),
    eeRechitCollection                = cms.InputTag("ecalRecHit","EcalRecHitsEE","RECO"),
    pfRechitCollection                = cms.InputTag("particleFlowRecHitECAL","","RECO"),
    pfClusterCollection               = cms.InputTag("particleFlowEGamma","EBEEClusters","RECO"),
    ebSuperClusterCollection          = cms.InputTag("particleFlowSuperClusterECAL","particleFlowSuperClusterECALBarrel","RECO"), 
    eeSuperClusterCollection          = cms.InputTag("particleFlowSuperClusterECAL","particleFlowSuperClusterECALEndcapWithPreshower","RECO"), 
    
    doCompression                     = cms.bool(True),  #do the compression of floats
    nBits                             = cms.int32(12),   #nbits for float compression (<=23)
    saveSimhits                       = cms.bool(True),  #save simHits information
    saveRechits                       = cms.bool(True),  #save recHits information
    savePFRechits                     = cms.bool(True),  #save pfRecHits information
    savePFCluster                     = cms.bool(True),  #save pfClusters information
    saveSuperCluster                  = cms.bool(False), #save superClusters information
    useEnergyRegression               = cms.bool(False), #save corrected energy
    motherID                          = cms.int32(22),   #save only caloParticles with this pdgId 
    #motherID                          = cms.int32(0),   #save only caloParticles with this pdgId 
)
