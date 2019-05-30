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
    
    useRechits                        = cms.bool(True), #save recHits infornmation
    usePFRechits                      = cms.bool(True), #save pfRecHits infornmation
    usePFCluster                      = cms.bool(True), #save pfClusters infornmation
    useSuperCluster                   = cms.bool(True), #save superClusters infornmation
    useEnergyRegression               = cms.bool(False),#save corrected energy
    motherID                          = cms.int32(22)   #save only caloParticles with this pdgId 
    #motherID                          = cms.int32(0)   #save only caloParticles with this pdgId 
)
