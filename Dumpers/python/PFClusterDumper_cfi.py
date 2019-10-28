import FWCore.ParameterSet.Config as cms

pfclusterdumper = cms.EDAnalyzer("PFClusterDumper",

    genParticleCollection             = cms.InputTag("genParticles",""),
    caloParticleCollection            = cms.InputTag("mix","MergedCaloTruth"),
    ebRechitCollection                = cms.InputTag("ecalRecHit","EcalRecHitsEB","RECO"),
    eeRechitCollection                = cms.InputTag("ecalRecHit","EcalRecHitsEE","RECO"),
    pfRechitCollection                = cms.InputTag("particleFlowRecHitECAL","","RECO"),
    pfClusterCollection               = cms.InputTag("particleFlowClusterECAL","","RECO"),
    
    doCompression                     = cms.bool(True),  #do the compression of floats
    nBits                             = cms.int32(12),   #nbits for float compression (<=23)

    saveHitsPosition                  = cms.bool(True),  #save hits information 
    useES                             = cms.bool(False), #use ES information in position computation
    genID                             = cms.vint32(22,11),  #save only caloParticles with this pdgId 
    #genID                            = cms.vdouble(0),  #save only caloParticles with this pdgId 
    saveScores                        = cms.bool(False) #save genParticles and CaloParticles matching scores
)
