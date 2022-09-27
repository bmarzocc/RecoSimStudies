import FWCore.ParameterSet.Config as cms

import Geometry.CaloEventSetup.caloTowerConstituents_cfi 
CaloTowerConstituentsMapBuilder = cms.ESProducer("CaloTowerConstituentsMapBuilder",
   MapFile = cms.untracked.string('Geometry/CaloTopology/data/CaloTowerEEGeometric.map.gz')
   )

rechitdumper = cms.EDAnalyzer("RechitDumper",

    rhoCollection                   = cms.InputTag("fixedGridRhoAll","","RECO"),
    pileupSummary                   = cms.InputTag("addPileupInfo","","RECO"),
    vertexCollection                = cms.InputTag("offlinePrimaryVertices","","RECO"), 
    ecalEBRechitCollection          = cms.InputTag("ecalRecHit","EcalRecHitsEB","RECO"),
    ecalEERechitCollection          = cms.InputTag("ecalRecHit","EcalRecHitsEE","RECO"),
    ecalESRechitCollection          = cms.InputTag("ecalPreshowerRecHit","EcalRecHitsES","RECO"),
    hcalEBHERechitCollection        = cms.InputTag("hbhereco","","RECO"),
    
    isMC                            = cms.bool(False),  #isMC
    doCompression                   = cms.bool(True),  #do the compression of floats
    nBits                           = cms.int32(23),   #nbits for float compression (<=23)
    
    saveEB                          = cms.bool(False),  #save EB rechits  
    saveEE                          = cms.bool(True),  #save EE rechits   
)
