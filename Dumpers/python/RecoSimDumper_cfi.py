import FWCore.ParameterSet.Config as cms

import Geometry.CaloEventSetup.caloTowerConstituents_cfi 
CaloTowerConstituentsMapBuilder = cms.ESProducer("CaloTowerConstituentsMapBuilder",
   MapFile = cms.untracked.string('Geometry/CaloTopology/data/CaloTowerEEGeometric.map.gz')
   )

recosimdumper = cms.EDAnalyzer("RecoSimDumper",

    rhoCollection                   = cms.InputTag("fixedGridRhoAll","","ECALClustering"),
    pileupSummary                   = cms.InputTag("addPileupInfo",""),
    #vertexCollection                = cms.InputTag("offlinePrimaryVertices","","ECALClustering"),
    vertexCollection                = cms.InputTag("offlineSlimmedPrimaryVertices","","ECALClustering"), 
    genParticleCollection           = cms.InputTag("genParticles",""),
    caloParticleCollection          = cms.InputTag("signalCaloParticles",""),
    puCaloParticleCollection        = cms.InputTag("reducedCaloParticlesPU",""),
    ootpuCaloParticleCollection     = cms.InputTag("reducedCaloParticlesOOTPU",""),
    ebRechitCollection              = cms.InputTag("ecalRecHit","EcalRecHitsEB","ECALClustering"),
    eeRechitCollection              = cms.InputTag("ecalRecHit","EcalRecHitsEE","ECALClustering"),
    pfRechitCollection              = cms.InputTag("particleFlowRecHitECAL","","ECALClustering"),
    pfClusterCollection             = cms.InputTag("particleFlowClusterECAL","","ECALClustering"),
    ebSuperClusterCollection        = cms.InputTag("particleFlowSuperClusterECAL","particleFlowSuperClusterECALBarrel","ECALClustering"), 
    eeSuperClusterCollection        = cms.InputTag("particleFlowSuperClusterECAL","particleFlowSuperClusterECALEndcapWithPreshower","ECALClustering"), 
    gsfElectronCollection           = cms.InputTag("gedGsfElectrons","","ECALClustering"), 
    gedPhotonCollection             = cms.InputTag("gedPhotons","","ECALClustering"), 
    patElectronCollection           = cms.InputTag("slimmedElectrons","","ECALClustering"), 
    patPhotonCollection             = cms.InputTag("slimmedPhotons","","ECALClustering"), 
    patJetCollection                = cms.InputTag("slimmedJets","","ECALClustering"), 
    patMETCollection                = cms.InputTag("slimmedMETs","","ECALClustering"), 
    ebRetunedSuperClusterCollection = cms.InputTag("particleFlowSuperClusterECALNewParams","particleFlowSuperClusterECALBarrelMustacheNewParams","ECALClustering"), 
    eeRetunedSuperClusterCollection = cms.InputTag("particleFlowSuperClusterECALNewParams","particleFlowSuperClusterECALEndcapWithPreshowerMustacheNewParams","ECALClustering"),
    ebDeepSuperClusterCollection    = cms.InputTag("particleFlowDeepSuperClusterECAL","particleFlowDeepSuperClusterECALBarrel","ECALClustering"), 
    eeDeepSuperClusterCollection    = cms.InputTag("particleFlowDeepSuperClusterECAL","particleFlowDeepSuperClusterECALEndcapWithPreshower","ECALClustering"),
    
    isMC                            = cms.bool(True),  #isMC
    doCompression                   = cms.bool(True),  #do the compression of floats
    nBits                           = cms.int32(23),   #nbits for float compression (<=23)
    
    #MC-only info (turned off if isMC == False)
    saveGenParticles                = cms.bool(True),  #save genParticles information   
    saveCaloParticles               = cms.bool(False),  #save caloParticles information
    saveCaloParticlesPU             = cms.bool(False),  #save PU caloParticles information
    saveCaloParticlesOOTPU          = cms.bool(False),  #save OOT PU caloParticles information
    subtractSignalCalo              = cms.bool(False),  #subtract signal caloParticle to PU caloParticle
    saveSimhits                     = cms.bool(False), #save simHits information
    saveSimhitsPU                   = cms.bool(False), #save simHits of PU information

    #Standard info
    saveRechits                     = cms.bool(False), #save recHits information
    savePFRechits                   = cms.bool(False), #save pfRecHits information
    savePFCluster                   = cms.bool(True),  #save pfClusters information
    savePFClusterhits               = cms.bool(False), #save pfClustershits information
    saveShowerShapes                = cms.bool(True),  #save showerShapes information
    saveSuperCluster                = cms.bool(True),  #save superClusters information
    saveRetunedSC                   = cms.bool(False),  #save additional retunedSCs information (missing from the standard RECO)
    saveDeepSC                      = cms.bool(False),  #save additional deepSCs information (missing from the standard RECO)
    saveGsfElectrons                = cms.bool(False),  #save gedGsfElectrons information
    saveGedPhotons                  = cms.bool(False),  #save gedPhotons information 
    savePatPhotons                  = cms.bool(True),  #save patPhotons and patMET information
    savePatElectrons                = cms.bool(True),  #save patElectrons and patMET information
    savePatJets                     = cms.bool(True),  #save patJets and patMET information

    egmCutBasedElectronIDVeto       = cms.string('cutBasedElectronID-Fall17-94X-V2-veto'),  #cutBasedEleID veto
    egmCutBasedElectronIDloose      = cms.string('cutBasedElectronID-Fall17-94X-V2-loose'),  #cutBasedEleID loose  
    egmCutBasedElectronIDmedium     = cms.string('cutBasedElectronID-Fall17-94X-V2-medium'),  #cutBasedEleID medium 
    egmCutBasedElectronIDtight      = cms.string('cutBasedElectronID-Fall17-94X-V2-tight'),  #cutBasedEleID tight
    egmMVAElectronIDloose           = cms.string('mvaEleID-Fall17-iso-V2-wpLoose'),  #mvaEleID loose 
    egmMVAElectronIDmedium          = cms.string('mvaEleID-Fall17-iso-V2-wp90'),  #mvaEleID medium 
    egmMVAElectronIDtight           = cms.string('mvaEleID-Fall17-iso-V2-wp80'),  #mvaEleID tight  
    egmMVAElectronIDlooseNoIso      = cms.string('mvaEleID-Fall17-noIso-V2-wpLoose'),  #mvaEleIDNoIso loose 
    egmMVAElectronIDmediumNoIso     = cms.string('mvaEleID-Fall17-noIso-V2-wp90'),  #mvaEleIDNoIso medium 
    egmMVAElectronIDtightNoIso      = cms.string('mvaEleID-Fall17-noIso-V2-wp80'),  #mvaEleIDNoIso tight  
    heepElectronID                  = cms.string('heepElectronID-HEEPV70'),  #mvaEleIDNoIso tight
    egmCutBasedPhotonIDloose        = cms.string('cutBasedPhotonID-Fall17-94X-V2-loose'),  #cutBasedPhoID loose  
    egmCutBasedPhotonIDmedium       = cms.string('cutBasedPhotonID-Fall17-94X-V2-medium'),  #cutBasedPhoID medium 
    egmCutBasedPhotonIDtight        = cms.string('cutBasedPhotonID-Fall17-94X-V2-tight'),  #cutBasedPhoID tight 
    egmMVAPhotonIDmedium            = cms.string('mvaPhoID-RunIIFall17-v2-wp90'),  #mvaPhoID medium 
    egmMVAPhotonIDtight             = cms.string('mvaPhoID-RunIIFall17-v2-wp80'),  #mvaPhoID tight   
)
