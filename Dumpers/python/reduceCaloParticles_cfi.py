import FWCore.ParameterSet.Config as cms

signalCaloParticles = cms.EDProducer('ReducedCaloParticleProducer',
                               caloParticleTag=cms.InputTag("mix","MergedCaloTruth"),
                               simClusterTag=cms.InputTag("mix","MergedCaloTruth"),
                               ##Parameters    
                               saveSignal = cms.bool(True), 
                               savePU = cms.bool(False), 
                               saveOOTPU = cms.bool(False)                                              
                              )

reducedCaloParticlesPU = cms.EDProducer('ReducedCaloParticleProducer',
                               caloParticleTag=cms.InputTag("mix","MergedCaloTruth"),
                               simClusterTag=cms.InputTag("mix","MergedCaloTruth"),
                               ##Parameters    
                               saveSignal = cms.bool(False), 
                               savePU = cms.bool(True), 
                               saveOOTPU = cms.bool(False)                                              
                              )

reducedCaloParticlesOOTPU = cms.EDProducer('ReducedCaloParticleProducer',
                               caloParticleTag=cms.InputTag("mix","MergedCaloTruth"),
                               simClusterTag=cms.InputTag("mix","MergedCaloTruth"),
                               ##Parameters    
                               saveSignal = cms.bool(False), 
                               savePU = cms.bool(False), 
                               saveOOTPU = cms.bool(True)                                              
                              )

reducedCaloParticlesSequence = cms.Sequence(signalCaloParticles+reducedCaloParticlesPU+reducedCaloParticlesOOTPU)
