import FWCore.ParameterSet.Config as cms

gendumper = cms.EDAnalyzer("GenDumper",

    genParticleCollection           = cms.InputTag("genParticles","")

)
