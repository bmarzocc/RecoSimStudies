# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: SingleGammaPt35_pythia8_cfi --conditions auto:phase1_2017_realistic -n 10 --era Run2_2017 --eventcontent FEVTDEBUG --relval 9000,50 -s GEN,SIM --datatier GEN-SIM --beamspot Realistic25ns13TeVEarly2017Collision --geometry DB:Extended --fileout file:step1.root
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('SIM',eras.Run2_2017)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.GeometrySimDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic25ns13TeVEarly2017Collision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(5)
)

# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('SingleGammaPt35_pythia8_cfi nevts:5000'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.FEVTDEBUGoutput = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    ),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:step1.root'),
    outputCommands = process.FEVTDEBUGEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
process.XMLFromDBSource.label = cms.string("Extended")
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2017_realistic', '')

#process.generator = cms.EDFilter("Pythia8PtGun",
#    PGunParameters = cms.PSet(
#        AddAntiParticle = cms.bool(True),
#        MaxEta = cms.double(2.5),
#        MaxPhi = cms.double(3.14159265359),
#        MaxPt = cms.double(35.01),
#        MinEta = cms.double(-2.5),
#        MinPhi = cms.double(-3.14159265359),
#        MinPt = cms.double(34.99),
#        ParticleID = cms.vint32(22)
#    ),
#    PythiaParameters = cms.PSet(
#        parameterSets = cms.vstring()
#    ),
#    Verbosity = cms.untracked.int32(0),
#    firstRun = cms.untracked.uint32(1),
#    psethack = cms.string('single gamma pt 35')
#)

process.generator = cms.EDProducer("CloseByParticleGunProducer",
   PGunParameters = cms.PSet(PartID = cms.vint32(22, 22),
   NParticles = cms.int32(1),
   EnMin = cms.double(1.),   # in GeV
   EnMax = cms.double(100.),
   RMin = cms.double(123.8), # in cm
   RMax = cms.double(123.8),
   ZMin = cms.double(-304.5),    # in cm
   ZMax = cms.double(304.5),
   Delta = cms.double(300),  # in cm  -> phi1-phi2 = Delta/R # for NParticles=1 irrelevant
   Pointing = cms.bool(True),# otherwise showers parallel/perpendicular to beam axis
   Overlapping = cms.bool(False),
   RandomShoot = cms.bool(False),
   MaxPhi = cms.double(3.14159265359),
   MinPhi = cms.double(-3.14159265359),
   MaxEta = cms.double(0.), # dummy, it is not used
   MinEta = cms.double(0.), # dummy, it is not used
   ),
  Verbosity = cms.untracked.int32(1),
  psethack = cms.string('two particles close to EB'),
  AddAntiParticle = cms.bool(False),
  firstRun = cms.untracked.uint32(1)
)



# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGoutput_step = cms.EndPath(process.FEVTDEBUGoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.endjob_step,process.FEVTDEBUGoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path).insert(0, process.generator)


#Setup FWK for multithreaded 
process.options.numberOfThreads=cms.untracked.uint32(8) 
process.options.numberOfStreams=cms.untracked.uint32(0) 
process.options.numberOfConcurrentLuminosityBlocks=cms.untracked.uint32(1) 


# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
