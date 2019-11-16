# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: SingleGammaPt35_pythia8_cfi --conditions auto:phase1_2017_realistic -n 10 --era Run2_2017 --eventcontent FEVTDEBUG --relval 9000,50 -s GEN,SIM --datatier GEN-SIM --beamspot Realistic25ns13TeVEarly2017Collision --geometry DB:Extended --fileout file:step1.root
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing('standard')
options.register('jobid',
                 1,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                "jobid")
options.register('seed1',
                 0,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                "Seed for generation")
options.register('seed2',
                 0,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                "Seed for generation")
options.register('seed3',
                 0,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                "Seed for generation")
options.register('seed4',
                 0,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                "Seed for generation")

options.parseArguments()
print options

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
    input = cms.untracked.int32(options.maxEvents)
)

# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('SingleGammaPt35_pythia8_cfi nevts:10'),
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

#Other statements
import SimGeneral.Configuration.ThrowAndSetRandomRun as ThrowAndSetRandomRun
ThrowAndSetRandomRun.throwAndSetRandomRun(process.source,[(options.jobid,1)])

process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",

    externalLHEProducer = cms.PSet(
       initialSeed = cms.untracked.uint32(options.seed1),
       engineName = cms.untracked.string('HepJamesRandom')
    ),
    generator = cms.PSet(
       initialSeed = cms.untracked.uint32(options.seed2),
       engineName = cms.untracked.string('HepJamesRandom')
    ),
    VtxSmeared = cms.PSet(
       initialSeed = cms.untracked.uint32(options.seed3),
       engineName = cms.untracked.string('HepJamesRandom')
    ),
    g4SimHits = cms.PSet(
       initialSeed = cms.untracked.uint32(options.seed4),
       engineName = cms.untracked.string('HepJamesRandom')
    )

)


# Additional output definition

# Other statements
process.XMLFromDBSource.label = cms.string("Extended")
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_mc2017_realistic_v6', '')

process.generator = cms.EDProducer("ManyParticleFlatEtGunProducer",
                                   PGunParameters = cms.PSet(
                                      #PartID = cms.vint32(22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22),
                                      PartID = cms.vint32(11, -11),
                                      #PtMin = cms.vdouble(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), 
                                      PtMin = cms.vdouble(1, 1, 1, 1),  
                                      #PtMax = cms.vdouble(100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100),
                                      PtMax = cms.vdouble(100, 100, 100, 100), 
                                      #PhiMin = cms.vdouble(-3.14159265359, -3.14159265359, -3.14159265359, -3.14159265359, -3.14159265359, -3.14159265359, -3.14159265359, -3.14159265359, -3.14159265359, -3.14159265359, -3.14159265359, -3.14159265359, -3.14159265359, -3.14159265359, -3.14159265359, -3.14159265359, -3.14159265359, -3.14159265359, -3.14159265359, -3.14159265359),
                                      PhiMin = cms.vdouble(-3.14159265359, -3.14159265359, -3.14159265359, -3.14159265359),
                                      #PhiMax = cms.vdouble(3.14159265359, 3.14159265359, 3.14159265359, 3.14159265359, 3.14159265359, 3.14159265359, 3.14159265359, 3.14159265359, 3.14159265359, 3.14159265359, 3.14159265359, 3.14159265359, 3.14159265359, 3.14159265359, 3.14159265359, 3.14159265359, 3.14159265359, 3.14159265359, 3.14159265359, 3.14159265359), 
                                      PhiMax = cms.vdouble(3.14159265359, 3.14159265359, 3.14159265359, 3.14159265359),
                                      #EtaMin = cms.vdouble(0., 0., 0., 0., 0., -1.479, -1.479, -1.479, -1.479, -1.479, 1.479, 1.479, 1.479, 1.479, 1.479, -3., -3, -3., -3., -3.), 
                                      EtaMin = cms.vdouble(-3., -1.479, 0., 1.479),    
                                      #EtaMax = cms.vdouble(1.479, 1.479, 1.479, 1.479, 1.479, 0., 0., 0., 0., 0., 3., 3., 3., 3., 3., -1.479, -1.479, -1.479, -1.479, -1.479),
                                      EtaMax = cms.vdouble(-1.479, 0, 1.479, 3),   
                                      MaxPhi = cms.double(0.), #not used
                                      MinPhi = cms.double(0.), #not used
                                      MaxEta = cms.double(0.), #not used
                                      MinEta = cms.double(0.), #not used
                                   ),
                                   Verbosity = cms.untracked.int32(0),
                                   AddAntiParticle = cms.bool(False),
                                   psethack = cms.string('20 random photons'),
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


# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
