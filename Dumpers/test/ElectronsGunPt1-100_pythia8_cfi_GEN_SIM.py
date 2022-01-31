# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: Configuration/GenProduction/python/EGM-Run3Summer21GS-00041-fragment.py --python_filename EGM-Run3Summer21GS-00041_1_cfg.py --eventcontent RAWSIM --customise Configuration/DataProcessing/Utils.addMonitoring --datatier GEN-SIM --fileout file:EGM-Run3Summer21GS-00041.root --conditions 120X_mcRun3_2021_realistic_v5 --beamspot Run3RoundOptics25ns13TeVLowSigmaZ --customise_commands process.g4SimHits.Physics.G4GeneralProcess = cms.bool(False) --step GEN,SIM --geometry DB:Extended --era Run3
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing('standard')
options.register('jobid',
                 1,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                "jobid")
options.register('seed1',
                 1,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                "Seed for generation")
options.register('seed2',
                 2,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                "Seed for generation")
options.register('seed3',
                 3,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                "Seed for generation")
options.register('seed4',
                 4,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                "Seed for generation")
options.register('seed5',
                 5,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                "Seed for generation")
options.register('outputFile',
                 'file:step1.root',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                "outputFile")

options.parseArguments()
#print options

from Configuration.Eras.Era_Run3_cff import Run3

process = cms.Process('SIM',Run3)

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
process.load('IOMC.EventVertexGenerators.VtxSmearedRun3RoundOptics25ns13TeVLowSigmaZ_cfi')
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
    FailPath = cms.untracked.vstring(),
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring(),
    SkipEvent = cms.untracked.vstring(),
    allowUnscheduled = cms.obsolete.untracked.bool,
    canDeleteEarly = cms.untracked.vstring(),
    deleteNonConsumedUnscheduledModules = cms.untracked.bool(True),
    dumpOptions = cms.untracked.bool(False),
    emptyRunLumiMode = cms.obsolete.untracked.string,
    eventSetup = cms.untracked.PSet(
        forceNumberOfConcurrentIOVs = cms.untracked.PSet(
            allowAnyLabel_=cms.required.untracked.uint32
        ),
        numberOfConcurrentIOVs = cms.untracked.uint32(0)
    ),
    fileMode = cms.untracked.string('FULLMERGE'),
    forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
    makeTriggerResults = cms.obsolete.untracked.bool,
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(0),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(1),
    printDependencies = cms.untracked.bool(False),
    sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(False)
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('Electron ParticleGun'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.1 $')
)

# Output definition

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    ),
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(1),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM'),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(20971520),
    fileName = cms.untracked.string(options.outputFile),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
import SimGeneral.Configuration.ThrowAndSetRandomRun as ThrowAndSetRandomRun
ThrowAndSetRandomRun.throwAndSetRandomRun(process.source,[(options.jobid,1)])

process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",

    LHCTransport = cms.PSet(
       initialSeed = cms.untracked.uint32(options.seed1),
       engineName = cms.untracked.string('TRandom3')
    ),
    externalLHEProducer = cms.PSet(
       initialSeed = cms.untracked.uint32(options.seed2),
       engineName = cms.untracked.string('HepJamesRandom')
    ),
    generator = cms.PSet(
       initialSeed = cms.untracked.uint32(options.seed3),
       engineName = cms.untracked.string('HepJamesRandom')
    ),
    VtxSmeared = cms.PSet(
       initialSeed = cms.untracked.uint32(options.seed4),
       engineName = cms.untracked.string('HepJamesRandom')
    ),
    g4SimHits = cms.PSet(
       initialSeed = cms.untracked.uint32(options.seed5),
       engineName = cms.untracked.string('HepJamesRandom')
    )

)

if hasattr(process, "XMLFromDBSource"): process.XMLFromDBSource.label="Extended"
if hasattr(process, "DDDetectorESProducerFromDB"): process.DDDetectorESProducerFromDB.label="Extended"
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '120X_mcRun3_2021_realistic_v5', '')

process.generator = cms.EDProducer("ManyParticleFlatEtGunProducer",
                                   PGunParameters = cms.PSet(
                                      PartID = cms.vint32(11,-11),
                                      PtMin = cms.vdouble(1, 1, 1, 1),  
                                      PtMax = cms.vdouble(100, 100, 100, 100), 
                                      PhiMin = cms.vdouble(-3.14159265359, -3.14159265359, -3.14159265359, -3.14159265359),
                                      PhiMax = cms.vdouble(3.14159265359, 3.14159265359, 3.14159265359, 3.14159265359),
                                      EtaMin = cms.vdouble(-3., -1.479, 0., 1.479),    
                                      EtaMax = cms.vdouble(-1.479, 0, 1.479, 3),   
                                      MaxPhi = cms.double(0.), #not used
                                      MinPhi = cms.double(0.), #not used
                                      MaxEta = cms.double(0.), #not used
                                      MinEta = cms.double(0.), #not used
                                   ),
                                   Verbosity = cms.untracked.int32(0),
                                   AddAntiParticle = cms.bool(False),
                                   psethack = cms.string('4 random electrons'),
                                   firstRun = cms.untracked.uint32(1)
)

process.ProductionFilterSequence = cms.Sequence(process.generator)

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step = cms.EndPath(process.RAWSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.endjob_step,process.RAWSIMoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path).insert(0, process.ProductionFilterSequence)

# customisation of the process.

# Automatic addition of the customisation function from Configuration.DataProcessing.Utils
from Configuration.DataProcessing.Utils import addMonitoring 

#call to customisation function addMonitoring imported from Configuration.DataProcessing.Utils
process = addMonitoring(process)

# End of customisation functions


# Customisation from command line

process.g4SimHits.Physics.G4GeneralProcess = cms.bool(False)
# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
