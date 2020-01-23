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

process.gj_filter = cms.EDFilter("PythiaFilterGammaGamma",
    AcceptPrompts = cms.bool(True),
    EnergyCut = cms.double(1.0),
    EtaElThr = cms.double(2.8),
    EtaGammaThr = cms.double(2.8),
    EtaMaxCandidate = cms.double(3.0),
    EtaSeedThr = cms.double(2.8),
    EtaTkThr = cms.double(2.2),
    InvMassMax = cms.double(14000.0),
    InvMassMin = cms.double(80.0),
    NTkConeMax = cms.int32(2),
    NTkConeSum = cms.int32(4),
    PromptPtThreshold = cms.double(15.0),
    PtElThr = cms.double(2.0),
    PtGammaThr = cms.double(0.0),
    PtMinCandidate1 = cms.double(15.0),
    PtMinCandidate2 = cms.double(15.0),
    PtSeedThr = cms.double(5.0),
    PtTkThr = cms.double(1.6),
    dEtaSeedMax = cms.double(0.12),
    dPhiSeedMax = cms.double(0.2),
    dRNarrowCone = cms.double(0.02),
    dRSeedMax = cms.double(0.0),
    dRTkMax = cms.double(0.2)
)


process.generator = cms.EDFilter("Pythia8GeneratorFilter",
    PythiaParameters = cms.PSet(
        parameterSets = cms.vstring(
            'pythia8CommonSettings', 
            'pythia8CP5Settings', 
            'processParameters'
        ),
        processParameters = cms.vstring(
            'PromptPhoton:qg2qgamma = on       ! prompt photon production', 
            'PromptPhoton:qqbar2ggamma = on    ! prompt photon production', 
            'PromptPhoton:gg2ggamma = on       ! prompt photon production', 
            'PhaseSpace:pTHatMin = 20.         ! minimum pt hat for hard interactions', 
            'PhaseSpace:pTHatMax = 40.          ! maximum pt hat for hard interactions'
        ),
        pythia8CP5Settings = cms.vstring(
            'Tune:pp 14', 
            'Tune:ee 7', 
            'MultipartonInteractions:ecmPow=0.03344', 
            'MultipartonInteractions:bProfile=2', 
            'MultipartonInteractions:pT0Ref=1.41', 
            'MultipartonInteractions:coreRadius=0.7634', 
            'MultipartonInteractions:coreFraction=0.63', 
            'ColourReconnection:range=5.176', 
            'SigmaTotal:zeroAXB=off', 
            'SpaceShower:alphaSorder=2', 
            'SpaceShower:alphaSvalue=0.118', 
            'SigmaProcess:alphaSvalue=0.118', 
            'SigmaProcess:alphaSorder=2', 
            'MultipartonInteractions:alphaSvalue=0.118', 
            'MultipartonInteractions:alphaSorder=2', 
            'TimeShower:alphaSorder=2', 
            'TimeShower:alphaSvalue=0.118', 
            'SigmaTotal:mode = 0', 
            'SigmaTotal:sigmaEl = 21.89', 
            'SigmaTotal:sigmaTot = 100.309', 
            'PDF:pSet=LHAPDF6:NNPDF31_nnlo_as_0118'
        ),
        pythia8CommonSettings = cms.vstring(
            'Tune:preferLHAPDF = 2', 
            'Main:timesAllowErrors = 10000', 
            'Check:epTolErr = 0.01', 
            'Beams:setProductionScalesFromLHEF = off', 
            'SLHA:keepSM = on', 
            'SLHA:minMassSM = 1000.', 
            'ParticleDecays:limitTau0 = on', 
            'ParticleDecays:tau0Max = 10', 
            'ParticleDecays:allowPhotonRadiation = on'
        )
    ),
    comEnergy = cms.double(14000.0),
    crossSection = cms.untracked.double(1.0),
    filterEfficiency = cms.untracked.double(1.0),
    maxEventsToPrint = cms.untracked.int32(0),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    pythiaPylistVerbosity = cms.untracked.int32(0)
)

process.ProductionFilterSequence = cms.Sequence(process.generator+process.gj_filter)

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
