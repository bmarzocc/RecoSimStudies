# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step1 --mc --eventcontent PREMIXRAW --datatier GEN-SIM-RAW --conditions 126X_mcRun3_2023_forPU65_v1 --step DIGI,DATAMIX,L1,DIGI2RAW,HLT:2022v15 --procModifiers premix_stage2,siPixelQualityRawToDigi --nThreads 4 --geometry DB:Extended --datamix PreMix --era Run3_2023 --filein file:step1.root --fileout file:step2.root --procModifiers premix_stage2 --customise Validation/RecoParticleFlow/customize_pfanalysis.customize_ecalclustering_caloparticle
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing('standard')
options.register('nThreads',
                 1,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "nThreads")
options.register('outputName',
                 'output.root',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "outputName")
options.register('pileupName',
                 'pileup.py',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "pileupName")

options.parseArguments()
print(options)

from Configuration.Eras.Era_Run3_2023_cff import Run3_2023
from Configuration.ProcessModifiers.premix_stage2_cff import premix_stage2
from Configuration.ProcessModifiers.siPixelQualityRawToDigi_cff import siPixelQualityRawToDigi
from Configuration.ProcessModifiers.premix_stage2_cff import premix_stage2

process = cms.Process('HLT',Run3_2023,premix_stage2,siPixelQualityRawToDigi,premix_stage2)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.DigiDM_cff')
process.load('Configuration.StandardSequences.DataMixerPreMix_cff')
process.load('Configuration.StandardSequences.SimL1EmulatorDM_cff')
process.load('Configuration.StandardSequences.DigiToRawDM_cff')
process.load('HLTrigger.Configuration.HLT_2022v15_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

# Fake empty source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(
    FailPath = cms.untracked.vstring(),
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring(),
    SkipEvent = cms.untracked.vstring(),
    accelerators = cms.untracked.vstring('*'),
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
    holdsReferencesToDeleteEarly = cms.untracked.VPSet(),
    makeTriggerResults = cms.obsolete.untracked.bool,
    modulesToIgnoreForDeleteEarly = cms.untracked.vstring(),
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(0),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    #numberOfThreads = cms.untracked.uint32(1),
    numberOfThreads = cms.untracked.uint32(options.nThreads),
    printDependencies = cms.untracked.bool(False),
    sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(False)
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step1 nevts:1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.PREMIXRAWoutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-RAW'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:'+options.outputName),
    outputCommands = process.PREMIXRAWEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '126X_mcRun3_2023_forPU65_v4', '')

# Path and EndPath definitions
process.digitisation_step = cms.Path(process.pdigi)
process.datamixing_step = cms.Path(process.pdatamix)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.digi2raw_step = cms.Path(process.DigiToRaw)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.PREMIXRAWoutput_step = cms.EndPath(process.PREMIXRAWoutput)

# Schedule definition
# process.schedule imported from cff in HLTrigger.Configuration
process.schedule.insert(0, process.digitisation_step)
process.schedule.insert(1, process.datamixing_step)
process.schedule.insert(2, process.L1simulation_step)
process.schedule.insert(3, process.digi2raw_step)
process.schedule.extend([process.endjob_step,process.PREMIXRAWoutput_step])
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

#Setup FWK for multithreaded
process.options.numberOfThreads = options.nThreads
process.options.numberOfStreams = 0

# customisation of the process.
process.hltParticleFlowClusterHBHE.seedFinder.thresholdsByDetector[0].seedingThreshold=[0.6,0.5,0.5,0.5]
process.hltParticleFlowClusterHBHE.initialClusteringStep.thresholdsByDetector[0].gatheringThreshold=[0.4,0.3,0.3,0.3]
process.hltParticleFlowClusterHBHE.pfClusterBuilder.recHitEnergyNorms[0].recHitEnergyNorm=[0.4,0.3,0.3,0.3]
process.hltParticleFlowClusterHBHE.pfClusterBuilder.positionCalc.logWeightDenominatorByDetector[0].logWeightDenominator=[0.4,0.3,0.3,0.3]
process.hltParticleFlowClusterHBHE.pfClusterBuilder.allCellsPositionCalc.logWeightDenominatorByDetector[0].logWeightDenominator=[0.4,0.3,0.3,0.3]
process.hltParticleFlowClusterHCAL.pfClusterBuilder.allCellsPositionCalc.logWeightDenominatorByDetector[0].logWeightDenominator=[0.4,0.3,0.3,0.3]
process.hltParticleFlowRecHitHBHE.producers[0].qualityTests[0].cuts[0].threshold=[0.4,0.3,0.3,0.3]
process.hltEgammaHoverE.eThresHB=[0.4,0.3,0.3,0.3]
process.hltEgammaHoverEUnseeded.eThresHB=[0.4,0.3,0.3,0.3]
process.hltEgammaHToverET.eThresHB=[0.4,0.3,0.3,0.3]
process.hltFixedGridRhoFastjetECALMFForMuons.eThresHB=[0.4,0.3,0.3,0.3]
process.hltFixedGridRhoFastjetAllCaloForMuons.eThresHB=[0.4,0.3,0.3,0.3]
process.hltFixedGridRhoFastjetHCAL.eThresHB=[0.4,0.3,0.3,0.3]
process.hltParticleFlowClusterHBHECPUOnly.seedFinder.thresholdsByDetector[0].seedingThreshold=[0.6,0.5,0.5,0.5]
process.hltParticleFlowClusterHBHECPUOnly.initialClusteringStep.thresholdsByDetector[0].gatheringThreshold=[0.4,0.3,0.3,0.3]
process.hltParticleFlowClusterHBHECPUOnly.pfClusterBuilder.recHitEnergyNorms[0].recHitEnergyNorm=[0.4,0.3,0.3,0.3]
process.hltParticleFlowClusterHBHECPUOnly.pfClusterBuilder.positionCalc.logWeightDenominatorByDetector[0].logWeightDenominator=[0.4,0.3,0.3,0.3]
process.hltParticleFlowClusterHBHECPUOnly.pfClusterBuilder.allCellsPositionCalc.logWeightDenominatorByDetector[0].logWeightDenominator=[0.4,0.3,0.3,0.3]
process.hltParticleFlowClusterHCALCPUOnly.pfClusterBuilder.allCellsPositionCalc.logWeightDenominatorByDetector[0].logWeightDenominator=[0.4,0.3,0.3,0.3]
process.hltParticleFlowRecHitHBHECPUOnly.producers[0].qualityTests[0].cuts[0].threshold=[0.4,0.3,0.3,0.3]
process.hltTowerMakerForAll.HBThreshold1=0.4
process.hltTowerMakerForAll.HBThreshold2=0.3
process.hltTowerMakerForAll.HBThreshold=0.3

# Automatic addition of the customisation function from Validation.RecoParticleFlow.customize_pfanalysis
from Validation.RecoParticleFlow.customize_pfanalysis import customize_ecalclustering_caloparticle 

#call to customisation function customize_ecalclustering_caloparticle imported from Validation.RecoParticleFlow.customize_pfanalysis
process = customize_ecalclustering_caloparticle(process)

# Automatic addition of the customisation function from HLTrigger.Configuration.customizeHLTforMC
from HLTrigger.Configuration.customizeHLTforMC import customizeHLTforMC 

#call to customisation function customizeHLTforMC imported from HLTrigger.Configuration.customizeHLTforMC
process = customizeHLTforMC(process)

# End of customisation functions


# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
