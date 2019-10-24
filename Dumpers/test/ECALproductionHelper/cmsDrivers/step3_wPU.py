# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step3 --conditions auto:phase1_2017_realistic -n 10 --era Run2_2017 --eventcontent RECOSIM,MINIAODSIM,DQM --runUnscheduled -s RAW2DIGI,L1Reco,RECO,RECOSIM,EI,PAT,VALIDATION:@standardValidation+@miniAODValidation,DQM:@standardDQM+@ExtraHLT+@miniAODDQM --datatier GEN-SIM-RECO,MINIAODSIM,DQMIO --geometry DB:Extended --filein file:step2.root --fileout file:step3.root
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing ('analysis')
# define the defaults here, changed from command line
options.maxEvents = -1 # -1 means all events, maxEvents considers the total over files considered
# add costum parameters
options.register ("pfrhMult",
                  1.0, # default value
                  VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.varType.float,          # string, int, or float
                  "multiplier of noise used for PFRH thresholds")
options.register ("seedMult",
                  3.0, # default value
                  VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.varType.float,          # string, int, or float
                  "multiplier of noise used for seeding threshold")
options.register('doRef',
                 0,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.int,
                "Use reference values for seeding,gathering,pfrechits")
options.register('nThr',
                 1,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.int,
                "Number of threads")
                  

options.parseArguments()

import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Run3_cff import Run3
from Configuration.ProcessModifiers.premix_stage2_cff import premix_stage2

process = cms.Process('RECO',Run3,premix_stage2)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.RecoSim_cff')
process.load('CommonTools.ParticleFlow.EITopPAG_cff')
process.load('PhysicsTools.PatAlgos.slimming.metFilterPaths_cff')
process.load('Configuration.StandardSequences.PATMC_cff')
process.load('Configuration.StandardSequences.Validation_cff')
#process.load('DQMServices.Core.DQMStoreNonLegacy_cff')
#process.load('DQMOffline.Configuration.DQMOfflineMC_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)


# Message Logger settings
#process.MessageLogger = cms.Service("MessageLogger",
#         destinations = cms.untracked.vstring(
#              'detailedInfo',
#              'critical',
#              #'cerr'
#         ),
#         critical = cms.untracked.PSet(threshold = cms.untracked.string('ERROR')),
#         detailedInfo = cms.untracked.PSet(threshold = cms.untracked.string('INFO')),
#        #cerr = cms.untracked.PSet(threshold  = cms.untracked.string('WARNING')),
#)


# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:step2.root'),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step3 nevts:5000'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-RECO'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:step3.root'),
    outputCommands = process.RECOSIMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)
process.RECOSIMoutput.outputCommands.extend(['keep *_mix_MergedCaloTruth_*',
                                            #'keep *PCaloHit*_g4SimHits_EcalHitsE*_*',
                                             'keep *_particleFlowRecHitECAL_*_*'])

# Additional output definition

# Other statements
process.mix.playback = True
process.mix.digitizers = cms.PSet()
for a in process.aliases: delattr(process, a)
process.RandomNumberGeneratorService.restoreStateLabel=cms.untracked.string("randomEngineStateProducer")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_mcRun3_2021_realistic_v3', '')

# override a global tag with the conditions from external module
from CalibCalorimetry.EcalTrivialCondModules.EcalTrivialCondRetriever_cfi import *
process.myCond = EcalTrivialConditionRetriever.clone()
# prefer these conditions over the globalTag's ones
process.es_prefer = cms.ESPrefer("EcalTrivialConditionRetriever","myCond")

### set all conditions producers to false except those I am interested in
if options.doRef == 0:
  process.myCond.producedEcalPFRecHitThresholds = cms.untracked.bool(True)
  process.myCond.EcalPFRecHitThresholdNSigmas = cms.untracked.double(options.pfrhMult/2.0)
  process.myCond.EcalPFRecHitThresholdNSigmasHEta = cms.untracked.double(options.pfrhMult/3.0)
  process.myCond.PFRecHitFile = cms.untracked.string("./data/noise/PFRecHitThresholds_EB.txt")
  process.myCond.PFRecHitFileEE = cms.untracked.string("./data/noise/PFRecHitThresholds_EE.txt")
else: # use the reference values
  process.myCond.producedEcalPFRecHitThresholds = cms.untracked.bool(False)

if options.doRef == 0:
  process.myCond.producedEcalPFSeedingThresholds = cms.untracked.bool(True)
  process.myCond.EcalPFSeedingThresholdNSigmas = cms.untracked.double(options.seedMult/2.0) # PFRHs files are at 2sigma of the noise for |eta|<2.5
  process.myCond.EcalPFSeedingThresholdNSigmasHEta = cms.untracked.double(options.seedMult/3.0) #                3sigma of the noise for |eta|>2.5
  process.myCond.PFSeedingFile = cms.untracked.string("./data/noise/PFRecHitThresholds_EB.txt")
  process.myCond.PFSeedingFileEE = cms.untracked.string("./data/noise/PFRecHitThresholds_EE.txt")
else: # use the reference values
  process.myCond.producedEcalPFSeedingThresholds = cms.untracked.bool(True)
  process.myCond.EcalPFSeedingThresholdNSigmas = cms.untracked.double(1.0) 
  process.myCond.EcalPFSeedingThresholdNSigmasHEta = cms.untracked.double(1.0) 
  process.myCond.PFSeedingFile = cms.untracked.string("./data/noise/fixed_SeedingThresholds_EB.txt")
  process.myCond.PFSeedingFileEE = cms.untracked.string("./data/noise/fixed_SeedingThresholds_EE.txt")

process.myCond.producedEcalPedestals = cms.untracked.bool(False)
process.myCond.producedEcalWeights = cms.untracked.bool(False)
process.myCond.producedEcalGainRatios = cms.untracked.bool(False)
process.myCond.producedEcalADCToGeVConstant = cms.untracked.bool(False)
process.myCond.producedEcalMappingElectronics = cms.untracked.bool(False)
process.myCond.producedEcalTimeOffsetConstant = cms.untracked.bool(False)
process.myCond.producedEcalLinearCorrections = cms.untracked.bool(False)
process.myCond.producedEcalIntercalibConstants = cms.untracked.bool(False)
process.myCond.producedEcalIntercalibConstantsMC = cms.untracked.bool(False)
process.myCond.producedEcalIntercalibErrors = cms.untracked.bool(False)
process.myCond.producedEcalTimeCalibConstants = cms.untracked.bool(False)
process.myCond.producedEcalTimeCalibErrors = cms.untracked.bool(False)
process.myCond.producedEcalSimPulseShape = cms.untracked.bool(False)
process.myCond.producedEcalChannelStatus = cms.untracked.bool(False)
process.myCond.producedEcalDQMChannelStatus = cms.untracked.bool(False)
process.myCond.producedEcalDCSTowerStatus = cms.untracked.bool(False)
process.myCond.producedEcalDAQTowerStatus = cms.untracked.bool(False)
process.myCond.producedEcalDQMTowerStatus = cms.untracked.bool(False)
process.myCond.producedEcalTrgChannelStatus = cms.untracked.bool(False)
process.myCond.producedEcalAlignmentEB = cms.untracked.bool(False)
process.myCond.producedEcalAlignmentEE = cms.untracked.bool(False)
process.myCond.producedEcalAlignmentEE = cms.untracked.bool(False)
process.myCond.producedEcalSampleMask = cms.untracked.bool(False)
# additional booleans to excplicitly set to False
process.myCond.producedEcalTimeBiasCorrections = cms.untracked.bool(False)
process.myCond.producedEcalSamplesCorrelation = cms.untracked.bool(False)
process.myCond.producedEcalLaserCorrection = cms.untracked.bool(False)
process.myCond.producedEcalClusterLocalContCorrParameters = cms.untracked.bool(False)
process.myCond.producedEcalClusterCrackCorrParameters = cms.untracked.bool(False)
process.myCond.producedEcalClusterEnergyCorrectionParameters = cms.untracked.bool(False)
process.myCond.producedEcalClusterEnergyUncertaintyParameters = cms.untracked.bool(False)
process.myCond.producedEcalClusterEnergyCorrectionObjectSpecificParameters = cms.untracked.bool(False)

########## end override

# Path and EndPath definitions 
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.recosim_step = cms.Path(process.recosim)
process.eventinterpretaion_step = cms.Path(process.EIsequence)
process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)

process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.recosim_step,
                            process.eventinterpretaion_step,  
                            process.RECOSIMoutput_step)
process.schedule.associate(process.patTask)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)


#Setup FWK for multithreaded 
process.options.numberOfThreads=cms.untracked.uint32(options.nThr)
process.options.numberOfStreams=cms.untracked.uint32(0)
process.options.numberOfConcurrentLuminosityBlocks=cms.untracked.uint32(1)


# customisation of the process.
if options.doRef == 0:
  process.particleFlowClusterECALUncorrected.initialClusteringStep.thresholdsByDetector = cms.VPSet(
        cms.PSet( detector = cms.string("ECAL_BARREL"),
               gatheringThreshold = cms.double(0.0),
               gatheringThresholdPt = cms.double(0.0)
               ),
        cms.PSet( detector = cms.string("ECAL_ENDCAP"),
               gatheringThreshold = cms.double(0.0),
               gatheringThresholdPt = cms.double(0.0)
               )
  )
else: # use reference values
  process.particleFlowClusterECALUncorrected.initialClusteringStep.thresholdsByDetector = cms.VPSet(
        cms.PSet( detector = cms.string("ECAL_BARREL"),
               gatheringThreshold = cms.double(0.08),
               gatheringThresholdPt = cms.double(0.0)
               ),
        cms.PSet( detector = cms.string("ECAL_ENDCAP"),
               gatheringThreshold = cms.double(0.30),
               gatheringThresholdPt = cms.double(0.0)
               )
  )

# Automatic addition of the customisation function from SimGeneral.MixingModule.fullMixCustomize_cff
from SimGeneral.MixingModule.fullMixCustomize_cff import setCrossingFrameOn 

#call to customisation function setCrossingFrameOn imported from SimGeneral.MixingModule.fullMixCustomize_cff
process = setCrossingFrameOn(process)

# End of customisation functions
#do not add changes to your config after this point (unless you know what you are doing)
from FWCore.ParameterSet.Utilities import convertToUnscheduled
process=convertToUnscheduled(process)

# customisation of the process.

# Automatic addition of the customisation function from PhysicsTools.PatAlgos.slimming.miniAOD_tools
from PhysicsTools.PatAlgos.slimming.miniAOD_tools import miniAOD_customizeAllMC 

#call to customisation function miniAOD_customizeAllMC imported from PhysicsTools.PatAlgos.slimming.miniAOD_tools
process = miniAOD_customizeAllMC(process)

# End of customisation functions

# Customisation from command line

process.PixelCPEGenericESProducer.IrradiationBiasCorrection = True
#Have logErrorHarvester wait for the same EDProducers to finish as those providing data for the OutputModule
from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
