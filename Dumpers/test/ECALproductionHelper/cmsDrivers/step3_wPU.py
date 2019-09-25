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
options.register ("seedMult",
                  3.0, # default value
                  VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.varType.float,          # string, int, or float
                  "multiplier of noise used for seeding threshold")
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

process.MINIAODSIMoutput = cms.OutputModule("PoolOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(4),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('MINIAODSIM'),
        filterName = cms.untracked.string('')
    ),
    dropMetaData = cms.untracked.string('ALL'),
    eventAutoFlushCompressedSize = cms.untracked.int32(-900),
    fastCloning = cms.untracked.bool(False),
    fileName = cms.untracked.string('file:step3_inMINIAODSIM.root'),
    outputCommands = process.MINIAODSIMEventContent.outputCommands,
    overrideBranchesSplitLevel = cms.untracked.VPSet(
        cms.untracked.PSet(
            branch = cms.untracked.string('patPackedCandidates_packedPFCandidates__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('recoGenParticles_prunedGenParticles__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('patTriggerObjectStandAlones_slimmedPatTrigger__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('patPackedGenParticles_packedGenParticles__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('patJets_slimmedJets__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('recoVertexs_offlineSlimmedPrimaryVertices__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('recoCaloClusters_reducedEgamma_reducedESClusters_*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('EcalRecHitsSorted_reducedEgamma_reducedEBRecHits_*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('EcalRecHitsSorted_reducedEgamma_reducedEERecHits_*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('recoGenJets_slimmedGenJets__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('patJets_slimmedJetsPuppi__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('EcalRecHitsSorted_reducedEgamma_reducedESRecHits_*'),
            splitLevel = cms.untracked.int32(99)
        )
    ),
    overrideInputFileSplitLevels = cms.untracked.bool(True),
    splitLevel = cms.untracked.int32(0)
)

# process.DQMoutput = cms.OutputModule("DQMRootOutputModule",
#     dataset = cms.untracked.PSet(
#         dataTier = cms.untracked.string('DQMIO'),
#         filterName = cms.untracked.string('')
#     ),
#     fileName = cms.untracked.string('file:step3_inDQM.root'),
#     outputCommands = process.DQMEventContent.outputCommands,
#     splitLevel = cms.untracked.int32(0)
# )

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
process.myCond.producedEcalPFRecHitThresholds = cms.untracked.bool(True)
process.myCond.EcalPFRecHitThresholdNSigmas = cms.untracked.double(1.0)
process.myCond.EcalPFRecHitThresholdNSigmasHEta = cms.untracked.double(1.0)
process.myCond.PFRecHitFile = cms.untracked.string("./data/noise/PFRecHitThresholds_EB.txt")
process.myCond.PFRecHitFileEE = cms.untracked.string("./data/noise/PFRecHitThresholds_EE.txt")

process.myCond.producedEcalPFSeedingThresholds = cms.untracked.bool(True)
process.myCond.EcalPFSeedingThresholdNSigmas = cms.untracked.double(options.seedMult/2.0) # PFRHs files are at 2sigma of the noise for |eta|<2.5
process.myCond.EcalPFSeedingThresholdNSigmasHEta = cms.untracked.double(options.seedMult/3.0) #                3sigma of the noise for |eta|>2.5
process.myCond.PFSeedingFile = cms.untracked.string("./data/noise/PFRecHitThresholds_EB.txt")
process.myCond.PFSeedingFileEE = cms.untracked.string("./data/noise/PFRecHitThresholds_EE.txt")

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
########## end override

# Path and EndPath definitions 
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.recosim_step = cms.Path(process.recosim)
process.eventinterpretaion_step = cms.Path(process.EIsequence)
process.Flag_trackingFailureFilter = cms.Path(process.goodVertices+process.trackingFailureFilter)
process.Flag_goodVertices = cms.Path(process.primaryVertexFilter)
# process.Flag_CSCTightHaloFilter = cms.Path(process.CSCTightHaloFilter)
# process.Flag_trkPOGFilters = cms.Path(process.trkPOGFilters)
# process.Flag_HcalStripHaloFilter = cms.Path(process.HcalStripHaloFilter)
# process.Flag_trkPOG_logErrorTooManyClusters = cms.Path(~process.logErrorTooManyClusters)
process.Flag_EcalDeadCellTriggerPrimitiveFilter = cms.Path(process.EcalDeadCellTriggerPrimitiveFilter)
process.Flag_ecalLaserCorrFilter = cms.Path(process.ecalLaserCorrFilter)
#process.Flag_globalSuperTightHalo2016Filter = cms.Path(process.globalSuperTightHalo2016Filter)
process.Flag_eeBadScFilter = cms.Path(process.eeBadScFilter)
process.Flag_METFilters = cms.Path(process.metFilters)
process.Flag_chargedHadronTrackResolutionFilter = cms.Path(process.chargedHadronTrackResolutionFilter)
# process.Flag_globalTightHalo2016Filter = cms.Path(process.globalTightHalo2016Filter)
# process.Flag_CSCTightHaloTrkMuUnvetoFilter = cms.Path(process.CSCTightHaloTrkMuUnvetoFilter)
# process.Flag_HBHENoiseIsoFilter = cms.Path(process.HBHENoiseFilterResultProducer+process.HBHENoiseIsoFilter)
process.Flag_BadChargedCandidateSummer16Filter = cms.Path(process.BadChargedCandidateSummer16Filter)
# process.Flag_hcalLaserEventFilter = cms.Path(process.hcalLaserEventFilter)
# process.Flag_BadPFMuonFilter = cms.Path(process.BadPFMuonFilter)
process.Flag_ecalBadCalibFilter = cms.Path(process.ecalBadCalibFilter)
# process.Flag_HBHENoiseFilter = cms.Path(process.HBHENoiseFilterResultProducer+process.HBHENoiseFilter)
# process.Flag_trkPOG_toomanystripclus53X = cms.Path(~process.toomanystripclus53X)
process.Flag_EcalDeadCellBoundaryEnergyFilter = cms.Path(process.EcalDeadCellBoundaryEnergyFilter)
process.Flag_BadChargedCandidateFilter = cms.Path(process.BadChargedCandidateFilter)
# process.Flag_trkPOG_manystripclus53X = cms.Path(~process.manystripclus53X)
# process.Flag_BadPFMuonSummer16Filter = cms.Path(process.BadPFMuonSummer16Filter)
# process.Flag_muonBadTrackFilter = cms.Path(process.muonBadTrackFilter)
# process.Flag_CSCTightHalo2015Filter = cms.Path(process.CSCTightHalo2015Filter)
# process.prevalidation_step = cms.Path(process.prevalidation)
# process.prevalidation_step1 = cms.Path(process.prevalidationMiniAOD)
# process.validation_step = cms.EndPath(process.validation)
# process.validation_step1 = cms.EndPath(process.validationMiniAOD)
# process.dqmoffline_step = cms.EndPath(process.DQMOffline)
# process.dqmoffline_1_step = cms.EndPath(process.DQMOfflineExtraHLT)
# process.dqmoffline_2_step = cms.EndPath(process.DQMOfflineMiniAOD)
# process.dqmofflineOnPAT_step = cms.EndPath(process.PostDQMOffline)
# process.dqmofflineOnPAT_1_step = cms.EndPath(process.PostDQMOffline)
# process.dqmofflineOnPAT_2_step = cms.EndPath(process.PostDQMOfflineMiniAOD)
process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)
process.MINIAODSIMoutput_step = cms.EndPath(process.MINIAODSIMoutput)
# process.DQMoutput_step = cms.EndPath(process.DQMoutput)

# process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.recosim_step,process.eventinterpretaion_step,process.Flag_HBHENoiseFilter,process.Flag_HBHENoiseIsoFilter,process.Flag_CSCTightHaloFilter,process.Flag_CSCTightHaloTrkMuUnvetoFilter,process.Flag_CSCTightHalo2015Filter,process.Flag_globalTightHalo2016Filter,process.Flag_globalSuperTightHalo2016Filter,process.Flag_HcalStripHaloFilter,process.Flag_hcalLaserEventFilter,process.Flag_EcalDeadCellTriggerPrimitiveFilter,process.Flag_EcalDeadCellBoundaryEnergyFilter,process.Flag_ecalBadCalibFilter,process.Flag_goodVertices,process.Flag_eeBadScFilter,process.Flag_ecalLaserCorrFilter,process.Flag_trkPOGFilters,process.Flag_chargedHadronTrackResolutionFilter,process.Flag_muonBadTrackFilter,process.Flag_BadChargedCandidateFilter,process.Flag_BadPFMuonFilter,process.Flag_BadChargedCandidateSummer16Filter,process.Flag_BadPFMuonSummer16Filter,process.Flag_trkPOG_manystripclus53X,process.Flag_trkPOG_toomanystripclus53X,process.Flag_trkPOG_logErrorTooManyClusters,process.Flag_METFilters,process.prevalidation_step,process.prevalidation_step1,process.validation_step,process.validation_step1,\
#     process.dqmoffline_step,process.dqmoffline_1_step,process.dqmoffline_2_step,process.dqmofflineOnPAT_step,process.dqmofflineOnPAT_1_step,process.dqmofflineOnPAT_2_step,\
#         process.RECOSIMoutput_step,process.MINIAODSIMoutput_step,\
#             process.DQMoutput_step)
# process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.recosim_step,
#                             process.eventinterpretaion_step,process.Flag_HBHENoiseFilter,process.Flag_HBHENoiseIsoFilter,  
#                             process.Flag_CSCTightHaloFilter,process.Flag_CSCTightHaloTrkMuUnvetoFilter,
#                             process.Flag_CSCTightHalo2015Filter,process.Flag_globalTightHalo2016Filter,
#                             process.Flag_globalSuperTightHalo2016Filter,process.Flag_HcalStripHaloFilter,
#                             process.Flag_hcalLaserEventFilter,process.Flag_EcalDeadCellTriggerPrimitiveFilter,
#                             process.Flag_EcalDeadCellBoundaryEnergyFilter,process.Flag_ecalBadCalibFilter,
#                             process.Flag_goodVertices,process.Flag_eeBadScFilter,process.Flag_ecalLaserCorrFilter,
#                             process.Flag_trkPOGFilters,process.Flag_chargedHadronTrackResolutionFilter,
#                             process.Flag_muonBadTrackFilter,process.Flag_BadChargedCandidateFilter,
#                             process.Flag_BadPFMuonFilter,process.Flag_BadChargedCandidateSummer16Filter,
#                             process.Flag_BadPFMuonSummer16Filter,process.Flag_trkPOG_manystripclus53X,
#                             process.Flag_trkPOG_toomanystripclus53X,process.Flag_trkPOG_logErrorTooManyClusters,
#                             process.Flag_METFilters,process.prevalidation_step,process.prevalidation_step1,
#                             process.validation_step,process.validation_step1,process.RECOSIMoutput_step,
#                             process.MINIAODSIMoutput_step)

# Remove because of errors in condor
#process.pfTausBaseSequence.remove(pfTausProducerSansRefs)

process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.recosim_step,
                            process.eventinterpretaion_step,  
                            process.Flag_EcalDeadCellTriggerPrimitiveFilter,
                            process.Flag_EcalDeadCellBoundaryEnergyFilter,process.Flag_ecalBadCalibFilter,
                            process.Flag_goodVertices,process.Flag_eeBadScFilter,process.Flag_ecalLaserCorrFilter,
                            process.Flag_BadChargedCandidateFilter,
                            process.Flag_METFilters,
                            process.RECOSIMoutput_step,
                            process.MINIAODSIMoutput_step)
process.schedule.associate(process.patTask)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)


#Setup FWK for multithreaded 
process.options.numberOfThreads=cms.untracked.uint32(options.nThr)
process.options.numberOfStreams=cms.untracked.uint32(0)
process.options.numberOfConcurrentLuminosityBlocks=cms.untracked.uint32(1)


# customisation of the process.

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
