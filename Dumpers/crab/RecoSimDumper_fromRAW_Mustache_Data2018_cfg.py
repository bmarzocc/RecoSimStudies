# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: -s RAW2DIGI,L1Reco,RECO,PAT --conditions 123X_dataRun2_v4 --data --era Run2_2018 --python_filename=RecoSimDumper_fromRAW_Mustache_Data2018_cfg.py --filein=/store/data/Run2022D/EGamma/RAW-RECO/ZElectron-PromptReco-v2/000/357/734/00000/0002e651-145c-4b78-8746-41629601d33f.root
import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing('standard')
options.register('inputFile',
                 'file:step2.root',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "intputFile")
options.register('outputFile',
                 'file:step3.root',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "outputFile")
                
options.parseArguments()


from Configuration.Eras.Era_Run2_2018_cff import Run2_2018

process = cms.Process('ECALClustering',Run2_2018)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')
process.load('PhysicsTools.PatAlgos.slimming.metFilterPaths_cff')
process.load('Configuration.StandardSequences.PAT_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFile),
    secondaryFileNames = cms.untracked.vstring()
)

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
    makeTriggerResults = cms.obsolete.untracked.bool,
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(1),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(2),
    printDependencies = cms.untracked.bool(False),
    sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(False)
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('-s nevts:1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string(''),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('step3.root'),
    outputCommands = process.RECOSIMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '123X_dataRun2_v4', '')

process.mySCReg = cms.ESSource("PoolDBESSource",
     toGet = cms.VPSet(
       cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("pfscecal_EBCorrection_offline_v2"),
         tag = cms.string("pfscecal_EBCorrection_offline_v2_2022ElectronsMustacheSC"),
         connect = cms.string("sqlite_file:scReg_2022ElectronsMustacheSC.db")),
       cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("pfscecal_EECorrection_offline_v2"),
         tag = cms.string("pfscecal_EECorrection_offline_v2_2022ElectronsMustacheSC"),
         connect = cms.string("sqlite_file:scReg_2022ElectronsMustacheSC.db")),
      cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("pfscecal_EBUncertainty_offline_v2"),
         tag = cms.string("pfscecal_EBUncertainty_offline_v2_2022ElectronsMustacheSC"),
         connect = cms.string("sqlite_file:scReg_2022ElectronsMustacheSC.db")),
      cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("pfscecal_EEUncertainty_offline_v2"),
         tag = cms.string("pfscecal_EEUncertainty_offline_v2_2022ElectronsMustacheSC"),
         connect = cms.string("sqlite_file:scReg_2022ElectronsMustacheSC.db")),
     )
)
process.es_prefer_scReg = cms.ESPrefer("PoolDBESSource","mySCReg")

process.myEleReg = cms.ESSource("PoolDBESSource",
     toGet = cms.VPSet(
      cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("electron_eb_ecalOnly_1To300_0p2To2_mean"),
         tag = cms.string("electron_eb_ecalOnly_1To300_0p2To2_mean_2022ElectronsMustacheSC"),
         connect = cms.string("sqlite_file:EleReg_2022ElectronsMustacheSC.db")),
      cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("electron_ee_ecalOnly_1To300_0p2To2_mean"),
         tag = cms.string("electron_ee_ecalOnly_1To300_0p2To2_mean_2022ElectronsMustacheSC"),
         connect = cms.string("sqlite_file:EleReg_2022ElectronsMustacheSC.db")),
      cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("electron_eb_ecalOnly_1To300_0p0002To0p5_sigma"),
         tag = cms.string("electron_eb_ecalOnly_1To300_0p0002To0p5_sigma_2022ElectronsMustacheSC"),
         connect = cms.string("sqlite_file:EleReg_2022ElectronsMustacheSC.db")),
      cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("electron_ee_ecalOnly_1To300_0p0002To0p5_sigma"),
         tag = cms.string("electron_ee_ecalOnly_1To300_0p0002To0p5_sigma_2022ElectronsMustacheSC"),
         connect = cms.string("sqlite_file:EleReg_2022ElectronsMustacheSC.db")),
      cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("electron_eb_ecalTrk_1To300_0p2To2_mean"),
         tag = cms.string("electron_eb_ecalTrk_1To300_0p2To2_mean_2022ElectronsMustacheSC"),
         connect = cms.string("sqlite_file:EleReg_2022ElectronsMustacheSC.db")),
     cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("electron_ee_ecalTrk_1To300_0p2To2_mean"),
         tag = cms.string("electron_ee_ecalTrk_1To300_0p2To2_mean_2022ElectronsMustacheSC"),
         connect = cms.string("sqlite_file:EleReg_2022ElectronsMustacheSC.db")),
     cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("electron_eb_ecalTrk_1To300_0p0002To0p5_sigma"),
         tag = cms.string("electron_eb_ecalTrk_1To300_0p0002To0p5_sigma_2022ElectronsMustacheSC"),
         connect = cms.string("sqlite_file:EleReg_2022ElectronsMustacheSC.db")),
     cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("electron_ee_ecalTrk_1To300_0p0002To0p5_sigma"),
         tag = cms.string("electron_ee_ecalTrk_1To300_0p0002To0p5_sigma_2022ElectronsMustacheSC"),
         connect = cms.string("sqlite_file:EleReg_2022ElectronsMustacheSC.db")),
     )
)
process.es_prefer_eleReg = cms.ESPrefer("PoolDBESSource","myEleReg")

process.TFileService = cms.Service("TFileService",
                                   #fileName = cms.string(options.outputFile)
                                   fileName = cms.string('output.root')
)

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.Flag_BadChargedCandidateFilter = cms.Path(process.BadChargedCandidateFilter)
process.Flag_BadChargedCandidateSummer16Filter = cms.Path(process.BadChargedCandidateSummer16Filter)
process.Flag_BadPFMuonDzFilter = cms.Path(process.BadPFMuonDzFilter)
process.Flag_BadPFMuonFilter = cms.Path(process.BadPFMuonFilter)
process.Flag_BadPFMuonSummer16Filter = cms.Path(process.BadPFMuonSummer16Filter)
process.Flag_CSCTightHalo2015Filter = cms.Path(process.CSCTightHalo2015Filter)
process.Flag_CSCTightHaloFilter = cms.Path(process.CSCTightHaloFilter)
process.Flag_CSCTightHaloTrkMuUnvetoFilter = cms.Path(process.CSCTightHaloTrkMuUnvetoFilter)
process.Flag_EcalDeadCellBoundaryEnergyFilter = cms.Path(process.EcalDeadCellBoundaryEnergyFilter)
process.Flag_EcalDeadCellTriggerPrimitiveFilter = cms.Path(process.EcalDeadCellTriggerPrimitiveFilter)
process.Flag_HBHENoiseFilter = cms.Path(process.HBHENoiseFilterResultProducer+process.HBHENoiseFilter)
process.Flag_HBHENoiseIsoFilter = cms.Path(process.HBHENoiseFilterResultProducer+process.HBHENoiseIsoFilter)
process.Flag_HcalStripHaloFilter = cms.Path(process.HcalStripHaloFilter)
process.Flag_METFilters = cms.Path(process.metFilters)
process.Flag_chargedHadronTrackResolutionFilter = cms.Path(process.chargedHadronTrackResolutionFilter)
process.Flag_ecalBadCalibFilter = cms.Path(process.ecalBadCalibFilter)
process.Flag_ecalLaserCorrFilter = cms.Path(process.ecalLaserCorrFilter)
process.Flag_eeBadScFilter = cms.Path(process.eeBadScFilter)
process.Flag_globalSuperTightHalo2016Filter = cms.Path(process.globalSuperTightHalo2016Filter)
process.Flag_globalTightHalo2016Filter = cms.Path(process.globalTightHalo2016Filter)
process.Flag_goodVertices = cms.Path(process.primaryVertexFilter)
process.Flag_hcalLaserEventFilter = cms.Path(process.hcalLaserEventFilter)
process.Flag_hfNoisyHitsFilter = cms.Path(process.hfNoisyHitsFilter)
process.Flag_muonBadTrackFilter = cms.Path(process.muonBadTrackFilter)
process.Flag_trackingFailureFilter = cms.Path(process.goodVertices+process.trackingFailureFilter)
process.Flag_trkPOGFilters = cms.Path(process.trkPOGFilters)
process.Flag_trkPOG_logErrorTooManyClusters = cms.Path(~process.logErrorTooManyClusters)
process.Flag_trkPOG_manystripclus53X = cms.Path(~process.manystripclus53X)
process.Flag_trkPOG_toomanystripclus53X = cms.Path(~process.toomanystripclus53X)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)

process.hcalDigis.saveQIE10DataNSamples = cms.untracked.vint32( 6)
process.hcalDigis.saveQIE10DataTags = cms.untracked.vstring( "QIE10Data" )

process.load("PhysicsTools.PatAlgos.patSequences_cff")

process.cleanPatTaus.preselection = cms.string('tauID("decayModeFinding") > 0.5 & tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits") > 0.5 & tauID("againstMuonTight3") > 0.5 ')
process.patMETs.addGenMET = cms.bool(False)

# filters
process.AllEvents                      = cms.EDProducer("EventCountProducer")
process.FilterL1FilterEvents           = cms.EDProducer("EventCountProducer")
process.FilterGoodVertexFilterEvents   = cms.EDProducer("EventCountProducer")
process.FilterNoScrapingFilterEvents   = cms.EDProducer("EventCountProducer")
process.FilterElectronFilterEvents     = cms.EDProducer("EventCountProducer")
process.FilterReRECOEvents             = cms.EDProducer("EventCountProducer")
process.FilterPatDefaultSequenceEvents = cms.EDProducer("EventCountProducer")

# filter on PhysDeclared bit
process.skimming = cms.EDFilter(
    "PhysDecl",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False),
    HLTriggerResults = cms.InputTag("TriggerResults","","HLT")
    )

# filter on bit = and (40 || 41) and !(bit36 || bit37 || bit38 || bit39)
process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')
process.load('HLTrigger/HLTfilters/hltLevel1GTSeed_cfi')
process.hltLevel1GTSeed.L1TechTriggerSeeding = cms.bool(True)
#process.hltLevel1GTSeed.L1SeedsLogicalExpression = cms.string('0 AND (40 OR 41) AND NOT (36 OR 37 OR 38 OR 39)')

from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.hltHighLevel = hltHighLevel.clone(HLTPaths = cms.vstring("HLT_Ele*","HLT_DoubleEle*") )

process.goodPrimaryVertices = cms.EDFilter("VertexSelector",
  src = cms.InputTag("offlineSlimmedPrimaryVertices"),     # "offlinePrimaryVertices"
  cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"),
  filter = cms.bool(True),
)

# filter on primary vertex
process.primaryVertexFilter = cms.EDFilter(
    "GoodVertexFilter",
    vertexCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),     # "offlinePrimaryVertices"
    minimumNDOF = cms.uint32(4) ,
    maxAbsZ = cms.double(24),
    maxd0 = cms.double(2)
)

# FilterOutScraping
process.noscraping = cms.EDFilter(
    "FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False),
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.25)
)

# select events with at least one gsf electron
process.highetele = cms.EDFilter(
    "PATElectronSelector",   #"GsfElectronSelector",
    src = cms.InputTag("slimmedElectrons"),  # "gedGsfElectrons"
    #"GsfElectronSelector",
    #src = cms.InputTag("gedGsfElectrons"),  # -> new!
    cut = cms.string("superCluster().get().energy()*sin(theta())> 0 ")
)

process.highetFilter = cms.EDFilter(
    "CandViewCountFilter",
    src = cms.InputTag("highetele"),
    minNumber = cms.uint32(1)
)

process.load('Calibration.EcalAlCaRecoProducers.WZElectronSkims_cff')
process.load('RecoSimStudies.Dumpers.RecoSimDumper_cfi')
process.dumper_step = cms.Path(
     process.ZeeSkimFilterSeq
    *process.AllEvents  
    #*process.hltLevel1GTSeed 
    *process.hltHighLevel
    *process.goodPrimaryVertices
    *process.highetele
    *process.highetFilter
    *process.FilterReRECOEvents   
    *process.FilterPatDefaultSequenceEvents 
    *process.recosimdumper
)

from RecoEgamma.EgammaElectronProducers.gsfElectrons_cfi import *
process.ecalDrivenGsfElectrons.EleDNNPFid= dict(
        modelsFiles = [
            "RecoEgamma/ElectronIdentification/data/Ele_PFID_dnn/Run3Summer21_120X/lowpT/lowpT_modelDNN.pb",
            "RecoEgamma/ElectronIdentification/data/Ele_PFID_dnn/Run3Summer21_120X/EB_highpT/barrel_highpT_modelDNN.pb",
            "RecoEgamma/ElectronIdentification/data/Ele_PFID_dnn/Run3Summer21_120X/EE_highpT/endcap_highpT_modelDNN.pb"
        ],
        scalersFiles = [
            "RecoEgamma/ElectronIdentification/data/Ele_PFID_dnn/Run3Summer21_120X/lowpT/lowpT_scaler.txt",
            "RecoEgamma/ElectronIdentification/data/Ele_PFID_dnn/Run3Summer21_120X/EB_highpT/barrel_highpT_scaler.txt",
            "RecoEgamma/ElectronIdentification/data/Ele_PFID_dnn/Run3Summer21_120X/EE_highpT/endcap_highpT_scaler.txt"
        ]
)
process.ecalDrivenGsfElectrons.EleDNNPFid.enabled = True

from RecoEgamma.EgammaPhotonProducers.gedPhotons_cfi import *
process.gedPhotons.PhotonDNNPFid = dict(
        modelsFiles = [ 
            "RecoEgamma/PhotonIdentification/data/Photon_PFID_dnn/Run3Summer21_120X/EB/barrel_modelDNN.pb",
            "RecoEgamma/PhotonIdentification/data/Photon_PFID_dnn/Run3Summer21_120X/EE/endcap_modelDNN.pb"
        ],
        scalersFiles = [
            "RecoEgamma/PhotonIdentification/data/Photon_PFID_dnn/Run3Summer21_120X/EB/barrel_scaler.txt",
            "RecoEgamma/PhotonIdentification/data/Photon_PFID_dnn/Run3Summer21_120X/EE/endcap_scaler.txt" 
        ]
)
process.gedPhotons.PhotonDNNPFid.enabled = True  

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.Flag_HBHENoiseFilter,process.Flag_HBHENoiseIsoFilter,process.Flag_CSCTightHaloFilter,process.Flag_CSCTightHaloTrkMuUnvetoFilter,process.Flag_CSCTightHalo2015Filter,process.Flag_globalTightHalo2016Filter,process.Flag_globalSuperTightHalo2016Filter,process.Flag_HcalStripHaloFilter,process.Flag_hcalLaserEventFilter,process.Flag_EcalDeadCellTriggerPrimitiveFilter,process.Flag_EcalDeadCellBoundaryEnergyFilter,process.Flag_ecalBadCalibFilter,process.Flag_goodVertices,process.Flag_eeBadScFilter,process.Flag_ecalLaserCorrFilter,process.Flag_trkPOGFilters,process.Flag_chargedHadronTrackResolutionFilter,process.Flag_muonBadTrackFilter,process.Flag_BadChargedCandidateFilter,process.Flag_BadPFMuonFilter,process.Flag_BadPFMuonDzFilter,process.Flag_hfNoisyHitsFilter,process.Flag_BadChargedCandidateSummer16Filter,process.Flag_BadPFMuonSummer16Filter,process.Flag_trkPOG_manystripclus53X,process.Flag_trkPOG_toomanystripclus53X,process.Flag_trkPOG_logErrorTooManyClusters,process.Flag_METFilters,process.dumper_step,process.endjob_step)
process.schedule.associate(process.patTask)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)


# customisation of the process.

# Automatic addition of the customisation function from PhysicsTools.PatAlgos.slimming.miniAOD_tools
from PhysicsTools.PatAlgos.slimming.miniAOD_tools import miniAOD_customizeAllData 

#call to customisation function miniAOD_customizeAllData imported from PhysicsTools.PatAlgos.slimming.miniAOD_tools
process = miniAOD_customizeAllData(process)

# End of customisation functions

# Customisation from command line

#Have logErrorHarvester wait for the same EDProducers to finish as those providing data for the OutputModule
from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
