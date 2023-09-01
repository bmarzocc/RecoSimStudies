from CRABClient.UserUtilities import config

config = config()

## General settings
config.General.requestName = 'FourElectronsGunPt1-500_pythia8_PremixRun3_bsRealistic2022_13p6TeV_126X_mcRun3_2023_forPU65_v4'
config.General.transferOutputs = True
config.General.transferLogs = True
## PrivateMC type with a fake miniAOD step to circunvent crab requests (official data-tier for PrivateMC)
config.JobType.pluginName  = 'PrivateMC'
config.JobType.psetName    = 'step2_DIGI_DATAMIX_L1_DIGI2RAW_HLT_PREMIX_Run3_2023_fake.py'
config.JobType.pyCfgParams = ['nThreads=4','outputName=output.root']
## To be executed on node with Arguments
config.JobType.scriptExe   = 'scriptExe.sh'
config.JobType.scriptArgs  = ['nEvents=100','nThreads=4','outputName=output.root']
config.JobType.inputFiles  = ['scriptExe.sh','ElectronsGunPt1-500_pythia8_cfi_GEN_SIM.py','GammasGunPt1-500_pythia8_cfi_GEN_SIM.py','pileup.py','step2_DIGI_DATAMIX_L1_DIGI2RAW_HLT_PREMIX_Run3_2023.py']
## Output file to be collected
config.JobType.outputFiles = ["output.root"]
config.JobType.disableAutomaticOutputCollection = True
## Memory, cores, cmssw
config.JobType.allowUndistributedCMSSW = True
config.JobType.maxMemoryMB = 5500
config.JobType.numCores    = 4

config.JobType.sendPythonFolder = True
## Data
config.Data.splitting   = 'EventBased'
config.Data.unitsPerJob = 100
config.Data.totalUnits  = 1000000

config.Data.outputPrimaryDataset = 'FourElectronsGunPt1-500_pythia8_PremixRun3_bsRealistic2022_13p6TeV_126X_mcRun3_2023_forPU65_v4'
config.Data.publication          = True
config.Data.publishDBS           = 'phys03'
#config.Data.ignoreLocality       = True
config.Data.outputDatasetTag     = 'GEN-SIM-RAW-CALO'
config.Data.outLFNDirBase        =  '/store/group/dpg_ecal/Clustering/'
config.Data.publishWithGroupName = True


## Site
config.Site.storageSite         = 'T2_CH_CSCS'
#config.Site.whitelist           = ['T2_CH_CSCS']
