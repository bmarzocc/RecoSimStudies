# CRAB3 config template for flashgg
# More options available on the twiki :
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookCRAB3Tutorial

from WMCore.Configuration import Configuration

config = Configuration()

config.section_('General')
config.General.requestName       = 'FourElectronsGunPt1-500_Mustache_thres235fb'
config.General.transferLogs      = True
config.General.transferOutputs   = True

config.section_('JobType')
config.JobType.pluginName        = 'Analysis'

# Name of the CMSSW configuration file
config.JobType.psetName          = 'RecoSimDumper_fromRAW_Mustache_StdMC_noise235_thres235fb_cfg.py'
config.JobType.priority          = 30
config.JobType.maxMemoryMB       = 3000
config.JobType.numCores          = 2

config.section_('Data')
# This string determines the primary dataset of the newly-produced outputs.
config.Data.inputDataset         = '/FourElectronsGunPt1-500_pythia8_PremixRun3_13p6TeV_122X_mcRun3_2021_realistic_v9_235fbNoise/dpg_ecal-GEN-SIM-RAW-CALO-6cb02bfc714e0c227af683c336cc381d/USER'
#config.Data.inputDataset         = '/FourGammasGunPt1-500_pythia8_PremixRun3_13p6TeV_122X_mcRun3_2021_realistic_v9_235fbNoise/dpg_ecal-GEN-SIM-RAW-CALO-6cb02bfc714e0c227af683c336cc381d/USER'
config.JobType.outputFiles       = ["output.root"]
config.JobType.inputFiles        = ["/afs/cern.ch/work/b/bmarzocc/Clustering_stdSeeding/CMSSW_13_0_0_pre4/src/RecoSimStudies/Dumpers/data/regressions/EleReg_13_0_0_235fbNoise_ElectronsDeepSCAlgoA_thres235fb_34sigma.db","/afs/cern.ch/work/b/bmarzocc/Clustering_stdSeeding/CMSSW_13_0_0_pre4/src/RecoSimStudies/Dumpers/data/regressions/EleReg_13_0_0_235fbNoise_ElectronsDeepSCAlgoA_thresUL18_2e3sigma.db","/afs/cern.ch/work/b/bmarzocc/Clustering_stdSeeding/CMSSW_13_0_0_pre4/src/RecoSimStudies/Dumpers/data/regressions/EleReg_13_0_0_235fbNoise_ElectronsMustache_thres235fb_34sigma.db","/afs/cern.ch/work/b/bmarzocc/Clustering_stdSeeding/CMSSW_13_0_0_pre4/src/RecoSimStudies/Dumpers/data/regressions/EleReg_13_0_0_235fbNoise_ElectronsMustache_thresUL18_2e3sigma.db","/afs/cern.ch/work/b/bmarzocc/Clustering_stdSeeding/CMSSW_13_0_0_pre4/src/RecoSimStudies/Dumpers/data/regressions/SCReg_13_0_0_235fbNoise_ElectronsDeepSCAlgoA_thres235fb_34sigma.db","/afs/cern.ch/work/b/bmarzocc/Clustering_stdSeeding/CMSSW_13_0_0_pre4/src/RecoSimStudies/Dumpers/data/regressions/SCReg_13_0_0_235fbNoise_ElectronsDeepSCAlgoA_thresUL18_2e3sigma.db","/afs/cern.ch/work/b/bmarzocc/Clustering_stdSeeding/CMSSW_13_0_0_pre4/src/RecoSimStudies/Dumpers/data/regressions/SCReg_13_0_0_235fbNoise_ElectronsMustache_thres235fb_34sigma.db","/afs/cern.ch/work/b/bmarzocc/Clustering_stdSeeding/CMSSW_13_0_0_pre4/src/RecoSimStudies/Dumpers/data/regressions/SCReg_13_0_0_235fbNoise_ElectronsMustache_thresUL18_2e3sigma.db"]  
config.Data.inputDBS             = 'global'   
config.Data.inputDBS             = 'phys03'   
config.Data.splitting            = 'FileBased'
config.Data.unitsPerJob          = 1
config.Data.publication          = False
config.Data.ignoreLocality       = True
config.Data.outputDatasetTag     = 'RECO_Mustache_thres235fb'

# This string is used to construct the output dataset name
#config.Data.publishDataName = 'CRAB3-tutorial'
config.Data.outLFNDirBase        =  '/store/group/dpg_ecal/alca_ecalcalib/bmarzocc/Clustering/'

config.section_('Site')
# Where the output files will be transmitted to
config.Site.storageSite         = 'T2_CH_CERN'
config.Site.whitelist           = ['T2_CH_CERN']


## config.Data.allowNonValidInputDataset=True
