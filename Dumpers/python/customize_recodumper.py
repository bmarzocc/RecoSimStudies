import FWCore.ParameterSet.Config as cms

#Produce cfg for running RecoSimDumper on RAW with Mustache using nonStdMC: cmsDriver.py --python_filename RecoSimDumper_fromRAW_Mustache_NonStdMC_cfg.py --eventcontent MINIAODSIM --datatier MINIAODSIM --fileout file:step3.root --conditions 125X_mcRun3_2022_realistic_v4 --step RAW2DIGI,L1Reco,RECO,RECOSIM,PAT --geometry DB:Extended --filein root://cms-xrd-global.cern.ch//store/group/dpg_ecal/Clustering/FourElectronsGunPt1-100_pythia8_StdMix_Flat55To75_14TeV_112X_mcRun3_2021_realistic_v15/GEN-SIM-RAW-Reduced/210519_101349/0003/cluster_job3235_step2_3245.root --era Run3,ctpps_2018 --no_exec --mc --customise RecoSimStudies/Dumpers/customize_recodumper.customize_baseline,RecoSimStudies/Dumpers/customize_recodumper.customize_dumperStepMC,RecoSimStudies/Dumpers/customize_recodumper.customize_linkedObjectsStep,RecoSimStudies/Dumpers/customize_recodumper.customize_addConditions,RecoSimStudies/Dumpers/customize_recodumper.customize_nonStdMC --processName ECALClustering --customise_commands 'from RecoEgamma.EgammaElectronProducers.gsfElectrons_cfi import * \n from RecoEgamma.EgammaPhotonProducers.gedPhotons_cfi import * \n from SimMuon.MCTruth.muonSimClassificationByHits_cff import *' -n 10

#Produce cfg for running RecoSimDumper on RAW with Mustache using StdMC: cmsDriver.py --python_filename RecoSimDumper_fromRAW_Mustache_StdMC_cfg.py --eventcontent MINIAODSIM --datatier MINIAODSIM --fileout file:step3.root --conditions 125X_mcRun3_2022_realistic_v4 --step RAW2DIGI,L1Reco,RECO,RECOSIM,PAT --geometry DB:Extended --filein root://cms-xrd-global.cern.ch//store/group/dpg_ecal/Clustering/FourElectronsGunPt1-500_pythia8_PremixRun3_13p6TeV_122X_mcRun3_2021_realistic_v9_235fbNoise/GEN-SIM-RAW-CALO/221026_074347/0001/cluster_job99_step2_1778.root --era Run3 --no_exec --mc --customise RecoSimStudies/Dumpers/customize_recodumper.customize_baseline,RecoSimStudies/Dumpers/customize_recodumper.customize_dumperStepMC,RecoSimStudies/Dumpers/customize_recodumper.customize_linkedObjectsStep,RecoSimStudies/Dumpers/customize_recodumper.customize_addConditions --processName ECALClustering --customise_commands 'from RecoEgamma.EgammaElectronProducers.gsfElectrons_cfi import * \n from RecoEgamma.EgammaPhotonProducers.gedPhotons_cfi import *' -n 10

#Produce cfg for running RecoSimDumper on RAW with DeepSC_algoA using StdMC: cmsDriver.py --python_filename RecoSimDumper_fromRAW_DeepSC_algoA_cfg.py --eventcontent MINIAODSIM --datatier MINIAODSIM --fileout file:step3.root --conditions 125X_mcRun3_2022_realistic_v4 --step RAW2DIGI,L1Reco,RECO,RECOSIM,PAT --geometry DB:Extended --filein root://cms-xrd-global.cern.ch//store/group/dpg_ecal/Clustering/FourElectronsGunPt1-500_pythia8_PremixRun3_13p6TeV_122X_mcRun3_2021_realistic_v9_235fbNoise/GEN-SIM-RAW-CALO/221026_074347/0001/cluster_job99_step2_1778.root --era Run3 --no_exec --mc --customise RecoSimStudies/Dumpers/customize_recodumper.customize_baseline,RecoSimStudies/Dumpers/customize_recodumper.customize_dumperStepMC,RecoSimStudies/Dumpers/customize_recodumper.customize_linkedObjectsStep,RecoSimStudies/Dumpers/customize_recodumper.customize_addConditions --processName ECALClustering --customise_commands 'from RecoEgamma.EgammaElectronProducers.gsfElectrons_cfi import * \n from RecoEgamma.EgammaPhotonProducers.gedPhotons_cfi import *' --procModifier ecal_deepsc -n 10

#Produce cfg for running RecoSimDumper on RAW with DeepSC_algoB using StdMC: cmsDriver.py --python_filename RecoSimDumper_fromRAW_DeepSC_algoB_cfg.py --eventcontent MINIAODSIM --datatier MINIAODSIM --fileout file:step3.root --conditions 125X_mcRun3_2022_realistic_v4 --step RAW2DIGI,L1Reco,RECO,RECOSIM,PAT --geometry DB:Extended --filein root://cms-xrd-global.cern.ch//store/group/dpg_ecal/Clustering/FourElectronsGunPt1-500_pythia8_PremixRun3_13p6TeV_122X_mcRun3_2021_realistic_v9_235fbNoise/GEN-SIM-RAW-CALO/221026_074347/0001/cluster_job99_step2_1778.root --era Run3 --no_exec --mc --customise RecoSimStudies/Dumpers/customize_recodumper.customize_baseline,RecoSimStudies/Dumpers/customize_recodumper.customize_dumperStepMC,RecoSimStudies/Dumpers/customize_recodumper.customize_linkedObjectsStep,RecoSimStudies/Dumpers/customize_recodumper.customize_addConditions,RecoSimStudies/Dumpers/customize_recodumper.customize_deepSCAlgoB --processName ECALClustering --customise_commands 'from RecoEgamma.EgammaElectronProducers.gsfElectrons_cfi import * \n from RecoEgamma.EgammaPhotonProducers.gedPhotons_cfi import * ' --procModifier ecal_deepsc -n 10

#Produce cfg for running RecoSimDumper on RAW with DeepSC_algoC using StdMC: cmsDriver.py --python_filename RecoSimDumper_fromRAW_DeepSC_algoC_cfg.py --eventcontent MINIAODSIM --datatier MINIAODSIM --fileout file:step3.root --conditions 125X_mcRun3_2022_realistic_v4 --step RAW2DIGI,L1Reco,RECO,RECOSIM,PAT --geometry DB:Extended --filein root://cms-xrd-global.cern.ch//store/group/dpg_ecal/Clustering/FourElectronsGunPt1-500_pythia8_PremixRun3_13p6TeV_122X_mcRun3_2021_realistic_v9_235fbNoise/GEN-SIM-RAW-CALO/221026_074347/0001/cluster_job99_step2_1778.root --era Run3 --no_exec --mc --customise RecoSimStudies/Dumpers/customize_recodumper.customize_baseline,RecoSimStudies/Dumpers/customize_recodumper.customize_dumperStepMC,RecoSimStudies/Dumpers/customize_recodumper.customize_linkedObjectsStep,RecoSimStudies/Dumpers/customize_recodumper.customize_addConditions,RecoSimStudies/Dumpers/customize_recodumper.customize_deepSCAlgoC --processName ECALClustering --customise_commands 'from RecoEgamma.EgammaElectronProducers.gsfElectrons_cfi import * \n from RecoEgamma.EgammaPhotonProducers.gedPhotons_cfi import * ' --procModifier ecal_deepsc -n 10

#Produce cfg for running RecoSimDumper on RAW with DeepSC_algoD using StdMC: cmsDriver.py --python_filename RecoSimDumper_fromRAW_DeepSC_algoD_cfg.py --eventcontent MINIAODSIM --datatier MINIAODSIM --fileout file:step3.root --conditions 125X_mcRun3_2022_realistic_v4 --step RAW2DIGI,L1Reco,RECO,RECOSIM,PAT --geometry DB:Extended --filein root://cms-xrd-global.cern.ch//store/group/dpg_ecal/Clustering/FourElectronsGunPt1-500_pythia8_PremixRun3_13p6TeV_122X_mcRun3_2021_realistic_v9_235fbNoise/GEN-SIM-RAW-CALO/221026_074347/0001/cluster_job99_step2_1778.root --era Run3 --no_exec --mc --customise RecoSimStudies/Dumpers/customize_recodumper.customize_baseline,RecoSimStudies/Dumpers/customize_recodumper.customize_dumperStepMC,RecoSimStudies/Dumpers/customize_recodumper.customize_linkedObjectsStep,RecoSimStudies/Dumpers/customize_recodumper.customize_addConditions,RecoSimStudies/Dumpers/customize_recodumper.customize_deepSCAlgoD --processName ECALClustering --customise_commands 'from RecoEgamma.EgammaElectronProducers.gsfElectrons_cfi import * \n from RecoEgamma.EgammaPhotonProducers.gedPhotons_cfi import * ' --procModifier ecal_deepsc -n 10

#Produce cfg for running RecoSimDumper on RAW with Mustache using Data2018: cmsDriver.py --python_filename RecoSimDumper_fromRAW_Mustache_Data2018_cfg.py --eventcontent MINIAOD --datatier MINIAOD --fileout file:step3.root --conditions 124X_dataRun2_v2 --step RAW2DIGI,L1Reco,RECO,PAT --geometry DB:Extended --filein root://cms-xrd-global.cern.ch//store/data/Run2018C/EGamma/RAW/v1/000/320/026/00000/8CF35F1E-678D-E811-84C0-FA163E8CC774.root --era Run2_2018 --no_exec --data --customise RecoSimStudies/Dumpers/customize_recodumper.customize_baseline,RecoSimStudies/Dumpers/customize_recodumper.customize_dumperStepData,RecoSimStudies/Dumpers/customize_recodumper.customize_linkedObjectsStep --processName ECALClustering --customise_commands 'from RecoEgamma.EgammaElectronProducers.gsfElectrons_cfi import * \n from RecoEgamma.EgammaPhotonProducers.gedPhotons_cfi import *' -n 10

#Produce cfg for running RecoSimDumper on RAW with Mustache using Data2018 with ZeeSkim: cmsDriver.py --python_filename RecoSimDumper_fromRAW_Mustache_Data2018_cfg.py --eventcontent MINIAOD --datatier MINIAOD --fileout file:step3.root --conditions 124X_dataRun2_v2 --step RAW2DIGI,L1Reco,RECO,PAT --geometry DB:Extended --filein root://cms-xrd-global.cern.ch//store/data/Run2018C/EGamma/RAW/v1/000/320/026/00000/8CF35F1E-678D-E811-84C0-FA163E8CC774.root --era Run2_2018 --no_exec --data --customise RecoSimStudies/Dumpers/customize_recodumper.customize_baseline,RecoSimStudies/Dumpers/customize_recodumper.customize_zeeFilterDumperData,RecoSimStudies/Dumpers/customize_recodumper.customize_linkedObjectsStep --processName ECALClustering --customise_commands 'from RecoEgamma.EgammaElectronProducers.gsfElectrons_cfi import * \n from RecoEgamma.EgammaPhotonProducers.gedPhotons_cfi import *' -n 10

#Produce cfg for running RecoSimDumper on RAW with DeepSC_algoA using Data2018: cmsDriver.py --python_filename RecoSimDumper_fromRAW_DeepSC_algoA_Data2018_cfg.py --eventcontent MINIAOD --datatier MINIAOD --fileout file:step3.root --conditions 124X_dataRun2_v2 --step RAW2DIGI,L1Reco,RECO,PAT --geometry DB:Extended --filein root://cms-xrd-global.cern.ch//store/data/Run2018C/EGamma/RAW/v1/000/320/026/00000/8CF35F1E-678D-E811-84C0-FA163E8CC774.root --era Run2_2018 --no_exec --data --customise RecoSimStudies/Dumpers/customize_recodumper.customize_baseline,RecoSimStudies/Dumpers/customize_recodumper.customize_dumperStepData,RecoSimStudies/Dumpers/customize_recodumper.customize_linkedObjectsStep --processName ECALClustering --customise_commands 'from RecoEgamma.EgammaElectronProducers.gsfElectrons_cfi import * \n from RecoEgamma.EgammaPhotonProducers.gedPhotons_cfi import *'  --procModifier ecal_deepsc -n 10

#Produce cfg for running RecoSimDumper on RAW with DeepSC_algoB using Data2018: cmsDriver.py --python_filename RecoSimDumper_fromRAW_DeepSC_algoB_Data2018_cfg.py --eventcontent MINIAOD --datatier MINIAOD --fileout file:step3.root --conditions 124X_dataRun2_v2 --step RAW2DIGI,L1Reco,RECO,PAT --geometry DB:Extended --filein root://cms-xrd-global.cern.ch//store/data/Run2018C/EGamma/RAW/v1/000/320/026/00000/8CF35F1E-678D-E811-84C0-FA163E8CC774.root --era Run2_2018 --no_exec --data --customise RecoSimStudies/Dumpers/customize_recodumper.customize_baseline,RecoSimStudies/Dumpers/customize_recodumper.customize_dumperStepData,RecoSimStudies/Dumpers/customize_recodumper.customize_linkedObjectsStep,RecoSimStudies/Dumpers/customize_recodumper.customize_deepSCAlgoB --processName ECALClustering --customise_commands 'from RecoEgamma.EgammaElectronProducers.gsfElectrons_cfi import * \n from RecoEgamma.EgammaPhotonProducers.gedPhotons_cfi import *'  --procModifier ecal_deepsc -n 10

#Produce cfg for running RecoSimDumper on RAW with Mustache using Data2022: cmsDriver.py --python_filename RecoSimDumper_fromRAW_Mustache_Data2022_cfg.py --eventcontent MINIAOD --datatier MINIAOD --fileout file:step3.root --conditions 125X_dataRun3_relval_v4 --step RAW2DIGI,L1Reco,RECO,PAT --geometry DB:Extended --filein root://cms-xrd-global.cern.ch//store/data/Run2022F/EGamma/RAW-RECO/ZElectron-PromptReco-v1/000/360/390/00000/00d2f8c6-db74-477a-b106-e723005f62ac.root --era Run3 --no_exec --data --customise RecoSimStudies/Dumpers/customize_recodumper.customize_baseline,RecoSimStudies/Dumpers/customize_recodumper.customize_dumperStepData,RecoSimStudies/Dumpers/customize_recodumper.customize_linkedObjectsStep --processName ECALClustering --customise_commands 'from RecoEgamma.EgammaElectronProducers.gsfElectrons_cfi import * \n from RecoEgamma.EgammaPhotonProducers.gedPhotons_cfi import *' -n 10

#Produce cfg for running RecoSimDumper on RAW with Mustache using Data2022 with ZeeSkim: cmsDriver.py --python_filename RecoSimDumper_fromRAW_Mustache_Data2022_cfg.py --eventcontent MINIAOD --datatier MINIAOD --fileout file:step3.root --conditions 125X_dataRun3_relval_v4 --step RAW2DIGI,L1Reco,RECO,PAT --geometry DB:Extended --filein root://cms-xrd-global.cern.ch//store/data/Run2022F/EGamma/RAW-RECO/ZElectron-PromptReco-v1/000/360/390/00000/00d2f8c6-db74-477a-b106-e723005f62ac.root --era Run3 --no_exec --data --customise RecoSimStudies/Dumpers/customize_recodumper.customize_baseline,RecoSimStudies/Dumpers/customize_recodumper.customize_zeeFilterDumperData,RecoSimStudies/Dumpers/customize_recodumper.customize_linkedObjectsStep --processName ECALClustering --customise_commands 'from RecoEgamma.EgammaElectronProducers.gsfElectrons_cfi import * \n from RecoEgamma.EgammaPhotonProducers.gedPhotons_cfi import *' -n 10

#Produce cfg for running RecoSimDumper on RAW with DeepSC_algoA using Data2022: cmsDriver.py --python_filename RecoSimDumper_fromRAW_DeepSC_algoA_Data2022_cfg.py --eventcontent MINIAOD --datatier MINIAOD --fileout file:step3.root --conditions 125X_dataRun3_relval_v4 --step RAW2DIGI,L1Reco,RECO,PAT --geometry DB:Extended --filein root://cms-xrd-global.cern.ch//store/data/Run2022F/EGamma/RAW-RECO/ZElectron-PromptReco-v1/000/360/390/00000/00d2f8c6-db74-477a-b106-e723005f62ac.root --era Run3 --no_exec --data --customise RecoSimStudies/Dumpers/customize_recodumper.customize_baseline,RecoSimStudies/Dumpers/customize_recodumper.customize_dumperStepData,RecoSimStudies/Dumpers/customize_recodumper.customize_linkedObjectsStep --processName ECALClustering --customise_commands 'from RecoEgamma.EgammaElectronProducers.gsfElectrons_cfi import * \n from RecoEgamma.EgammaPhotonProducers.gedPhotons_cfi import *'  --procModifier ecal_deepsc -n 10

#Produce cfg for running RecoSimDumper on RAW with DeepSC_algoB using Data2022: cmsDriver.py --python_filename RecoSimDumper_fromRAW_DeepSC_algoB_Data2022_cfg.py --eventcontent MINIAOD --datatier MINIAOD --fileout file:step3.root --conditions 125X_dataRun3_relval_v4 --step RAW2DIGI,L1Reco,RECO,PAT --geometry DB:Extended --filein root://cms-xrd-global.cern.ch//store/data/Run2022F/EGamma/RAW-RECO/ZElectron-PromptReco-v1/000/360/390/00000/00d2f8c6-db74-477a-b106-e723005f62ac.root --era Run3 --no_exec --data --customise RecoSimStudies/Dumpers/customize_recodumper.customize_baseline,RecoSimStudies/Dumpers/customize_recodumper.customize_dumperStepData,RecoSimStudies/Dumpers/customize_recodumper.customize_linkedObjectsStep,RecoSimStudies/Dumpers/customize_recodumper.customize_deepSCAlgoB --processName ECALClustering --customise_commands 'from RecoEgamma.EgammaElectronProducers.gsfElectrons_cfi import * \n from RecoEgamma.EgammaPhotonProducers.gedPhotons_cfi import *'  --procModifier ecal_deepsc -n 10


def customize_baseline(process):

 '''
 process.mySCReg = cms.ESSource("PoolDBESSource",
     toGet = cms.VPSet(
       cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("pfscecal_EBCorrection_offline_v2"),
         tag = cms.string("pfscecal_EBCorrection_offline_v2_2022GammasMustacheSC"),
         connect = cms.string("sqlite_file:scReg_2022GammasMustacheSC.db")),
       cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("pfscecal_EECorrection_offline_v2"),
         tag = cms.string("pfscecal_EECorrection_offline_v2_2022GammasMustacheSC"),
         connect = cms.string("sqlite_file:scReg_2022GammasMustacheSC.db")),
      cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("pfscecal_EBUncertainty_offline_v2"),
         tag = cms.string("pfscecal_EBUncertainty_offline_v2_2022GammasMustacheSC"),
         connect = cms.string("sqlite_file:scReg_2022GammasMustacheSC.db")),
      cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("pfscecal_EEUncertainty_offline_v2"),
         tag = cms.string("pfscecal_EEUncertainty_offline_v2_2022GammasMustacheSC"),
         connect = cms.string("sqlite_file:scReg_2022GammasMustacheSC.db")),
     )
 )
 '''

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

 process.myPhoReg = cms.ESSource("PoolDBESSource",
     toGet = cms.VPSet(
      cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("photon_eb_ecalOnly_5To300_0p2To2_mean"),
         tag = cms.string("photon_eb_ecalOnly_5To300_0p2To2_mean_2022PhotonsMustacheSC"),
         connect = cms.string("sqlite_file:phoReg_2022PhotonsMustacheSC.db")),
      cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("photon_ee_ecalOnly_5To300_0p2To2_mean"),
         tag = cms.string("photon_ee_ecalOnly_5To300_0p2To2_mean_2022PhotonsMustacheSC"),
         connect = cms.string("sqlite_file:phoReg_2022PhotonsMustacheSC.db")),
      cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("photon_eb_ecalOnly_5To300_0p0002To0p5_sigma"),
         tag = cms.string("photon_eb_ecalOnly_5To300_0p0002To0p5_sigma_2022PhotonsMustacheSC"),
         connect = cms.string("sqlite_file:phoReg_2022PhotonsMustacheSC.db")),
      cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("photon_ee_ecalOnly_5To300_0p0002To0p5_sigma"),
         tag = cms.string("photon_ee_ecalOnly_5To300_0p0002To0p5_sigma_2022PhotonsMustacheSC"),
         connect = cms.string("sqlite_file:phoReg_2022PhotonsMustacheSC.db")),
     ) 
 )
 process.es_prefer_phoReg = cms.ESPrefer("PoolDBESSource","myPhoReg")

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

 #update electron_pfDNN
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

 #update photon_pfDNN
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

 return process

def customize_addConditions(process):

 process.myTPGLinear = cms.ESSource("PoolDBESSource",
     connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
     toGet = cms.VPSet(
         cms.PSet(
             record = cms.string('EcalTPGLinearizationConstRcd'),
             tag = cms.string('EcalTPGLinearizationConst_UL_2018_mc_EE_BTCP_116_SIC_1')
         )
     )
 )
 process.es_prefer_tpg = cms.ESPrefer("PoolDBESSource","myTPGLinear")

 process.myLaser = cms.ESSource("PoolDBESSource",
     connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
     toGet = cms.VPSet(
         cms.PSet(
             record = cms.string('EcalLaserAPDPNRatiosRcd'),
             tag = cms.string('EcalLaserAPDPNRatios_UL_2018_mc_3sigma_v2')
         )
     )
 )
 process.es_prefer_laser = cms.ESPrefer("PoolDBESSource","myLaser")

 process.myNoise = cms.ESSource("PoolDBESSource",
     connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
     toGet = cms.VPSet(
         cms.PSet(
             record = cms.string('EcalPedestalsRcd'),
             tag = cms.string('EcalPedestals_mid2021_235fb_mc')
         )
     )
 )
 process.es_prefer_noise = cms.ESPrefer("PoolDBESSource","myNoise")

 process.myPFRechitThres = cms.ESSource("PoolDBESSource",
     connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
     toGet = cms.VPSet(
         cms.PSet(
             record = cms.string('EcalPFRecHitThresholdsRcd'),
             tag = cms.string('EcalPFRecHitThresholds_UL_2018_2e3sig')
         )
     )
 )
 process.es_prefer_pfRechitThres = cms.ESPrefer("PoolDBESSource","myPFRechitThres")

 return process

def customize_addIdealICs(process):

 process.myICs = cms.ESSource("PoolDBESSource",
     connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
     toGet = cms.VPSet(
         cms.PSet(
             record = cms.string('EcalIntercalibConstantsRcd'),
             tag = cms.string('EcalIntercalibConstants_MC_Digi_2018')
         )
     )
 )
 process.es_prefer_icReco = cms.ESPrefer("PoolDBESSource","myICs")
 
 return process

def customize_nonStdMC(process):

 process.prunedTpClusterProducer.throwOnMissingCollections = cms.bool(False)
 process.prunedTrackMCMatch.throwOnMissingTPCollection = cms.bool(False)
 process.MINIAODSIMoutput.outputCommands.extend(['drop *_muonSimClassifier_*_*'])
 process.muonSimClassificationByHitsTask.remove(process.muonSimClassifier)

 return process

def customize_deepSCAlgoA(process):

 process.particleFlowSuperClusterECAL.deepSuperClusterConfig.collectionStrategy =  'Cascade'
 return process

def customize_deepSCAlgoB(process):

 process.particleFlowSuperClusterECAL.deepSuperClusterConfig.collectionStrategy =  'CollectAndMerge'
 return process

def customize_deepSCAlgoC(process):

 process.particleFlowSuperClusterECAL.deepSuperClusterConfig.collectionStrategy =  'SeedsFirst'
 return process

def customize_deepSCAlgoD(process):

 process.particleFlowSuperClusterECAL.deepSuperClusterConfig.collectionStrategy =  'CascadeHighest'
 return process

def customize_linkedObjectsStep(process):

 process.linkedObjects = cms.EDProducer("PATObjectCrossLinker",
    jets=cms.InputTag("slimmedJets"),
    muons=cms.InputTag("slimmedMuons"),
    electrons=cms.InputTag("slimmedElectrons"),
    lowPtElectrons=cms.InputTag("slimmedLowPtElectrons"),
    taus=cms.InputTag("slimmedTaus"),
    photons=cms.InputTag("slimmedPhotons"),
 )
 process.linkedObjects_step = cms.Path(process.linkedObjects)
 process.schedule += cms.Schedule(process.linkedObjects_step)

 return process

def customize_dumperStepMC(process):

 #dumper step
 process.TFileService = cms.Service("TFileService",
                                   #fileName = cms.string(options.outputFile)
                                   fileName = cms.string('output.root')
 )
 process.load('RecoSimStudies.Dumpers.RecoSimDumperMC_cfi')
 process.dumper_step = cms.Path(process.recosimdumper)

 process.schedule += cms.Schedule(process.dumper_step)

 return process

def customize_dumperStepData(process):
 
 #dumper step
 process.TFileService = cms.Service("TFileService",
                                   #fileName = cms.string(options.outputFile)
                                   fileName = cms.string('output.root')
 )
 process.load('RecoSimStudies.Dumpers.RecoSimDumperData_cfi')
 process.dumper_step = cms.Path(process.recosimdumper)

 process.schedule += cms.Schedule(process.dumper_step)

 return process

def customize_zeeFilterDumperMC(process):

 #dumper step
 process.TFileService = cms.Service("TFileService",
                                   #fileName = cms.string(options.outputFile)
                                   fileName = cms.string('output.root')
 )
 process.load("PhysicsTools.PatAlgos.patSequences_cff")

 process.cleanPatTaus.preselection = cms.string('tauID("decayModeFinding") > 0.5 & tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits") > 0.5 & tauID("againstMuonTight3") > 0.5')
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
 process.load('RecoSimStudies.Dumpers.RecoSimDumperMC_cfi')
 process.filter_step = cms.Path(
    #process.AllEvents
    process.ZeeSkimFilterSeq
    *process.AllEvents  
    #*process.hltLevel1GTSeed ---> Don't use it
    *process.hltHighLevel
    *process.goodPrimaryVertices
    *process.highetele
    *process.highetFilter
    *process.FilterReRECOEvents   
    *process.FilterPatDefaultSequenceEvents 
    *process.recosimdumper
 )


 process.schedule += cms.Schedule(process.filter_step)

 return process

def customize_zeeFilterDumperData(process):

 #dumper step
 process.TFileService = cms.Service("TFileService",
                                   #fileName = cms.string(options.outputFile)
                                   fileName = cms.string('output.root')
 )
 process.load("PhysicsTools.PatAlgos.patSequences_cff")

 process.cleanPatTaus.preselection = cms.string('tauID("decayModeFinding") > 0.5 & tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits") > 0.5 & tauID("againstMuonTight3") > 0.5')
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
 process.load('RecoSimStudies.Dumpers.RecoSimDumperData_cfi')
 process.filter_step = cms.Path(
    #process.AllEvents
    process.ZeeSkimFilterSeq
    *process.AllEvents  
    #*process.hltLevel1GTSeed ---> Don't use it
    *process.hltHighLevel
    *process.goodPrimaryVertices
    *process.highetele
    *process.highetFilter
    *process.FilterReRECOEvents   
    *process.FilterPatDefaultSequenceEvents 
    *process.recosimdumper
 )


 process.schedule += cms.Schedule(process.filter_step)

 return process




 
