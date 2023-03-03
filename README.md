# RecoSimStudies

1) Install:

    * scram project CMSSW_13_0_0_pre4
    * cd CMSSW_13_0_0_pre4/src/
    * cmsenv
    * git cms-init
    * git cms-checkout-topic valsdav:deepsc-opt-zeropad-13_0_0pre4
    * cd RecoEcal/EgammaClusterProducers/data
    * git clone -b deepsc-pfthres https://github.com/valsdav/RecoEcal-EgammaClusterProducers/
    * cd -
    * git clone https://github.com/bmarzocc/RecoSimStudies
    * cd RecoSimStudies
    * git checkout CMSSW_13_0_0
    * cd -
    * scram b -j 5

2) Produce RECO:

    * cd RecoSimStudies/Dumpers/test/
    * cmsRun step3_RAW2DIGI_L1Reco_RECO_RECOSIM_EI_PAT_Run3_2021.py

3) Run general dumper (GenParticle, CaloParticle, PFcluster, superCluster, pat::objects infos) on a RECO sample (produced in the previous steps):
    
    Turn on the collections that you want to save in RecoSimDumper_cfg.py, and select the proper data/MC GT

    * cd RecoSimStudies/Dumpers/
    * cmsRun python/RecoSimDumper_cfg.py

4) Run on condor RECO+Dumper (set properly the voms-key in the cfg before):

    Turn on the collections that you want to save in RecoSimDumper_cfg.py, and select the proper data/MC GT in RecoSimDumper_cfi.py

    Run only special 2021-MC RAW samples produced for these studies

    * cd RecoSimStudies/Dumpers/condor/
    * module load lxbatch/tzero
    * voms-proxy-init --voms cms --valid 168:00
    * . proxy.sh x509up_u35923
    * python condor_production_recodumper.py -i /eos/cms/store/group/dpg_ecal/alca_ecalcalib/bmarzocc/Clustering/FourJetsGunPt1-100_EMEnriched_pythia8_PU_Flat55To75_14TeV_112X_mcRun3_2021_realistic_v15_GEN-SIM-RAW -o /eos/cms/store/group/dpg_ecal/alca_ecalcalib/bmarzocc/Clustering/FourJetsGunPt1-100_EMEnriched_pythia8_PU_Flat55To75_14TeV_112X_mcRun3_2021_realistic_v16_Dumper -c /afs/cern.ch/work/b/bmarzocc/Clustering_stdSeeding/CMSSW_11_2_2_patch1 -d RecoSimDumper -e EOS -q tomorrow -p x509up_u35923 #if input sample on EOS
    * python condor_production_recodumper.py -D FourJetsGunPt1-100_EMEnriched_StdMix_Flat55To75_14TeV_112X_mcRun3_2021_realistic_v15_GEN-SIM-RAW/dpg_ecal-GEN-SIM-RAW-993992ae0f1e290312fe9cf793f811bc/USER -o /eos/cms/store/group/dpg_ecal/alca_ecalcalib/bmarzocc/Clustering/FourJetsGunPt1-100_EMEnriched_pythia8_PU_Flat55To75_14TeV_112X_mcRun3_2021_realistic_v16_Dumper -c /afs/cern.ch/work/b/bmarzocc/Clustering_stdSeeding/CMSSW_11_2_2_patch1 -d RecoSimDumper -e EOS -q tomorrow -p x509up_u35923 #if input sample on DAS 

    Run JPsi-Dumper on 2022 data RAW: 

    * python condor_production_jpsidumper.py -D ParkingDoubleElectronLowMass0/Run2022G-v1/RAW -o /eos/cms/store/group/dpg_ecal/alca_ecalcalib/bmarzocc/Clustering/ParkingDoubleElectronLowMass0_Run2022G-v1_DeepSC_algoA_jpsiDumper -c /afs/cern.ch/work/b/bmarzocc/Clustering_stdSeeding/DeepSC_JPsi_12_5_0 -e EOS -q tomorrow -p x509up_u35923 -s T0_CH_CERN_Disk 

5) Run on crab RECO+Dumper (set properly the voms-key in the cfg before):
    
    Turn on the collections that you want to save in python/RecoSimDumperData_cfi.py or python/RecoSimDumperMC_cfi.py, and put the proper GT in python/RecoSimDumperData_cfg.py or python/RecoSimDumperMC_cfg.py 

    The configs to be run on CRAB can be produced using cmsDriver.py and python/customize_recodumper.py (add new customizations if needed):

    Produce cfg for running RecoSimDumper on RAW with Mustache using nonStdMC: 

    * cmsDriver.py --python_filename RecoSimDumper_fromRAW_Mustache_NonStdMC_cfg.py --eventcontent MINIAODSIM --datatier MINIAODSIM --fileout file:step3.root --conditions 125X_mcRun3_2022_realistic_v4 --step RAW2DIGI,L1Reco,RECO,RECOSIM,PAT --geometry DB:Extended --filein root://cms-xrd-global.cern.ch//store/group/dpg_ecal/Clustering/FourElectronsGunPt1-100_pythia8_StdMix_Flat55To75_14TeV_112X_mcRun3_2021_realistic_v15/GEN-SIM-RAW-Reduced/210519_101349/0003/cluster_job3235_step2_3245.root --era Run3,ctpps_2018 --no_exec --mc --customise RecoSimStudies/Dumpers/customize_recodumper.customize_mustache,RecoSimStudies/Dumpers/customize_recodumper.customize_baseline,RecoSimStudies/Dumpers/customize_recodumper.customize_dumperStepMC,RecoSimStudies/Dumpers/customize_recodumper.customize_linkedObjectsStep,RecoSimStudies/Dumpers/customize_recodumper.customize_addConditions,RecoSimStudies/Dumpers/customize_recodumper.customize_nonStdMC --processName ECALClustering --customise_commands 'from RecoEgamma.EgammaElectronProducers.gsfElectrons_cfi import * \n from RecoEgamma.EgammaPhotonProducers.gedPhotons_cfi import * \n from SimMuon.MCTruth.muonSimClassificationByHits_cff import *' -n 10

    Produce cfg for running RecoSimDumper on RAW with Mustache using StdMC: 

    * cmsDriver.py --python_filename RecoSimDumper_fromRAW_Mustache_StdMC_cfg.py --eventcontent MINIAODSIM --datatier MINIAODSIM --fileout file:step3.root --conditions 125X_mcRun3_2022_realistic_v4 --step RAW2DIGI,L1Reco,RECO,RECOSIM,PAT --geometry DB:Extended --filein root://cms-xrd-global.cern.ch//store/group/dpg_ecal/Clustering/FourElectronsGunPt1-500_pythia8_PremixRun3_13p6TeV_122X_mcRun3_2021_realistic_v9_235fbNoise/GEN-SIM-RAW-CALO/221026_074347/0001/cluster_job99_step2_1778.root --era Run3 --no_exec --mc --customise RecoSimStudies/Dumpers/customize_recodumper.customize_mustache,RecoSimStudies/Dumpers/customize_recodumper.customize_baseline,RecoSimStudies/Dumpers/customize_recodumper.customize_dumperStepMC,RecoSimStudies/Dumpers/customize_recodumper.customize_linkedObjectsStep,RecoSimStudies/Dumpers/customize_recodumper.customize_addConditions --processName ECALClustering --customise_commands 'from RecoEgamma.EgammaElectronProducers.gsfElectrons_cfi import * \n from RecoEgamma.EgammaPhotonProducers.gedPhotons_cfi import *' -n 10

    Produce cfg for running RecoSimDumper on RAW with DeepSC_algoA using StdMC: 

    * cmsDriver.py --python_filename RecoSimDumper_fromRAW_DeepSC_algoA_cfg.py --eventcontent MINIAODSIM --datatier MINIAODSIM --fileout file:step3.root --conditions 125X_mcRun3_2022_realistic_v4 --step RAW2DIGI,L1Reco,RECO,RECOSIM,PAT --geometry DB:Extended --filein root://cms-xrd-global.cern.ch//store/group/dpg_ecal/Clustering/FourElectronsGunPt1-500_pythia8_PremixRun3_13p6TeV_122X_mcRun3_2021_realistic_v9_235fbNoise/GEN-SIM-RAW-CALO/221026_074347/0001/cluster_job99_step2_1778.root --era Run3 --no_exec --mc --customise RecoSimStudies/Dumpers/customize_recodumper.customize_baseline,RecoSimStudies/Dumpers/customize_recodumper.customize_deepSCAlgoA,RecoSimStudies/Dumpers/customize_recodumper.customize_dumperStepMC,RecoSimStudies/Dumpers/customize_recodumper.customize_linkedObjectsStep,RecoSimStudies/Dumpers/customize_recodumper.customize_addConditions --processName ECALClustering --customise_commands 'from RecoEgamma.EgammaElectronProducers.gsfElectrons_cfi import * \n from RecoEgamma.EgammaPhotonProducers.gedPhotons_cfi import *' --procModifier ecal_deepsc -n 10

    Produce cfg for running RecoSimDumper on RAW with DeepSC_algoB using StdMC: 

    * cmsDriver.py --python_filename RecoSimDumper_fromRAW_DeepSC_algoB_cfg.py --eventcontent MINIAODSIM --datatier MINIAODSIM --fileout file:step3.root --conditions 125X_mcRun3_2022_realistic_v4 --step RAW2DIGI,L1Reco,RECO,RECOSIM,PAT --geometry DB:Extended --filein root://cms-xrd-global.cern.ch//store/group/dpg_ecal/Clustering/FourElectronsGunPt1-500_pythia8_PremixRun3_13p6TeV_122X_mcRun3_2021_realistic_v9_235fbNoise/GEN-SIM-RAW-CALO/221026_074347/0001/cluster_job99_step2_1778.root --era Run3 --no_exec --mc --customise RecoSimStudies/Dumpers/customize_recodumper.customize_baseline,RecoSimStudies/Dumpers/customize_recodumper.customize_dumperStepMC,RecoSimStudies/Dumpers/customize_recodumper.customize_linkedObjectsStep,RecoSimStudies/Dumpers/customize_recodumper.customize_addConditions,RecoSimStudies/Dumpers/customize_recodumper.customize_deepSCAlgoB --processName ECALClustering --customise_commands 'from RecoEgamma.EgammaElectronProducers.gsfElectrons_cfi import * \n from RecoEgamma.EgammaPhotonProducers.gedPhotons_cfi import * ' --procModifier ecal_deepsc -n 10

    Produce cfg for running RecoSimDumper on RAW with DeepSC_algoC using StdMC: 

    * cmsDriver.py --python_filename RecoSimDumper_fromRAW_DeepSC_algoC_cfg.py --eventcontent MINIAODSIM --datatier MINIAODSIM --fileout file:step3.root --conditions 125X_mcRun3_2022_realistic_v4 --step RAW2DIGI,L1Reco,RECO,RECOSIM,PAT --geometry DB:Extended --filein root://cms-xrd-global.cern.ch//store/group/dpg_ecal/Clustering/FourElectronsGunPt1-500_pythia8_PremixRun3_13p6TeV_122X_mcRun3_2021_realistic_v9_235fbNoise/GEN-SIM-RAW-CALO/221026_074347/0001/cluster_job99_step2_1778.root --era Run3 --no_exec --mc --customise RecoSimStudies/Dumpers/customize_recodumper.customize_baseline,RecoSimStudies/Dumpers/customize_recodumper.customize_dumperStepMC,RecoSimStudies/Dumpers/customize_recodumper.customize_linkedObjectsStep,RecoSimStudies/Dumpers/customize_recodumper.customize_addConditions,RecoSimStudies/Dumpers/customize_recodumper.customize_deepSCAlgoC --processName ECALClustering --customise_commands 'from RecoEgamma.EgammaElectronProducers.gsfElectrons_cfi import * \n from RecoEgamma.EgammaPhotonProducers.gedPhotons_cfi import * ' --procModifier ecal_deepsc -n 10

    Produce cfg for running RecoSimDumper on RAW with DeepSC_algoD using StdMC: 

    * cmsDriver.py --python_filename RecoSimDumper_fromRAW_DeepSC_algoD_cfg.py --eventcontent MINIAODSIM --datatier MINIAODSIM --fileout file:step3.root --conditions 125X_mcRun3_2022_realistic_v4 --step RAW2DIGI,L1Reco,RECO,RECOSIM,PAT --geometry DB:Extended --filein root://cms-xrd-global.cern.ch//store/group/dpg_ecal/Clustering/FourElectronsGunPt1-500_pythia8_PremixRun3_13p6TeV_122X_mcRun3_2021_realistic_v9_235fbNoise/GEN-SIM-RAW-CALO/221026_074347/0001/cluster_job99_step2_1778.root --era Run3 --no_exec --mc --customise RecoSimStudies/Dumpers/customize_recodumper.customize_baseline,RecoSimStudies/Dumpers/customize_recodumper.customize_dumperStepMC,RecoSimStudies/Dumpers/customize_recodumper.customize_linkedObjectsStep,RecoSimStudies/Dumpers/customize_recodumper.customize_addConditions,RecoSimStudies/Dumpers/customize_recodumper.customize_deepSCAlgoD --processName ECALClustering --customise_commands 'from RecoEgamma.EgammaElectronProducers.gsfElectrons_cfi import * \n from RecoEgamma.EgammaPhotonProducers.gedPhotons_cfi import * ' --procModifier ecal_deepsc -n 10

    Produce cfg for running RecoSimDumper on RAW with Mustache using Data2018: 

    * cmsDriver.py --python_filename RecoSimDumper_fromRAW_Mustache_Data2018_cfg.py --eventcontent MINIAOD --datatier MINIAOD --fileout file:step3.root --conditions 124X_dataRun2_v2 --step RAW2DIGI,L1Reco,RECO,PAT --geometry DB:Extended --filein root://cms-xrd-global.cern.ch//store/data/Run2018C/EGamma/RAW/v1/000/320/026/00000/8CF35F1E-678D-E811-84C0-FA163E8CC774.root --era Run2_2018 --no_exec --data --customise RecoSimStudies/Dumpers/customize_recodumper.customize_mustache,RecoSimStudies/Dumpers/customize_recodumper.customize_baseline,RecoSimStudies/Dumpers/customize_recodumper.customize_dumperStepData,RecoSimStudies/Dumpers/customize_recodumper.customize_linkedObjectsStep --processName ECALClustering --customise_commands 'from RecoEgamma.EgammaElectronProducers.gsfElectrons_cfi import * \n from RecoEgamma.EgammaPhotonProducers.gedPhotons_cfi import *' -n 10

    Produce cfg for running RecoSimDumper on RAW with Mustache using Data2018 with ZeeSkim: 

    * cmsDriver.py --python_filename RecoSimDumper_fromRAW_Mustache_Data2018_cfg.py --eventcontent MINIAOD --datatier MINIAOD --fileout file:step3.root --conditions 124X_dataRun2_v2 --step RAW2DIGI,L1Reco,RECO,PAT --geometry DB:Extended --filein root://cms-xrd-global.cern.ch//store/data/Run2018C/EGamma/RAW/v1/000/320/026/00000/8CF35F1E-678D-E811-84C0-FA163E8CC774.root --era Run2_2018 --no_exec --data --customise RecoSimStudies/Dumpers/customize_recodumper.customize_mustache,RecoSimStudies/Dumpers/customize_recodumper.customize_baseline,RecoSimStudies/Dumpers/customize_recodumper.customize_zeeFilterDumperData,RecoSimStudies/Dumpers/customize_recodumper.customize_linkedObjectsStep --processName ECALClustering --customise_commands 'from RecoEgamma.EgammaElectronProducers.gsfElectrons_cfi import * \n from RecoEgamma.EgammaPhotonProducers.gedPhotons_cfi import *' -n 10

    Produce cfg for running RecoSimDumper on RAW with DeepSC_algoA using Data2018: 

    * cmsDriver.py --python_filename RecoSimDumper_fromRAW_DeepSC_algoA_Data2018_cfg.py --eventcontent MINIAOD --datatier MINIAOD --fileout file:step3.root --conditions 124X_dataRun2_v2 --step RAW2DIGI,L1Reco,RECO,PAT --geometry DB:Extended --filein root://cms-xrd-global.cern.ch//store/data/Run2018C/EGamma/RAW/v1/000/320/026/00000/8CF35F1E-678D-E811-84C0-FA163E8CC774.root --era Run2_2018 --no_exec --data --customise RecoSimStudies/Dumpers/customize_recodumper.customize_baseline,RecoSimStudies/Dumpers/customize_recodumper.customize_deepSCAlgoA,RecoSimStudies/Dumpers/customize_recodumper.customize_dumperStepData,RecoSimStudies/Dumpers/customize_recodumper.customize_linkedObjectsStep --processName ECALClustering --customise_commands 'from RecoEgamma.EgammaElectronProducers.gsfElectrons_cfi import * \n from RecoEgamma.EgammaPhotonProducers.gedPhotons_cfi import *'  --procModifier ecal_deepsc -n 10

    Produce cfg for running RecoSimDumper on RAW with DeepSC_algoB using Data2018: 

    * cmsDriver.py --python_filename RecoSimDumper_fromRAW_DeepSC_algoB_Data2018_cfg.py --eventcontent MINIAOD --datatier MINIAOD --fileout file:step3.root --conditions 124X_dataRun2_v2 --step RAW2DIGI,L1Reco,RECO,PAT --geometry DB:Extended --filein root://cms-xrd-global.cern.ch//store/data/Run2018C/EGamma/RAW/v1/000/320/026/00000/8CF35F1E-678D-E811-84C0-FA163E8CC774.root --era Run2_2018 --no_exec --data --customise RecoSimStudies/Dumpers/customize_recodumper.customize_baseline,RecoSimStudies/Dumpers/customize_recodumper.customize_dumperStepData,RecoSimStudies/Dumpers/customize_recodumper.customize_linkedObjectsStep,RecoSimStudies/Dumpers/customize_recodumper.customize_deepSCAlgoB --processName ECALClustering --customise_commands 'from RecoEgamma.EgammaElectronProducers.gsfElectrons_cfi import * \n from RecoEgamma.EgammaPhotonProducers.gedPhotons_cfi import *'  --procModifier ecal_deepsc -n 10

    Produce cfg for running RecoSimDumper on RAW with Mustache using Data2022: 

    * cmsDriver.py --python_filename RecoSimDumper_fromRAW_Mustache_Data2022_cfg.py --eventcontent MINIAOD --datatier MINIAOD --fileout file:step3.root --conditions 125X_dataRun3_relval_v4 --step RAW2DIGI,L1Reco,RECO,PAT --geometry DB:Extended --filein root://cms-xrd-global.cern.ch//store/data/Run2022F/EGamma/RAW-RECO/ZElectron-PromptReco-v1/000/360/390/00000/00d2f8c6-db74-477a-b106-e723005f62ac.root --era Run3 --no_exec --data --customise RecoSimStudies/Dumpers/customize_recodumper.customize_mustache,RecoSimStudies/Dumpers/customize_recodumper.customize_baseline,RecoSimStudies/Dumpers/customize_recodumper.customize_dumperStepData,RecoSimStudies/Dumpers/customize_recodumper.customize_linkedObjectsStep --processName ECALClustering --customise_commands 'from RecoEgamma.EgammaElectronProducers.gsfElectrons_cfi import * \n from RecoEgamma.EgammaPhotonProducers.gedPhotons_cfi import *' -n 10

    Produce cfg for running RecoSimDumper on RAW with Mustache using Data2022 with ZeeSkim: 

    * cmsDriver.py --python_filename RecoSimDumper_fromRAW_Mustache_Data2022_cfg.py --eventcontent MINIAOD --datatier MINIAOD --fileout file:step3.root --conditions 125X_dataRun3_relval_v4 --step RAW2DIGI,L1Reco,RECO,PAT --geometry DB:Extended --filein root://cms-xrd-global.cern.ch//store/data/Run2022F/EGamma/RAW-RECO/ZElectron-PromptReco-v1/000/360/390/00000/00d2f8c6-db74-477a-b106-e723005f62ac.root --era Run3 --no_exec --data --customise RecoSimStudies/Dumpers/customize_recodumper.customize_mustache,RecoSimStudies/Dumpers/customize_recodumper.customize_baseline,RecoSimStudies/Dumpers/customize_recodumper.customize_zeeFilterDumperData,RecoSimStudies/Dumpers/customize_recodumper.customize_linkedObjectsStep --processName ECALClustering --customise_commands 'from RecoEgamma.EgammaElectronProducers.gsfElectrons_cfi import * \n from RecoEgamma.EgammaPhotonProducers.gedPhotons_cfi import *' -n 10

    Produce cfg for running RecoSimDumper on RAW with DeepSC_algoA using Data2022: 

    * cmsDriver.py --python_filename RecoSimDumper_fromRAW_DeepSC_algoA_Data2022_cfg.py --eventcontent MINIAOD --datatier MINIAOD --fileout file:step3.root --conditions 125X_dataRun3_relval_v4 --step RAW2DIGI,L1Reco,RECO,PAT --geometry DB:Extended --filein root://cms-xrd-global.cern.ch//store/data/Run2022F/EGamma/RAW-RECO/ZElectron-PromptReco-v1/000/360/390/00000/00d2f8c6-db74-477a-b106-e723005f62ac.root --era Run3 --no_exec --data --customise RecoSimStudies/Dumpers/customize_recodumper.customize_baseline,RecoSimStudies/Dumpers/customize_recodumper.customize_deepSCAlgoA,RecoSimStudies/Dumpers/customize_recodumper.customize_dumperStepData,RecoSimStudies/Dumpers/customize_recodumper.customize_linkedObjectsStep --processName ECALClustering --customise_commands 'from RecoEgamma.EgammaElectronProducers.gsfElectrons_cfi import * \n from RecoEgamma.EgammaPhotonProducers.gedPhotons_cfi import *'  --procModifier ecal_deepsc -n 10

    Produce cfg for running RecoSimDumper on RAW with DeepSC_algoB using Data2022: 

    * cmsDriver.py --python_filename RecoSimDumper_fromRAW_DeepSC_algoB_Data2022_cfg.py --eventcontent MINIAOD --datatier MINIAOD --fileout file:step3.root --conditions 125X_dataRun3_relval_v4 --step RAW2DIGI,L1Reco,RECO,PAT --geometry DB:Extended --filein root://cms-xrd-global.cern.ch//store/data/Run2022F/EGamma/RAW-RECO/ZElectron-PromptReco-v1/000/360/390/00000/00d2f8c6-db74-477a-b106-e723005f62ac.root --era Run3 --no_exec --data --customise RecoSimStudies/Dumpers/customize_recodumper.customize_baseline,RecoSimStudies/Dumpers/customize_recodumper.customize_dumperStepData,RecoSimStudies/Dumpers/customize_recodumper.customize_linkedObjectsStep,RecoSimStudies/Dumpers/customize_recodumper.customize_deepSCAlgoB --processName ECALClustering --customise_commands 'from RecoEgamma.EgammaElectronProducers.gsfElectrons_cfi import * \n from RecoEgamma.EgammaPhotonProducers.gedPhotons_cfi import *'  --procModifier ecal_deepsc -n 10

    Produce cfg for producing standard MINIAOD from RAW with Mustache using Data2022: 

    * cmsDriver.py --python_filename MiniAOD_fromRAW_Mustache_Data2022_cfg.py --eventcontent MINIAOD --datatier MINIAOD --fileout file:step3.root --conditions 125X_dataRun3_relval_v4 --step RAW2DIGI,L1Reco,RECO,PAT --geometry DB:Extended --filein root://cms-xrd-global.cern.ch//store/data/Run2022F/EGamma/RAW-RECO/ZElectron-PromptReco-v1/000/360/390/00000/00d2f8c6-db74-477a-b106-e723005f62ac.root --era Run3 --no_exec --data --customise RecoSimStudies/Dumpers/customize_recodumper.customize_mustache,RecoSimStudies/Dumpers/customize_recodumper.customize_baseline --processName ECALClustering --customise_commands 'from RecoEgamma.EgammaElectronProducers.gsfElectrons_cfi import * \n from RecoEgamma.EgammaPhotonProducers.gedPhotons_cfi import *' -n 10

    Produce cfg for producing standard MINIAOD from RAW with DeepSC_algoA using Data2022: 

    * cmsDriver.py --python_filename MiniAOD_fromRAW_DeepSC_algoA_Data2022_cfg.py --eventcontent MINIAOD --datatier MINIAOD --fileout file:step3.root --conditions 125X_dataRun3_relval_v4 --step RAW2DIGI,L1Reco,RECO,PAT --geometry DB:Extended --filein root://cms-xrd-global.cern.ch//store/data/Run2022F/EGamma/RAW-RECO/ZElectron-PromptReco-v1/000/360/390/00000/00d2f8c6-db74-477a-b106-e723005f62ac.root --era Run3 --no_exec --data --customise RecoSimStudies/Dumpers/customize_recodumper.customize_baseline,RecoSimStudies/Dumpers/customize_recodumper.customize_deepSCAlgoA --processName ECALClustering --customise_commands 'from RecoEgamma.EgammaElectronProducers.gsfElectrons_cfi import * \n from RecoEgamma.EgammaPhotonProducers.gedPhotons_cfi import *'  --procModifier ecal_deepsc -n 10

    Prepare the proper CRAB cfg options

    * cd RecoSimStudies/Dumpers/crab/
    * crab submit crab_cfg_reco.py
    
