# RecoSimStudies

1) Install:

    * scram project CMSSW_12_5_0
    * cd CMSSW_12_5_0/src/
    * cmsenv
    * git cms-init
    * git clone https://github.com/bmarzocc/RecoSimStudies
    * cd RecoSimStudies
    * git checkout CMSSW_12_5_0
    * cd -
    * scram b -j 5

2) Produce RECO:

    * cd RecoSimStudies/Dumpers/test/
    * cmsRun step3_RAW2DIGI_L1Reco_RECO_RECOSIM_EI_PAT_Run3_2021.py

3) Run general dumper (GenParticle, CaloParticle, PFcluster, superCluster, pat::objects infos) on a RECO sample (produced in the previous steps):
    
    #Turn on the collections that you want to save in RecoSimDumper_cfg.py, and select the proper data/MC GT
    * cd RecoSimStudies/Dumpers/
    * cmsRun python/RecoSimDumper_cfg.py

4) Run on condor RECO+Dumper (set properly the voms-key in the cfg before):

    #Turn on the collections that you want to save in RecoSimDumper_cfg.py, and select the proper data/MC GT in RecoSimDumper_cfi.py
    #Run only special 2021-MC RAW samples produced for these studies
    * cd RecoSimStudies/Dumpers/condor/
    * module load lxbatch/tzero
    * voms-proxy-init --voms cms --valid 168:00
    * . proxy.sh x509up_u35923
    * python condor_production_recodumper.py -i /eos/cms/store/group/dpg_ecal/alca_ecalcalib/bmarzocc/Clustering/FourJetsGunPt1-100_EMEnriched_pythia8_PU_Flat55To75_14TeV_112X_mcRun3_2021_realistic_v15_GEN-SIM-RAW -o /eos/cms/store/group/dpg_ecal/alca_ecalcalib/bmarzocc/Clustering/FourJetsGunPt1-100_EMEnriched_pythia8_PU_Flat55To75_14TeV_112X_mcRun3_2021_realistic_v16_Dumper -c /afs/cern.ch/work/b/bmarzocc/Clustering_stdSeeding/CMSSW_11_2_2_patch1 -d RecoSimDumper -e EOS -q tomorrow -p x509up_u35923 #if input sample on EOS
    * python condor_production_recodumper.py -D FourJetsGunPt1-100_EMEnriched_StdMix_Flat55To75_14TeV_112X_mcRun3_2021_realistic_v15_GEN-SIM-RAW/dpg_ecal-GEN-SIM-RAW-993992ae0f1e290312fe9cf793f811bc/USER -o /eos/cms/store/group/dpg_ecal/alca_ecalcalib/bmarzocc/Clustering/FourJetsGunPt1-100_EMEnriched_pythia8_PU_Flat55To75_14TeV_112X_mcRun3_2021_realistic_v16_Dumper -c /afs/cern.ch/work/b/bmarzocc/Clustering_stdSeeding/CMSSW_11_2_2_patch1 -d RecoSimDumper -e EOS -q tomorrow -p x509up_u35923 #if input sample on DAS 

5) Run on crab RECO+Dumper (set properly the voms-key in the cfg before):
    
    #Turn on the collections that you want to save in RecoSimDumper_cfg.py, and select the proper data/MC GT
    #To run on MC use: RecoSimDumper_fromRAW_*_cfg.py
    #To run on 2018 Data use: RecoSimDumper_fromRAW_*_Data2018_cfg.py
    #Select the filters you want in process.dumper_step while running on data
    #Prepare the proper CRAB cfg options
    * cd RecoSimStudies/Dumpers/crab/
    * crab submit crab_cfg_reco.py
    
