# RecoSimStudies

1) Install:

    * scram project CMSSW_12_3_0_pre5
    * cd CMSSW_12_3_0_pre5/src/
    * cmsenv
    * git cms-init
    * git cms-merge-topic valsdav:GraphSC_CMSSW_12_3_0_pre5 #if you produce DeepSC in RECO step
    * git clone https://github.com/bmarzocc/RecoSimStudies
    * cd RecoSimStudies
    * git checkout CMSSW_12_3_0
    * cd -
    * scram b -j 5

2) Produce RECO:

    * cd RecoSimStudies/Dumpers/test/
    * cmsRun step3_RAW2DIGI_L1Reco_RECO_RECOSIM_EI_PAT_Run3_2021.py

3) Run general dumper (per crystal, CaloParticle, PFcluster, superCluster infos) on a RECO sample (produced in the previous steps):
    
    * cd RecoSimStudies/Dumpers/
    * cmsRun python/RecoSimDumper_cfg.py

4) Run on condor (set properly the voms-key in the cfg before):

    * cd RecoSimStudies/Dumpers/condor/
    * module load lxbatch/tzero
    * voms-proxy-init --voms cms --valid 168:00
    * . proxy.sh x509up_u35923
    * python condor_production_recodumper.py -i /eos/cms/store/group/dpg_ecal/alca_ecalcalib/bmarzocc/Clustering/FourJetsGunPt1-100_EMEnriched_pythia8_PU_Flat55To75_14TeV_112X_mcRun3_2021_realistic_v15_GEN-SIM-RAW -o /eos/cms/store/group/dpg_ecal/alca_ecalcalib/bmarzocc/Clustering/FourJetsGunPt1-100_EMEnriched_pythia8_PU_Flat55To75_14TeV_112X_mcRun3_2021_realistic_v16_Dumper -c /afs/cern.ch/work/b/bmarzocc/Clustering_stdSeeding/CMSSW_11_2_2_patch1 -d RecoSimDumper -e EOS -q tomorrow -p x509up_u35923 #if input sample on EOS
    * python condor_production_recodumper.py -D FourJetsGunPt1-100_EMEnriched_StdMix_Flat55To75_14TeV_112X_mcRun3_2021_realistic_v15_GEN-SIM-RAW/dpg_ecal-GEN-SIM-RAW-993992ae0f1e290312fe9cf793f811bc/USER -o /eos/cms/store/group/dpg_ecal/alca_ecalcalib/bmarzocc/Clustering/FourJetsGunPt1-100_EMEnriched_pythia8_PU_Flat55To75_14TeV_112X_mcRun3_2021_realistic_v16_Dumper -c /afs/cern.ch/work/b/bmarzocc/Clustering_stdSeeding/CMSSW_11_2_2_patch1 -d RecoSimDumper -e EOS -q tomorrow -p x509up_u35923 #if input sample on DAS 


