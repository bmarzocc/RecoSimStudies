# RecoSimStudies

1) Install:

    * scram project CMSSW_11_2_2_patch1
    * cd CMSSW_11_2_2_patch1/src/
    * cmsenv
    * git cms-init
    * git cms-merge-topic bmarzocc:ParticleGuns_CMSSW_11_2_2_patch1 #if you produce particle-guns
    * git cms-merge-topic bmarzocc:CaloParticles_CMSSW_11_2_2_patch1 #if you produce RAW samples with caloParticles
    * git clone https://github.com/bmarzocc/RecoSimStudies
    * cd RecoSimStudies
    * git checkout CMSSW_11_2_2_patch1
    * cd -
    * scram b -j 5

2) Produce GEN-SIM (ParticleGuns and QCD):
    
    * cd RecoSimStudies/Dumpers/test/
    * cmsRun GammasGunPt1-100_pythia8_cfi_GEN_SIM.py maxEvents=10 #Photon gun
    * cmsRun ElectronsGunPt1-100_pythia8_cfi_GEN_SIM.py maxEvents=10 #Electron gun
    * cmsRun JetsGunPt1-100_EMEnriched_pythia8_cfi_GEN_SIM.py maxEvents=10 #Jet gun with EMEnriched filter
    * cmsRun JetsGunPt1-100_pythia8_cfi_GEN_SIM.py maxEvents=10 #Jet gun
    * cmsRun QCD_Pt-15to7000_TuneCUETP8M1_Flat_14TeV-pythia8_cfi_GEN_SIM.py maxEvents=10 #QCD

3) Produce DIGI-RAW (Standard Mixing):
    
    * cd RecoSimStudies/Dumpers/test/
    * cmsRun step2_DIGI_L1_DIGI2RAW_HLT_PU_Run3_2021.py

4) Produce RECO:

    * cd RecoSimStudies/Dumpers/test/
    * cmsRun step3_RAW2DIGI_L1Reco_RECO_RECOSIM_EI_PAT_Run3_2021.py

5) Run general dumper (per crystal, CaloParticle, PFcluster, superCluster infos) on a RECO sample (produced in the previous steps):
    
    * cd RecoSimStudies/Dumpers/
    * cmsRun python/RecoSimDumper_cfg.py

6) Run on condor (set properly the voms-key in the cfg before):

    * cd RecoSimStudies/Dumpers/condor/
    * voms-proxy-init --voms cms --valid 168:00
    * . proxy.sh x509up_u35923
    * python condor_production_raw.py -n 100000 -s 10 -o /eos/cms/store/group/dpg_ecal/alca_ecalcalib/bmarzocc/Clustering/FourJetsGunPt1-100_EMEnriched_pythia8_PU_Flat55To75_14TeV_112X_mcRun3_2021_realistic_v15_GEN-SIM-RAW/ -c /afs/cern.ch/work/b/bmarzocc/Clustering_stdSeeding/CMSSW_11_2_2_patch1 -e EOS -q tomorrow -g Jets -S 0 -p x509up_u35923
    * python condor_production_recodumper.py -i /eos/cms/store/group/dpg_ecal/alca_ecalcalib/bmarzocc/Clustering/FourJetsGunPt1-100_EMEnriched_pythia8_PU_Flat55To75_14TeV_112X_mcRun3_2021_realistic_v15_GEN-SIM-RAW -o /eos/cms/store/group/dpg_ecal/alca_ecalcalib/bmarzocc/Clustering/FourJetsGunPt1-100_EMEnriched_pythia8_PU_Flat55To75_14TeV_112X_mcRun3_2021_realistic_v16_Dumper -c /afs/cern.ch/work/b/bmarzocc/Clustering_stdSeeding/CMSSW_11_2_2_patch1 -d RecoSimDumper -e EOS -q tomorrow #if input sample on EOS
    * python condor_production_recodumper.py -D FourJetsGunPt1-100_EMEnriched_StdMix_Flat55To75_14TeV_112X_mcRun3_2021_realistic_v15_GEN-SIM-RAW/dpg_ecal-GEN-SIM-RAW-993992ae0f1e290312fe9cf793f811bc/USER -o /eos/cms/store/group/dpg_ecal/alca_ecalcalib/bmarzocc/Clustering/FourJetsGunPt1-100_EMEnriched_pythia8_PU_Flat55To75_14TeV_112X_mcRun3_2021_realistic_v16_Dumper -c /afs/cern.ch/work/b/bmarzocc/Clustering_stdSeeding/CMSSW_11_2_2_patch1 -d RecoSimDumper -e EOS -q tomorrow #if input sample on DAS 


