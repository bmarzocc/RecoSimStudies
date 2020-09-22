# RecoSimStudies

1) Install in CMSSW_10_6_X:

    * scram project CMSSW_10_6_11
    * cd CMSSW_10_6_11/src/
    * cmsenv
    * git cms-init
    * git cms-merge-topic bmarzocc:PR_CaloParticles
    * git cms-merge-topic bmarzocc:PR_ParticleGuns
    * git cms-merge-topic bmarzocc:10_6_11_DeepSC_noParticleFlow #if you want to reco the deepSCs too 
    * git clone https://github.com/bmarzocc/RecoSimStudies
    * cd RecoSimStudies
    * git checkout 10_6_X
    * scram b -j 5

2) Produce GEN-SIM of standard PhotonGun:
    
    * cd RecoSimStudies/Dumpers/test/
    * cmsRun SingleGammaPt35_pythia8_cfi_GEN_SIM.py

3) Produce GEN-SIM of PhotonGun in front of ECAL:
    
    * cd RecoSimStudies/Dumpers/test/
    * cmsRun DoubleGammaE50_CloseEcal_cfi_GEN_SIM.py

4) Produce DIGI-RAW:
    
    In case of studies including PCaloHits uncomment line #80 to keep the PCaloHit collection. 

    * cd RecoSimStudies/Dumpers/test/
    * cmsRun step2_DIGI_L1_DIGI2RAW_HLT.py

5) Produce RECO:

    In case of studies including PCaloHits uncomment line #65 to keep the PCaloHit collection. 
    
    * cd RecoSimStudies/Dumpers/test/
    * cmsRun step3_RAW2DIGI_L1Reco_RECO_RECOSIM_EI_PAT_VALIDATION_DQM.py

6) Produce TwentyPhotons Run3_2021 Recos with condor:

    * cd RecoSimStudies/Dumpers/test/
    * python condor_production.py  -o /eos/cms/store/group/dpg_ecal/alca_ecalcalib/bmarzocc/Clustering/TwentyGammasGunPt1-100_pythia8_withPU_withTracker_Run3_2021/ -c /afs/cern.ch/work/b/bmarzocc/Clustering/CMSSW_10_6_4/ -q tomorrow -n 100000 -s 100 -e cms

7) Run general dumper (per crystal, PFcluster, superCluster infos) on a RECO sample (produced in the previous steps):
    
    * cd RecoSimStudies/Dumpers/
    * cmsRun python/RecoSimDumper_cfg.py


