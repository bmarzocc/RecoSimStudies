# RecoSimStudies

1) Install:

    * scram project CMSSW_10_6_0
    * cd CMSSW_10_6_0/src/
    * cmsenv
    * git cms-init
    * git cms-merge-topic bmarzocc:RecoSimStudies 
    * git clone https://github.com/bmarzocc/RecoSimStudies
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

6) Run the dumper on a RECO sample (produced in the previous steps):
    
    * cd RecoSimStudies/Dumpers/
    * cmsRun python/RecoSimDumper_cfg.py

7) Dumper tree info, a vector of caloParticles per event is saved with:
    
    * genParticles information;
    * caloParticle information;
    * In std::vector<float> information of hits of caloParticle simClusters;
    * In std::vector<float> energy of recHits associated with simClusters hit (-1. if the recHit is missing);
    * In std::vector<float> energy of pfRecHits associated with simClusters hit (-1. if the pfRecHit is missing);
    * In std::vector<float> energy of PFCluster hits associated with simClusters hit (-1. if the PFCluster hit is missing);
    * In std::vector<float> energy of SuperCluster hits associated with simClusters hit (-1. if the SuperCluster hit is missing);
    * PFClusters information;
    * SuperClusters information;
    * In std::map<int,int> mapping between std::vector<float> of simCluster hits and std::vector<float> of PFClusters
    * In std::map<int,int> mapping between std::vector<float> of superCluster hits and std::vector<float> of SuperClusters

