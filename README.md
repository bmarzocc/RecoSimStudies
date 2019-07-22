# RecoSimStudies

In this repository, you find all the necessary codes for the production of samples, from the generation to the dumping.

## Installation

    * cmsrel CMSSW_10_6_0
If you get an error, make sure that the remote machine on which you are working on is new enough to be compatible with the CMSSW_10_6_0 release. At the moment of writing, this release only works for machines with SL7 architecture at least, and one has typically to ask for a t3ui07 account to the PSI-T3 administrators.

    * cd CMSSW_10_6_0/src/
    * cmsenv
    * git cms-init
    * git cms-merge-topic bmarzocc:RecoSimStudies 
    * git clone git@github.com:pfclustering/RecoSimStudies.git
    * scram b -j 5

## Generation
### GEN_SIM production
```    
cd RecoSimStudies/Dumpers/test/
```

Note: the following instructions will allow you to generate single photon events in front of ECAL.

Choose the parameters you want in the "User's decision board" in launch_step1.sh. You may also want to change the output directories SERESULTDIR and TOPWORKDIR.

Finally, launch this bash script following the instructions written directly in the file. For instance, to run the file on the batch, do
```
sbatch -p wn -o logs/step1.out -e logs/step1.err -q long.q --ntasks=8 launch_step1.sh
```

### DIGI-RAW production

```                         
cd RecoSimStudies/Dumpers/test/
```

Choose the parameters you want in the "User's decision board" in launch_step2.sh. You may also want to change the output directories SERESULTDIR and TOPWORKDIR.

Finally, launch this bash script following the instructions written directly in the file. For instance, to run the file on the batch, do
```
sbatch -p wn -o logs/step2.out -e logs/step2.err -q long.q --ntasks=8 launch_step2.sh
```

### RECO production

```                         
cd RecoSimStudies/Dumpers/test/
```

Choose the parameters you want in the "User's decision board" in launch_step3.sh. You may also want to change the output directories SERESULTDIR and TOPWORKDIR.

Finally, launch this bash script following the instructions written directly in the file. For instance, to run the file on the batch, do
```
sbatch -p wn -o logs/step3.out -e logs/step3.err -q long.q --ntasks=8 launch_step3.sh
```


    
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

