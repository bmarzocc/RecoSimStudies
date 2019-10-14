# RecoSimStudies

In this repository, you find all the necessary codes for the production of samples, from the generation to the dumping.


## Installation

First installation:

    * cmsrel CMSSW_10_6_1_patch1
If you get an error, make sure that the remote machine on which you are working on is new enough to be compatible with the CMSSW_10_6_1_patch1 release. At the moment of writing, this release only works for machines with SL7 architecture at least, and one has typically to ask for a t3ui07 account to the PSI-T3 administrators.

    * cd CMSSW_10_6_1_patch1/src/
    * cmsenv
    * git cms-init
    * git cms-merge-topic bmarzocc:RecoSimStudies 
    * git clone git@github.com:pfclustering/RecoSimStudies.git
    * scram b -j 5


In case you want to interact with the Storage Element, don't forget to set up your proxy:
```    
voms-proxy-init --voms cms --valid 186:00
```

## Workflow

### Development of ```cmssw```
Development of ```cmssw``` by members of pfclustering team are done via fork
One creates his/her own branch locally, pushes it to his/her own fork and then opens pull request for https://github.com/bmarzocc/cmssw/tree/RecoSimStudies

More information and tricks on how to work with cmssw and github here: http://cms-sw.github.io/faq.html

    * cd CMSSW_X_Y_Z/src
    * git cms-merge-topic bmarzocc:PR_CaloParticles
    * git cms-merge-topic bmarzocc:PR_EcalPFSeedingThresholds
    * git cms-merge-topic mgratti:bmarzocc/PR_ParticleGuns
    * git remote add my-cmssw git@github.com:mgratti/cmssw.git # only first time
    * git checkout -b RecoSimStudies-reco-mg # this is an example
    * git cms-addpkg CalibCalorimetry/EcalTrivialCondModules # this is an exmaple
    * developments (git add bla.cpp, git commit -m "bla") 
    * git push my-cmssw RecoSimStudies-reco-mg
    * open pull request to relevant topic branch under bmarzocc repo

### Development of ```RecoSimStudies```
Development of ```RecoSimStudies``` by members of pfclustering team happens within the pfclustering fork; 
each member has his/her own branch where to develop the new features. When development is over, he/she opens a pull request,
(ideally another member checks) and merges with master branch.

After master is in sync, developments of bmarzocc/RecoSimStudies are fetched via a pull request (from web page) with a brief comment about the changes.

## Generation
For all steps of generation until reco files 
```
cd Dumpers/test/ECALproductionHelper
```
See available options:
```
python prodHelper.py --help
```
Example commands in ```Dumpers/test/ECALproductionHelper/README.md```

### Dumper
```                         
cd RecoSimStudies/Dumpers/python/
```

Run the dumper on a RECO sample. Example of commands are given in Cfg_RecoSimDumper_cfg.py. For instance, use

```
cmsRun Cfg_RecoSimDumper_cfg.py outputFile=../test/outputfiles/dumpedFiles/dumped_singlePhoton_5k_EB.root inputFiles=file:../test/outputfiles/singlePhoton_5k_EB/step3.root
```


