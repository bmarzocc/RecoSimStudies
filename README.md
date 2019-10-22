# RecoSimStudies

This is a customised version of https://github.com/bmarzocc/RecoSimStudies
In this repository, you find all the necessary codes for the production of samples, from the generation to the dumping.

## Installation

First installation:
```
cmsrel CMSSW_10_6_1_patch1
```
If you get an error, make sure that the remote machine on which you are working on is new enough to be compatible with the CMSSW_10_6_1_patch1 release. At the moment of writing, this release only works for machines with SL7 architecture at least, and one has typically to ask for a t3ui07 account to the PSI-T3 administrators.

```
cd CMSSW_10_6_1_patch1/src/
cmsenv
git cms-init
git branch -b base
git cms-merge-topic bmarzocc:PR_CaloParticles
git cms-merge-topic bmarzocc:PR_EcalPFSeedingThresholds
git cms-merge-topic bmarzocc:PR_ParticleGuns
git clone git@github.com:pfclustering/RecoSimStudies.git
scram b -j 5
```

In case you want to interact with the Storage Element, don't forget to set up your proxy:
```    
voms-proxy-init --voms cms --valid 186:00
```

## Workflow

### Development of ```cmssw```
Developments Development of ```cmssw``` by members of pfclustering team are done via fork of cmssw.
Currently there are three topics that can be changed (see above), under bmarzocc repo.
Changes to a given topic are pushed first to own fork, and then PR is done https://github.com/bmarzocc/cmssw/tree/<TOPIC_BRANCH>

*IMPORTANT NOTE* B. is still using an old release! So if you try to push, will get also release changes in the PR

The developments should be tested in the full (meaning three topics) configuration, but only the commits relevant to a given topic should be pushed
to the relevant topic, with the following workflow:

* within the same area you usually work on, create new local branch with following convention and do developments:
```
cd CMSSW_X_Y_Z/src
git checkout -b mg-PR_<TOPIC>
git cms-merge-topic bmarzocc:<TOPIC>
git remote add my-cmssw git@github.com:mgratti/cmssw.git # only first time
```
* make the relevant changes (which should already have been tested)
```
git add bla.cpp
git commit -m "bla" 
```
* push to remote branch, with same name convention:
```
git push my-cmssw mg-PR_<TOPIC>
```
* pull request to relevant topic branch under bmarzocc repo
* once PR has been accepted, go back to `base` branch, and DELETE both local and remote branches
```
git checkout base
git branch -d mg-PR_<topic>
git push my-cmssw --delete mg-PR_<topic>
```
* re-do all the relevant merge-topic 

More information and tricks on how to work with cmssw and github here: http://cms-sw.github.io/faq.html


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


