#!/bin/bash                                                                                                                                                                                       
set -e
BASE=$PWD
RELEASE_BASE=$CMSSW_BASE

export SCRAM_ARCH=slc7_amd64_gcc10
source /cvmfs/cms.cern.ch/cmsset_default.sh

echo "setting up CMSSW environment"
cd $RELEASE_BASE
eval `scram runtime -sh`
cd $BASE

echo "cmsRun -e -j FrameworkJobReport.xml gen_step.py jobNum="$1" "$2" "$3" outputName=genStep.root" $5
#cmsRun -e -j FrameworkJobReport.xml GammasGunPt1-500_pythia8_cfi_GEN_SIM.py jobNum=$1 $2 $3 outputName=gensimStep.root $5
cmsRun -e -j FrameworkJobReport.xml ElectronsGunPt1-500_pythia8_cfi_GEN_SIM.py jobNum=$1 $2 $3 outputName=gensimStep.root $5

echo "cmsRun -e -j FrameworkJobReport.xml sim_digi_premix_step.py "$3" inputName=genStep.root outputName=output.root"
cmsRun -e -j FrameworkJobReport.xml step2_DIGI_DATAMIX_L1_DIGI2RAW_HLT_PREMIX_Run3_2023.py $3 inputName=gensimStep.root outputName=output.root
rm gensimStep*.root
