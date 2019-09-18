#!/bin/bash

#NOTE: files with same path and name in the Storage Element are now overwritten by default

###############################################################
#How to launch this script:
#if you want to run it locally: 
#
# source launch_step1.sh
#
#if you want it to run it with slurm, two possibilities:
#
#1. either run it with the wn partition (will use by default the processor t3wn38, 2.6GHz, 16cores) 
#
# sbatch -p wn -o logs/step1_EB.out -e logs/step1_EB.err --job-name=step1_EB --ntasks=8 launch_step1.sh
# sbatch -p wn -o logs/step1_EE.out -e logs/step1_EE.err --job-name=step1_EE --ntasks=8 launch_step1.sh
#
#2. or use the gpu ressources
#
# sbatch --account=gpu_gres --partition=gpu --gres=gpu:2 --time=2-23:59 --job-name=step1_EB -o logs/step1_EB.out -e logs/step1_EB.err launch_step1.sh 
# sbatch --account=gpu_gres --partition=gpu --gres=gpu:2 --time=2-23:59 --job-name=step1_EE -o logs/step1_EE.out -e logs/step1_EE.err launch_step1.sh
#
# Add nodes: --nodes=4 (max for wn) --nodes=2 (max for gpu)
###############################################################



###############################################################
#                  User's decision board                      #

#Do you want to launch the production for EE or EB
#(choose one at a time)
doEB=true
doEEP=false
doEEM=false

#Do you want to store the output file in your work are or in the 
#storage element? (choose one at a time)
saveWork=false
saveSE=true

#Choose name of the directory
DIRNAME="singlePhoton_closeECAL_0to100GeV_150k_test2"


#Choose the number of events that you want to generate
#Please enter an EVEN number
#NEVENTS=150000
NEVENTS=50

#Choose the energy range of the photon gun
ETMIN=0.
ETMAX=100.
###############################################################


#Geometry configuration
if [ "$doEB" = true ] && [ "$doEEP" = false ] && [ "$doEEM" = false ] ; then
   RMIN=123.8
   RMAX=123.8
   ZMIN=-304.5
   ZMAX=304.5
fi
if [ "$doEEP" = true ] && [ "$doEEM" = false ] && [ "$doEB" = false ] ; then
   RMIN=31.6
   RMAX=171.1
   ZMIN=317.0
   ZMAX=317.0
fi
if [ "$doEEM" = true ] && [ "$doEEP" = false ] && [ "$doEB" = false ] ; then
   RMIN=31.6
   RMAX=171.1
   ZMIN=-317.0
   ZMAX=-317.0
fi


#in case of EE, half of the total number of events is produces in EEM, the other half in EEP
if [ "$doEE" = true ] && [ "$doEB" = false ] ; then
   NEVENTS=$((NEVENTS/1))
   echo $NEVENTS
fi

#File configuration
JOBOPFILENAME="step1_SingleGammaPt35_pythia8_cfi_GEN_SIM.py"

if [ "$doEB" = true ] && [ "$doEEP" = false ] && [ "$doEEP" = false ] ; then
   DIRNAME=$DIRNAME"_EB"
fi

if [ "$doEEP" = true ] && [ "$doEEM" = false ] && [ "$doEB" = false ] ; then
   DIRNAME=$DIRNAME"_EEP" 
fi

if [ "$doEEM" = true ] && [ "$doEEP" = false ] && [ "$doEB" = false ] ; then
   DIRNAME=$DIRNAME"_EEM" 
fi


# Job configuration
FILENAME="step1.root"

if [ "$saveSE" = true ] && [ "$saveWork" = false ] ; then
   SERESULTDIR="/pnfs/psi.ch/cms/trivcat/store/user/"$USER"/EcalProd/"$DIRNAME 
fi

if [ "$saveWork" = true ] && [ "$saveSE" = false ] ; then
   SERESULTDIR="/t3home/"$USER"/CMSSW_10_6_0/src/RecoSimStudies/Dumpers/test/outputfiles/"$DIRNAME
fi

STARTDIR=`pwd`
TOPWORKDIR="/scratch/"$USER
JOBDIR="gen_"$SERESULTDIR
WORKDIR=$TOPWORKDIR/$JOBDIR
SEPREFIX="root://t3dcachedb.psi.ch:1094/"


# Job instructions
source $VO_CMS_SW_DIR/cmsset_default.sh
shopt -s expand_aliases

echo ""
echo "Going to set up cms environment"
cd $STARTDIR
cmsenv

echo ""
echo "Going to create work dir"
mkdir -p $WORKDIR

echo "workdir: "
echo $WORKDIR

echo ""
echo "Going to create the output dir"
echo "May give an error if the directory already exists, which can be safely ignored"
if [ "$saveSE" = true ] && [ "$saveWork" = false ] ; then
   xrdfs t3dcachedb03.psi.ch mkdir $SERESULTDIR 
fi
if [ "$saveWork" = true ] && [ "$saveSE" = false ] ; then
   mkdir -p $SERESULTDIR
fi


echo ""
echo "Going to copy cms driver"

cp $JOBOPFILENAME $WORKDIR/$JOBOPFILENAME
cd $WORKDIR

echo ""
echo "Going to run"

DATE_START=`date +%s`
   cmsRun $JOBOPFILENAME maxEvents=$NEVENTS etmin=$ETMIN etmax=$ETMAX rmin=$RMIN rmax=$RMAX zmin=$ZMIN zmax=$ZMAX
DATE_END=`date +%s`

echo ""
echo "Finished running"

echo ""
echo "Content of current directory"
pwd
ls -al

echo ""
echo "Going to copy the output to the output directory"

if [ "$saveSE" = true ] && [ "$saveWork" = false ] ; then
   xrdcp -f $FILENAME $SEPREFIX/$SERESULTDIR/$FILENAME
fi

if [ "$saveWork" = true ] && [ "$saveSE" = false ] ; then
  cp $FILENAME $SERESULTDIR/$FILENAME
fi

echo ""
echo "Cleaning up $WORKDIR"
rm -rf $WORKDIR

RUNTIME=$((DATE_END-DATE_START))
echo "Wallclock running time: $RUNTIME s"

cd $STARTDIR


