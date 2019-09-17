#!/bin/bash

#NOTE: files with same path and name in the Storage Element are now overwritten by default

###############################################################
#How to PU_launch this script:
#if you want to run it locally: 
#
# source PU_launch_step3.sh
#
#if you want it to run it with slurm, two possibilities:
#
#1. either run it with the wn partition (will use by default the processor t3wn38, 2.6GHz, 16cores) 
#
# sbatch -p wn -o logs/step3_EB.out -e logs/step3_EB.err --job-name=step3_EB --ntasks=8 PU_launch_step3.sh
# sbatch -p wn -o logs/step3_EE.out -e logs/step3_EE.err --job-name=step3_EE --ntasks=8 PU_launch_step3.sh
#
#2. or use the gpu ressources
#
# sbatch --account=gpu_gres --partition=gpu --gres=gpu:2 --time=2-23:59 --job-name=step3_EB -o logs/step3_EB.out -e logs/step3_EB.err PU_launch_step3.sh 
# sbatch --account=gpu_gres --partition=gpu --gres=gpu:2 --time=2-23:59 --job-name=step3_EE -o logs/step3_EE.out -e logs/step3_EE.err PU_launch_step3.sh
#
# Add nodes: --nodes=4 (max for wn) --nodes=2 (max for gpu)
###############################################################


###############################################################
#                  User's decision board                      #

#Do you want to PU_launch the production for EE or EB
#(choose one at a time)
doEB=true
doEEP=false
doEEM=false

#Do you want to store the output file in your work are or in the 
#storage element? (choose one at a time)
saveWork=false
saveSE=true

#Choose name of the directory
DIRNAME="singlePhoton_closeECAL_0to100GeV_1K_PU"


###############################################################

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
JOBOPFILENAME="PU_step3_RAW2DIGI_L1Reco_RECO_RECOSIM_EI_PAT_VALIDATION_DQM.py"
FILENAME="step3.root"
INFILENAME="step2.root"

if [ "$saveSE" = true ] && [ "$saveWork" = false ] ; then
   SERESULTDIR="/pnfs/psi.ch/cms/trivcat/store/user/"$USER"/EcalProd/"$DIRNAME
fi
if [ "$saveWork" = true ] && [ "$saveSE" = false ] ; then
   SERESULTDIR="/t3home/"$USER"/CMSSW_10_6_0/src/RecoSimStudies/Dumpers/test/outputfiles/"$DIRNAME
fi



STARTDIR=`pwd`
TOPWORKDIR="/scratch/"$USER/
JOBDIR="gen_"$SERESULTDIR
WORKDIR=$TOPWORKDIR/test/$JOBDIR
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
echo "Going to copy cms driver and files with noise information per crystal"
cp $JOBOPFILENAME $WORKDIR/$JOBOPFILENAME
cp -r ../data/ $WORKDIR/../data

#in case we decide to clean the scratch area
echo ""
echo "Going to copy input file"
echo "May give an error if the file still exists in the scratch , but can be safely ignored"
if [ "$saveSE" = true ] && [ "$saveWork" = false ] ; then
   xrdcp $SEPREFIX/$SERESULTDIR/$INFILENAME $WORKDIR/$INFILENAME  
fi
if [ "$saveWork" = true ] && [ "$saveSE" = false ] ; then
   cp $SERESULTDIR/$INFILENAME $WORKDIR/$INFILENAME
fi

cd $WORKDIR

echo ""
echo "Going to run"

DATE_START=`date +%s`
cmsRun $JOBOPFILENAME 
DATE_END=`date +%s`

echo ""
echo "Finished running"

echo ""
echo "Content of current directory"
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


