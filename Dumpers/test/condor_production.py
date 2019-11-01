import sys
import os
import argparse
import random
from math import *

with open("command.txt", "w") as of:
    of.write(" ".join(["python"]+sys.argv))

'''
This scripts runs hadd on single crystal files to 
group them in strips reading a DOF file
'''
parser = argparse.ArgumentParser()

#parser.add_argument("-f", "--files", type=str, help="input file", required=True)
parser.add_argument("-n", "--nevents", type=int, help="N events for each eta and energy", required=True)
parser.add_argument("-s", "--split", type=int, help="Divide N events in S jobs", default=10)
parser.add_argument("-o", "--outputdir", type=str, help="Outputdir", required=True)
parser.add_argument("-c", "--cmssw", type=str, help="Absolute path to CMSSW release", required=True)
parser.add_argument("-q", "--queue", type=str, help="Condor queue", default="longlunch", required=True)
parser.add_argument("-e", "--eos", type=str, default="user", help="EOS instance user/cms", required=False)
parser.add_argument("--redo", action="store_true", default=False, help="Redo all files")
args = parser.parse_args()


# Prepare condor jobs
condor = '''executable              = run_script.sh
output                  = output/strips.$(ClusterId).$(ProcId).out
error                   = error/strips.$(ClusterId).$(ProcId).err
log                     = log/strips.$(ClusterId).log
transfer_input_files    = run_script.sh

+JobFlavour             = "{queue}"
queue arguments from arguments.txt
'''

condor = condor.replace("{queue}", args.queue)
user = os.environ["USER"]

script = '''#!/bin/sh -e

export X509_USER_PROXY=/afs/cern.ch/user/{user1}/{user}/x509up_u35923
voms-proxy-info

cp -r {cmssw_loc} .
cd {cmssw_file}/src

echo -e "evaluate"
eval `scramv1 ru -sh`
export HOME='/afs/cern.ch/user/{user1}/{user}'

JOBID=$1;
OUTPUTFILE=$2;
NEVENTS=$3;
SEED1=$4
SEED2=$5
SEED3=$6
SEED4=$7
SEED5=$8
SEED6=$9

cd RecoSimStudies/Dumpers/test

echo -e "cmsRun..";
echo -e ">>> STEP1";
cmsRun GammasGunPt1-100_pythia8_cfi_GEN_SIM.py jobid=$JOBID  maxEvents=$NEVENTS seed1=$SEED1 seed2=$SEED2 seed3=$SEED3 seed4=$SEED4;

echo -e ">>> STEP2";
cmsRun step2_DIGI_DATAMIX_L1_DIGI2RAW_HLT_PreMix_Run3_2021.py jobid=$JOBID seed1=$SEED5 seed2=$SEED6

echo -e ">>> STEP3";
cmsRun step3_RAW2DIGI_L1Reco_RECO_RECOSIM_EI_Run3_2021.py

xrdcp --nopbar step3.root root://eos{eosinstance}.cern.ch/${OUTPUTFILE}_RECO_$JOBID.root;

echo -e "DONE";
'''

script = script.replace("{eosinstance}", args.eos)
script = script.replace("{user1}", user[:1])
script = script.replace("{user}", user)
cmssw_file = args.cmssw.split("/")[-2]
script = script.replace("{cmssw_loc}", args.cmssw)
script = script.replace("{cmssw_file}", cmssw_file)

arguments= []
if not os.path.exists(args.outputdir):
    os.makedirs(args.outputdir)

outputfiles = [args.outputdir +"/"+f for f in os.listdir(args.outputdir)]

jobid = 0
njobs = args.nevents // args.split

for ijob in range(njobs):
            jobid +=1
            outputfile = args.outputdir + "/cluster_job{}".format(jobid)
        
            if not args.redo and outputfile+".root" in outputfiles:
                continue
            arguments.append("{} {} {} {} {} {} {} {} {}".format(
                jobid,outputfile,args.split,
                random.randint(1,10000),random.randint(1,10000),random.randint(1,10000),
                random.randint(1,10000),random.randint(1,10000),random.randint(1,10000)))
        

print("Njobs: ", len(arguments))
    
with open("condor_job.txt", "w") as cnd_out:
    cnd_out.write(condor)

with open("arguments.txt", "w") as args:
    args.write("\n".join(arguments))

with open("run_script.sh", "w") as rs:
    rs.write(script)

#os.system("condor_submit condor_job.txt")




