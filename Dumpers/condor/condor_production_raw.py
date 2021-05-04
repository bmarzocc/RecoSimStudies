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
parser.add_argument("-g", "--gun", type=str, help="GEN-SIM cfg", required=True)
parser.add_argument("-S", "--seed", type=int, help="Seed offset", default=0)
parser.add_argument("-p", "--proxy", type=str, help="Proxy key", required=False)
parser.add_argument("--redo", action="store_true", default=False, help="Redo all files")
args = parser.parse_args()


# Prepare condor jobs
condor = '''executable              = run_script.sh
output                  = output/strips.$(ClusterId).$(ProcId).out
error                   = error/strips.$(ClusterId).$(ProcId).err
log                     = log/strips.$(ClusterId).log
transfer_input_files    = run_script.sh
on_exit_remove          = (ExitBySignal == False) && (ExitCode == 0)
periodic_release        = (NumJobStarts < 3) && ((CurrentTime - EnteredCurrentStatus) > (60*60))

+JobFlavour             = "{queue}"
+AccountingGroup = "group_u_CMS.CAF.ALCA"
queue arguments from arguments.txt
'''

condor = condor.replace("{queue}", args.queue)
user = os.environ["USER"]

script = '''#!/bin/sh -e

export X509_USER_PROXY=DIR/PROXY
voms-proxy-info

cp -r {cmssw_loc} ./
cd {cmssw_file}/src

echo -e "evaluate"
eval `scramv1 ru -sh`
export HOME='/afs/cern.ch/user/{user1}/{user}'

JOBID=$1; shift; 
OUTPUTFILE=$1;
NEVENTS=$2;
SEED1=$3
SEED2=$4
SEED3=$5
SEED4=$6
SEED5=$7
SEED6=$8
SEED7=$9


cd RecoSimStudies/Dumpers/test

echo -e "cmsRun..";
echo -e ">>> STEP1";

cmsRun GEN_SIM_CFG jobid=$JOBID  maxEvents=$NEVENTS seed1=$SEED1 seed2=$SEED2 seed3=$SEED3 seed4=$SEED4 seed5=$SEED5 

echo -e ">>> STEP2";
cmsRun step2_DIGI_L1_DIGI2RAW_HLT_PU_Run3_2021.py jobid=$JOBID seed1=$SEED6 seed2=$SEED7 

xrdcp --nopbar step2.root root://eoscms.cern.ch/${OUTPUTFILE}_step2.root;

echo -e "DONE";
'''

script = script.replace("{eosinstance}", args.eos)
script = script.replace("{user1}", user[:1])
script = script.replace("{user}", user)
cmssw_file = args.cmssw.split("/")[-1]
script = script.replace("{cmssw_loc}", args.cmssw)
script = script.replace("{cmssw_file}", cmssw_file)
script = script.replace("DIR", os.getcwd())
script = script.replace("PROXY", args.proxy)

if args.gun == "Photons": script = script.replace("GEN_SIM_CFG", "GammasGunPt1-100_pythia8_cfi_GEN_SIM.py")
elif args.gun == "Electrons": script = script.replace("GEN_SIM_CFG", "ElectronsGunPt1-100_pythia8_cfi_GEN_SIM.py")
elif args.gun == "Jets": script = script.replace("GEN_SIM_CFG", "JetsGunPt1-100_EMEnriched_pythia8_cfi_GEN_SIM.py")
elif args.gun == "QCD": script = script.replace("GEN_SIM_CFG", "QCD_Pt-15to7000_TuneCUETP8M1_Flat_14TeV-pythia8_cfi_GEN_SIM.py")
else: 
  print "Wrong GUN option, please use: 'Photons' or 'Electrons' or 'Jets' or 'QCD'"
  exit() 

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

            arguments.append("{} {} {} {} {} {} {} {} {} {}".format(
                jobid,outputfile,args.split, 
                jobid+args.seed+1,jobid+args.seed+2, 
                jobid+args.seed+3,jobid+args.seed+4,
                jobid+args.seed+5,jobid+args.seed+6,
                jobid+args.seed+7))

print("Njobs: ", len(arguments))
    
with open("condor_job.txt", "w") as cnd_out:
    cnd_out.write(condor)

with open("arguments.txt", "w") as args:
    args.write("\n".join(arguments))

with open("run_script.sh", "w") as rs:
    rs.write(script)

#os.system("condor_submit condor_job.txt")




