import sys
import os
import argparse
import random
import commands
from math import *

with open("command.txt", "w") as of:
    of.write(" ".join(["python"]+sys.argv))

'''
This scripts runs hadd on single crystal files to 
group them in strips reading a DOF file
'''
parser = argparse.ArgumentParser()

#parser.add_argument("-f", "--files", type=str, help="input file", required=True)
parser.add_argument("-i", "--inputdir", type=str, default="", help="Inputdir", required=False)
parser.add_argument("-D", "--das", type=str, default="", help="input DAS dataset", required=False)
parser.add_argument("-o", "--outputdir", type=str, help="Outputdir", required=True)
parser.add_argument("-c", "--cmssw", type=str, help="CMSSW tar", required=True)
parser.add_argument("-d", "---dumper", type=str, help="Dumper to run", required=True, default="RecoSimDumper")
parser.add_argument("-q", "--queue", type=str, help="Condor queue", default="longlunch", required=True)
parser.add_argument("-e", "--eos", type=str, default="user", help="EOS instance user/cms", required=False)
parser.add_argument("--redo", action="store_true", default=False, help="Redo all files")
args = parser.parse_args()

if (args.inputdir=="" and args.das=="") or (args.inputdir!="" and args.das!=""):
  print "ERROR: Give either inputdir, either DAS option"
  exit()  

# Prepare condor jobs
condor = '''executable              = run_script.sh
output                  = output/strips.$(ClusterId).$(ProcId).out
error                   = error/strips.$(ClusterId).$(ProcId).err
log                     = log/strips.$(ClusterId).log
transfer_input_files    = run_script.sh
on_exit_remove          = (ExitBySignal == False) && (ExitCode == 0)
periodic_release        = (NumJobStarts < 3) && ((CurrentTime - EnteredCurrentStatus) > (60*60))

+JobFlavour             = "{queue}"
+AccountingGroup        = "group_u_CMS.CAF.ALCA" 
queue arguments from arguments.txt

'''

condor = condor.replace("{queue}", args.queue)
user = os.environ["USER"]

script = '''#!/bin/sh -e

cp -r {cmssw_loc} ./
cd {cmssw_file}/src

echo -e "evaluate"
eval `scramv1 ru -sh`
export HOME='/afs/cern.ch/user/{user1}/{user}'

JOBID=$1;  
INPUTFILE=$2;
OUTPUTFILE=$3;

cd RecoSimStudies/Dumpers/test/
rm *.root

echo -e "cmsRun..";

echo -e ">>> STEP3";
cmsRun step3_RAW2DIGI_L1Reco_RECO_RECOSIM_EI_PAT_Run3_2021.py inputFile=${INPUTFILE}

echo -e "Running dumper.."

cd ..
cmsRun python/{dumper}_cfg.py inputFile=test/step3.root outputFile=output_${JOBID}.root


echo -e "Copying result to: $OUTPUTFILE";
xrdcp -f --nopbar  output_${JOBID}.root root://eoscms.cern.ch/${OUTPUTFILE};

echo -e "DONE";
'''

script = script.replace("{eosinstance}", args.eos)
script = script.replace("{user1}", user[:1])
script = script.replace("{user}", user)
cmssw_file = args.cmssw.split("/")[-1]
script = script.replace("{cmssw_loc}", args.cmssw)
script = script.replace("{cmssw_file}", cmssw_file)
script = script.replace("{dumper}", args.dumper)

arguments= []
if not os.path.exists(args.outputdir):
  os.makedirs(args.outputdir)

outputfiles = [args.outputdir +"/"+f for f in os.listdir(args.outputdir)]

inputfiles = ""
if args.inputdir!="":
  inputfiles = [ f for f in os.listdir(args.inputdir)]
if args.das!="":
  query = "dasgoclient --query='file dataset=/"+args.das+" instance=prod/phys03'"
  status, output = commands.getstatusoutput(query) 
  #print output
  inputfiles = output.split()

jobid = 0
for ifile in inputfiles:
    jobid +=1
    if args.inputdir!="": 
      inputfile = args.inputdir + "/" + ifile
      outputfile = args.outputdir + "/" + ifile[:-5] + "_output.root"
    else: 
      inputfile = "root://cms-xrd-global.cern.ch/" + ifile
      outputfile = args.outputdir + "/" + ifile.split('/')[-1]
    

    if not args.redo and outputfile in outputfiles:
        continue
    
    arguments.append("{} {} {}".format(jobid,inputfile,outputfile))

print("Njobs: ", len(arguments))
    
with open("condor_job.txt", "w") as cnd_out:
    cnd_out.write(condor)

with open("arguments.txt", "w") as args:
    args.write("\n".join(arguments))

with open("run_script.sh", "w") as rs:
    rs.write(script)

#os.system("condor_submit condor_job.txt")




