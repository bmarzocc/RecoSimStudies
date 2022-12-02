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
parser.add_argument("-q", "--queue", type=str, help="Condor queue", default="longlunch", required=True)
parser.add_argument("-e", "--eos", type=str, default="user", help="EOS instance user/cms", required=False)
parser.add_argument("-p", "--proxy", type=str, help="Proxy key", required=False)
parser.add_argument("-s", "--site", type=str, help="DAS site", required=False)
parser.add_argument("--redo", action="store_true", default=False, help="Redo all files")
args = parser.parse_args()

if (args.inputdir=="" and args.das=="") or (args.inputdir!="" and args.das!=""):
  print "ERROR: Give either inputdir, either DAS option"
  exit()  

if not os.path.isdir('error'): os.mkdir('error') 
if not os.path.isdir('output'): os.mkdir('output') 
if not os.path.isdir('log'): os.mkdir('log') 

command = "tar -czvf cmssw.tar DIR"
command = command.replace("DIR", args.cmssw)
os.system("rm -rf cmssw.tar")
os.system(command)

# Prepare condor jobs
condor = '''executable              = run_script.sh
output                  = output/strips.$(ClusterId).$(ProcId).out
error                   = error/strips.$(ClusterId).$(ProcId).err
log                     = log/strips.$(ClusterId).log
transfer_input_files    = run_script.sh
transfer_output_files   = ""
on_exit_remove          = (ExitBySignal == False) && (ExitCode == 0)
periodic_release        = (NumJobStarts < 3) && ((CurrentTime - EnteredCurrentStatus) > (60*60))

+JobFlavour             = "{queue}"
+AccountingGroup        = "group_u_CMS.CAF.ALCA" 
queue arguments from arguments.txt

'''

condor = condor.replace("{queue}", args.queue)
user = os.environ["USER"]

script = '''#!/bin/sh -e
export X509_USER_PROXY=DIR/PROXYKEY
voms-proxy-info
cp -r DIR/cmssw.tar ./
tar -xzvf cmssw.tar
cd {cmssw_file}/src
echo -e "evaluate"
eval `scramv1 ru -sh`
export HOME='/afs/cern.ch/user/{user1}/{user}'

JOBID=$1;  
INPUTFILE=$2;
OUTPUTFILE=$3;

cd RecoSimStudies/Dumpers/crab/

echo -e "cmsRun..";

echo -e ">>> STEP3";
cmsRun MiniAOD_fromRaw_Run3_rereco_DeepSC_algoA_Data2022_cfg.py inputFile=${INPUTFILE} outputFile=../../../output.root

echo -e "Running JPsi-dumper.."

cd ../../..;
cd Bmmm/Analysis/test/;
python3 inspector_kee_analysis.py --inputFiles ../../../output.root --filename output_${JOBID}

echo -e "Copying result to: $OUTPUTFILE";
xrdcp -f --nopbar  output_${JOBID}.root root://eoscms.cern.ch/${OUTPUTFILE};

echo -e "DONE";
'''

script = script.replace("{eosinstance}", args.eos)
script = script.replace("{user1}", user[:1])
script = script.replace("{user}", user)
cmssw_file = str(args.cmssw)[1:]
script = script.replace("{cmssw_loc}", args.cmssw)
script = script.replace("{cmssw_file}", cmssw_file)
script = script.replace("DIR", os.getcwd())
script = script.replace("PROXYKEY", args.proxy)

arguments= []
if not os.path.exists(args.outputdir):
  os.makedirs(args.outputdir)

outputfiles = [args.outputdir +"/"+f for f in os.listdir(args.outputdir)]

inputfiles = ""
if args.inputdir!="":
  inputfiles = [ f for f in os.listdir(args.inputdir)]
if args.das!="":
  print "Site:",args.site 
  if args.site=="":
    query = "dasgoclient --query='file dataset=/"+args.das+"'"
  else: 
    query = "dasgoclient --query='file dataset=/"+args.das+" site="+args.site+"'"
    print "QUERY:",query 
  #query = "dasgoclient --query='file dataset=/"+args.das+" site=T0_CH_CERN_Disk'"
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




