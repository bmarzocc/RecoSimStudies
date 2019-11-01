# https://docs.google.com/spreadsheets/d/1rJ4FppHxdD6TQpCHnnnX0QD0Y0d4aTEN6ibgTVn7Vyg/edit#gid=0

import sys
import os
import subprocess


def getOptions():
  from argparse import ArgumentParser

  parser = ArgumentParser(description='Production helper for Clustering studies', add_help=True)

  parser.add_argument('-v','--ver', type=str, dest='ver', help='version of production, e.g. V00_v00', default='V00_v00')
  #parser.add_argument('-r','--rel', type=str, dest='rel', help='cmssw release', default='10_6_1_patch1')
  #parser.add_argument(     '--gt',      type=str, dest='gt', help='global tag', default='')

  parser.add_argument('-n','--nevts', type=int, dest='nevts', help='total number of events to be generated', default=10)
  parser.add_argument('-c','--ch', type=str, dest='ch', help='channel, e.g. photon', default='photon', choices=['photon'])
  parser.add_argument('--etmax', type=int, dest='etmax', help='max Et (GeV)', default=100)
  parser.add_argument('--etmin', type=float, dest='etmin', help='min Et (GeV)', default=1)
  parser.add_argument('--doflatenergy', dest='doflatenergy', help='generate flat in energy, otherwise in pt', action='store_true', default=False)
  parser.add_argument('--npart', type=int, dest='npart', help='number of particles to generate per event for closeEcal configuration, specify only if you want to override the default', default=None)
  parser.add_argument('-g','--geo',type=str, dest='geo', help='detector configuration: wTk, noTk, closeEcal', default='closeEcal', choices=['wTk', 'noTk', 'closeEcal'])
  parser.add_argument('-d','--det', type=str, dest='det', help='sub-detector: EB, EE or all', default='EB', choices=['EB', 'EE', 'all'])

  parser.add_argument('--pu', type=str, dest='pu', help='PU configuration', default='noPU', choices=['noPU', 'wPU'])

  parser.add_argument('--pfrhmult', type=float, dest='pfrhmult', help='how many sigma of the noise to use for PFRH thresholds', default=1.)
  parser.add_argument('--seedmult', type=float, dest='seedmult', help='how many sigma of the noise to use for seeding thresholds', default=3.)
  parser.add_argument('--doref', dest='doref', help='use reference values for pfrh,seed, thresholds', action='store_true', default=False)
  parser.add_argument('--doreco', dest='doreco', help='do only step 3 (reconstruction) starting from an existing production', action='store_true', default=False)
  parser.add_argument('--custominput', type=str, dest='custominput', help='full path of the input file that you want to use for the reconstruction (also SE is supported)', default=None)
  parser.add_argument('--custominputdir', type=str, dest='custominputdir', help='full path of the input directory that you want to use for the reconstruction (also SE is supported)', default=None)
  parser.add_argument('--doold', dest='doold', help='use old_ version of the scripts', action='store_true', default=False)

  #parser.add_argument('--dostep3only', dest='dostep3only', help='do only step 3', action='store_true', default=False)
  parser.add_argument('--doshort', dest='doshort', help='set 2 hours as wall clock time instead of 1 day', action='store_true', default=False)
  parser.add_argument('--domedium', dest='domedium', help='set 2 days as wall clock time instead of 1 day', action='store_true', default=False)
  parser.add_argument('--dolong', dest='dolong', help='set 3 days as wall clock time instead of 1 day', action='store_true', default=False)
  parser.add_argument('--docustomtime', dest='docustomtime', help='set custom time', action='store_true', default=False)
  parser.add_argument('--time', type=str, dest='time', help='allowed time for each job of each step in hours, example 01,02,01 for 1 hour for step1, 2 for step2, 2 for step3', default='2,2,5')
  parser.add_argument('--domultithread', dest='domultithread', help='run multithreaded', action='store_true', default=False)
  parser.add_argument('--domultijob', dest='domultijob', help='run several separate jobs', action='store_true', default=False)
  parser.add_argument('--njobs', type=int, dest='njobs', help='number of parallel jobs to submit', default=10)
  parser.add_argument('--dosavehome', dest='dosavehome', help='save in home, otherwise save to SE', action='store_true', default=False)
  parser.add_argument('--doskipdumper', dest='doskipdumper', help='do not run the dumper at the end', action='store_true', default=False)
  
  return parser.parse_args()


if __name__ == "__main__":

  opt = getOptions()

  ##############################
  # job configurations
  #############################
  user = os.environ["USER"]
  evar = 'E' if opt.doflatenergy else 'Et'
  etRange='{}{}to{}GeV'.format(evar,opt.etmin,opt.etmax)
  prodLabel='{c}_{e}_{g}_{d}_{pu}_pfrh{pf}_seed{s}_{v}_n{n}'.format(c=opt.ch,e=etRange,g=opt.geo,d=opt.det,pu=opt.pu,pf=opt.pfrhmult,s=opt.seedmult,v=opt.ver,n=opt.nevts)
  if opt.doref: # change to ref label
    prodLabel='{c}_{e}_{g}_{d}_{pu}_pfrh{pf}_seed{s}_{v}_n{n}'.format(c=opt.ch,e=etRange,g=opt.geo,d=opt.det,pu=opt.pu,pf='Ref',s='Ref',v=opt.ver,n=opt.nevts)
  doref = 1 if opt.doref else 0
  doflatenergy = 1 if opt.doflatenergy else 0
  nthr = 8 if opt.domultithread else 1
  njobs = opt.njobs if opt.domultijob else 1
  if opt.domultijob and opt.njobs <= 1: raise RuntimeError('when running multiple jobs, the number of parallel jobs should be larger than 1')
  if opt.domultijob and opt.nevts % opt.njobs != 0: raise RuntimeError('cannot split events in njobs evenly, please change njobs / nevts')
  nevtsjob = opt.nevts if not opt.domultijob else opt.nevts/opt.njobs
  nevtspremixfile = 600 # current number of events in each premixed file
  npremixfiles = nevtsjob / nevtspremixfile + 1
  if opt.doreco and opt.custominput == None and opt.custominputdir == None: raise RuntimeError('you must supply the custom input, when running with doreco activated')
  if opt.doreco and not opt.domultijob and not os.path.isfile(opt.custominput): raise RuntimeError('custominput {} not found').format(opt.custominput)
  if opt.doreco and opt.domultijob and not os.path.isdir(opt.custominputdir): raise RuntimeError('custominputdir {} not found').format(opt.custominputdir)
  if opt.doreco and opt.domultijob: print 'A gentle reminder that if you are running with doreco and domultijob activated, you should use the same job splitting that was used for the original production'
  #if opt.doreco and opt.domultijob: raise RuntimeError('combination not supported, currently cannot re-reco from job run over multiple files')
  if opt.docustomtime: 
    time1=opt.time.split(',')[0]
    time2=opt.time.split(',')[1]
    time3=opt.time.split(',')[2]
    times = [time1,time2,time3]
    sbatch_times = map(lambda x: '--time=0-{}:00'.format(x), times)
  else:
    if opt.domedium:
      time = '--time=1-23:59'
    elif opt.dolong:
      time = '--time=2-23:59'
    elif opt.doshort:
      time = '--time=0-02:00'
    else:
      time = '--time=1-00:00'
    sbatch_times = [time, time, time]

  ##############################
  # create production directory and logs directory within
  #############################
  prodDir = './{}'.format(prodLabel)
  command = 'mkdir -p {}'.format(prodDir)
  if not os.path.isdir(prodDir):
    os.system(command)
  #else: raise RuntimeError('directory {} already present, not going to overwrite'.format(prodDir))

  command = 'mkdir {}/logs'.format(prodDir)
  os.system(command)

  ############################
  # copy the relevant cmsDriver to prod directory
  ############################
  ## find the names first
  step1_driverName = 'step1_{c}_{g}.py'.format(c=opt.ch,g=opt.geo)
  step2_driverName = 'step2_{pu}.py'.format(pu=opt.pu)
  step3_driverName = 'step3_{pu}.py'.format(pu=opt.pu)
  drivers = [step1_driverName, step2_driverName, step3_driverName]
  if opt.doold: drivers=map(lambda x : 'old_' + x, drivers)
  target_drivers = ['step1.py', 'step2.py', 'step3.py']
  infiles  = ['', 'step1_nj{nj}.root', 'step2_nj{nj}.root']
  infiles_loc = ['', 'step1.root', 'step2.root']
  outfiles = ['step1_nj{nj}.root', 'step2_nj{nj}.root', 'step3_nj{nj}.root']
  outfiles_loc = ['step1.root', 'step2.root', 'step3.root']

  ## copy them to dir
  for i,idriver in enumerate(drivers):
    if opt.doreco and i!=2: continue # skip everything that is not related to step3
    if not os.path.isfile('cmsDrivers/{idr}'.format(idr=idriver)):
      raise RuntimeError('cmsDriver {idr} not found, please check naming'.format(idr=idriver))
    command = 'cp cmsDrivers/{idr} {d}/{td}'.format(idr=idriver,d=prodDir,td=target_drivers[i])
    os.system(command) 

  ############################
  # write the cmsRun commands for all steps
  ############################
  ## step1
  if opt.geo == 'closeEcal':
    if opt.det == 'EB':
      rmin = 123.8
      rmax = 123.8
      zmin = -304.5
      zmax = 304.5
      npart = 10
    else:
      rmin = 31.6
      rmax = 171.1
      zmin = 317.0
      zmax = 317.0
      npart = 10
  
    if opt.npart!=None:
      npart = opt.npart

    step1_cmsRun = 'cmsRun {jo} maxEvents={n} etmin={etmin} etmax={etmax} rmin={r1} rmax={r2} zmin={z1} zmax={z2} np={np} nThr={nt} doFlatEnergy={dfe}'.format(
                    jo=target_drivers[0], n=nevtsjob, etmin=opt.etmin, etmax=opt.etmax, r1=rmin, r2=rmax, z1=zmin, z2=zmax, np=npart, nt=nthr, dfe=doflatenergy)
    step1_cmsRun_add = 'seedOffset={nj}' # format at a later stage
  else:
    raise RuntimeError('this option is not currently supported')
  ## other steps  
  step2_cmsRun = 'cmsRun {jo} nThr={nt} nPremixFiles={npf}'.format(jo=target_drivers[1], nt=nthr, npf=npremixfiles)
  step2_cmsRun_add = 'randomizePremix=True' if opt.domultijob else ''
  step3_cmsRun = 'cmsRun {jo} pfrhMult={pfrhm} seedMult={sm} nThr={nt} doRef={dr}'.format(jo=target_drivers[2], pfrhm=opt.pfrhmult, sm=opt.seedmult, nt=nthr, dr=doref)
  cmsRuns = [step1_cmsRun, step2_cmsRun, step3_cmsRun]
  cmsRuns_add = [step1_cmsRun_add, step2_cmsRun_add, '']
  ############################
  # write the launching scripts
  ############################
  for i,idriver in enumerate(drivers):

    if opt.doreco and i!=2: continue # skip everything that is not related to step3

    for nj in range(0,njobs):
  
      ### configurations for the template script
      outputDir = '"/pnfs/psi.ch/cms/trivcat/store/user/"$USER"/EcalProd/"' 
      if opt.dosavehome: outputDir = '`pwd`/../' 

      mkdiroutput_command = 'xrdfs t3dcachedb03.psi.ch mkdir $SERESULTDIR'
      if opt.dosavehome: mkdiroutput_command = 'mkdir -p $SERESULTDIR'

      cpinput_command = ''
      if infiles[i]!='':
        cpinput_command = 'xrdcp $SEPREFIX/$SERESULTDIR/{infile} $WORKDIR/{infile_loc}'.format(infile=infiles[i].format(nj=nj),infile_loc=infiles_loc[i])
        if opt.dosavehome: cpinput_command = 'cp $SERESULTDIR/{infile} $WORKDIR/{infile_loc}'.format(infile=infiles[i].format(nj=nj),infile_loc=infiles_loc[i])
        if opt.doreco: 
          custominput = opt.custominput if not opt.domultijob else opt.custominputdir + '/step2_nj{nj}.root'.format(nj=nj) # could have used also infiles[i].format(nj=nj)
          if opt.dosavehome: cpinput_command = 'cp {ci} $WORKDIR/{infile_loc}'.format(ci=custominput,infile_loc=infiles_loc[i])
          else:              cpinput_command = 'xrdcp $SEPREFIX/{ci} $WORKDIR/{infile_loc}'.format(ci=custominput,infile_loc=infiles_loc[i])

      cpoutput_command = 'xrdcp -f {outfile_loc} $SEPREFIX/$SERESULTDIR/{outfile}'.format(outfile_loc=outfiles_loc[i],outfile=outfiles[i].format(nj=nj))
      if opt.dosavehome: cpoutput_command = 'cp {outfile_loc} $SERESULTDIR/{outfile}'.format(outfile_loc=outfiles_loc[i],outfile=outfiles[i].format(nj=nj))
       
      cpaux_command = ''
      if 'step3' in idriver:
        cpaux_command = 'cp -r $CMSSW_BASE/src/RecoSimStudies/Dumpers/data $WORKDIR'

      ### define a template script  
      template = [
      '#!/bin/bash',
      '',

      #### variables
      'DIRNAME="{ind}"',
      'STARTDIR=`pwd`',
      'TOPWORKDIR="/scratch/"$USER/',
      'JOBDIR="gen_"$SLURM_JOB_ID',
      'WORKDIR=$TOPWORKDIR/$JOBDIR',
      'SEPREFIX="root://t3dcachedb.psi.ch:1094/"',
      'SERESULTDIR={od}/$DIRNAME',
      'JOBOPFILENAME="{jo}"',
      '',

      #### environment
      'source $VO_CMS_SW_DIR/cmsset_default.sh',
      'shopt -s expand_aliases',
      'echo ""',
      'echo "Going to set up cms environment"',
      'cd $STARTDIR',
      'cmsenv',
      'echo ""',
      '',

      #### workdir
      'echo "Going to create work dir"',
      'mkdir -p $WORKDIR',
      'echo "workdir: "',
      'echo $WORKDIR',
      'echo ""',
      '',

      #### outputdir
      'echo "Going to create the output dir"',
      'echo "May give an error if the directory already exists, which can be safely ignored"',
      '{mkdir}',
      'echo ""',
      '',

      #### copy driver and other aux files 
      'echo "Going to copy cms driver"',
      'cp $JOBOPFILENAME $WORKDIR/$JOBOPFILENAME',
      '{cpaux}',
      'echo ""',
      '',

      #### copy input file 
      'echo "Going to copy input file if needed"',
      '{cpin}',
      'echo ""',
      '',

      #### run
      'cd $WORKDIR',
      'echo ""',
      '',
      'echo "Going to run"',
      'DATE_START=`date +%s`',
      '{cmsRun}',
      'DATE_END=`date +%s`',
      'echo ""',
      '',
      'echo "Finished running"',
      'echo "Content of current directory"',
      'ls -al',

      #### copy back output
      'echo "Going to copy the output to the output directory"',
      '{cpout}',
      '',
      'echo ""',

      #### clean and go 
      'echo "Cleaning up $WORKDIR"',
      'rm -rf $WORKDIR',
      'RUNTIME=$((DATE_END-DATE_START))',
      'echo "Wallclock running time: $RUNTIME s"',
      'cd $STARTDIR',

      ] 
      template = '\n'.join(template)
      template = template.format(ind=prodLabel,od=outputDir,jo=target_drivers[i],mkdir=mkdiroutput_command,
                                 #cpin=cpinput_command,cpout=cpoutput_command,cmsRun=cmsRuns[i]+' '+cmsRuns_add[i].format(nj=nj+nthr+1),cpaux=cpaux_command) 
                                 cpin=cpinput_command,cpout=cpoutput_command,cmsRun=cmsRuns[i]+' '+cmsRuns_add[i].format(nj=nj+1),cpaux=cpaux_command) 

      launcherFile = '{}/launch_{}.sh'.format(prodDir,outfiles[i].format(nj=nj).split('.root')[0])
      with open(launcherFile, 'w') as f:
        f.write(template)

  #############################
  # write the template script to run the dumper
  # one template for the full task
  ############################
  if not opt.doskipdumper:
 
    # the template should run a python script to get the sample list
    # and then run the actual cmsRun command for the dumper
    # all done locally
    template_dumper = [
      '#!/bin/bash',
      '',
      'STARTDIR=$PWD/../',
      #### environment
      'source $VO_CMS_SW_DIR/cmsset_default.sh',
      'shopt -s expand_aliases',
      'echo ""',
      'echo "Going to set up cms environment"',
      'cd $STARTDIR',
      'cmsenv',
      'echo ""',
      '',
      ### running part
      'echo "Going to run postProdHelper.py"',
      'DATE_START=`date +%s`',
      'python postProdHelper.py --pl {pl} --user {u}'.format(pl=prodLabel, u=user),
      'echo ""',
      '',
      'echo "Going to run the Dumper"',
      'cmsRun ../../python/Cfg_RecoSimDumper_cfg.py outputFile=../../test/outputfiles/{pl}.root inputFiles_load=../../data/samples/{pl}.txt'.format(pl=prodLabel),
      'DATE_END=`date +%s`',
      'echo ""',
      '',
      ### print job time
      'echo "Finished running"',
      'RUNTIME=$((DATE_END-DATE_START))',
      'echo "Wallclock running time: $RUNTIME s"',
    ]
    
    template_dumper = '\n'.join(template_dumper)

    launcherFile_dumper = '{}/launch_dumper.sh'.format(prodDir)
    with open(launcherFile_dumper, 'w') as f:
      f.write(template_dumper)

  #############################
  # finally write the submitter
  ############################# 
  submitter_template = []

  dumper_dependencies = ''

  for nj in range(0,njobs):

    sbatch_command_step1 = 'jid1_nj{nj}=$(sbatch -p wn --account=cn-test -o logs/step1_nj{nj}.log -e logs/step1_nj{nj}.log --job-name=step1_{pl} {t} --ntasks={nt} launch_step1_nj{nj}.sh)'.format(nj=nj,pl=prodLabel,t=sbatch_times[0],nt=nthr)

    sbatch_command_step2 = 'jid2_nj{nj}=$(sbatch -p wn --account=cn-test -o logs/step2_nj{nj}.log -e logs/step2_nj{nj}.log --job-name=step2_{pl} {t} --ntasks={nt} --dependency=afterany:$jid1_nj{nj} launch_step2_nj{nj}.sh)'.format(nj=nj,pl=prodLabel,t=sbatch_times[1],nt=nthr)

    sbatch_command_step3 = 'jid3_nj{nj}=$(sbatch -p wn --account=cn-test -o logs/step3_nj{nj}.log -e logs/step3_nj{nj}.log --job-name=step3_{pl} {t} --ntasks={nt} --dependency=afterany:$jid2_nj{nj} launch_step3_nj{nj}.sh)'.format(nj=nj,pl=prodLabel,t=sbatch_times[2],nt=nthr)
    if opt.doreco: # strip the dependency away
      sbatch_command_step3 = 'jid3_nj{nj}=$(sbatch -p wn --account=cn-test -o logs/step3_nj{nj}.log -e logs/step3_nj{nj}.log --job-name=step3_{pl} {t} --ntasks={nt}  launch_step3_nj{nj}.sh)'.format(nj=nj,pl=prodLabel,t=sbatch_times[2],nt=nthr)

    dumper_dependencies += ':$jid3_nj{nj}'.format(nj=nj)

    if not opt.doreco:
      submitter_template.append(sbatch_command_step1)
      submitter_template.append('echo "$jid1_nj%i"' % nj)
      submitter_template.append('jid1_nj%i=${jid1_nj%i#"Submitted batch job "}' % (nj,nj))
      submitter_template.append(sbatch_command_step2)
      submitter_template.append('echo "$jid2_nj%i"' % nj)
      submitter_template.append('jid2_nj%i=${jid2_nj%i#"Submitted batch job "}' % (nj,nj))

    submitter_template.append(sbatch_command_step3)
    submitter_template.append('echo "$jid3_nj%i"' % nj)
    submitter_template.append('jid3_nj%i=${jid3_nj%i#"Submitted batch job "}' % (nj,nj))

  # add the dumper part
  if not opt.doskipdumper:
    sbatch_command_dumper = 'jid_d=$(sbatch -p wn -o logs/dumper.log -e logs/dumper.log --job-name=dumper_{pl} {t} --ntasks=1 --dependency=afterany{dd} launch_dumper.sh)'.format(pl=prodLabel,t='--time=0-00:59',dd=dumper_dependencies)
    submitter_template.append(sbatch_command_dumper)
    submitter_template.append('echo "$jid_d"')
    submitter_template.append('jid_d=${jid_d#"Submitted batch job "}')
  
  submitter_template = '\n\n'.join(submitter_template)

  submitterFile = '{}/submit.sh'.format(prodDir)
  with open(submitterFile, 'w') as f:
    f.write(submitter_template)


