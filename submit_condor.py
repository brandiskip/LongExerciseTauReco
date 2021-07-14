# ! /bin/env python

import os
import subprocess
import datetime
from argparse import ArgumentParser
import pdb
import math
#from samples import samples
from glob import glob
from pdb import set_trace

parser = ArgumentParser()
parser.add_argument("analyzer", help = "which analyser to run", default = 'read_taus.py' )
parser.add_argument("sample", help = "sample", default = 'SMS' )
# parser.add_argument("sample", help = "sample", nargs = '+', choices = 'SMS', default = 'SMS' )
parser.add_argument("-n"  , "--njobs"  , dest = "njobs"  ,  type = int, help = "tot number of input files to be read. All = -1" , default = -1                            )
parser.add_argument("-d"  , "--outdir" , dest = "outdir" ,  help = "output dir"                                     , default = "ntuples" )
parser.add_argument("-a"  , "--addtag" , dest = "addtag" ,  help = "add tag to output dir"                          , default = "ntuples" )
parser.add_argument("-t"  , "--test"   , dest = "test"   ,  help = "do not submit to queue"                        , default = False, action='store_true')
# parser.add_argument("--print"          , dest = "printN" ,  help = "print infos"                                   , default = False, action='store_true'      )
# parser.add_argument("-S"  , "--start"  , dest = "start"  ,  help = "choose starting file"                           , default =  0                            )
# parser.add_argument("-c"  , "--chan"   , dest = "channel",  help = "LMNR, Psi"                                  , default = "LMNR"                 )
# parser.add_argument("-g"  , "--dogen"  , dest = "dogen"  ,  help = "produce gen ntuples "                       , default=False, action='store_true'      )
# parser.add_argument("-f"  , "--flavour",                    help = "job flavour (sets running time) https://indico.cern.ch/event/731021/contributions/3013463/attachments/1656036/2651022/18-05-24_HTCondor_CMG.pdf", default =  'microcentury', choices = ['espresso', 'microcentury', 'longlunch', 'workday', 'tomorrow', 'testmatch', 'nextweek'])
args = parser.parse_args()


script_loc = os.path.realpath(args.analyzer)
# channel = args.channel
# gen_flag = '--dogen' if (args.dogen) else ''
# isMC    = samples[args.samples[0]]['isMC']
# mc_flag = '--mc' if isMC==False else ''

base_out = '%s/%s' %(args.outdir, args.addtag)
# for sample_name in args.samples:
#     sample = samples[sample_name]
    ## local out folder for logs
#     base_out = f'{args.outdir}/{channel}/{sample_name}_{args.addtag}'
os.makedirs('%s/scripts'%base_out)
os.makedirs('%s/outCondor'%base_out)

##  output folder for root files
full_eos_out = '{eos_out_folder}/{base_out}/'.format(eos_out_folder = os.getcwd(), base_out = base_out)


## retrieve filelist from das and calculate n jobs to be submitted
# import subprocess
# process = subprocess.Popen(['dasgoclient', '--query=file dataset=%s'%ds_name, '--format=list'], \
#                            stdout=subprocess.PIPE,\
#                            universal_newlines=True)
# files = process.stdout.readlines()
# njobs = len(files)
#

sample_dict = {
  'SMS'             : '/SMS-TStauStau_ctau-0p01to10_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-GridpackScan_102X_upgrade2018_realistic_v15-v1/MINIAODSIM',
  'SMS_mstau90'     : '/SMS-TStauStau_ctau-0p01to10_mStau-90_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-GridpackScan_102X_upgrade2018_realistic_v15-v1/MINIAODSIM',
  'SMS_mstau250'    : '/SMS-TStauStau_ctau-0p01to10_mStau-250to500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-GridpackScan_102X_upgrade2018_realistic_v15-v1/MINIAODSIM',
  'SMS_mlsp50'      : '/SMS-TStauStau_ctau-0p01to10_mLSP-50to100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-GridpackScan_102X_upgrade2018_realistic_v15-v1/MINIAODSIM',
  'HNL_M_10'        : '/HeavyNeutrino_trilepton_M-10_V-0_00108_tau_massiveAndCKM_LO/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM',
  'HNL_M_5'         : '/HeavyNeutrino_trilepton_M-5_V-0_00639_tau_massiveAndCKM_LO/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM',
  'HNL_M_10_ReMINI' : 'mysample',
  'gmsb'            : '/Staus_M_500_0p01mm_13TeV_2018MC/jalimena-Reco-2f1667a4ab974bdf4cb2916f291c3603/USER'

}

ds_name = sample_dict[args.sample]

# njobs = 0
filelistname = '%s_filelist.txt'%args.sample

# instance --instance=cms_dbs_ph_analysis_03

# from subprocess import Popen, PIPE, STDOUT

if not os.path.isfile(filelistname):
    with open(filelistname, 'w') as f:
        if 'gmsb' not in args.sample:
            process = subprocess.Popen(['dasgoclient', '--query=file dataset=%s'%ds_name, '--format=list'], stdout=f)
            process.wait()
        else:
            print 'gmsb -> does not work'
            process = subprocess.Popen(['dasgoclient', '--query="file dataset=%s'%ds_name, '--format=list', 'instance=prod/phys03"'], stdout=f)
            process.wait()

            pdb.set_trace()

# pdb.set_trace()
# files = open(filelistname, 'r').readlines()
# with open(filelistname, 'r') as f:
# njobs = len(open(filelistname, 'r').readlines())
njobs=len(open(filelistname, 'r').readlines())
os.system('cp {fname} {base_out}'.format(fname=filelistname, base_out=base_out))

# files = process.stdout.readlines()
# njobs = len(files)


# import fnmatch
# njobs = int(args.njobs   )
# if args.njobs == -1:
#     flist = glob(samples[args.samples[0]]['path']+'/000*/*.root')
#     njobs = len(flist)

## find missing files from crab
# job_n = []
# for i,f in enumerate(flist):  job_n.append(int(f.split('_')[-1].split('.')[-2]))
# job_n.sort()

# pdb.set_trace()
# missing =  sorted(set(range(job_n[0], job_n[-1])) - set(job_n))
# print 'missing files from crab:', missing
# njobs = njobs + len(missing) + 1 ## adding 1 since condor start from 0 and crab_0.root does not exists
print 'n jobs to be submitted: ', njobs

bname = os.path.realpath('%s/scripts/script_condor.sh'%base_out)

## is nstart!=0 -> instead of process ID maybe use X+processID?
    ## write script_flat.sh script
with open(bname, 'w') as batch:
    batch.write('''#!/bin/tcsh
source  /cvmfs/cms.cern.ch/cmsset_default.csh
pushd /afs/cern.ch/user/b/bskipwor/work/CMSSW_10_6_25/src/RecoTauTag/RecoTau/test/LongExerciseTauReco
eval `scram runtime -csh`
setenv PYTHONPATH `which python`
setenv PYTHONHOME `scram tool info python | grep PYTHON_BASE | sed 's/PYTHON_BASE=//'`
setenv KRB5CCNAME /afs/cern.ch/user/b/bskipwor/tmp/krb5cc_`id -u bskipwor`
popd

echo "python {script_loc} --sample {thesample} --file $1 "
time python {script_loc}  --sample {thesample} --file $1
mv tau*.root {full_eos_out}
'''
.format(script_loc   = script_loc,
        full_eos_out = full_eos_out,
        thesample    = args.sample
        )
)
subprocess.call(['chmod', '+x', bname])


## write the cfg for condor submission condor_multiple_readnano.cfg
with open('%s/condor_sub.cfg'%base_out, 'w') as cfg:
    cfg.write('''Universe = vanilla
Executable = {bname}
use_x509userproxy = True
transfer_input_files = {filelistname}
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
getenv = True
requirements = (OpSysAndVer =?= "CentOS7")
Log    = {base_out}/outCondor/condor_job_$(Process).log
Output = {base_out}/outCondor/condor_job_$(Process).out
Error  = {base_out}/outCondor/condor_job_$(Process).err
Arguments = $(Process) {sample}
+JobFlavour = "longlunch"
Queue {njobs}'''.format( bname = bname,
                    base_out = base_out,
                    outdir = args.outdir,
                    sample = args.sample,
                    filelistname = filelistname,
#                      era = args.era,
                    njobs = njobs )
)
    # submit to the queue
    print('condor_submit {base_out}/condor_sub.cfg'.format(base_out=base_out))
    if not args.test:
        os.system("condor_submit {base_out}/condor_sub.cfg".format(base_out=base_out))


# +JobFlavour = "{flavour}"
