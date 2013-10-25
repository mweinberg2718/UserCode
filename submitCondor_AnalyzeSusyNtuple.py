#!/usr/bin/env python

import sys
import optparse
import commands
import os
import glob
import time

#######################
# Get options
#######################

parser = optparse.OptionParser("usage: %prog [options]\
<input directory> \n")
parser.add_option ('--m', dest='mod', type='string',
                   default = 'mod',
                   help="name modifier")

parser.add_option ('--prefix', dest='prefix', type='string',
                   default = 'NONE',
                   help="directory containing input susy ntuples")
parser.add_option ('--indir', dest='indir', type='string',
                   default = './',
                   help="directory containing input susy ntuples")
parser.add_option ('--infile', dest='infile', type='string',
                   default = 'USEDIR',
                   help="directory containing input susy ntuples")
parser.add_option ('--cfg', dest='cfgBase', type='string',
                   default = 'condor_filelist.template',
                   help="crab config base")
parser.add_option ('--o', dest='baseOutdir', type='string',
                   default = './',
                   help="base out directory")
parser.add_option ('--ds', dest='dataset', type='string',
                   default = 'test',
                   help="base out directory")
parser.add_option ('--njobs', dest='njobs', type='int',
                   default = '-1',
                   help="number of jobs, default = 1 job per file")
parser.add_option ('--nevts', dest='nevts', type='int',
                   default = '-1',
                   help="number of events")
parser.add_option ('--test', action="store_true",
                   dest="test", default=False,
                   help="Just testing")
parser.add_option ('--lumi', dest='lumi', type='string',
                   default = '1.0',
                   help="luminosity used to scale events")
parser.add_option ('--xsec', dest='xsec', type='string',
                   default = '-1.0',
                   help="cross section used to scale MC events (= -1.0 for data)")
parser.add_option ('--json', dest='json', type='string',
                   default = '',
                   help="JSON file used to select events")
parser.add_option ('--jec', dest='jec', type='string',
                   default = '0.0',
                   help="JEC corrections (1.0 for 1 sigma up, -1.0 for 1 sigma down)")

options, args = parser.parse_args()

prefix = options.prefix
dataset= options.dataset
mod = options.mod
test = options.test
baseOutdir = options.baseOutdir
cfgBase = options.cfgBase
njobs = options.njobs
indir = options.indir
infile = options.infile
lumi = options.lumi
xsec = options.xsec
json = options.json
jec = options.jec

cwd = os.getcwd()

if not os.path.isdir(baseOutdir) :
    os.system("mkdir -p "+baseOutdir)
    print "Making directory %s." % baseOutdir
#elif force:
#    timelist = time.localtime()
#    backup_code = "_backup"
#    for ind in range(1,6):
#        backup_code += "_"+str(timelist[ind])
#    os.system("mv "+outdir+" "+outdir+backup_code)
#    os.system("mkdir "+outdir)
else:
    print "Output directory %s already exists.  Exiting." % baseOutdir
    sys.exit()


if infile == "USEDIR":
    files0 = glob.glob(indir+"/*.root")
else:
    files0 = glob.glob(infile)
    njobs = 1

if prefix != "NONE":
    files = []
    for file in files0:
        files.append(prefix+file)
else:
    files = files0

nfiles = len(files)
# If -1, make 1 job per file:
if njobs == -1: njobs = nfiles

# Split files into jobs as evenly as possible:

baseFilesPerJob = nfiles/njobs
extraFiles = nfiles%njobs

numberOfFiles = {}
for ijob in range(njobs):
    numberOfFiles[ijob] = baseFilesPerJob

for ijob in range(extraFiles):
    numberOfFiles[ijob] += 1


files_for_jobs = {}
IFILE = 0
for ijob in range(njobs):
    nfiles_ijob = numberOfFiles[ijob]
    files_for_jobs[ijob] = []
    for ifile in range(nfiles_ijob):
        files_for_jobs[ijob].append( files [IFILE] )
        IFILE += 1

if IFILE != nfiles:
    print "Mismatch.", IFILE, nfiles-1,"  Exiting."
    sys.exit()

for ijob in range(njobs):
    jobid = str(ijob)
    jobdir = baseOutdir+"/condor_"+jobid

    if not os.path.isdir(jobdir): os.system("mkdir "+jobdir)
    else:
        print jobdir,"exists.  Exiting.";sys.exit()

    filelistname = "filelist_"+jobid
    filelist = open(jobdir+"/"+filelistname, 'w')
    config = cfgBase+"_"+jobid
    path_config = jobdir+"/"+config
    nfiles = len(files_for_jobs[ijob])

    for ifile in files_for_jobs[ijob]:
        filelist.write(ifile+"\n")

    print "Creating job %s with %s files." % (jobid, nfiles)


    commandList = []
    commandList.append("cp "+cfgBase+" "+path_config)
    commandList.append("replace FILELIST "+filelistname+" -- "+path_config)
    commandList.append("replace JOBDIR "+jobdir+" -- "+path_config)
    commandList.append("replace CWD "+cwd+" -- "+path_config)
    commandList.append("replace JOBID "+jobid+" -- "+path_config)

    if json == '':
        commandList.append('replace JSON  " " -- '+path_config)
    else:
        commandList.append('replace JSON  '+json+' -- '+path_config)

    commandList.append("replace DATASET "+dataset+" -- "+path_config)
    #commandList.append("replace OUTDIR "+outdir+" -- "+path_config)
    commandList.append("replace LUMI "+lumi+" -- "+path_config)
    commandList.append("replace XSEC "+xsec+" -- "+path_config)
    commandList.append("replace JEC "+jec+" -- "+path_config)

    for command in commandList:
        os.system(command)
    if not test:
        os.chdir(jobdir)
        os.system("condor_submit "+config)
        os.chdir(cwd)





