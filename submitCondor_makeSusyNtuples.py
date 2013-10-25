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
parser.add_option ('--f', action="store_true",
                   dest="force", default=False,
                   help="Force overwrite.")
parser.add_option ('--m', dest='mod', type='string',
                   default = 'mod',
                   help="name mode")
parser.add_option ('--intype', dest='intype', type='string',
                   default = '-1',
                   help="eos, file, or dcache")
parser.add_option ('--indir', dest='indir', type='string',
                   default = 'None',
                   help="input directory")
parser.add_option ('--infilelist', dest='infilelist', type='string',
                   default = 'None',
                   help="list of input infiles")
parser.add_option ('--outfile', dest='outfile', type='string',
                   default = 'susyEvents',
                   help="output file name")
parser.add_option ('--cfg', dest='cfgBase', type='string',
                   default = 'condor.template',
                   help="crab config base")
parser.add_option ('--pycfg', dest='pycfgBase', type='string',
                   default = 'runOverAOD_btag_TEMPLATE.py',
                   help="CMSSW config base")
parser.add_option ('--script', dest='scriptBase', type='string',
                   default = 'runJob_TEMPLATE.csh',
                   help="condor executable")
parser.add_option ('--store', dest='store', type='string',
                   default = '-1',
                   help="storage type: eos or dcache")
parser.add_option ('--eos_area', dest='eos_area', type='string',
                   default = 'jhirsch',
                   help="the script prepends /eos/uscms/store/user/")
parser.add_option ('--dcache_area', dest='dcache_area', type='string',
                   default = 'lpcsusystealth',
                   help="the script prepends /pnfs/cms/WAX/11/store/user/")
parser.add_option ('--area', dest='area', type='string',
                   default = '-1',
                   help="extension of path within eos or dcache (not necessary)")
parser.add_option ('--jobdir_area', dest='jobdir_area', type='string',
                   default = '-1',
                   help="prepended to jobdir (for collecting all jobdirs in a single dir, for instance)")
parser.add_option ('--o', dest='baseOutdir', type='string',
                   default = '-1',
                   help="base output directory (must be unique - dedicated to this run")
parser.add_option ('--njobs', dest='njobs', type='int',
                   default = '-1',
                   help="number of jobs")
parser.add_option ('--nevts', dest='nevts', type='int',
                   default = '-1',
                   help="number of events")
parser.add_option ('--test', action="store_true",
                   dest="test", default=False,
                   help="Just testing")

options, args = parser.parse_args()


mod = options.mod
test = options.test
baseOutdir = options.baseOutdir
cfgBase = options.cfgBase
pycfgBase = options.pycfgBase
scriptBase = options.scriptBase
store = options.store
nevts = options.nevts
njobs = options.njobs
force = options.force
indir = options.indir
infilelist = options.infilelist
intype = options.intype
outfile = options.outfile
area = options.area
jobdir_area = options.jobdir_area
eos_area = options.eos_area
dcache_area = options.dcache_area

createdDirs = []

# WORKING HERE
if store == "dcache":
    cmd = commands.getoutput("voms-proxy-info")
    # Check if voms-proxy-info exists
    if len(cmd.split("\n"))==1 and cmd.split(":")[2] == " command not found": # note leading space " command ..."
        print 'Cannot find "voms-proxy-info" command.'
        print 'Please set up grid access -- maybe with "source /uscmst1/prod/grid/gLite_SL5.csh"'
        print 'Then get proxy with "voms-proxy-init"'
        print 'Exiting.'
        sys.exit()
    elif len(cmd.split("\n"))==3 and cmd.split("\n")[1] == "Couldn't find a valid proxy.":
        print 'Please get proxy with "voms-proxy-init"'
        print 'Exiting.'
        sys.exit()

    lines = cmd.split("\n")
    for line in lines:
        sline = [x.strip() for x in line.split(":")]
        if sline[0] == "path": proxy = sline[1]

    proxyname = proxy.split("/")[2]
    home = commands.getoutput("echo $HOME")
    if not os.path.isfile(home+"/"+proxyname):
        print "Please copy %s to %s" % (proxy, home+"/"+proxyname)
        print 'Also consider renewing proxy with "voms-proxy-init"'
        print "Exiting."
        sys.exit()

    proxy = home+"/"+proxyname
    
# Make sure things are set:
if baseOutdir == "-1":
    print "Please specify base output name for naming jobdir and output_dir.  Exiting."
    sys.exit()
if intype == "-1":
    print "Please specify input type = eos, dcache, or file.  Exiting."
    sys.exit()    
if store == "-1":
    print "Please specify output storage = eos or dcache."
    print "Please edit this script to specify area (for eos and dcache).  Exiting."
    sys.exit()

# Make list of input files and check for existence:
if indir != "None" and infilelist != "None":
    print "Please specify only --indir or --infilelist, not both.  Exiting."
    sys.exit()
elif indir == "None" and infilelist == "None":
    print "Please specify --indir or --infilelist.  Exiting."
    sys.exit()
elif indir != "None":
    files = glob.glob(indir+"/*.root")
elif infilelist != "None":
    files = open(infilelist, "r").readlines()
    files = [file.strip() for file in files]


if len(files) == 0:
    print "No input files found."
    print "Check inputfile list, or indir should be something like:"
    print "/pnfs/cms/WAX/11/store/user/XXXXX"
    print "/eos/uscms/store/user/XXXXX"
    print "/uscms/home/XXXXX"
    print "Exiting."
    sys.exit()

indir_first = files[0].split("/")[1]

# Make sure input type matches location
if indir_first != "pnfs" and indir_first != "eos" and indir_first != "uscms":
    print 'indir must start with \"/pnfs\" or \"/eos\" or \"/uscms\".  Exiting.'
    sys.exit()
elif indir_first == "pnfs" and intype != "dcache":
    print "For "+indir+", intype should be dcache, not "+intype
elif indir_first == "eos" and intype != "eos":
    print "For "+indir+", intype should be eos, not "+intype
elif indir_first == "uscms" and intype != "file":
    print "For "+indir+", intype should be file, not "+intype
    

# Check and make output directories
# 1) First check that base eos/dcache area exists.
# 2) Then check that final output directory (including full path)
# does NOT exist, and make output directory.

if store == "eos":
    output_dir = "/eos/uscms/store/user/"+eos_area
elif store == "dcache":
    output_dir = "/pnfs/cms/WAX/11/store/user/"+dcache_area

if not os.path.isdir(output_dir):
    print output_dir+" does not exist.  Exiting."
    sys.exit()

output_dir += "/"+area
output_dir += "/"+baseOutdir
if os.path.isdir(output_dir):
    print "Output directory %s already exists.  Exiting." % output_dir
    sys.exit()

os.system("mkdir -p "+output_dir)
createdDirs.append(output_dir)


# Make job directories for submitting to condor
if area != "-1":
    basejobdir = "cond_"+area+"_"+baseOutdir
else:
    basejobdir = "cond_"+baseOutdir

if jobdir_area != "-1":
    basejobdir = jobdir_area+"/"+basejobdir
    
if not os.path.isdir(basejobdir) :
    os.system("mkdir -p "+basejobdir)
    createdDirs.append(basejobdir)

else:
    print basejobdir+" already exists.  Exiting."
    sys.exit()


# Loop over input files starting one job per file
for infile in files:
    print " "
    print "======================================================"
    print "======================================================"
    print "Processing "+infile
    # Get jobid
    sinfile = infile.split("/")
    lastind = len(sinfile)-1
    jobid   = sinfile[lastind].split(".root")[0]

    # Concoct name for output jobdir
    cwd = os.getcwd()
    jobdir = basejobdir+"/"+jobid

    # Create output jobdir as necessary accounting for force and mod options
    if not os.path.isdir(jobdir) :
        os.system("mkdir "+jobdir)
        createdDirs.append(jobdir)

        print "Making directory %s." % jobdir
    elif force:
        timelist = time.localtime()
        backup_code = "_backup"
        for ind in range(1,6):
            backup_code += "_"+str(timelist[ind])
        os.system("mv "+jobdir+" "+jobdir+backup_code)
        os.system("mkdir "+jobdir)
        createdDirs.append(jobdir)
    #elif mod != "mod":
    # No longer valid
    #    print "Directory %s already exists." % jobdir
    #    jobdir = cwd+"/"+baseOutdir+"/"+jobid+"_"+mod
    #    if not os.path.isdir(jobdir) :
    #        os.system("mkdir "+jobdir)
    #        print "Making directory %s." % jobdir
    #    else:
    #        print "Directory %s already exists.  Exiting." % jobdir
    #        sys.exit()
    else:
        print "Directory %s already exists. Exiting." % jobdir
        sys.exit()

    # Make names for necessary files:
    outfile = jobid+"_susyNtuple.root"
    pycfg   = "runOverAOD_btag_"+jobid+".py"
    config  = "condor."+jobid
    script  = "runJob_"+jobid+".csh"

    print "Creating job %s." % jobid


    commandList = []
    # Make cmsRun python config in JOBDIR/
    # INFILE (with full dcap/eos path)
    # OUTFILE (with no path)
    commandList.append("cp "+pycfgBase+" "+jobdir+"/"+pycfg)  # cp runOverAOD_btag_TEMPLATE.py JOBDIR/runOverAOD_btag_JOBID.py
    if intype == "dcache":
        commandList.append("replace INFILE dcap://"+infile+" -- "+jobdir+"/"+pycfg)
    else:
        commandList.append("replace INFILE file:"+infile+" -- "+jobdir+"/"+pycfg)
        
    commandList.append("replace OUTFILE "+outfile+" -- "+jobdir+"/"+pycfg)

    # Make condor config
    # JOBDIR/PYCFG
    # JOBDIR/SCRIPT
    # PYCFG
    commandList.append("cp "+cfgBase+" "+jobdir+"/"+config)  # cp condor.template JOBDIR/condor.JOBID
    commandList.append("replace JOBDIR "+jobdir+" -- "+jobdir+"/"+config)
    commandList.append("replace PYCFG "+pycfg+" -- "+jobdir+"/"+config)
    commandList.append("replace SCRIPT "+script+" -- "+jobdir+"/"+config)
    
    # Make submit script
    # PYCFG
    # CP_COMMAND

    # EOS
    # cp file.root /eos/uscms/store/user/jhirsch/file.root
    
    # dCache
    # /opt/d-cache/srm/bin/srmcp -debug=true -2 "file:////${_CONDOR_SCRATCH_DIR}/file.root" "srm://cmssrm.fnal.gov:8443/11/store/user/lpcsusystealth/file.root"

    if store == "eos":
        cp_command = "cp "+outfile+" "+"/eos/uscms/store/user/"+eos_area+"/"+area+"/"+baseOutdir+"/"+outfile
    elif store == "dcache":
        output_area = "srm://cmssrm.fnal.gov:8443/11/store/user/"+dcache_area+"/"+area+"/"+baseOutdir
        cp_command = '/opt/d-cache/srm/bin/srmcp -debug=true -2 "file:////${_CONDOR_SCRATCH_DIR}/'+outfile+'" "'+output_area+"/"+outfile+'"'

    
    commandList.append("cp "+scriptBase+" "+jobdir+"/"+script)  # cp condor.template JOBDIR/condor.JOBID
    commandList.append("replace PYCFG "+pycfg+" -- "+jobdir+"/"+script)
    commandList.append("replace CP_COMMAND '%s' -- %s/%s" % (cp_command, jobdir, script))
    if store == "dcache":
        commandList.append('replace LOCPROXY "'+proxy+'" -- '+jobdir+'/'+script) 
    commandList.append('replace OUTFILE "'+outfile+'" -- '+jobdir+'/'+script) 
    commandList.append('replace OUTPUTDIR "'+output_dir+'" -- '+jobdir+'/'+script) 
    commandList.append('replace JOBDIR "'+jobdir+'" -- '+jobdir+'/'+script) 
    
    for command in commandList:
        os.system(command)
    if not test:
        os.chdir(jobdir)
        os.system("condor_submit "+config)
        os.chdir(cwd)
        
print "Created directories:"
for dir in createdDirs:
    print dir
    
        


