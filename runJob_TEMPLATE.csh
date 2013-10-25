#!/bin/csh -f

setenv X509_USER_PROXY LOCPROXY
setenv SRM_PATH /opt/d-cache/srm

echo "Changing directories to ${_CONDOR_SCRATCH_DIR}"
cd ${_CONDOR_SCRATCH_DIR}
echo "Starting job: cmsRun PYCFG"
cmsRun PYCFG

CP_COMMAND
if (-e OUTPUTDIR/OUTFILE) then
    rm -f OUTFILE
else if (-e OUTFILE) then
    set size = `ls -lh OUTFILE | awk '{print $5}'`
    echo "Copy to storage failed, but output file seems to exist."   
    echo "Copying OUTFILE back to JOBDIR."
    echo "Size of OUTFILE is ${size} bytes."
else 
    echo "File OUTFILE not produced."
endif

# EOS
# cp file.root /eos/uscms/store/user/jhirsch/file.root

# dCache
# /opt/d-cache/srm/bin/srmcp -debug=true -2 "file:////${_CONDOR_SCRATCH_DIR}/file.root" "srm://cmssrm.fnal.gov:8443/11/store/user/lpcsusystealth/file.root"



