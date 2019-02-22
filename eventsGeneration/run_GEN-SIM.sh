#!/bin/bash
export XRD_NETWORKSTACK=IPv4
export CMSSWVER="CMSSW_9_4_6_patch1"
export CMSSWDIR="/afs/cern.ch/work/j/jniedzie/private/CMSSW_9_4_6_patch1/src/"
export SCRAM_ARCH="slc6_amd64_gcc630"

source /afs/cern.ch/cms/cmsset_default.sh
cd ${CMSSWDIR}
eval `scramv1 runtime -sh`
edmPluginRefresh -p ../lib/$SCRAM_ARCH

## Execute job and retrieve the outputs
echo "Job running on `hostname` at `date`"

cmsRun chargino300GeV_ctau10cm_GEN-SIM_Run2_2017_cfg.py


