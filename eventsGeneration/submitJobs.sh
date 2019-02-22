#!/bin/bash

number_of_jobs=$1
number_of_events_per_job=$2

echo "usage: ./submitJobs.sh number_of_jobs number_of_events_per_job"

# create condor configs and generation scripts
for iJob in $(seq 1 $number_of_jobs);
do
FILE="scripts/run_batch_${iJob}.sub"
/bin/cat <<EOM >$FILE
executable  = scripts/run_GEN-SIM_${iJob}.sh
arguments   = \$(ClusterID) \$(ProcId)
output      = output/GEN-SIM.\$(ClusterId).\$(ProcId).out
error       = error/GEN-SIM.\$(ClusterId).\$(ProcId).err
log         = log/GEN-SIM.\$(ClusterId).log
transfer_input_files   = scripts/chargino300GeV_ctau10cm_GEN-SIM_${iJob}.py
+JobFlavour     = "nextweek"
queue
EOM

cmsDriver.py Configuration/GenProduction/python/ThirteenTeV/AMSB_chargino/AMSB_chargino300GeV_ctau10cm_NoFilter_13TeV.py \
--fileout file:generatedEvents/chargino300GeV_ctau10cm_GEN-SIM_${iJob}.root \
--step GEN,SIM \
--mc \
--datatier GEN-SIM \
--beamspot Realistic25ns13TeVEarly2017Collision \
--conditions auto:phase1_2017_realistic \
--eventcontent RAWSIM \
--era Run2_2017 \
--python_filename scripts/chargino300GeV_ctau10cm_GEN-SIM_${iJob}.py \
--customise Configuration/DataProcessing/Utils.addMonitoring,SimG4Core/CustomPhysics/Exotica_HSCP_SIM_cfi,SimG4Core/CustomPhysics/GenPlusSimParticles_cfi.customizeProduce,SimG4Core/CustomPhysics/GenPlusSimParticles_cfi.customizeKeep,Configuration/GenProduction/RandomSeed_cfi.customizeRandomSeed \
--no_exec \
-n ${number_of_events_per_job}

FILE2="scripts/run_GEN-SIM_${iJob}.sh"
/bin/cat <<EOM >$FILE2
#!/bin/bash
export XRD_NETWORKSTACK=IPv4
export CMSSWVER="CMSSW_9_4_6_patch1"
export SCRAM_ARCH="slc6_amd64_gcc630"

source /afs/cern.ch/cms/cmsset_default.sh
eval \`scramv1 runtime -sh\`
edmPluginRefresh -p ../lib/\$SCRAM_ARCH

## Execute job and retrieve the outputs
echo "Job running on `hostname` at `date`"

cmsRun scripts/chargino300GeV_ctau10cm_GEN-SIM_${iJob}.py
EOM

condor_submit scripts/run_batch_${iJob}.sub

done;


