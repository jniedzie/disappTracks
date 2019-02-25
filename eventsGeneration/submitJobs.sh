#!/bin/bash

number_of_jobs=-1
number_of_events_per_job=-1
output_path="unknown"

echo "------------------------------------------------------------------------------------------------------------------------"
echo "usage: ./submitJobs.sh -j number_of_jobs -e number_of_events_per_job -s { GEN-SIM, GEN-SIM-RAW } -i input_path -o output_path"
echo "Output path is expected to be an EOS area and the directory should already exist."
echo "------------------------------------------------------------------------------------------------------------------------"

while getopts ":j:e:s:i:o:" opt; do
	case $opt in
		j) number_of_jobs=$OPTARG
		;;
		e) number_of_events_per_job=$OPTARG
		;;
		s) run_step=$OPTARG
		;;
		i) input_path=$OPTARG
		;;
		o) output_path=$OPTARG
		;;
		\?) echo "Invalid option -$OPTARG" >&2
		;;
	esac
done

if [ $number_of_jobs = -1 ]; then
	echo "Number of jobs not specified. Use -j option."
	exit
fi

if [ $run_step = "GEN-SIM" ] && [ $number_of_events_per_job = -1 ]; then
	echo "Number of events per jobs not specified. Use -e option."
	exit
fi

if [ "$output_path" = "unknown" ]; then
	echo "Output path not specified. Use -o option."
	exit
fi

mkdir -p scripts

for iJob in $(seq 1 $number_of_jobs); do
    echo "Creating job ${iJob}"
  
    if [ "$run_step" = "GEN-SIM" ]; then
        echo "Running GEN-SIM step"
        # Create HTCondor config file
        FILE="scripts/run_batch_${iJob}.sub"

        #----------------------------------------------------------------------------------------
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
        #----------------------------------------------------------------------------------------


        # Create generation script using cmsDriver
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

        # create script that will be run by HTCondor
        FILE2="scripts/run_GEN-SIM_${iJob}.sh"

        #----------------------------------------------------------------------------------------
        /bin/cat <<EOM >$FILE2
#!/bin/bash
export XRD_NETWORKSTACK=IPv4
export CMSSWVER="CMSSW_9_4_6_patch1"
export SCRAM_ARCH="slc6_amd64_gcc630"

cd `pwd`
cmsenv

source /afs/cern.ch/cms/cmsset_default.sh
eval \`scramv1 runtime -sh\`
edmPluginRefresh -p ../lib/\$SCRAM_ARCH

## Execute job and retrieve the outputs
echo "Job running on `hostname` at `date`"

cmsRun scripts/chargino300GeV_ctau10cm_GEN-SIM_${iJob}.py
xrdcp -N -v generatedEvents/chargino300GeV_ctau10cm_GEN-SIM_${iJob}.root root://eoscms.cern.ch/$output_path/
rm generatedEvents/chargino300GeV_ctau10cm_GEN-SIM_${iJob}.root
EOM
        #----------------------------------------------------------------------------------------

        # submit the job
        condor_submit scripts/run_batch_${iJob}.sub
    elif [ "$run_step" = "GEN-SIM-RAW" ]; then
        echo "Running GEN-SIM-RAW step"
        # Create HTCondor config file
        FILE="scripts/run_batch_GEN-SIM-RAW_${iJob}.sub"

        #----------------------------------------------------------------------------------------
        /bin/cat <<EOM >$FILE
executable  = scripts/run_GEN-SIM-RAW_${iJob}.sh
arguments   = \$(ClusterID) \$(ProcId)
output      = output/GEN-SIM.\$(ClusterId).\$(ProcId).out
error       = error/GEN-SIM.\$(ClusterId).\$(ProcId).err
log         = log/GEN-SIM.\$(ClusterId).log
transfer_input_files   = scripts/chargino300GeV_ctau10cm_GEN-SIM-RAW_${iJob}.py
+JobFlavour     = "nextweek"
queue
EOM
        #----------------------------------------------------------------------------------------


        # Create generation script using cmsDriver
        cmsDriver.py \
        --step DIGI,L1,DIGI2RAW,HLT \
        --datatier GEN-SIM-RAW \
        --conditions auto:phase1_2017_realistic \
        --eventcontent RAWSIM \
        --era Run2_2017 \
        --filein file:$input_path/chargino300GeV_ctau10cm_GEN-SIM_${iJob}.root \
        --fileout file:generatedEvents/chargino300GeV_ctau10cm_GEN-SIM-RAW_${iJob}.root \
        --python_filename scripts/chargino300GeV_ctau10cm_GEN-SIM-RAW_${iJob}.py \
        --customise Configuration/DataProcessing/Utils.addMonitoring,SimG4Core/CustomPhysics/GenPlusSimParticles_cfi.customizeProduce,SimG4Core/CustomPhysics/GenPlusSimParticles_cfi.customizeKeep \
        -n -1 \
        --no_exec

        # create script that will be run by HTCondor
        FILE2="scripts/run_GEN-SIM-RAW_${iJob}.sh"

        #----------------------------------------------------------------------------------------
        /bin/cat <<EOM >$FILE2
#!/bin/bash
export XRD_NETWORKSTACK=IPv4
export CMSSWVER="CMSSW_9_4_6_patch1"
export SCRAM_ARCH="slc6_amd64_gcc630"

cd `pwd`
cmsenv

source /afs/cern.ch/cms/cmsset_default.sh
eval \`scramv1 runtime -sh\`
edmPluginRefresh -p ../lib/\$SCRAM_ARCH

## Execute job and retrieve the outputs
echo "Job running on `hostname` at `date`"

cmsRun scripts/chargino300GeV_ctau10cm_GEN-SIM-RAW_${iJob}.py
xrdcp -N -v generatedEvents/chargino300GeV_ctau10cm_GEN-SIM-RAW_${iJob}.root root://eoscms.cern.ch/$output_path/
rm generatedEvents/chargino300GeV_ctau10cm_GEN-SIM-RAW_${iJob}.root
EOM
        #----------------------------------------------------------------------------------------

        # submit the job
        condor_submit scripts/run_batch_GEN-SIM-RAW_${iJob}.sub
    else
      echo "Unknown step passed: ${run_step}"
    fi

done;


