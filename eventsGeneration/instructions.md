### Chargino samples generation instructions

1. **Setup CMSSW**

* go to the working directory and initialize CMSSW:

```
cmsrel CMSSW_9_4_6_patch1
cd CMSSW_9_4_6_patch1
cmsenv
```

2. **Install additional necessary stuff**

* For your convenience, all required commands are packed in a bash script. Download it to the current location ($WORK_DIR/CMSSW_9_4_6_patch1/):

```
wget https://raw.githubusercontent.com/jniedzie/disappTracks/master/eventsGeneration/install.sh
chmod 777 install.sh
```

* Run the script and go get a coffee:

```
./install.sh
```

3. **Generate GEN-SIM events**

* once everything is set up, we need to prepare a script for samples generation. Go th the $WORK_DIR/CMSSW_9_4_6_patch1/src/:

`cd src`

* run this command (**modify the last parameter to set desired number of events!**):

```
cmsDriver.py Configuration/GenProduction/python/ThirteenTeV/AMSB_chargino/AMSB_chargino300GeV_ctau10cm_NoFilter_13TeV.py \
--fileout file:chargino300GeV_ctau10cm_GEN-SIM.root \
--step GEN,SIM \
--mc \
--datatier GEN-SIM \
--beamspot Realistic25ns13TeVEarly2017Collision \
--conditions auto:phase1_2017_realistic \
--eventcontent RAWSIM \
--era Run2_2017 \
--python_filename chargino300GeV_ctau10cm_GEN-SIM.py \
--customise Configuration/DataProcessing/Utils.addMonitoring,SimG4Core/CustomPhysics/Exotica_HSCP_SIM_cfi,SimG4Core/CustomPhysics/GenPlusSimParticles_cfi.customizeProduce,SimG4Core/CustomPhysics/GenPlusSimParticles_cfi.customizeKeep,Configuration/GenProduction/RandomSeed_cfi.customizeRandomSeed \
--no_exec \
-n 10
```

* then, you can just run the script:

`cmsRun chargino300GeV_ctau10cm_GEN-SIM.py`

* If everything goes fine, you should see `chargino300GeV_ctau10cm_GEN-SIM.root` file containing generated events.
