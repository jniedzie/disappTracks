### Chargino samples generation instructions

1. Setup CMSSW

* go to the working directory and initialize CMSSW:

```
cmsrel CMSSW_9_4_6_patch1
cd CMSSW_9_4_6_patch1
cmsenv
```

2. Install additional necessary stuff

* For your convenience, all required commands are packed in a bash script. Download it to the current location ($WORK_DIR/CMSSW_9_4_6_patch1/):

```
wget -c https://raw.githubusercontent.com/jniedzie/disappTracks/master/eventsGeneration/install.sh
chmod 777 install.sh
```

* Run the script and go get a coffee (a small one):

`./install.sh`

3. Generate GEN-SIM events

* once everything is set up, we need to prepare a script for samples generation. Go th the $WORK_DIR/CMSSW_9_4_6_patch1/src/:

`cd src`

* run this scripts with proper arguments and go get a coffee (a big one if you scheduled many jobs):

`./submitJobs.sh number_of_jobs number_of_events_per_job`

* If everything goes fine, after many hours you should see files like `chargino300GeV_ctau10cm_GEN-SIM_0.root` containing generated events in `generatedEvents` directory. To check if jobs are running, use `condor_q` command.
