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
wget -N https://github.com/jniedzie/disappTracks/raw/master/eventsGeneration/install.sh
chmod 777 install.sh
```

* Run the script and go get a coffee (a small one):

`./install.sh`

3. **Generate GEN-SIM events**

* once everything is set up, go to the $WORK_DIR/CMSSW_9_4_6_patch1/src/:

`cd src`

* run this scripts with proper arguments and go get a coffee (a big one if you scheduled many jobs):

`./submitJobs.sh number_of_jobs number_of_events_per_job`

* If everything goes fine, after many hours you should see files like `chargino300GeV_ctau10cm_GEN-SIM_0.root` containing generated events in `generatedEvents` directory. To check if jobs are running, use `condor_q` command.


4. **Some random info that may be useful**

* if you have some problems with jobs submission, try adding something like this to your ~/.bashrc file:
`export X509_USER_PROXY=/afs/cern.ch/user/a/aalibaba/x509up_u12345`
changing the path to point to an existing file in your user directory
* make sure that you are in BASH shell
* sometimes, for some reason, wget doesn't manage to override previously downloaded files. If you want to update things, it better to remove old files (from: $WORK_DIR/CMSSW_9_4_6_patch1/ directory):

```
rm install.sh
rm src/Configuration/GenProduction/python/ThirteenTeV/AMSB_chargino/AMSB_chargino300GeV_ctau10cm_NoFilter_13TeV.py
rm src/Configuration/GenProduction/python/RandomSeed_cfi.py
rm -fr src/DisappTrks/SignalMC/data/geant4
rm src/submitJobs.sh
```
