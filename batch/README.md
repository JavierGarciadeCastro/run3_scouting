## Run3 Info
*(Last updated: 31 Jan 2025)*

The RAW data are skimmed and/or slimmed first, as in Run 2. The skim requires at least 2 Scouting Muons and at least 1 Scouting Vertex for data.

* Do `. install_cmssw.sh`. This sets up a CMSSW environment with some EDFilters and some EDAnalyzers. (Note that `Scouting/NtupleMaker/test/producer_Run3.py` is the relevant PSet for slim/skimming)

##### crab commands (used for initial test on signal and data)

* Submit crab jobs
  * `source /cvmfs/cms.cern.ch/crab3/crab.sh` for crab commands
  * submit with

```bash
#For data/PFMonitor/public mc
python3 crabcfg_run3.py 2022B
python3 crabcfg_run3.py 2022C
python3 crabcfg_run3.py 2022D
python3 crabcfg_run3.py 2022E
python3 crabcfg_run3.py 2022F
python3 crabcfg_run3.py 2022G
python3 crabcfg_run3.py 2022B PFMonitor
python3 crabcfg_run3.py 2022C PFMonitor
python3 crabcfg_run3.py 2022D PFMonitor
python3 crabcfg_run3.py 2022E PFMonitor
python3 crabcfg_run3.py 2022F PFMonitor
python3 crabcfg_run3.py 2022G PFMonitor
python3 crabcfg_run3.py 2022X /InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/Run3Summer22DR-Pilot_124X_mcRun3_2022_realistic_v12-v4/AODSIM
python3 crabcfg_run3.py 2022X /HTo2ZdTo2mu2x_MZd-10_Epsilon-1e-06_TuneCP5_13p6TeV_madgraph-pythia8/Run3Summer22EEDRPremix-124X_mcRun3_2022_realistic_postEE_v1-v2/AODSIM
python3 crabcfg_run3.py 2022X /HTo2ZdTo2mu2x_MZd-10_Epsilon-1e-07_TuneCP5_13p6TeV_madgraph-pythia8/Run3Summer22EEDRPremix-124X_mcRun3_2022_realistic_postEE_v1-v2/AODSIM
python3 crabcfg_run3.py 2022X /HTo2ZdTo2mu2x_MZd-10_Epsilon-3e-08_TuneCP5_13p6TeV_madgraph-pythia8/Run3Summer22EEDRPremix-124X_mcRun3_2022_realistic_postEE_v1-v2/AODSIM
python3 crabcfg_run3.py 2022X /HTo2ZdTo2mu2x_MZd-10_Epsilon-5e-07_TuneCP5_13p6TeV_madgraph-pythia8/Run3Summer22EEDRPremix-124X_mcRun3_2022_realistic_postEE_v1-v2/AODSIM

#For local MC
python3 crabcfg_run3_mc.py /ceph/cms/store/user/isuarez/ProjectMetis/DarkShower_ScenarioA_default_Run3Summer22GS_v0p30_AODSIM_v0p30
python3 crabcfg_run3_mc.py /ceph/cms/store/user/isuarez/ProjectMetis/DarkShower_ScenarioA_default_Run3Summer22GS_v1p3_AODSIM_v1p3
python3 crabcfg_run3_mc.py /ceph/cms/store/user/isuarez/ProjectMetis/DarkShower_ScenarioA_default_Run3Summer22GS_v1p4_AODSIM_v1p4
python3 crabcfg_run3_mc.py /ceph/cms/store/user/isuarez/ProjectMetis/DarkShower_ScenarioA_default_Run3Summer22GS_v1p5_AODSIM_v1p5
python3 crabcfg_run3_mc.py /ceph/cms/store/user/jthakral/ProjectMetis/DarkShower_ScenarioB_default_Run3Summer22GS_v0p32_AODSIM_v0p32
python3 crabcfg_run3_mc.py /ceph/cms/store/user/jthakral/ProjectMetis/DarkShower_ScenarioC_default_Run3Summer22GS_v0p34_AODSIM_v0p34
python3 crabcfg_run3_mc.py ///path to mc signal sample, but must be on T2 ucsd for this cfg to work///

# For grid of central MC (will automatically access the published datasets for central signals, to be updated as signals are produced)
. install_cmssw.sh [2022central or 2023central]
python3 crabcfg_run3_centralmc.py [era = 2022, 2022postEE, 2023, 2023BPix] [model = HTo2ZdTo2mu2x]

#2023 data/PFMonitor - must launch from separate env under cms 13_1_X
python3 crabcfg_run3.py 2023B
python3 crabcfg_run3.py 2023C  triggerV10
python3 crabcfg_run3.py 2023C
python3 crabcfg_run3.py 2023D
python3 crabcfg_run3.py 2023B PFMonitor
python3 crabcfg_run3.py 2023C PFMonitor triggerV10
python3 crabcfg_run3.py 2023C PFMonitor
python3 crabcfg_run3.py 2023D PFMonitor
```

### Run on central signal MC:

Should setup the corresponding CMSSW version for each campaign e.g.
```bash
. install_cmssw.sh 2022central
. install_cmssw.sh 2023central
```

and then run:
```bash
python3 crabcfg_run3_centralmc.py [era] [signal]
```

### Test run (local) on 2024 central signal MC:
```bash
cd batch
. install_cmssw.sh 2024central # This should create a CMSSW_15_0_2 rel (if it already existed before the previous step, delete it before doing .install)
cmsRun Scouting/NtupleMaker/test/producer_Run3.py data=False era=2024  inputs=/store/mc/RunIII2024Summer24MiniAODv6/GluGluHToDarkShowers-ScenarioA_Par-ctau-1p0-mA-0p25-mpi-1_TuneCP5_13p6TeV_powheg-pythia8/MINIAODSIM/150X_mcRun3_2024_realistic_v2-v2/2560000/883f003e-5660-4d9d-be26-5099dad5aed8.root nevents=100 output=/tmp/skimmer_test_signal.root
```

### Full run (crab) on 2024 central signal MC:
```bash
cd batch
. install_cmssw.sh 2024central # This should create a CMSSW_15_0_2 rel (if it already existed before the previous step, delete it before doing .install)
python3 crabcfg_run3_centralmc.py 2024 DQCD_sig # Dark shower signals, reading from data/datasets_dqcd_2024_signal.txt
python3 crabcfg_run3_centralmc.py 2024 DQCD_bkg # QCD and MinBias background, reading from data/datasets_dqcd_2024_background.txt
python3 crabcfg_run3_centralmc.py 2024 HTo2ZdTo2mu2x # HAHM signal, reading from data/datasets_hahm_2024_signal.txt
```

**Note 1: Right now setup for central HAHM, central BtoPhi and part of private DQCD samples are available**

**Note 2: Since BtoPhi samples, output could be published in DAS so it can be accessed out UAF through xrd, this is used for DQCD samples produced outside**

**Note 3: If running out of UAF and/or out of SnT Scouting group, crab details should be modifed to re-direct the output to a site with writting permission**

To prepare the files that are used to run in condor the looper and the filler you just need to run:
```
python3 setupCentralSignals_remote.py
```
**Note:** You will have to select inside the model you want to run on.

Then, the following files will be created:
```
centralDatasets.txt                           # input for looper
runScoutingLooper_onCondor_CentralSignal.sub  # input for looper
runScoutingHistos_onCondor_CentralSignal.sub  # input for filler
```


## Run2 Info

RAW ntuples are generated/slimmed first, then these are made into babies.

### (Slimmed) RAW ntuples

#### Data

RAW scouting data is already available, but it stores huge FED and track info, so we first skim it with CRAB to events with
at least 2 Scouting Muons and at least 1 Scouting Vertex, and we add some extra information, while preserving the EDM structure.
So in principle, one can run the same CMSSW code on these new ntuples to add/compute other stuff.
* Do `. install_cmssw.sh`. This sets up a CMSSW environment with some EDFilters and some EDAnalyzers. (Note that `Scouting/NtupleMaker/test/producer.py` is the relevant PSet for slim/skimming)

##### Manual way
* Submit crab jobs
  * `source /cvmfs/cms.cern.ch/crab3/crab.sh` for crab commands
  * submit with
```bash
# A-D
crab submit -c crabcfg.py General.requestName="skim_2018A_v4" Data.inputDataset="/ScoutingCaloMuon/Run2018A-v1/RAW" ;
crab submit -c crabcfg.py General.requestName="skim_2018B_v4" Data.inputDataset="/ScoutingCaloMuon/Run2018B-v1/RAW" ;
crab submit -c crabcfg.py General.requestName="skim_2018C_v4" Data.inputDataset="/ScoutingCaloMuon/Run2018C-v1/RAW" ;
crab submit -c crabcfg.py General.requestName="skim_2018D_v4" Data.inputDataset="/ScoutingCaloMuon/Run2018D-v1/RAW" ;

# C again, but just the unblinded set
crab submit -c crabcfg.py General.requestName="skim_2018D_v4_unblind1fb" Data.inputDataset="/ScoutingCaloMuon/Run2018D-v1/RAW" Data.lumiMask="data/unblind_2018C_1fb_JSON.txt" Data.unitsPerJob=2000000;
```
  * Spam `crab status -c crab/<requestname>`, `crab resubmit -c ...`

##### Nicer way
* Submit crab jobs
  * `cmsenv` in the checked out release directory
  * `source /cvmfs/cms.cern.ch/crab3/crab.sh` for crab commands
  * Edit `multicrab.py` and `python multicrab.py` to submit
  * It also gives you a nice monitoring page (see the bottom of the script)

#### MC
* First submit jobs to generate and make RAW/scouting-tier data in the [generation folder](../generation/)
* Once those are done, proceed to the section on making babies below.


### Flattened babies (non-EDM)

* Run `python babymaker.py -h` to see all the options. A list of files can be specified. Note that the `-a` option will retain all events:
no skim of >=1 DV and >=2 Muon is performed (useful for MC samples where acceptance needs to be calculated).

* Clone [ProjectMetis](https://github.com/aminnj/ProjectMetis/) and source its environment
* Run the babymaker on a file locally to test
* Edit `submit_baby_jobs.py` to have the right paths (search for "hadoop") to the crab output, including the request names, output tags, and skimming options (search for "skim")
* `./make_tar.sh` to make the tarball for the jobs
* `python submit_baby_jobs.py`
