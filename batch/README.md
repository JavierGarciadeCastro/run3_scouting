## How to run the skimmer (by year)

Run from `lxplus8`, in the `batch` folder. All the releases are built for el8.

## 2025

### Data

Setup:

```bash
. install_cmssw.sh 2025   # CMSSW_15_0_15_patch1
source /cvmfs/cms.cern.ch/crab3/crab.sh
voms-proxy-init -voms cms -rfc -valid 192:00
```

Local test:

```bash
cmsRun Scouting/NtupleMaker/test/producer_Run3.py data=True era=2025G-ScoutingPF0 nevents=2000 output=test.root inputs=root://xrootd-cms.infn.it//store/data/Run2025G/ScoutingPF0/HLTSCOUT/v1/000/398/459/00000/893ee3a7-6176-4b82-8845-5d182eea1bf5.root
```

Full run (CRAB, if the local test runs correctly, with the same setup as the local test):

```bash
python3 crabcfg_run3.py 2025B ScoutingPFRun3   # B to D
python3 crabcfg_run3.py 2025E ScoutingPF0      # from E
python3 crabcfg_run3.py 2025E ScoutingPF1
```

### Signal

Not available yet

## 2024

### Data

Setup:

```bash
. install_cmssw.sh 2024   # CMSSW_14_0_18
source /cvmfs/cms.cern.ch/crab3/crab.sh
voms-proxy-init -voms cms -rfc -valid 192:00
```

Local test:

```bash
# 2024C
cmsRun Scouting/NtupleMaker/test/producer_Run3.py data=True era=2024C-ScoutingPFRun3 nevents=2000 output=test_2024C.root inputs=root://xrootd-cms.infn.it//store/data/Run2024C/ScoutingPFRun3/HLTSCOUT/v1/000/380/115/00000/ead74af7-66f7-4421-b11a-3d0f1be00f5b.root

# 2024D
cmsRun Scouting/NtupleMaker/test/producer_Run3.py data=True era=2024D-ScoutingPFRun3 nevents=2000 output=test_2024D.root inputs=root://cmsxrootd.fnal.gov//store/data/Run2024D/ScoutingPFRun3/HLTSCOUT/v1/000/380/360/00001/7999a7bf-8328-4922-861a-1c8286c5826e.root

# 2024E
cmsRun Scouting/NtupleMaker/test/producer_Run3.py data=True era=2024E-ScoutingPFRun3 nevents=2000 output=test_2024E.root inputs=root://xrootd-cms.infn.it//store/data/Run2024E/ScoutingPFRun3/HLTSCOUT/v1/000/381/384/00000/c12a3482-a47a-4154-864a-5b5dc26393e2.root

# 2024F
cmsRun Scouting/NtupleMaker/test/producer_Run3.py data=True era=2024F-ScoutingPFRun3 nevents=2000 output=test_2024F.root inputs=root://xrootd-cms.infn.it//store/data/Run2024F/ScoutingPFRun3/HLTSCOUT/v1/000/382/258/00000/c8d0d817-e6d2-445e-bcbb-0ce613fd75c5.root

# 2024G
cmsRun Scouting/NtupleMaker/test/producer_Run3.py data=True era=2024G-ScoutingPFRun3 nevents=2000 output=test_2024G.root inputs=root://xrootd-cms.infn.it//store/data/Run2024G/ScoutingPFRun3/HLTSCOUT/v1/000/385/764/00000/7d47b0c7-408a-4185-b4eb-6bbdcabebc74.root

# 2024H
cmsRun Scouting/NtupleMaker/test/producer_Run3.py data=True era=2024H-ScoutingPFRun3 nevents=2000 output=test_2024H.root inputs=root://xrootd-cms.infn.it//store/data/Run2024H/ScoutingPFRun3/HLTSCOUT/v1/000/385/986/00001/c2c19093-764c-4fc1-a1cd-5bec888248f1.root

# 2024I
cmsRun Scouting/NtupleMaker/test/producer_Run3.py data=True era=2024I-ScoutingPFRun3 nevents=2000 output=test_2024I.root inputs=root://xrootd-cms.infn.it//store/data/Run2024I/ScoutingPFRun3/HLTSCOUT/v1/000/386/604/00001/d146630e-0151-4b4d-9571-5c7cf98879e8.root
```

Full run (CRAB):

```bash
python3 crabcfg_run3.py 2024C ScoutingPFRun3
python3 crabcfg_run3.py 2024D ScoutingPFRun3
python3 crabcfg_run3.py 2024E ScoutingPFRun3
python3 crabcfg_run3.py 2024F ScoutingPFRun3
python3 crabcfg_run3.py 2024G ScoutingPFRun3
python3 crabcfg_run3.py 2024H ScoutingPFRun3
python3 crabcfg_run3.py 2024I ScoutingPFRun3
```

### Central signal and background

Setup:

```bash
. install_cmssw.sh 2024central   # CMSSW_15_0_2
source /cvmfs/cms.cern.ch/crab3/crab.sh
voms-proxy-init -voms cms -rfc -valid 192:00
```

Local test:

```bash
# Dark shower signal
cmsRun Scouting/NtupleMaker/test/producer_Run3.py data=False era=2024 nevents=100 output=test.root inputs=/store/mc/RunIII2024Summer24MiniAODv6/GluGluHToDarkShowers-ScenarioA_Par-ctau-1p0-mA-0p25-mpi-1_TuneCP5_13p6TeV_powheg-pythia8/MINIAODSIM/150X_mcRun3_2024_realistic_v2-v2/2560000/883f003e-5660-4d9d-be26-5099dad5aed8.root

# HAHM signal
cmsRun Scouting/NtupleMaker/test/producer_Run3.py data=False era=2024 nevents=100 output=test.root inputs=/store/mc/RunIII2024Summer24MiniAODv6/HTo2ZdTo2Mu2X_Par-MZd-5p0-Ctau-10mm_TuneCP5_13p6TeV_madgraph-pythia8/MINIAODSIM/150X_mcRun3_2024_realistic_v2-v2/2540000/b68d87fb-241a-4048-8490-d4b29347afe1.root
```

Full run (CRAB):

```bash
python3 crabcfg_run3_centralmc.py 2024 minBias_2024          # MinBias
python3 crab_throttle.py DQCD_sig --match ScenarioA --sleep 1800   # dark shower signal
python3 crabcfg_run3_centralmc.py 2024 HTo2ZdTo2mu2x   # HAHM signal, all the samples at once
```

| Model | Samples | Dataset list |
|---|---|---|
| `minBias_2024` | MinBias | `data/datasets_dqcd_2024_background.txt` |
| `DQCD_sig` | dark shower signal (without `--match`: all the scenarios) | `data/datasets_dqcd_2024_signal.txt` |
| `HTo2ZdTo2mu2x` | HAHM signal | `data/datasets_hahm_2024_signal.txt` |

### Private dark shower signal

Setup:

```bash
. install_cmssw.sh 2024private   # CMSSW_14_0_21
source /cvmfs/cms.cern.ch/crab3/crab.sh
voms-proxy-init -voms cms -rfc -valid 192:00
```

Local test:

```bash
cmsRun Scouting/NtupleMaker/test/producer_Run3.py data=False era=2024 isMiniAOD=False nevents=1000 output=test.root inputs=<file>
```

Full run (CRAB):

```bash
python3 crab_throttle.py DQCD_sig_private --sleep 1800
```

This reads `data/datasets_dqcd_2024_signal_private.txt`, and works as the throttle for the central samples.

## 2023

### Data

Setup:

```bash
. install_cmssw.sh 2023   # CMSSW_13_0_10
source /cvmfs/cms.cern.ch/crab3/crab.sh
voms-proxy-init -voms cms -rfc -valid 192:00
```

Full run (CRAB):

```bash
python3 crabcfg_run3.py 2023B                          # B, C, D
python3 crabcfg_run3.py 2023C triggerV10
python3 crabcfg_run3.py 2023B PFMonitor                # B, C, D
python3 crabcfg_run3.py 2023C PFMonitor triggerV10
```

### Central signal

Setup:

```bash
. install_cmssw.sh 2023central   # CMSSW_13_0_14
source /cvmfs/cms.cern.ch/crab3/crab.sh
voms-proxy-init -voms cms -rfc -valid 192:00
```

Full run (CRAB):

```bash
python3 crabcfg_run3_centralmc.py 2023 <model>
python3 crabcfg_run3_centralmc.py 2023BPix <model>
```

| Model | Samples | Dataset list |
|---|---|---|
| `HTo2ZdTo2mu2x` | HAHM signal | in the script |
| `BToPhi`, `BToPhiExtra` | BToPhi signal | in the script |
| `DQCD_sig` | dark shower signal | `data/datasets_dqcd_2023_1000mm.txt`, `data/datasets_dqcd_2023BPix_1000mm.txt` |

## 2022

### Data

Setup:

```bash
. install_cmssw.sh   # CMSSW_12_6_0_patch1
source /cvmfs/cms.cern.ch/crab3/crab.sh
voms-proxy-init -voms cms -rfc -valid 192:00
```

Full run (CRAB):

```bash
python3 crabcfg_run3.py 2022B                          # B to G
python3 crabcfg_run3.py 2022B PFMonitor                # B to G
```

### Central signal

Setup:

```bash
. install_cmssw.sh 2022central   # CMSSW_12_4_16
source /cvmfs/cms.cern.ch/crab3/crab.sh
voms-proxy-init -voms cms -rfc -valid 192:00
```

Full run (CRAB):

```bash
python3 crabcfg_run3_centralmc.py 2022 <model>
python3 crabcfg_run3_centralmc.py 2022postEE <model>
```

| Model | Samples | Dataset list |
|---|---|---|
| `HTo2ZdTo2mu2x` | HAHM signal | in the script |
| `BToPhi`, `BToPhiExtra` | BToPhi signal | in the script |
| `DQCD_sig` | dark shower signal | `data/datasets_dqcd_2022_1000mm.txt`, `data/datasets_dqcd_2022postEE_1000mm.txt` |

### Local MC and single datasets

```bash
python3 crabcfg_run3_mc.py <path to the sample at T2_US_UCSD>
python3 crabcfg_run3.py 2022X <dataset>
```

For example:

```bash
python3 crabcfg_run3_mc.py /ceph/cms/store/user/isuarez/ProjectMetis/DarkShower_ScenarioA_default_Run3Summer22GS_v0p30_AODSIM_v0p30
python3 crabcfg_run3.py 2022X /InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/Run3Summer22DR-Pilot_124X_mcRun3_2022_realistic_v12-v4/AODSIM
```

BToPhi and private DQCD outputs produced outside UAF can be published in DAS and read with xrootd.

## Recovering part of a data era

To process only some lumisections, with automatic splitting, set these variables and use the same command:

| Variable | Meaning |
|---|---|
| `CRAB_RECOVERY2` | JSON with the lumis to process (turns this mode on) |
| `CRAB_REQSUFFIX` | suffix of the task name, different for each submission |
| `CRAB_MINUTES` | target job length in minutes (default 720) |
| `CRAB_SITES` | comma-separated list of sites where the jobs can run |
| `CRAB_OUTDIR` | output directory, e.g. `/store/group/Run3Scouting/RAWScouting_2p0_recovery` (default: the usual one) |

```bash
export CRAB_RECOVERY2=data/<lumis>.json CRAB_REQSUFFIX=_recovery1 CRAB_MINUTES=1440
export CRAB_OUTDIR=/store/group/Run3Scouting/RAWScouting_2p0_recovery
python3 crabcfg_run3.py 2024F ScoutingPFRun3
unset CRAB_RECOVERY2 CRAB_REQSUFFIX CRAB_MINUTES CRAB_OUTDIR
```

* Find where the data are on disk with `dasgoclient -query="site dataset=/ScoutingPFRun3/Run2024F-v1/HLTSCOUT"` and choose sites close to them. Exit code 8028 means that the input could not be opened.
* Submitting the same lumis twice, and also automatic splitting, can leave the same lumi in several files. Remove the duplicates before running the looper.

## Check the skim

Run a local test before submitting, and check that the L1 bits are filled:

```python
import ROOT
f = ROOT.TFile.Open("test.root")
t = f.Get("Events")
fired = t.GetEntries("Sum$(bools_triggerMaker_l1result_SLIM.obj && doubles_triggerMaker_l1prescale_SLIM.obj==1)>0")
print(fired / t.GetEntries())
```

The fraction of events with an unprescaled L1 seed must not be 0: in 2024 it is about 0.87 in the data F file above, and 0.3 in a private signal file. Repeat the check on the CRAB output, for every sample and every task.

| Input | L1 input tag |
|---|---|
| Data and private MC | `gtStage2Digis` (unpacked in the skim) |
| Central MC | `gtStage2Digis::RECO` |

## Inputs for the looper and filler

Select the model inside the script and run:

```bash
python3 setupCentralSignals_remote.py
```

It creates:

| File | Used by |
|---|---|
| `centralDatasets.txt` | looper |
| `runScoutingLooper_onCondor_CentralSignal.sub` | looper |
| `runScoutingHistos_onCondor_CentralSignal.sub` | filler |
