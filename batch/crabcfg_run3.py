from CRABClient.UserUtilities import config #, getUsernameFromSiteDB
from CRABAPI.RawCommand import crabCommand

# https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3ConfigurationFile
config = config()

import sys

era = sys.argv[1]

year=0
if ("2022") in era:
    year=2022
elif ("2023") in era:
    year=2023
elif ("2024") in era:
    year=2024
elif ("2025") in era:
    year=2025
else:
    quit()

data=False
config.JobType.pyCfgParams=["era={}".format(era),"data=False",]
if era in ["2022B", "2022C", "2022D", "2022E", "2022F", "2022G", "2023B", "2023C", "2023D", "2024C", "2024D", "2024E", "2024F", "2024G", "2024H", "2024I"]:
    print("Era identified")
    data=True
if era in ["2025B", "2025C", "2025D", "2025E", "2025F", "2025G"]:
    print("Era identified (2025)")
    data=True


extra=""
if (len(sys.argv)>2):
    if (data and "PFMonitor" in sys.argv[2] and len(sys.argv)>3 and "triggerV10" in sys.argv[3]): #triggerV10
        extra=sys.argv[2]+"_"+sys.argv[3]
        config.Data.inputDataset = '/ScoutingPFMonitor/Run{}-v1/RAW'.format(era)
        config.JobType.pyCfgParams=["era={}".format(era+"-triggerV10"),"data=True","monitor=True"]  
    elif (data and "PFMonitor" not in sys.argv[2] and (year==2022 or year==2023)): # Data in 2022 and 2023
        print("A")
        extra=sys.argv[2]
        config.Data.inputDataset = '/ScoutingPFRun3/Run{}-v1/RAW'.format(era)
        config.JobType.pyCfgParams=["era={}".format(era+"-"+extra),"data=True",]
    elif (data and "PFMonitor" not in sys.argv[2] and (year==2024 or year==2025)): # Data in 2024 and 2025
        print("B")
        extra=sys.argv[2]
        config.Data.inputDataset = '/ScoutingPFRun3/Run{}-v1/HLTSCOUT'.format(era)
        config.JobType.pyCfgParams=["era={}".format(era+"-"+extra),"data=True",]
    elif(data): #PFMonitor
        extra=sys.argv[2]
        config.Data.inputDataset = '/ScoutingPFMonitor/Run{}-v1/RAW'.format(era)
        config.JobType.pyCfgParams=["era={}".format(era),"data=True","monitor=True"]
    else: #MC
        config.Data.inputDataset = sys.argv[2]
        extra = sys.argv[2].split("/")[1]
elif(data): #other data
    config.Data.inputDataset = '/ScoutingPFRun3/Run{}-v1/RAW'.format(era)
    config.JobType.pyCfgParams=["era={}".format(era),"data=True",]
else:
    quit()

if (data and year==2025 and len(sys.argv)>2 and "PFMonitor" not in sys.argv[2]):
    extra=sys.argv[2]
    config.Data.inputDataset = '/{}/Run{}-v1/HLTSCOUT'.format(extra, era)
    config.JobType.pyCfgParams=["era={}".format(era+"-"+extra),"data=True",]

ntuple_version = "2p0"

config.General.requestName = 'skim__{}_{}_{}'.format(
        era,
        extra,
        ntuple_version,
        )

import os
base = os.environ["CMSSW_BASE"]
config.General.workArea = base+'/..'

config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'Scouting/NtupleMaker/test/producer_Run3.py'

config.Data.splitting = 'EventAwareLumiBased'

if (data):
    config.Data.unitsPerJob = int(10e7/3)
else:
    config.Data.unitsPerJob = int(10e4)
if (data and (year==2025)):
    config.Data.unitsPerJob = 8000000


#something like this can be useful for limited disk availability
#config.Data.inputBlocks = [
#'/ScoutingPFRun3/Run2022A-v1/RAW#eb476276-22f6-47f0-8e0a-17245b460227',
#'/ScoutingPFRun3/Run2022A-v1/RAW#ff292f67-5ce4-4240-92da-949bda18ade9',
#'/ScoutingPFRun3/Run2022A-v1/RAW#ff2aeb80-96a6-4007-a019-ff75eed7a527',
#]

if (data and year==2023):
    config.Data.lumiMask = "data/Cert_Collisions2023_366442_370790_Golden.json"
    if(era=="2023C" and "triggerV10" in extra):
        config.Data.lumiMask = "data/Cert_Collisions2023_eraC_367095_368823_Golden_1.json"
    elif(era=="2023C"):
        config.Data.lumiMask = "data/Cert_Collisions2023_eraC_367095_368823_Golden_2.json"

if (data and year==2022):
   config.Data.lumiMask = "data/Cert_Collisions2022_355100_362760_Golden.json"

if (data and year==2024):
   config.Data.lumiMask = "data/Cert_Collisions2024_378981_386951_Golden.json"

if (data and year==2025):
   config.Data.lumiMask = "data/Cert_Collisions2025_391658_398903_Golden.json"


#edit the area and user name
config.Data.outLFNDirBase = '/store/group/Run3Scouting/RAWScouting_'+ntuple_version # DB no
config.Data.publication = False
config.Site.storageSite = "T2_US_UCSD"
if os.environ.get("CRAB_RECOVERY2"):
    config.General.requestName = config.General.requestName + os.environ.get("CRAB_REQSUFFIX", "_recovery2")
    config.Data.splitting = 'Automatic'
    config.Data.unitsPerJob = int(os.environ.get("CRAB_MINUTES", "720"))
    config.Data.lumiMask = os.environ["CRAB_RECOVERY2"]
    config.Data.outLFNDirBase = os.environ.get("CRAB_OUTDIR", config.Data.outLFNDirBase)
    config.JobType.maxMemoryMB = 2500
    config.Site.whitelist = ['T2_US_UCSD', 'T2_US_Caltech', 'T2_US_Wisconsin', 'T2_US_Purdue', 'T2_US_Nebraska', 'T2_US_MIT', 'T2_US_Florida', 'T2_US_Vanderbilt']
    config.Data.ignoreLocality = True
    if os.environ.get("CRAB_SITES"):
        config.Site.whitelist = os.environ["CRAB_SITES"].split(",")

us_sites = ['T2_US_UCSD', 'T2_US_Caltech', 'T2_US_Wisconsin', 'T2_US_Purdue', 'T2_US_Nebraska', 'T2_US_MIT', 'T2_US_Florida']
eu_sites = ['T2_CH_CERN', 'T2_DE_RWTH', 'T2_DE_DESY', 'T2_CH_CSCS', 'T2_BE_IIHE', 'T2_IT_Pisa', 'T2_IT_Legnaro', 'T2_ES_IFCA']
auto_eras = {
    "2024F": (381944, 383779, 3, us_sites),
    "2024G": (383780, 385813, 3, eu_sites),
    "2024H": (385814, 386408, 1, eu_sites),
    "2024I": (386409, 386951, 1, us_sites),
}

def submit():
    print(config)
    crabCommand('submit', config = config, dryrun = os.environ.get("CRAB_DRYRUN", "0")=="1")

if era in auto_eras and extra=="ScoutingPFRun3" and not os.environ.get("CRAB_RECOVERY2"):
    import json
    first, last, nparts, sites = auto_eras[era]
    lumis = {r: v for r, v in json.load(open(config.Data.lumiMask)).items() if first <= int(r) <= last}
    nlumis = {r: sum(y-x+1 for x, y in v) for r, v in lumis.items()}
    parts = [{} for i in range(nparts)]
    done = 0
    for r in sorted(lumis, key=int):
        parts[min(nparts-1, done*nparts//sum(nlumis.values()))][r] = lumis[r]
        done += nlumis[r]
    name = config.General.requestName
    config.Data.splitting = 'Automatic'
    config.Data.unitsPerJob = 1440
    config.Data.ignoreLocality = True
    config.Data.outLFNDirBase = '/store/group/Run3Scouting/RAWScouting_'+ntuple_version+'_L1seeds'
    config.JobType.maxMemoryMB = 2500
    config.Site.whitelist = os.environ["CRAB_SITES"].split(",") if os.environ.get("CRAB_SITES") else sites
    for i in range(nparts):
        label = "abc"[i]
        config.General.requestName = "{}_L1seeds_{}".format(name, label)
        if os.path.isdir(os.path.join(config.General.workArea, "crab_"+config.General.requestName)):
            print("{} was already submitted, skipping".format(config.General.requestName))
            continue
        config.Data.lumiMask = "data/L1seeds_{}_{}.json".format(era, label)
        with open(config.Data.lumiMask, "w") as f:
            json.dump(parts[i], f)
        submit()
else:
    submit()
