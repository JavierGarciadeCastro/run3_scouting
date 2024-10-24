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
else:
    quit()

data=False
config.JobType.pyCfgParams=["era={}".format(era),"data=False",]
if era in ["2022B", "2022C", "2022D", "2022E", "2022F", "2022G", "2023B", "2023C", "2023D"]:
    data=True

extra=""
if (len(sys.argv)>2):
    if (data and "PFMonitor" in sys.argv[2] and len(sys.argv)>3 and "triggerV10" in sys.argv[3]): #triggerV10
        extra=sys.argv[2]+"_"+sys.argv[3]
        config.Data.inputDataset = '/ScoutingPFMonitor/Run{}-v1/RAW'.format(era)
        config.JobType.pyCfgParams=["era={}".format(era+"-triggerV10"),"data=True","monitor=True"]  
    elif (data and "PFMonitor" not in sys.argv[2]): #triggerV10
        extra=sys.argv[2]
        config.Data.inputDataset = '/ScoutingPFRun3/Run{}-v1/RAW'.format(era)
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

ntuple_version = "5"

config.General.requestName = 'skim4__{}_{}_{}'.format(
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
    config.Data.unitsPerJob = int(10e6/3)
else:
    config.Data.unitsPerJob = int(10e4)

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

#edit the area and user name
config.Data.outLFNDirBase = '/store/group/Run3Scouting/RAWScouting_'+ntuple_version # DB no
config.Data.publication = False
config.Site.storageSite = "T2_US_UCSD"

print(config)
crabCommand('submit', config = config, dryrun = False) ## dryrun = True for local test
