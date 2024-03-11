from CRABClient.UserUtilities import config #, getUsernameFromSiteDB
from CRABAPI.RawCommand import crabCommand

# https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3ConfigurationFile
#config = config()

import sys

era = sys.argv[1] # [2022, 2022postEE, 2023 or 2023BPix]
signal = sys.argv[2] # [HTo2ZdTo2mu2x, BToPhi]

year=0
if ("2022") in era:
    year=2022
elif ("2023") in era:
    year=2023
else:
    quit()

# This is only MC, should not be used to run on data
data=False

# ntuple version defined now
ntuple_version = "5p0"

# Setup working environment
import os
base = os.environ["CMSSW_BASE"]
#if not os.path.exists(base + '/src/centralTasks'):
#
#    os.mkdir(base + '/src/centralTasks')

# List of datasets
config_list = []
extras = []
if (len(sys.argv)>2):
    mass_points = []
    config = config()
    config.General.workArea = base+'/../'
    config.General.transferLogs = True
    config.JobType.pluginName = 'Analysis'
    config.JobType.psetName = 'Scouting/NtupleMaker/test/producer_Run3.py'
    config.Data.splitting = 'EventAwareLumiBased'
    config.Data.unitsPerJob = int(10e4)
    config.Data.publication = False # By defailt but set to true below
    config.Site.storageSite = "T2_US_UCSD"
    if "HTo2ZdTo2mu2x" in sys.argv[2]:
        config.Data.outLFNDirBase = '/store/group/Run3Scouting/RAWScouting_HTo2ZdTo2mu2x_v'++ntuple_version # DB no
        config.Data.inputDBS = 'global'
        # Set the points to produce
        mass_points.append(['0p5', '1'])
        mass_points.append(['0p5', '10'])
        mass_points.append(['0p7', '1'])
        mass_points.append(['0p7', '10'])
        mass_points.append(['1p5', '1'])
        mass_points.append(['1p5', '10'])
        mass_points.append(['1p5', '100'])
        mass_points.append(['2p0', '1'])
        mass_points.append(['2p0', '10'])
        mass_points.append(['2p0', '100'])
        mass_points.append(['2p5', '1'])
        mass_points.append(['2p5', '10'])
        mass_points.append(['2p5', '100'])
        mass_points.append(['3p0', '1'])
        mass_points.append(['3p0', '10'])
        mass_points.append(['3p0', '100'])
        mass_points.append(['4p0', '1'])
        mass_points.append(['4p0', '10'])
        mass_points.append(['4p0', '100'])
        mass_points.append(['5p0', '1'])
        mass_points.append(['5p0', '10'])
        mass_points.append(['5p0', '100'])
        mass_points.append(['6p0', '1'])
        mass_points.append(['6p0', '10'])
        mass_points.append(['6p0', '100'])
        mass_points.append(['7p0', '1'])
        mass_points.append(['7p0', '10'])
        mass_points.append(['7p0', '100'])
        mass_points.append(['8p0', '1'])
        mass_points.append(['8p0', '10'])
        mass_points.append(['8p0', '100'])
        mass_points.append(['10p0', '1'])
        mass_points.append(['10p0', '10'])
        mass_points.append(['10p0', '100'])
        mass_points.append(['12p0', '1'])
        mass_points.append(['12p0', '10'])
        mass_points.append(['12p0', '100'])
        mass_points.append(['14p0', '1'])
        mass_points.append(['14p0', '10'])
        mass_points.append(['14p0', '100'])
        mass_points.append(['16p0', '1'])
        mass_points.append(['16p0', '10'])
        mass_points.append(['16p0', '100'])
        mass_points.append(['20p0', '1'])
        mass_points.append(['20p0', '10'])
        mass_points.append(['20p0', '100'])
        mass_points.append(['22p0', '1'])
        mass_points.append(['22p0', '10'])
        mass_points.append(['22p0', '100'])
        mass_points.append(['24p0', '1'])
        mass_points.append(['24p0', '10'])
        mass_points.append(['24p0', '100'])
        mass_points.append(['30p0', '1'])
        mass_points.append(['30p0', '10'])
        mass_points.append(['30p0', '100'])
        mass_points.append(['30p0', '1000'])
        mass_points.append(['34p0', '1'])
        mass_points.append(['34p0', '10'])
        mass_points.append(['34p0', '100'])
        mass_points.append(['34p0', '1000'])
        mass_points.append(['40p0', '1'])
        mass_points.append(['40p0', '10'])
        mass_points.append(['40p0', '100'])
        mass_points.append(['40p0', '1000'])
        mass_points.append(['44p0', '1'])
        mass_points.append(['44p0', '10'])
        mass_points.append(['44p0', '100'])
        mass_points.append(['44p0', '1000'])
        mass_points.append(['50p0', '1'])
        mass_points.append(['50p0', '10'])
        mass_points.append(['50p0', '100'])
        mass_points.append(['50p0', '1000'])
        for [m,t] in mass_points:
            config_list.append(config)
            config_list[-1].JobType.pyCfgParams=["era={}".format(era),"data=False",]
            if era=="2022":
                dataset_name = '/HTo2ZdTo2mu2x_MZd-{}_ctau-{}mm_TuneCP5_13p6TeV_madgraph-pythia8/Run3Summer22DRPremix-124X_mcRun3_2022_realistic_v12-v2/AODSIM'.format(m, t)
                config_list[-1].Data.inputDataset = dataset_name
            elif era=="2022postEE":
                dataset_name = '/HTo2ZdTo2mu2x_MZd-{}_ctau-{}mm_TuneCP5_13p6TeV_madgraph-pythia8/Run3Summer22EEDRPremix-124X_mcRun3_2022_realistic_postEE_v1-v2/AODSIM'.format(m, t)
                config_list[-1].Data.inputDataset = dataset_name
            elif era=="2023":
                dataset_name = '/HTo2ZdTo2mu2x_MZd-{}_ctau-{}mm_TuneCP5_13p6TeV_madgraph-pythia8/Run3Summer23DRPremix-130X_mcRun3_2023_realistic_v14-v1/AODSIM'.format(m, t)
                config_list[-1].Data.inputDataset = dataset_name
            elif era=="2023BPix":
                dataset_name = '/HTo2ZdTo2mu2x_MZd-{}_ctau-{}mm_TuneCP5_13p6TeV_madgraph-pythia8/Run3Summer23BPixDRPremix-130X_mcRun3_2023_realistic_postBPix_v2-v1/AODSIM'.format(m, t)
                config_list[-1].Data.inputDataset = dataset_name
            config_list[-1].General.requestName = 'centralSkim__{}_{}_m-{}_ctau-{}mm_{}'.format(signal, era, m, t, ntuple_version)
            print(config)
            crabCommand('submit', config = config, dryrun = False) ## dryrun = True for local test
    elif "BToPhi" in sys.argv[2]:
        # This setup is provisional as it is tested with private signal crab produced samples
        #   -> Will be replaced by central datasets when done
        config.Data.outLFNDirBase = '/store/group/Run3Scouting/RAWScouting_privBToPhi_v'+ntuple_version # DB no
        config.Data.inputDBS = 'phys03'
        config.Data.splitting = 'FileBased'
        config.Data.publication = True
        config.Data.unitsPerJob = int(100) # Increased to match 10 jobs per file aprox
        config.Data.outputDatasetTag = "private-Skim_{era}-v1".format(era=era)
        if era=="2022":
            mass_points.append(['BToPhi_MPhi-2p0_ctau-1mm', '/BToPhi_MPhi-2p0_ctau-1mm-pythia8/Run3Scouting-private-Run3Summer22-AODSIM-b87ef10f6cfee71a9c25d28c950fbc4d/USER'])
            mass_points.append(['BToPhi_MPhi-2p0_ctau-10mm', '/BToPhi_MPhi-2p0_ctau-10mm-pythia8/Run3Scouting-private-Run3Summer22-AODSIM-b87ef10f6cfee71a9c25d28c950fbc4d/USER'])
            mass_points.append(['BToPhi_MPhi-2p0_ctau-100mm', '/BToPhi_MPhi-2p0_ctau-100mm-pythia8/Run3Scouting-private-Run3Summer22-AODSIM-b87ef10f6cfee71a9c25d28c950fbc4d/USER'])
        elif era=="2022postEE":
            mass_points.append(['BToPhi_MPhi-2p0_ctau-1mm', '/BToPhi_MPhi-2p0_ctau-1mm-pythia8/Run3Scouting-private-Run3Summer22EE-AODSIM-59a22edf0600a784f6c900595d24e883/USER'])
            mass_points.append(['BToPhi_MPhi-2p0_ctau-10mm', '/BToPhi_MPhi-2p0_ctau-10mm-pythia8/Run3Scouting-private-Run3Summer22EE-AODSIM-59a22edf0600a784f6c900595d24e883/USER'])
            mass_points.append(['BToPhi_MPhi-2p0_ctau-100mm', '/BToPhi_MPhi-2p0_ctau-100mm-pythia8/Run3Scouting-private-Run3Summer22EE-AODSIM-59a22edf0600a784f6c900595d24e883/USER'])
        for [signal_name,dataset_name] in mass_points:
            config_list.append(config)
            config_list[-1].JobType.pyCfgParams=["era={}".format(era),"data=False",]
            config_list[-1].Data.inputDataset = dataset_name
            config_list[-1].General.requestName = 'centralSkim_{}_{}_{}'.format(signal_name, era, ntuple_version)
            print(config_list[-1])
            crabCommand('submit', config = config_list[-1], dryrun = False) ## dryrun = True for local test
            #print(config)
            #crabCommand('submit', config = config, dryrun = False) ## dryrun = True for local test
    #elif "[signal]" in sys.argv[2]: (<--- Add additional signals here)
    else:
        quit()
    
else:
    quit()



