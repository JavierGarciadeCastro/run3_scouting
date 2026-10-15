from CRABClient.UserUtilities import config #, getUsernameFromSiteDB
from CRABAPI.RawCommand import crabCommand
from http.client import HTTPException
from CRABClient.ClientExceptions import ClientException

# https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3ConfigurationFile
#config = config()

import sys

era = sys.argv[1] # [2022, 2022postEE, 2023 or 2023BPix]
signal = sys.argv[2] # [HTo2ZdTo2mu2x, BToPhi, ScenarioB1]
point = ''
if (len(sys.argv) > 3):
    point = sys.argv[3] # mass,ctau

year=0
if ("2022") in era:
    year=2022
elif ("2023") in era:
    year=2023
elif ("2024") in era:
    year=2024
else:
    quit()

# This is only MC, should not be used to run on data
data=False

# ntuple version defined now

#ntuple_version = "vhahm_8p0"
#ntuple_version = "vdqcd_final_8p0"
#ntuple_version = "btophi_extra_8p0"
ntuple_version = "_dqcd_2024"
#ntuple_version = "_minBias_2024"

# Setup working environment
import os
base = os.environ["CMSSW_BASE"]
print(base)
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
    #config.Data.splitting = 'EventAwareLumiBased'
    #config.Data.unitsPerJob = int(10e3)
    config.Data.splitting = 'FileBased'
    config.Data.unitsPerJob = int(10)
    config.Data.publication = False # By defailt but set to true below
    config.Site.storageSite = "T2_US_UCSD"
    #config.Site.blacklist = ["T2_UK_London_IC", "T2_UK_London_Brunel", "T1_RU_JINR", "T2_FR_GRIF", "T2_ES_CIEMAT"]
    config.Site.whitelist = ["T2_US_UCSD", "T2_US_Florida", "T2_US_Wisconsin", "T2_US_Nebraska", "T2_US_Purdue"]
    config.Data.ignoreLocality = True
    if "HTo2ZdTo2mu2x" in sys.argv[2]:
        config.Data.outLFNDirBase = "/store/group/Run3Scouting/RAWScouting_HTo2ZdTo2mu2x_" + str(year) # DB no
        config.Data.inputDBS = 'global'
        config.Data.publication = True
        if era == "2024":
            # Central HAHM 2024 MiniAOD signals, published to DBS global as
            # /HTo2ZdTo2Mu2X_Par-MZd-{mass}-Ctau-{ctau}mm_TuneCP5_13p6TeV_madgraph-pythia8/RunIII2024Summer24MiniAODv6-150X_mcRun3_2024_realistic_v2-v2/MINIAODSIM
            inputfile = 'data/datasets_hahm_2024_signal.txt'
            with open(inputfile,'r') as f:
                dataset_list = f.readlines()
            for dataset_name in dataset_list:
                dataset_name = dataset_name.strip()
                if not dataset_name or dataset_name[0] == '#':
                    continue
                config_list.append(config)
                config_list[-1].JobType.pyCfgParams=["era={}".format(era),"data=False",]
                config_list[-1].Data.inputDataset = dataset_name
                m = dataset_name.split('MZd-')[1].split('-Ctau')[0]
                t = dataset_name.split('-Ctau-')[1].split('mm_')[0]
                config_list[-1].General.requestName = 'centralSkim__{}_{}_m-{}_ctau-{}mm_{}'.format(signal, era, m, t, ntuple_version)
                print(config)
                try:
                    crabCommand('submit', config = config, dryrun = False) ## dryrun = True for local test
                except (HTTPException, ClientException) as e:
                    print('centralSkim__{}_{}_m-{}_ctau-{}mm_{} cant be launched! Skipping... ({})'.format(signal, era, m, t, ntuple_version, e))
        else:
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
            ## Extension
            #mass_points.append(['1p5', '1000'])
            mass_points.append(['2p0', '1000'])
            mass_points.append(['2p5', '1000'])
            mass_points.append(['3p0', '1000'])
            mass_points.append(['4p0', '1000'])
            mass_points.append(['5p0', '1000'])
            mass_points.append(['6p0', '1000'])
            mass_points.append(['7p0', '1000'])
            mass_points.append(['8p0', '1000'])
            mass_points.append(['10p0', '1000'])
            mass_points.append(['12p0', '1000'])
            mass_points.append(['14p0', '1000'])
            mass_points.append(['16p0', '1000'])
            mass_points.append(['20p0', '1000'])
            mass_points.append(['22p0', '1000'])
            mass_points.append(['24p0', '1000'])
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
                    #dataset_name = '/HTo2ZdTo2mu2x_MZd-{}_ctau-{}mm_TuneCP5_13p6TeV_madgraph-pythia8/Run3Summer23DRPremix-130X_mcRun3_2023_realistic_v15-v2/AODSIM'.format(m, t)
                    dataset_name = '/HTo2ZdTo2mu2x_MZd-{}_ctau-{}mm_TuneCP5_13p6TeV_madgraph-pythia8/Run3Summer23DRPremix-130X_mcRun3_2023_realistic_v14-v2/AODSIM'.format(m, t)
                    config_list[-1].Data.inputDataset = dataset_name
                elif era=="2023BPix":
                    #dataset_name = '/HTo2ZdTo2mu2x_MZd-{}_ctau-{}mm_TuneCP5_13p6TeV_madgraph-pythia8/Run3Summer23BPixDRPremix-130X_mcRun3_2023_realistic_postBPix_v6-v2/AODSIM'.format(m, t)
                    dataset_name = '/HTo2ZdTo2mu2x_MZd-{}_ctau-{}mm_TuneCP5_13p6TeV_madgraph-pythia8/Run3Summer23BPixDRPremix-130X_mcRun3_2023_realistic_postBPix_v2-v2/AODSIM'.format(m, t)
                    config_list[-1].Data.inputDataset = dataset_name
                config_list[-1].General.requestName = 'centralSkim__{}_{}_m-{}_ctau-{}mm_{}'.format(signal, era, m, t, ntuple_version)
                print(config)
                crabCommand('submit', config = config, dryrun = False) ## dryrun = True for local test
    elif "BToPhiExtra" in sys.argv[2]:
        config.Data.outLFNDirBase = "/store/group/Run3Scouting/RAWScouting_BToPhiExtra_" + str(year) + "_v"+ntuple_version 
        config.Data.publication = False
        config.Data.inputDBS = 'phys03'
        config.Data.splitting = 'FileBased'
        config.Data.unitsPerJob = int(30)
        # Set the points to produce
        if era=="2022":
            mass_points.append(['4p6', '0p1', '/BToPhi_MPhi-4p6_ctau-0p1mm_TuneCP5_13p6TeV_pythia8_ext/fernance-private-HLTRAW-2022-b722d1cf11a99a4476f09a94f34c768e/USER'])
            mass_points.append(['4p6', '1', '/BToPhi_MPhi-4p6_ctau-1mm_TuneCP5_13p6TeV_pythia8_ext/fernance-private-HLTRAW-2022-b722d1cf11a99a4476f09a94f34c768e/USER'])
            mass_points.append(['4p6', '10', '/BToPhi_MPhi-4p6_ctau-10mm_TuneCP5_13p6TeV_pythia8_ext/fernance-private-HLTRAW-2022-b722d1cf11a99a4476f09a94f34c768e/USER'])
            mass_points.append(['4p6', '100', '/BToPhi_MPhi-4p6_ctau-100mm_TuneCP5_13p6TeV_pythia8_ext/fernance-private-HLTRAW-2022-b722d1cf11a99a4476f09a94f34c768e/USER'])
            mass_points.append(['0p25', '1000', '/BToPhi_MPhi-0p25_ctau-1000mm_TuneCP5_13p6TeV_pythia8_ext3/fernance-private-HLTRAW-2022-b722d1cf11a99a4476f09a94f34c768e/USER'])
            mass_points.append(['0p3', '1000', '/BToPhi_MPhi-0p3_ctau-1000mm_TuneCP5_13p6TeV_pythia8_ext3/fernance-private-HLTRAW-2022-b722d1cf11a99a4476f09a94f34c768e/USER'])
            mass_points.append(['0p4', '1000', '/BToPhi_MPhi-0p4_ctau-1000mm_TuneCP5_13p6TeV_pythia8_ext3/fernance-private-HLTRAW-2022-b722d1cf11a99a4476f09a94f34c768e/USER'])
            mass_points.append(['0p5', '1000', '/BToPhi_MPhi-0p5_ctau-1000mm_TuneCP5_13p6TeV_pythia8_ext3/fernance-private-HLTRAW-2022-b722d1cf11a99a4476f09a94f34c768e/USER'])
            mass_points.append(['0p6', '1000', '/BToPhi_MPhi-0p6_ctau-1000mm_TuneCP5_13p6TeV_pythia8_ext3/fernance-private-HLTRAW-2022-b722d1cf11a99a4476f09a94f34c768e/USER'])
            mass_points.append(['0p7', '1000', '/BToPhi_MPhi-0p7_ctau-1000mm_TuneCP5_13p6TeV_pythia8_ext3/fernance-private-HLTRAW-2022-b722d1cf11a99a4476f09a94f34c768e/USER'])
            mass_points.append(['0p9', '1000', '/BToPhi_MPhi-0p9_ctau-1000mm_TuneCP5_13p6TeV_pythia8_ext3/fernance-private-HLTRAW-2022-b722d1cf11a99a4476f09a94f34c768e/USER'])
            mass_points.append(['1p25', '1000', '/BToPhi_MPhi-1p25_ctau-1000mm_TuneCP5_13p6TeV_pythia8_ext3/fernance-private-HLTRAW-2022-b722d1cf11a99a4476f09a94f34c768e/USER'])
            mass_points.append(['1p5', '1000', '/BToPhi_MPhi-1p5_ctau-1000mm_TuneCP5_13p6TeV_pythia8_ext3/fernance-private-HLTRAW-2022-b722d1cf11a99a4476f09a94f34c768e/USER'])
            mass_points.append(['2p0', '1000', '/BToPhi_MPhi-2p0_ctau-1000mm_TuneCP5_13p6TeV_pythia8_ext3/fernance-private-HLTRAW-2022-b722d1cf11a99a4476f09a94f34c768e/USER'])
            mass_points.append(['2p85', '1000', '/BToPhi_MPhi-2p85_ctau-1000mm_TuneCP5_13p6TeV_pythia8_ext3/fernance-private-HLTRAW-2022-b722d1cf11a99a4476f09a94f34c768e/USER'])
            mass_points.append(['3p35', '1000', '/BToPhi_MPhi-3p35_ctau-1000mm_TuneCP5_13p6TeV_pythia8_ext3/fernance-private-HLTRAW-2022-b722d1cf11a99a4476f09a94f34c768e/USER'])
            mass_points.append(['4p0', '1000', '/BToPhi_MPhi-4p0_ctau-1000mm_TuneCP5_13p6TeV_pythia8_ext3/fernance-private-HLTRAW-2022-b722d1cf11a99a4476f09a94f34c768e/USER'])
            mass_points.append(['4p6', '1000', '/BToPhi_MPhi-4p6_ctau-1000mm_TuneCP5_13p6TeV_pythia8_ext3/fernance-private-HLTRAW-2022-b722d1cf11a99a4476f09a94f34c768e/USER'])
        if era=="2022postEE":
            mass_points.append(['0p25', '1000', '/BToPhi_MPhi-0p25_ctau-1000mm_TuneCP5_13p6TeV_pythia8_ext3/fernance-private-HLTRAW-2022postEE-31b3ee15c0b04cb98bdf3666445e3e75/USER'])
            mass_points.append(['0p3', '1000', '/BToPhi_MPhi-0p3_ctau-1000mm_TuneCP5_13p6TeV_pythia8_ext3/fernance-private-HLTRAW-2022postEE-31b3ee15c0b04cb98bdf3666445e3e75/USER'])
            mass_points.append(['0p4', '1000', '/BToPhi_MPhi-0p4_ctau-1000mm_TuneCP5_13p6TeV_pythia8_ext3/fernance-private-HLTRAW-2022postEE-31b3ee15c0b04cb98bdf3666445e3e75/USER'])
            mass_points.append(['0p5', '1000', '/BToPhi_MPhi-0p5_ctau-1000mm_TuneCP5_13p6TeV_pythia8_ext3/fernance-private-HLTRAW-2022postEE-31b3ee15c0b04cb98bdf3666445e3e75/USER'])
            mass_points.append(['0p6', '1000', '/BToPhi_MPhi-0p6_ctau-1000mm_TuneCP5_13p6TeV_pythia8_ext3/fernance-private-HLTRAW-2022postEE-31b3ee15c0b04cb98bdf3666445e3e75/USER'])
            mass_points.append(['0p7', '1000', '/BToPhi_MPhi-0p7_ctau-1000mm_TuneCP5_13p6TeV_pythia8_ext3/fernance-private-HLTRAW-2022postEE-31b3ee15c0b04cb98bdf3666445e3e75/USER'])
            mass_points.append(['0p9', '1000', '/BToPhi_MPhi-0p9_ctau-1000mm_TuneCP5_13p6TeV_pythia8_ext3/fernance-private-HLTRAW-2022postEE-31b3ee15c0b04cb98bdf3666445e3e75/USER'])
            mass_points.append(['1p5', '1000', '/BToPhi_MPhi-1p5_ctau-1000mm_TuneCP5_13p6TeV_pythia8_ext3/fernance-private-HLTRAW-2022postEE-31b3ee15c0b04cb98bdf3666445e3e75/USER'])
            mass_points.append(['2p0', '1000', '/BToPhi_MPhi-2p0_ctau-1000mm_TuneCP5_13p6TeV_pythia8_ext3/fernance-private-HLTRAW-2022postEE-31b3ee15c0b04cb98bdf3666445e3e75/USER'])
            mass_points.append(['2p85', '1000', '/BToPhi_MPhi-2p85_ctau-1000mm_TuneCP5_13p6TeV_pythia8_ext3/fernance-private-HLTRAW-2022postEE-31b3ee15c0b04cb98bdf3666445e3e75/USER'])
            mass_points.append(['3p35', '1000', '/BToPhi_MPhi-3p35_ctau-1000mm_TuneCP5_13p6TeV_pythia8_ext3/fernance-private-HLTRAW-2022postEE-31b3ee15c0b04cb98bdf3666445e3e75/USER'])
            mass_points.append(['4p0', '1000', '/BToPhi_MPhi-4p0_ctau-1000mm_TuneCP5_13p6TeV_pythia8_ext3/fernance-private-HLTRAW-2022postEE-31b3ee15c0b04cb98bdf3666445e3e75/USER'])
            mass_points.append(['4p6', '1000', '/BToPhi_MPhi-4p6_ctau-1000mm_TuneCP5_13p6TeV_pythia8_ext3/fernance-private-HLTRAW-2022postEE-31b3ee15c0b04cb98bdf3666445e3e75/USER'])
            mass_points.append(['4p6', '0p1', '/BToPhi_MPhi-4p6_ctau-0p1mm_TuneCP5_13p6TeV_pythia8_ext/fernance-private-HLTRAW-2022postEE-31b3ee15c0b04cb98bdf3666445e3e75/USER'])
            mass_points.append(['4p6', '1', '/BToPhi_MPhi-4p6_ctau-1mm_TuneCP5_13p6TeV_pythia8_ext/fernance-private-HLTRAW-2022postEE-31b3ee15c0b04cb98bdf3666445e3e75/USER'])
            mass_points.append(['4p6', '10', '/BToPhi_MPhi-4p6_ctau-10mm_TuneCP5_13p6TeV_pythia8_ext/fernance-private-HLTRAW-2022postEE-31b3ee15c0b04cb98bdf3666445e3e75/USER'])
            mass_points.append(['4p6', '100', '/BToPhi_MPhi-4p6_ctau-100mm_TuneCP5_13p6TeV_pythia8_ext/fernance-private-HLTRAW-2022postEE-31b3ee15c0b04cb98bdf3666445e3e75/USER'])
        if era=="2023":
            mass_points.append(['0p25', '1000', '/BToPhi_MPhi-0p25_ctau-1000mm_TuneCP5_13p6TeV_pythia8_ext3/fernance-private-HLTRAW-2023-c017b2c35ae16f5766f4c67c30206b8e/USER'])
            mass_points.append(['0p3', '1000', '/BToPhi_MPhi-0p3_ctau-1000mm_TuneCP5_13p6TeV_pythia8_ext3/fernance-private-HLTRAW-2023-c017b2c35ae16f5766f4c67c30206b8e/USER'])
            mass_points.append(['0p4', '1000', '/BToPhi_MPhi-0p4_ctau-1000mm_TuneCP5_13p6TeV_pythia8_ext3/fernance-private-HLTRAW-2023-c017b2c35ae16f5766f4c67c30206b8e/USER'])
            mass_points.append(['0p5', '1000', '/BToPhi_MPhi-0p5_ctau-1000mm_TuneCP5_13p6TeV_pythia8_ext3/fernance-private-HLTRAW-2023-c017b2c35ae16f5766f4c67c30206b8e/USER'])
            mass_points.append(['0p6', '1000', '/BToPhi_MPhi-0p6_ctau-1000mm_TuneCP5_13p6TeV_pythia8_ext3/fernance-private-HLTRAW-2023-c017b2c35ae16f5766f4c67c30206b8e/USER'])
            mass_points.append(['0p7', '1000', '/BToPhi_MPhi-0p7_ctau-1000mm_TuneCP5_13p6TeV_pythia8_ext3/fernance-private-HLTRAW-2023-c017b2c35ae16f5766f4c67c30206b8e/USER'])
            mass_points.append(['0p9', '1000', '/BToPhi_MPhi-0p9_ctau-1000mm_TuneCP5_13p6TeV_pythia8_ext3/fernance-private-HLTRAW-2023-c017b2c35ae16f5766f4c67c30206b8e/USER'])
            mass_points.append(['1p25', '1000', '/BToPhi_MPhi-1p25_ctau-1000mm_TuneCP5_13p6TeV_pythia8_ext3/fernance-private-HLTRAW-2023-c017b2c35ae16f5766f4c67c30206b8e/USER'])
            mass_points.append(['1p5', '1000', '/BToPhi_MPhi-1p5_ctau-1000mm_TuneCP5_13p6TeV_pythia8_ext3/fernance-private-HLTRAW-2023-c017b2c35ae16f5766f4c67c30206b8e/USER'])
            mass_points.append(['2p0', '1000', '/BToPhi_MPhi-2p0_ctau-1000mm_TuneCP5_13p6TeV_pythia8_ext3/fernance-private-HLTRAW-2023-c017b2c35ae16f5766f4c67c30206b8e/USER'])
            mass_points.append(['2p85', '1000', '/BToPhi_MPhi-2p85_ctau-1000mm_TuneCP5_13p6TeV_pythia8_ext3/fernance-private-HLTRAW-2023-c017b2c35ae16f5766f4c67c30206b8e/USER'])
            mass_points.append(['3p35', '1000', '/BToPhi_MPhi-3p35_ctau-1000mm_TuneCP5_13p6TeV_pythia8_ext3/fernance-private-HLTRAW-2023-c017b2c35ae16f5766f4c67c30206b8e/USER'])
            mass_points.append(['4p0', '1000', '/BToPhi_MPhi-4p0_ctau-1000mm_TuneCP5_13p6TeV_pythia8_ext3/fernance-private-HLTRAW-2023-c017b2c35ae16f5766f4c67c30206b8e/USER'])
            mass_points.append(['4p6', '1000', '/BToPhi_MPhi-4p6_ctau-1000mm_TuneCP5_13p6TeV_pythia8_ext3/fernance-private-HLTRAW-2023-c017b2c35ae16f5766f4c67c30206b8e/USER'])
            mass_points.append(['4p6', '0p1', '/BToPhi_MPhi-4p6_ctau-0p1mm_TuneCP5_13p6TeV_pythia8_ext/fernance-private-HLTRAW-2023-c017b2c35ae16f5766f4c67c30206b8e/USER'])
            mass_points.append(['4p6', '1', '/BToPhi_MPhi-4p6_ctau-1mm_TuneCP5_13p6TeV_pythia8_ext/fernance-private-HLTRAW-2023-c017b2c35ae16f5766f4c67c30206b8e/USER'])
            mass_points.append(['4p6', '10', '/BToPhi_MPhi-4p6_ctau-10mm_TuneCP5_13p6TeV_pythia8_ext/fernance-private-HLTRAW-2023-c017b2c35ae16f5766f4c67c30206b8e/USER'])
            mass_points.append(['4p6', '100', '/BToPhi_MPhi-4p6_ctau-100mm_TuneCP5_13p6TeV_pythia8_ext/fernance-private-HLTRAW-2023-c017b2c35ae16f5766f4c67c30206b8e/USER'])
        if era=="2023BPix":
            mass_points.append(['0p25', '1000', '/BToPhi_MPhi-0p25_ctau-1000mm_TuneCP5_13p6TeV_pythia8_ext3/fernance-private-HLTRAW-2023BPix-8a4c70c5aaaf26ad44e675df103c838b/USER'])
            mass_points.append(['0p3', '1000', '/BToPhi_MPhi-0p3_ctau-1000mm_TuneCP5_13p6TeV_pythia8_ext3/fernance-private-HLTRAW-2023BPix-8a4c70c5aaaf26ad44e675df103c838b/USER'])
            mass_points.append(['0p4', '1000', '/BToPhi_MPhi-0p4_ctau-1000mm_TuneCP5_13p6TeV_pythia8_ext3/fernance-private-HLTRAW-2023BPix-8a4c70c5aaaf26ad44e675df103c838b/USER'])
            mass_points.append(['0p5', '1000', '/BToPhi_MPhi-0p5_ctau-1000mm_TuneCP5_13p6TeV_pythia8_ext3/fernance-private-HLTRAW-2023BPix-8a4c70c5aaaf26ad44e675df103c838b/USER'])
            mass_points.append(['0p6', '1000', '/BToPhi_MPhi-0p6_ctau-1000mm_TuneCP5_13p6TeV_pythia8_ext3/fernance-private-HLTRAW-2023BPix-8a4c70c5aaaf26ad44e675df103c838b/USER'])
            mass_points.append(['0p7', '1000', '/BToPhi_MPhi-0p7_ctau-1000mm_TuneCP5_13p6TeV_pythia8_ext3/fernance-private-HLTRAW-2023BPix-8a4c70c5aaaf26ad44e675df103c838b/USER'])
            mass_points.append(['0p9', '1000', '/BToPhi_MPhi-0p9_ctau-1000mm_TuneCP5_13p6TeV_pythia8_ext3/fernance-private-HLTRAW-2023BPix-8a4c70c5aaaf26ad44e675df103c838b/USER'])
            mass_points.append(['1p25', '1000', '/BToPhi_MPhi-1p25_ctau-1000mm_TuneCP5_13p6TeV_pythia8_ext3/fernance-private-HLTRAW-2023BPix-8a4c70c5aaaf26ad44e675df103c838b/USER'])
            mass_points.append(['1p5', '1000', '/BToPhi_MPhi-1p5_ctau-1000mm_TuneCP5_13p6TeV_pythia8_ext3/fernance-private-HLTRAW-2023BPix-8a4c70c5aaaf26ad44e675df103c838b/USER'])
            mass_points.append(['2p0', '1000', '/BToPhi_MPhi-2p0_ctau-1000mm_TuneCP5_13p6TeV_pythia8_ext3/fernance-private-HLTRAW-2023BPix-8a4c70c5aaaf26ad44e675df103c838b/USER'])
            mass_points.append(['2p85', '1000', '/BToPhi_MPhi-2p85_ctau-1000mm_TuneCP5_13p6TeV_pythia8_ext3/fernance-private-HLTRAW-2023BPix-8a4c70c5aaaf26ad44e675df103c838b/USER'])
            mass_points.append(['3p35', '1000', '/BToPhi_MPhi-3p35_ctau-1000mm_TuneCP5_13p6TeV_pythia8_ext3/fernance-private-HLTRAW-2023BPix-8a4c70c5aaaf26ad44e675df103c838b/USER'])
            mass_points.append(['4p0', '1000', '/BToPhi_MPhi-4p0_ctau-1000mm_TuneCP5_13p6TeV_pythia8_ext3/fernance-private-HLTRAW-2023BPix-8a4c70c5aaaf26ad44e675df103c838b/USER'])
            mass_points.append(['4p6', '1000', '/BToPhi_MPhi-4p6_ctau-1000mm_TuneCP5_13p6TeV_pythia8_ext3/fernance-private-HLTRAW-2023BPix-8a4c70c5aaaf26ad44e675df103c838b/USER'])
            mass_points.append(['4p6', '0p1', '/BToPhi_MPhi-4p6_ctau-0p1mm_TuneCP5_13p6TeV_pythia8_ext/fernance-private-HLTRAW-2023BPix-8a4c70c5aaaf26ad44e675df103c838b/USER'])
            mass_points.append(['4p6', '1', '/BToPhi_MPhi-4p6_ctau-1mm_TuneCP5_13p6TeV_pythia8_ext/fernance-private-HLTRAW-2023BPix-8a4c70c5aaaf26ad44e675df103c838b/USER'])
            mass_points.append(['4p6', '10', '/BToPhi_MPhi-4p6_ctau-10mm_TuneCP5_13p6TeV_pythia8_ext/fernance-private-HLTRAW-2023BPix-8a4c70c5aaaf26ad44e675df103c838b/USER'])
            mass_points.append(['4p6', '100', '/BToPhi_MPhi-4p6_ctau-100mm_TuneCP5_13p6TeV_pythia8_ext/fernance-private-HLTRAW-2023BPix-8a4c70c5aaaf26ad44e675df103c838b/USER'])
            #mass_points.append(['4p6', '1000', '/BToPhi_MPhi-4p6_ctau-1000mm_TuneCP5_13p6TeV_pythia8_ext/fernance-private-HLTRAW-2023BPix-8a4c70c5aaaf26ad44e675df103c838b/USER'])
        for [m,t,dataset_name] in mass_points:
            config_list.append(config)
            config_list[-1].JobType.pyCfgParams=["era={}".format(era),"data=False",]
            config_list[-1].Data.inputDataset = dataset_name
            config_list[-1].General.requestName = 'centralSkim__{}_{}_m-{}_ctau-{}mm_{}'.format(signal, era, m, t, ntuple_version)
            if os.path.exists('crab_centralSkim__BToPhi_{}_m-{}_ctau-{}mm_{}'.format(era, m, t, ntuple_version)):
                print('Skipping because it exists: ' + 'crab_centralSkim__BToPhi_{}_m-{}_ctau-{}mm_{}'.format(era, m, t, ntuple_version))
            else:
                print(config)
                crabCommand('submit', config = config, dryrun = False) ## dryrun = True for local test
    elif "BToPhi" in sys.argv[2]:
        config.Data.outLFNDirBase = "/store/group/Run3Scouting/RAWScouting_BToPhiExtra_" + str(year) + "_v"+ntuple_version 
        config.Data.inputDBS = 'global'
        config.Data.publication = True
        # Set the points to produce
        if point:
            mass_points.append(point.split(','))
        else:
            mass_points.append(['0p25', '0p0'])
            mass_points.append(['0p25', '0p1'])
            mass_points.append(['0p25', '100'])
            mass_points.append(['0p25', '10'])
            mass_points.append(['0p25', '1'])
            mass_points.append(['0p3', '0p0'])
            mass_points.append(['0p3', '0p1'])
            mass_points.append(['0p3', '100'])
            mass_points.append(['0p3', '10'])
            mass_points.append(['0p3', '1'])
            mass_points.append(['0p4', '0p0'])
            mass_points.append(['0p4', '0p1'])
            mass_points.append(['0p4', '100'])
            mass_points.append(['0p4', '10'])
            mass_points.append(['0p4', '1'])
            mass_points.append(['0p5', '0p0'])
            mass_points.append(['0p5', '0p1'])
            mass_points.append(['0p5', '100'])
            mass_points.append(['0p5', '10'])
            mass_points.append(['0p5', '1'])
            mass_points.append(['0p6', '0p0'])
            mass_points.append(['0p6', '0p1'])
            mass_points.append(['0p6', '100'])
            mass_points.append(['0p6', '10'])
            mass_points.append(['0p6', '1'])
            mass_points.append(['0p7', '0p0'])
            mass_points.append(['0p7', '0p1'])
            mass_points.append(['0p7', '100'])
            mass_points.append(['0p7', '10'])
            mass_points.append(['0p7', '1'])
            mass_points.append(['0p9', '0p0'])
            mass_points.append(['0p9', '0p1'])
            mass_points.append(['0p9', '100'])
            mass_points.append(['0p9', '10'])
            mass_points.append(['0p9', '1'])
            mass_points.append(['1p25', '0p0'])
            mass_points.append(['1p25', '0p1'])
            mass_points.append(['1p25', '100'])
            mass_points.append(['1p25', '10'])
            mass_points.append(['1p25', '1'])
            mass_points.append(['1p5', '0p0'])
            mass_points.append(['1p5', '0p1'])
            mass_points.append(['1p5', '100'])
            mass_points.append(['1p5', '10'])
            mass_points.append(['1p5', '1'])
            mass_points.append(['2p0', '0p0'])
            mass_points.append(['2p0', '0p1'])
            mass_points.append(['2p0', '100'])
            mass_points.append(['2p0', '10'])
            mass_points.append(['2p0', '1'])
            mass_points.append(['2p85', '0p0'])
            mass_points.append(['2p85', '0p1'])
            mass_points.append(['2p85', '100'])
            mass_points.append(['2p85', '10'])
            mass_points.append(['2p85', '1'])
            mass_points.append(['3p35', '0p0'])
            mass_points.append(['3p35', '0p1'])
            mass_points.append(['3p35', '100'])
            mass_points.append(['3p35', '10'])
            mass_points.append(['3p35', '1'])
            mass_points.append(['4p0', '0p0'])
            mass_points.append(['4p0', '0p1'])
            mass_points.append(['4p0', '100'])
            mass_points.append(['4p0', '10'])
            mass_points.append(['4p0', '1'])
            mass_points.append(['5p0', '0p0'])
            mass_points.append(['5p0', '0p1'])
            mass_points.append(['5p0', '100'])
            mass_points.append(['5p0', '10'])
            mass_points.append(['5p0', '1'])
        for [m,t] in mass_points:
            config_list.append(config)
            config_list[-1].JobType.pyCfgParams=["era={}".format(era),"data=False",]
            if era=="2022":
                dataset_name = '/BToPhi_MPhi-{}_ctau-{}mm_TuneCP5_13p6TeV_pythia8/Run3Summer22DRPremix-124X_mcRun3_2022_realistic_v12-v2/AODSIM'.format(m, t)
                config_list[-1].Data.inputDataset = dataset_name
            elif era=="2022postEE":
                dataset_name = '/BToPhi_MPhi-{}_ctau-{}mm_TuneCP5_13p6TeV_pythia8/Run3Summer22EEDRPremix-124X_mcRun3_2022_realistic_postEE_v1-v2/AODSIM'.format(m, t)
                config_list[-1].Data.inputDataset = dataset_name
            elif era=="2023":
                dataset_name = '/BToPhi_MPhi-{}_ctau-{}mm_TuneCP5_13p6TeV_pythia8/Run3Summer23DRPremix-130X_mcRun3_2023_realistic_v15-v2/AODSIM'.format(m, t)
                config_list[-1].Data.inputDataset = dataset_name
            elif era=="2023BPix":
                dataset_name = '/BToPhi_MPhi-{}_ctau-{}mm_TuneCP5_13p6TeV_pythia8/Run3Summer23BPixDRPremix-130X_mcRun3_2023_realistic_postBPix_v6-v2/AODSIM'.format(m, t)
                config_list[-1].Data.inputDataset = dataset_name
            config_list[-1].General.requestName = 'centralSkim__{}_{}_m-{}_ctau-{}mm_{}'.format(signal, era, m, t, ntuple_version)
            if os.path.exists('crab_centralSkim__BToPhi_{}_m-{}_ctau-{}mm_{}'.format(era, m, t, ntuple_version)):
                print('Skipping because it exists: ' + 'crab_centralSkim__BToPhi_{}_m-{}_ctau-{}mm_{}'.format(era, m, t, ntuple_version))
            else:
                print(config)
                crabCommand('submit', config = config, dryrun = False) ## dryrun = True for local test
    elif "DQCD_sig_private" in sys.argv[2]:
        # Private ScenarioA signals produced by garciaja, published to phys03/USER
        config.Data.outLFNDirBase = '/store/group/Run3Scouting/RAWScouting_DQCD_sig2024_private'
        config.Data.inputDBS = 'phys03'
        config.Data.splitting = 'FileBased'
        config.Data.unitsPerJob = int(10)
        config.Data.publication = True
        config.Data.outputDatasetTag = "private-Skim_final_{era}-v2".format(era=era)
        inputfile = 'data/datasets_dqcd_2024_signal_private.txt'
        with open(inputfile,'r') as f:
            dataset_list = f.readlines()
        for dataset_name in dataset_list:
            config_list.append(config)
            config_list[-1].JobType.pyCfgParams=["era={}".format(era),"data=False","isMiniAOD=False",]
            print(dataset_name)
            if dataset_name[0]=='#':
                continue
            config_list[-1].Data.inputDataset = dataset_name.strip()
            model_name = dataset_name.split('-')[1]
            t = dataset_name.split('ctau-')[1].split('-')[0]
            mA = dataset_name.split('mA-')[1].split('-')[0]
            mpi = dataset_name.split('mpi-')[1].split('/')[0]
            config_list[-1].General.requestName = 'centralSkim__{}_{}_mpi-{}_mA-{}_ctau-{}mm{}'.format(model_name, era, mpi, mA, t, ntuple_version)
            print(config)
            try:
                crabCommand('submit', config = config, dryrun = False)
            except:
                print('centralSkim__{}_{}_mpi-{}_mA-{}_ctau-{}mm_{} cant be launched! Skipping...'.format(model_name, era, mpi, mA, t, ntuple_version))
    elif "DQCD_sig" in sys.argv[2]:
        config.Data.outLFNDirBase = '/store/group/Run3Scouting/RAWScouting_DQCD_sig2024'
        #config.Data.inputDBS = 'phys03'
        config.Data.inputDBS = 'global'
        #config.Data.splitting = 'FileBased'
        config.Data.splitting = 'EventAwareLumiBased'
        config.Data.unitsPerJob = int(10e4)
        config.Data.publication = True
        config.Data.unitsPerJob = int(20) # Increased to match 10 jobs per file aprox
        config.Data.outputDatasetTag = "private-Skim_final_{era}-v2".format(era=era)
        if era=="2022":
            inputfile = 'data/datasets_dqcd_2022_1000mm.txt'
        if era=="2022postEE":
            inputfile = 'data/datasets_dqcd_2022postEE_1000mm.txt'
        if era=="2023":
            inputfile = 'data/datasets_dqcd_2023_1000mm.txt'
        if era=="2023BPix":
            inputfile = 'data/datasets_dqcd_2023BPix_1000mm.txt'
        if era=="2024":
            inputfile = 'data/datasets_dqcd_2024_signal.txt'
        with open(inputfile,'r') as f:
            dataset_list = f.readlines()
        for dataset_name in dataset_list:
            config_list.append(config)
            config_list[-1].JobType.pyCfgParams=["era={}".format(era),"data=False",]
            print(dataset_name)
            if dataset_name[0]=='#':
                continue
                
            if era == "2024":
                config_list[-1].Data.inputDataset = dataset_name.strip()
                model_name = dataset_name.split('-')[1]
                t = dataset_name.split('ctau-')[1].split('-')[0]
                mA = dataset_name.split('mA-')[1].split('-')[0]
                mpi = dataset_name.split('mpi-')[1].split('_')[0]
            else:
                config_list[-1].Data.inputDataset = dataset_name[:-1]
                model_name = 'S'+dataset_name.split('_')[0][2:]
                mpi = dataset_name.split('mpi_')[1].split('_')[0]
                mA = dataset_name.split('mA_')[1].split('_')[0]
                t = dataset_name.split('ctau_')[1].split('/')[0]
            config_list[-1].General.requestName = 'centralSkim__{}_{}_mpi-{}_mA-{}_ctau-{}mm{}'.format(model_name, era, mpi, mA, t, ntuple_version)
            config_list[-1].Data.inputDataset = dataset_name.strip()

            print(config)
            try:
                crabCommand('submit', config = config, dryrun = False) ## dryrun = True for local test
            except:
                print('centralSkim__{}_{}_mpi-{}_mA-{}_ctau-{}mm_{} cant be launched! Skipping...'.format(model_name, era, mpi, mA, t, ntuple_version))
    elif "ScenarioA" in sys.argv[2]:
        config.Data.outLFNDirBase = '/store/group/Run3Scouting/RAWScouting_privScenarioA_v'+ntuple_version # DB no
        config.Data.inputDBS = 'phys03'
        config.Data.splitting = 'FileBased'
        config.Data.publication = True
        config.Data.unitsPerJob = int(10) # Increased to match 10 jobs per file aprox
        config.Data.outputDatasetTag = "private-Skim_{era}-v2".format(era=era)
        mass_points.append(['1', '0p25', '0p1'])
        mass_points.append(['1', '0p25', '1p0'])
        mass_points.append(['1', '0p25', '10'])
        mass_points.append(['1', '0p25', '100'])
        mass_points.append(['1', '0p33', '0p1'])
        mass_points.append(['1', '0p33', '1p0'])
        mass_points.append(['1', '0p33', '10'])
        mass_points.append(['1', '0p33', '100'])
        mass_points.append(['1', '0p45', '0p1'])
        mass_points.append(['1', '0p45', '1p0'])
        mass_points.append(['1', '0p45', '10'])
        mass_points.append(['1', '0p45', '100'])
        mass_points.append(['2', '0p25', '0p1'])
        mass_points.append(['2', '0p25', '1p0'])
        mass_points.append(['2', '0p25', '10'])
        mass_points.append(['2', '0p25', '100'])
        mass_points.append(['2', '0p40', '0p1'])
        mass_points.append(['2', '0p40', '1p0'])
        mass_points.append(['2', '0p40', '10'])
        mass_points.append(['2', '0p40', '100'])
        mass_points.append(['2', '0p50', '0p1'])
        mass_points.append(['2', '0p50', '1p0'])
        mass_points.append(['2', '0p50', '10'])
        mass_points.append(['2', '0p50', '100'])
        mass_points.append(['2', '0p67', '0p1'])
        mass_points.append(['2', '0p67', '1p0'])
        mass_points.append(['2', '0p67', '10'])
        mass_points.append(['2', '0p67', '100'])
        mass_points.append(['2', '0p90', '0p1'])
        mass_points.append(['2', '0p90', '1p0'])
        mass_points.append(['2', '0p90', '10'])
        mass_points.append(['2', '0p90', '100'])
        mass_points.append(['4', '0p25', '0p1'])
        mass_points.append(['4', '0p25', '1p0'])
        mass_points.append(['4', '0p25', '10'])
        mass_points.append(['4', '0p25', '100'])
        mass_points.append(['4', '0p40', '0p1'])
        mass_points.append(['4', '0p40', '1p0'])
        mass_points.append(['4', '0p40', '10'])
        mass_points.append(['4', '0p40', '100'])
        mass_points.append(['4', '0p80', '0p1'])
        mass_points.append(['4', '0p80', '1p0'])
        mass_points.append(['4', '0p80', '10'])
        mass_points.append(['4', '0p80', '100'])
        mass_points.append(['4', '1p30', '0p1'])
        mass_points.append(['4', '1p30', '1p0'])
        mass_points.append(['4', '1p30', '10'])
        mass_points.append(['4', '1p30', '100'])
        mass_points.append(['4', '1p90', '0p1'])
        mass_points.append(['4', '1p90', '1p0'])
        mass_points.append(['4', '1p90', '10'])
        mass_points.append(['4', '1p90', '100'])
        mass_points.append(['5', '0p50', '0p1'])
        mass_points.append(['5', '0p50', '1p0'])
        mass_points.append(['5', '0p50', '10'])
        mass_points.append(['5', '0p50', '100'])
        mass_points.append(['5', '1p00', '0p1'])
        mass_points.append(['5', '1p00', '1p0'])
        mass_points.append(['5', '1p00', '10'])
        mass_points.append(['5', '1p00', '100'])
        mass_points.append(['5', '1p67', '0p1'])
        mass_points.append(['5', '1p67', '1p0'])
        mass_points.append(['5', '1p67', '10'])
        mass_points.append(['5', '1p67', '100'])
        mass_points.append(['5', '2p40', '0p1'])
        mass_points.append(['5', '2p40', '1p0'])
        mass_points.append(['5', '2p40', '10'])
        mass_points.append(['5', '2p40', '100'])
        mass_points.append(['10', '1p00', '0p1'])
        mass_points.append(['10', '1p00', '1p0'])
        mass_points.append(['10', '1p00', '10'])
        mass_points.append(['10', '1p00', '100'])
        mass_points.append(['10', '2p00', '0p1'])
        mass_points.append(['10', '2p00', '1p0'])
        mass_points.append(['10', '2p00', '10'])
        mass_points.append(['10', '2p00', '100'])
        mass_points.append(['10', '3p33', '0p1'])
        mass_points.append(['10', '3p33', '1p0'])
        mass_points.append(['10', '3p33', '10'])
        mass_points.append(['10', '3p33', '100'])
        mass_points.append(['10', '4p90', '0p1'])
        mass_points.append(['10', '4p90', '1p0'])
        mass_points.append(['10', '4p90', '10'])
        mass_points.append(['10', '4p90', '100'])
        for [mpi,mA,t] in mass_points:
            config_list.append(config)
            config_list[-1].JobType.pyCfgParams=["era={}".format(era),"data=False",]
            if era=="2022":
                dataset_name = '/scenarioA_mpi_{}_mA_{}_ctau_{}/jleonhol-AODSIM_2022_ext-bd8711905ed05b0226f084d42f06d7ac/USER'.format(mpi, mA, t)
                config_list[-1].Data.inputDataset = dataset_name
            elif era=="2022postEE":
                dataset_name = '/scenarioA_mpi_{}_mA_{}_ctau_{}/tafoyava-AODSIM_2022-6c55dbd4f99ed6c824b000a7a99348b1/USER'.format(mpi, mA, t)
                config_list[-1].Data.inputDataset = dataset_name
            # No 2023 eras for the moment, when adding them, you have to run . install_cmssw.sh 2023central first
            config_list[-1].General.requestName = 'centralSkim__{}_{}_mpi-{}_mA-{}_ctau-{}mm_{}'.format(signal, era, mpi, mA, t, ntuple_version)
            print(config)
            crabCommand('submit', config = config, dryrun = True) ## dryrun = True for local test
    elif "ScenarioB1" in sys.argv[2]:
        config.Data.outLFNDirBase = '/store/group/Run3Scouting/RAWScouting_privScenarioA_v'+ntuple_version # DB no
        config.Data.inputDBS = 'phys03'
        config.Data.splitting = 'FileBased'
        config.Data.publication = True
        config.Data.unitsPerJob = int(10) # Increased to match 10 jobs per file aprox
        config.Data.outputDatasetTag = "private-Skim_{era}-v2".format(era=era)
        mass_points.append(['1', '0p40', '0p1'])
        mass_points.append(['1', '0p40', '1p0'])
        mass_points.append(['1', '0p40', '10'])
        mass_points.append(['1', '0p40', '100'])
        mass_points.append(['2', '0p33', '0p1'])
        mass_points.append(['2', '0p33', '1p0'])
        mass_points.append(['2', '0p33', '10'])
        mass_points.append(['2', '0p33', '100'])
        mass_points.append(['2', '0p90', '0p1'])
        mass_points.append(['2', '0p90', '1p0'])
        mass_points.append(['2', '0p90', '10'])
        mass_points.append(['2', '0p90', '100'])
        mass_points.append(['4', '0p67', '0p1'])
        mass_points.append(['4', '0p67', '1p0'])
        mass_points.append(['4', '0p67', '10'])
        mass_points.append(['4', '0p67', '100'])
        mass_points.append(['4', '1p90', '0p1'])
        mass_points.append(['4', '1p90', '1p0'])
        mass_points.append(['4', '1p90', '10'])
        mass_points.append(['4', '1p90', '100'])
        mass_points.append(['5', '0p83', '0p1'])
        mass_points.append(['5', '0p83', '1p0'])
        mass_points.append(['5', '0p83', '10'])
        mass_points.append(['5', '0p83', '100'])
        mass_points.append(['5', '1p00', '0p1'])
        mass_points.append(['5', '1p00', '1p0'])
        mass_points.append(['5', '1p00', '10'])
        mass_points.append(['5', '1p00', '100'])
        mass_points.append(['5', '1p67', '0p1'])
        mass_points.append(['5', '1p67', '1p0'])
        mass_points.append(['5', '1p67', '10'])
        mass_points.append(['5', '1p67', '100'])
        mass_points.append(['5', '2p40', '0p1'])
        mass_points.append(['5', '2p40', '1p0'])
        mass_points.append(['5', '2p40', '10'])
        mass_points.append(['5', '2p40', '100'])
        for [mpi,mA,t] in mass_points:
            config_list.append(config)
            config_list[-1].JobType.pyCfgParams=["era={}".format(era),"data=False",]
            if era=="2022":
                os.system("""dasgoclient -query="/scenarioB1_mpi_{}_mA_{}_ctau_{}/jleonhol-AODSIM_2022_ext*bd8711905ed05b0226f084d42f06d7ac/USER instance=prod/phys03" > temp.txt""".format(mpi, mA, t))
                with open('temp.txt', 'r') as file:
                    dataset_name = file.readlines()[-1]
                #dataset_name = '/scenarioB1_mpi_{}_mA_{}_ctau_{}/jleonhol-AODSIM_2022_ext-bd8711905ed05b0226f084d42f06d7ac/USER'.format(mpi, mA, t)
                config_list[-1].Data.inputDataset = dataset_name
            elif era=="2022postEE":
                os.system("""dasgoclient -query="/scenarioB1_mpi_{}_mA_{}_ctau_{}/tafoyava-AODSIM_2022-6c55dbd4f99ed6c824b000a7a99348b1/USER instance=prod/phys03" > temp.txt""".format(mpi, mA, t))
                with open('temp.txt', 'r') as file:
                    dataset_name = file.readlines()[-1]
                #dataset_name = '/scenarioB1_mpi_{}_mA_{}_ctau_{}/tafoyava-AODSIM_2022-6c55dbd4f99ed6c824b000a7a99348b1/USER'.format(mpi, mA, t)
                config_list[-1].Data.inputDataset = dataset_name
            # No 2023 eras for the moment, when adding them, you have to run . install_cmssw.sh 2023central first
            config_list[-1].General.requestName = 'centralSkim__{}_{}_mpi-{}_mA-{}_ctau-{}mm_{}'.format(signal, era, mpi, mA, t, ntuple_version)
            print(config)
            #crabCommand('submit', config = config, dryrun = False) ## dryrun = True for local test

    elif "ScenarioB1priv" in sys.argv[2]:
        # This setup is provisional as it is tested with private signal crab produced samples
        #   -> Will be replaced by central datasets when done
        config.Data.outLFNDirBase = '/store/group/Run3Scouting/RAWScouting_privScenarioB1_v'+ntuple_version # DB no
        config.Data.inputDBS = 'phys03'
        config.Data.splitting = 'FileBased'
        config.Data.publication = True
        config.Data.unitsPerJob = int(100) # Increased to match 10 jobs per file aprox
        config.Data.outputDatasetTag = "private-Skim_{era}-v1".format(era=era)
        if era=="2022":
            mass_points.append(['ScenarioB1_mpi_4_mA_1p33_ctau_0p1', '/scenarioB1_mpi_4_mA_1p33_ctau_0p1/jleonhol-AODSIM_2022-bd8711905ed05b0226f084d42f06d7ac/USER'])
            mass_points.append(['ScenarioB1_mpi_4_mA_1p33_ctau_1p0', '/scenarioB1_mpi_4_mA_1p33_ctau_1p0/jleonhol-AODSIM_2022-bd8711905ed05b0226f084d42f06d7ac/USER'])
            mass_points.append(['ScenarioB1_mpi_4_mA_1p33_ctau_10', '/scenarioB1_mpi_4_mA_1p33_ctau_10/jleonhol-AODSIM_2022-bd8711905ed05b0226f084d42f06d7ac/USER'])
            mass_points.append(['ScenarioB1_mpi_4_mA_1p33_ctau_100', '/scenarioB1_mpi_4_mA_1p33_ctau_100/jleonhol-AODSIM_2022-bd8711905ed05b0226f084d42f06d7ac/USER'])
        elif era=="2022postEE":
            mass_points.append(['ScenarioB1_mpi_4_mA_1p33_ctau_0p1', '/scenarioB1_mpi_4_mA_1p33_ctau_0p1/jleonhol-AODSIM_2022-bd8711905ed05b0226f084d42f06d7ac/USER'])
            mass_points.append(['ScenarioB1_mpi_4_mA_1p33_ctau_1p0', '/scenarioB1_mpi_4_mA_1p33_ctau_1p0/jleonhol-AODSIM_2022-bd8711905ed05b0226f084d42f06d7ac/USER'])
            mass_points.append(['ScenarioB1_mpi_4_mA_1p33_ctau_10', '/scenarioB1_mpi_4_mA_1p33_ctau_10/jleonhol-AODSIM_2022-bd8711905ed05b0226f084d42f06d7ac/USER'])
            mass_points.append(['ScenarioB1_mpi_4_mA_1p33_ctau_100', '/scenarioB1_mpi_4_mA_1p33_ctau_100/jleonhol-AODSIM_2022-bd8711905ed05b0226f084d42f06d7ac/USER'])
        for [signal_name,dataset_name] in mass_points:
            config_list.append(config)
            config_list[-1].JobType.pyCfgParams=["era={}".format(era),"data=False",]
            config_list[-1].Data.inputDataset = dataset_name
            config_list[-1].General.requestName = 'centralSkim_{}_{}_{}'.format(signal_name, era, ntuple_version)
            print(config_list[-1])
            crabCommand('submit', config = config_list[-1], dryrun = False) ## dryrun = True for local test
            #print(config)
            #crabCommand('submit', config = config, dryrun = False) ## dryrun = True for local test
    elif "DQCD_bkg" in sys.argv[2]:
        if era != "2024":
            print("DQCD_bkg only supported for era=2024"); quit()
        config.Data.outLFNDirBase = '/store/group/Run3Scouting/RAWScouting_DQCD_bkg2024_v_dqcd'
        config.Data.inputDBS = 'global'
        config.Data.splitting = 'EventAwareLumiBased'
        config.Data.unitsPerJob = int(10e4)
        config.Data.publication = False
        inputfile = 'data/datasets_dqcd_2024_background.txt'
        with open(inputfile,'r') as f:
            dataset_list = f.readlines()
        for dataset_name in dataset_list:
            dataset_name = dataset_name.strip()
            if not dataset_name or dataset_name[0] == '#':
                continue
            config_list.append(config)
            config_list[-1].JobType.pyCfgParams = ["era={}".format(era), "data=False", "background=True"]
            config_list[-1].Data.inputDataset = dataset_name
            short_name = dataset_name.split('/')[1].split('_TuneCP5')[0]
            config_list[-1].General.requestName = 'centralSkim__{}_{}_{}'.format(short_name, era, ntuple_version)
            print(config)
            try:
                crabCommand('submit', config=config, dryrun=False)
            except (HTTPException, ClientException) as e:
                print('{} cant be launched! Skipping... ({})'.format(short_name, e))
    elif "minBias_2024" in sys.argv[2]:
        if era != "2024":
            print("minBias_2024 only supported for era=2024"); quit()
        ntuple_version_mb = "_minBias_2024"
        config.Data.outLFNDirBase = '/store/group/Run3Scouting/RAWScouting_DQCD_bkg2024_v_minBias'
        config.Data.inputDBS = 'global'
        config.Data.splitting = 'FileBased'
        config.Data.unitsPerJob = 10
        config.Data.publication = False
        inputfile = 'data/datasets_dqcd_2024_background.txt'
        with open(inputfile,'r') as f:
            dataset_list = f.readlines()
        for dataset_name in dataset_list:
            dataset_name = dataset_name.strip()
            if not dataset_name or dataset_name[0] == '#':
                continue
            config_list.append(config)
            config_list[-1].JobType.pyCfgParams = ["era={}".format(era), "data=False", "background=True"]
            config_list[-1].Data.inputDataset = dataset_name
            short_name = dataset_name.split('/')[1].split('_TuneCP5')[0]
            config_list[-1].General.requestName = 'centralSkim__{}_{}_{}'.format(short_name, era, ntuple_version_mb)
            print(config)
            try:
                crabCommand('submit', config=config, dryrun=False)
            except (HTTPException, ClientException) as e:
                print('{} cant be launched! Skipping... ({})'.format(short_name, e))
    #elif "[signal]" in sys.argv[2]: (<--- Add additional signals here)
    else:
        quit()
    
else:
    quit()



