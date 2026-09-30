#include <cstdlib>
#include <iostream>
#include <fstream>
#include "run3ScoutingLooper.C"

unsigned char2unsigned(const char *c) {
  char p = *c;
  unsigned res = 0;
  while (p) {
      res = res*10 + (p - '0');
      c++;
      p = *c;
  }
  return res;
}
int char2int(const char *c) {
  return (*c == '-') ? -char2unsigned(c+1) : char2unsigned(c);
}


std::vector<TString> getFiles(const std::string inputDir, const int startFile, const int nFiles, const bool isCondor, const bool fromCrab) {
  std::vector<TString> files;
  std::string fullInputDir;
  unsigned int iFile=0;
  if (fromCrab) {
    std::cout << "Trying to identify the files with dasgoclient" << std::endl;
    std::string command;
    command = "/cvmfs/cms.cern.ch/common/dasgoclient --query=\"file dataset=";
    command += inputDir; 
    command += " instance=prod/phys03\" -format string > infiles.txt";
    std::system(command.c_str());
    std::ifstream infiles("infiles.txt");
    std::string line;
    while(getline(infiles, line)) {
      if (iFile<startFile) {
	iFile++;
	continue;
      }
      if (!TString(line.c_str()).Contains(".root"))
	continue;
      else {
	std::cout << "...adding file: " << line << std::endl;
	files.push_back(TString("root://cmsxrootd.fnal.gov//"+line));
      }
      iFile++;
      if (iFile == startFile+nFiles)
	break;
    }
  } else if (!isCondor) {
    if (inputDir.rfind("/ceph/cms", 0) == 0)
      fullInputDir = inputDir;
    else
      fullInputDir = "/ceph/cms"+inputDir;
    const fs::path dir{fullInputDir};
  
    fs::recursive_directory_iterator it{dir}, end;
    for (; it != end; ++it) {
      if (it->is_directory()) {
	std::string dname = it->path().filename().string();
	if (dname=="log" || dname=="failed" || dname=="temp")
	  it.disable_recursion_pending();
	continue;
      }
      if (!TString(it->path().c_str()).Contains(".root"))
	continue;
      if (iFile<startFile) {
	iFile++;
	continue;
      }
      files.push_back(TString(it->path().c_str()));
      iFile++;
      if (iFile == startFile+nFiles)
	break;
    }
  }
  else {
    std::string temp_str(inputDir);
    if (inputDir.rfind("/ceph/cms", 0) == 0)
      fullInputDir = temp_str.replace(temp_str.find("/ceph/cms"),sizeof("/ceph/cms")-1,"");
    else
      fullInputDir = inputDir;
    std::string command;
    command = "python3 recursive_xrdls.py ";
    command += fullInputDir;
    command += " > infiles.txt";
    std::system(command.c_str());
    std::ifstream infiles("infiles.txt");
    std::string line;
    while(getline(infiles, line)) {
      if (iFile<startFile) {
	iFile++;
	continue;
      }
      if (!TString(line.c_str()).Contains(".root"))
	continue;
      else {
	files.push_back(TString("root://redirector.t2.ucsd.edu:1095/"+line));
      }
      iFile++;
      if (iFile == startFile+nFiles)
	break;
    }
  }
  std::system("rm -f infiles.txt");
  return files;
}

int main(int argc, char **argv) {
  // Arguments
  const char* outdir      = ( argc > 1  ? argv[1]            : "temp_data" );
  TString year            = ( argc > 2  ? argv[2]            : "2022" );
  TString sampleArg       = ( argc > 3  ? argv[3]            : "DataF" );
  int startFile           = ( argc > 4  ? char2int(argv[4])  : 0 );
  int nFiles              = ( argc > 5  ? char2int(argv[5])  : 1000000 ); // Large number as default to run over all files
  bool isCondor           = ( argc > 6  ? char2int(argv[6])  : 0 );
  bool fromCrab           = ( argc > 7  ? char2int(argv[7])  : 0 );
  std::vector<TString> files;
  TString process;
  // Sample list: Data
  if ( sampleArg=="DataB" && year=="2022") {
    files = getFiles("/ceph/cms/store/group/Run3Scouting/Run3ScoutingSamples/Nov-13-2023/Data/2022B/", startFile, nFiles, isCondor, fromCrab);
    process = "DataB";
  }
  if ( sampleArg=="DataC" && year=="2022") {
    files = getFiles("/ceph/cms/store/group/Run3Scouting/Run3ScoutingSamples/Nov-13-2023/Data/2022C/", startFile, nFiles, isCondor, fromCrab);
    process = "DataC";
  }
  if ( sampleArg=="DataD" && year=="2022") {
    files = getFiles("/ceph/cms/store/group/Run3Scouting/Run3ScoutingSamples/Nov-13-2023/Data/2022D/", startFile, nFiles, isCondor, fromCrab);
    process = "DataD";
  }
  if ( sampleArg=="DataE" && year=="2022") {
    files = getFiles("/ceph/cms/store/group/Run3Scouting/Run3ScoutingSamples/Nov-13-2023/Data/2022E/", startFile, nFiles, isCondor, fromCrab);
    process = "DataE";
  }
  if ( sampleArg=="DataF" && year=="2022") {
    files = getFiles("/ceph/cms/store/group/Run3Scouting/Run3ScoutingSamples/Nov-13-2023/Data/2022F/", startFile, nFiles, isCondor, fromCrab);
    process = "DataF";
  }
  if ( sampleArg=="DataG" && year=="2022") {
    files = getFiles("/ceph/cms/store/group/Run3Scouting/Run3ScoutingSamples/Nov-13-2023/Data/2022G/", startFile, nFiles, isCondor, fromCrab);
    process = "DataG";
  }
  if ( sampleArg=="DataB" && year=="2023") { 
    files = getFiles("/ceph/cms/store/group/Run3Scouting/Run3ScoutingSamples/May-24-2024/Data/2023B/", startFile, nFiles, isCondor, fromCrab);
    process = "DataB";
  }
  if ( sampleArg=="DataC-triggerV10" && year=="2023") {
    files = getFiles("/ceph/cms/store/group/Run3Scouting/Run3ScoutingSamples/May-24-2024/Data/2023C-triggerV10/", startFile, nFiles, isCondor, fromCrab);
    process = "DataC-triggerV10";
  }
  if ( sampleArg=="DataC" && year=="2023") {
    files = getFiles("/ceph/cms/store/group/Run3Scouting/Run3ScoutingSamples/May-24-2024/Data/2023C/", startFile, nFiles, isCondor, fromCrab);
    process = "DataC";
  }
  if ( sampleArg=="DataD" && year=="2023") {
    files = getFiles("/ceph/cms/store/group/Run3Scouting/Run3ScoutingSamples/May-24-2024/Data/2023D/", startFile, nFiles, isCondor, fromCrab);
    process = "DataD";
  }
  //
  // Sample list: Monte Carlo
  if ( sampleArg=="DileptonMinBias" && year=="2022") {
    files = getFiles("/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/ppradeep-crab_skim4__2022X_InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8_5_syst-cf0b66eebbd8a1efdeefe3d3090bd687/USER", startFile, nFiles, isCondor, true);
    process = "DileptonMinBias";
  }
  //
  // Sample list: Central Signal
  std::string line,sampleName,samplePath;
  std::string delimiter = ",";
  std::vector<std::string> centralLists = {"input/centralDatasets.txt", "input/centralDatasets_"+std::string(year.Data())+".txt"};
  for (auto const& centralList : centralLists) {
    std::ifstream centralFiles(centralList);
    while(getline(centralFiles, line)) {
      std::string s = line;
      sampleName = s.substr(0, s.find(delimiter));
      s.erase(0, s.find(delimiter) + delimiter.length());
      samplePath = s;
      if (sampleName==sampleArg) {
        files = getFiles(samplePath, startFile, nFiles, isCondor, fromCrab);
        process = sampleArg;
        break;
      }
    }
  }
  //
  // Sample list: Signal
  if ( sampleArg=="Signal_HTo2ZdTo2mu2x_MZd-2p0_ctau-1mm" && year=="2022") {
    files = getFiles("/ceph/cms/store/group/Run3Scouting/Run3ScoutingSamples/Nov-13-2023/Signal/HTo2ZdTo2mu2x_MZd-2p0_ctau-1mm_2022/", startFile, nFiles, isCondor, fromCrab); 
    process = "Signal_HTo2ZdTo2mu2x_MZd-2p0_ctau-1mm";
  }
  if ( sampleArg=="Signal_HTo2ZdTo2mu2x_MZd-2p0_ctau-10mm" && year=="2022") {
    files = getFiles("/ceph/cms/store/group/Run3Scouting/Run3ScoutingSamples/Nov-13-2023/Signal/HTo2ZdTo2mu2x_MZd-2p0_ctau-10mm_2022/", startFile, nFiles, isCondor, fromCrab); 
    process = "Signal_HTo2ZdTo2mu2x_MZd-2p0_ctau-10mm";
  }
  if ( sampleArg=="Signal_HTo2ZdTo2mu2x_MZd-2p0_ctau-100mm" && year=="2022") {
    files = getFiles("/ceph/cms/store/group/Run3Scouting/Run3ScoutingSamples/Nov-13-2023/Signal/HTo2ZdTo2mu2x_MZd-2p0_ctau-100mm_2022/", startFile, nFiles, isCondor, fromCrab); 
    process = "Signal_HTo2ZdTo2mu2x_MZd-2p0_ctau-100mm";
  }
  if ( sampleArg=="Signal_HTo2ZdTo2mu2x_MZd-7p0_ctau-1mm" && year=="2022") {
    files = getFiles("/ceph/cms/store/group/Run3Scouting/Run3ScoutingSamples/Nov-13-2023/Signal/HTo2ZdTo2mu2x_MZd-7p0_ctau-1mm_2022/", startFile, nFiles, isCondor, fromCrab); 
    process = "Signal_HTo2ZdTo2mu2x_MZd-7p0_ctau-1mm";
  }
  if ( sampleArg=="Signal_HTo2ZdTo2mu2x_MZd-7p0_ctau-10mm" && year=="2022") {
    files = getFiles("/ceph/cms/store/group/Run3Scouting/Run3ScoutingSamples/Nov-13-2023/Signal/HTo2ZdTo2mu2x_MZd-7p0_ctau-10mm_2022/", startFile, nFiles, isCondor, fromCrab); 
    process = "Signal_HTo2ZdTo2mu2x_MZd-7p0_ctau-10mm";
  }
  if ( sampleArg=="Signal_HTo2ZdTo2mu2x_MZd-7p0_ctau-100mm" && year=="2022") {
    files = getFiles("/ceph/cms/store/group/Run3Scouting/Run3ScoutingSamples/Nov-13-2023/Signal/HTo2ZdTo2mu2x_MZd-7p0_ctau-100mm_2022/", startFile, nFiles, isCondor, fromCrab); 
    process = "Signal_HTo2ZdTo2mu2x_MZd-7p0_ctau-100mm";
  }
  if ( sampleArg=="Signal_HTo2ZdTo2mu2x_MZd-50_Epsilon-2e-07_testL1" && year=="2022") {
    files = getFiles("/ceph/cms/store/group/Run3Scouting/Run3ScoutingSamples/Nov-13-2023/Signal/HTo2ZdTo2mu2x_MZd-50_Epsilon-2e-07_testL1/", startFile, nFiles, isCondor, fromCrab); 
    process = "Signal_HTo2ZdTo2mu2x_MZd-50_Epsilon-2e-07_testL1";
  }
  if ( sampleArg=="Signal_HTo2ZdTo2mu2x_MZd-50_Epsilon-6e-08_testL1" && year=="2022") {
    files = getFiles("/ceph/cms/store/group/Run3Scouting/Run3ScoutingSamples/Nov-13-2023/Signal/HTo2ZdTo2mu2x_MZd-50_Epsilon-6e-08_testL1/", startFile, nFiles, isCondor, fromCrab); 
    process = "Signal_HTo2ZdTo2mu2x_MZd-50_Epsilon-6e-08_testL1";
  }
  if ( sampleArg=="Signal_HTo2ZdTo2mu2x_MZd-50_Epsilon-1e-08_testL1" && year=="2022") {
    files = getFiles("/ceph/cms/store/group/Run3Scouting/Run3ScoutingSamples/Nov-13-2023/Signal/HTo2ZdTo2mu2x_MZd-50_Epsilon-1e-08_testL1/", startFile, nFiles, isCondor, fromCrab); 
    process = "Signal_HTo2ZdTo2mu2x_MZd-50_Epsilon-1e-08_testL1";
  }
  if ( sampleArg=="Signal_HTo2ZdTo2mu2x_MZd-50_Epsilon-4e-09_testL1" && year=="2022") {
    files = getFiles("/ceph/cms/store/group/Run3Scouting/Run3ScoutingSamples/Nov-13-2023/Signal/HTo2ZdTo2mu2x_MZd-50_Epsilon-4e-09_testL1/", startFile, nFiles, isCondor, fromCrab); 
    process = "Signal_HTo2ZdTo2mu2x_MZd-50_Epsilon-4e-09_testL1";
  }
  if ( sampleArg=="Signal_HTo2ZdTo2mu2x_MZd-60_Epsilon-1e-07_testL1" && year=="2022") {
    files = getFiles("/ceph/cms/store/group/Run3Scouting/Run3ScoutingSamples/Nov-13-2023/Signal/HTo2ZdTo2mu2x_MZd-60_Epsilon-1e-07_testL1/", startFile, nFiles, isCondor, fromCrab); 
    process = "Signal_HTo2ZdTo2mu2x_MZd-60_Epsilon-1e-07_testL1";
  }
  if ( sampleArg=="Signal_HTo2ZdTo2mu2x_MZd-60_Epsilon-4e-08_testL1" && year=="2022") {
    files = getFiles("/ceph/cms/store/group/Run3Scouting/Run3ScoutingSamples/Nov-13-2023/Signal/HTo2ZdTo2mu2x_MZd-60_Epsilon-4e-08_testL1/", startFile, nFiles, isCondor, fromCrab); 
    process = "Signal_HTo2ZdTo2mu2x_MZd-60_Epsilon-4e-08_testL1";
  }
  if ( sampleArg=="Signal_HTo2ZdTo2mu2x_MZd-60_Epsilon-7e-09_testL1" && year=="2022") {
    files = getFiles("/ceph/cms/store/group/Run3Scouting/Run3ScoutingSamples/Nov-13-2023/Signal/HTo2ZdTo2mu2x_MZd-60_Epsilon-7e-09_testL1/", startFile, nFiles, isCondor, fromCrab); 
    process = "Signal_HTo2ZdTo2mu2x_MZd-60_Epsilon-7e-09_testL1";
  }
  if ( sampleArg=="Signal_HTo2ZdTo2mu2x_MZd-60_Epsilon-2e-09_testL1" && year=="2022") {
    files = getFiles("/ceph/cms/store/group/Run3Scouting/Run3ScoutingSamples/Nov-13-2023/Signal/HTo2ZdTo2mu2x_MZd-60_Epsilon-2e-09_testL1/", startFile, nFiles, isCondor, fromCrab); 
    process = "Signal_HTo2ZdTo2mu2x_MZd-60_Epsilon-2e-09_testL1";
  }
  if ( sampleArg=="Signal_ScenB1_30_9p9_4p8_ctau_1mm" && year=="2022") {
    files = getFiles("/ceph/cms/store/group/Run3Scouting/Run3ScoutingSamples/Nov-13-2023/Signal/ScenB1_30_9p9_4p8_ctau_1mm_2022/", startFile, nFiles, isCondor, fromCrab); 
    process = "Signal_ScenB1_30_9p9_4p8_ctau_1mm";
  }
  if ( sampleArg=="Signal_ScenB1_30_9p9_4p8_ctau_10mm" && year=="2022") {
    files = getFiles("/ceph/cms/store/group/Run3Scouting/Run3ScoutingSamples/Nov-13-2023/Signal/ScenB1_30_9p9_4p8_ctau_10mm_2022/", startFile, nFiles, isCondor, fromCrab); 
    process = "Signal_ScenB1_30_9p9_4p8_ctau_10mm";
  }
  if ( sampleArg=="Signal_ScenB1_30_9p9_4p8_ctau_100mm" && year=="2022") {
    files = getFiles("/ceph/cms/store/group/Run3Scouting/Run3ScoutingSamples/Nov-13-2023/Signal/ScenB1_30_9p9_4p8_ctau_100mm_2022/", startFile, nFiles, isCondor, fromCrab); 
    process = "Signal_ScenB1_30_9p9_4p8_ctau_100mm";
  }
  //
  // Sample list: PF Monitor 2022
  if ( sampleArg=="MonDataTest" && year=="2022") {
    files = getFiles("/ceph/cms/store/group/Run3Scouting/Run3ScoutingSamples/Jan-9-2024/Data/Mon2022Test/", startFile, nFiles, isCondor, fromCrab);  // 1 files
    process = "MonDataTest";
  }
  if ( sampleArg=="MonDataB" && year=="2022") {
    files = getFiles("/ceph/cms/store/group/Run3Scouting/Run3ScoutingSamples/Jan-9-2024/Data/Mon2022B/", startFile, nFiles, isCondor, fromCrab);  // 5 files
    process = "MonDataB";
  }
  if ( sampleArg=="MonDataC" && year=="2022") {
    files = getFiles("/ceph/cms/store/group/Run3Scouting/Run3ScoutingSamples/Jan-9-2024/Data/Mon2022C/", startFile, nFiles, isCondor, fromCrab);  // 24 files
    process = "MonDataC";
  }
  if ( sampleArg=="MonDataD" && year=="2022") {
    files = getFiles("/ceph/cms/store/group/Run3Scouting/Run3ScoutingSamples/Jan-9-2024/Data/Mon2022D/", startFile, nFiles, isCondor, fromCrab);  // 5 files
    process = "MonDataD";
  }
  if ( sampleArg=="MonDataE" && year=="2022") {
    files = getFiles("/ceph/cms/store/group/Run3Scouting/Run3ScoutingSamples/Jan-9-2024/Data/Mon2022E/", startFile, nFiles, isCondor, fromCrab);  // 6 files
    process = "MonDataE";
  }
  if ( sampleArg=="MonDataF" && year=="2022") {
    files = getFiles("/ceph/cms/store/group/Run3Scouting/Run3ScoutingSamples/Jan-9-2024/Data/Mon2022F/", startFile, nFiles, isCondor, fromCrab);  // 13 files
    process = "MonDataF";
  }
  if ( sampleArg=="MonDataG" && year=="2022") {
    files = getFiles("/ceph/cms/store/group/Run3Scouting/Run3ScoutingSamples/Jan-9-2024/Data/Mon2022G/", startFile, nFiles, isCondor, fromCrab);  // 5 files
    process = "MonDataG";
  }
  if (sampleArg == "LocalTestSignal" && year == "2024") {
    files.push_back("/home/users/garciaja/fullRun3/CMSSW_15_0_2/src/run3_scouting/batch/skimmer_test_signal.root");
    process = "LocalTestSignal";
  }
  if (sampleArg == "LocalTestBkg" && year == "2024") {
    files.push_back("/home/users/garciaja/fullRun3/CMSSW_15_0_2/src/run3_scouting/batch/skimmer_test_bkg.root");
    process = "LocalTestBkg";
  }
  std::cout << "################################## \n";
  std::cout << "Number of files to process: " << files.size() << "\n";
  std::cout << "################################## \n";
  if (files.empty()) {
    std::cerr << "No files found for sample " << sampleArg << " year " << year << ". Exiting.\n";
    return 1;
  }
  run3ScoutingLooper(files, year, process, outdir, "_"+std::to_string(startFile)+"To"+std::to_string(startFile+nFiles-1));

  return 0;
}
