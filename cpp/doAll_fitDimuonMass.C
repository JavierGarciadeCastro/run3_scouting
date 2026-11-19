{
  gROOT->ProcessLine(".L ./cpp/fit_dimuon.C+");  // Macro that performs the fitting
  gROOT->ProcessLine(".L ./cpp/helper.C+");  // Helper with handles 

  bool useData = true;
  bool useSignalMC = true;
  bool mergeEras = true;
  bool writeWS = true;
  bool reweighting = true;
  bool doUpAndDownVariations = true;
  if (!useSignalMC)
    doUpAndDownVariations = false;
TString period = "2024"; // Either 2022 or 2023
TString model = "ScenarioA"; // Either HTo2ZdTo2mu2x : ScenarioB1 : ScenarioA : BToPhi
  bool bdtRegions = true;

  // Dir with the RooDataSets
  //TString inDir = "/ceph/cms/store/user/fernance/Run3ScoutingOutput/outputHistograms_Jul-02-2024_2022_SRsOnly"; // last 2022
  //TString inDir = "/ceph/cms/store/user/fernance/Run3ScoutingOutput/outputHistograms_Jun-14-2024_SRsOnly_2023"; // last 2023
  //TString inDir = "/ceph/cms/store/user/garciaja/Run3ScoutingOutput/BToPhi_allCuts"; // BToPhi
  //TString inDir = "/ceph/cms/store/user/fernance/Run3ScoutingOutput/outputHistograms_Sep-25-2024_RooDatasets_unblind"; // last 2022 (unblinded)
  //TString inDir = "/ceph/cms/store/user/fernance/Run3ScoutingOutput/outputHistograms_Nov-13-2024_ctauReweighting";
  //TString inDir = "/ceph/cms/store/user/fernance/Run3ScoutingOutput/outputHistograms_Dec-03-2024_2022_complete"; //este
  //TString inDir = "/ceph/cms/store/user/fernance/Run3ScoutingOutput/outputHistograms_Dec-03-2024_2023_complete"; // y este
  //TString inDir = "/ceph/cms/store/user/fernance/Run3ScoutingOutput/outputHistograms_Dec-03-2024_2022_complete_reweighting/"; // Re-weighting
  //TString inDir = "/ceph/cms/store/user/fernance/Run3ScoutingOutput/outputHistograms_Dec-08-2024_2023_complete_reweighting//"; // Re-weighting
  //TString inDir = "/ceph/cms/store/user/fernance/Run3ScoutingOutput/outputHistograms_Dec-18-2024_2023_reweightingValidation/"; // Re-weighting validation 2023
  //TString inDir = "/ceph/cms/store/user/fernance/Run3ScoutingOutput/outputHistograms_Feb-20-2025_2023"; //este
  //TString inDir = "/ceph/cms/store/user/fernance/Run3ScoutingOutput/outputHistograms_Feb-20-2025_2022";
  //TString inDir = "/ceph/cms/store/user/fernance/Run3ScoutingOutput/outputHistograms_Feb-25-2025_2023";
  //TString inDir = "/ceph/cms/store/user/fernance/Run3ScoutingOutput/outputHistograms_Feb-25-2025_2022";
  //TString inDir = "/ceph/cms/store/user/fernance/Run3ScoutingOutput/outputHistograms_Feb-25-2025_allYears"; // This is for HAHM
  TString inDir = "fitInputs/2024_pnn/dqcd_mratio_WP_tight";
  //TString inDir = "/ceph/cms/store/user/fernance/Run3ScoutingOutput/outputHistograms_Jun-05-2025_allEras_forDQCDAnalysis"; // This is for DQCD
  //TString inDir = "/ceph/cms/store/user/fernance/Run3ScoutingOutput/outputHistograms_Sep-24-2025_allEras_forDQCDAnalysisAndBPhi/"; // This is for DQCD
  //TString inDir = "/ceph/cms/store/user/fernance/Run3ScoutingOutput/outputHistograms_Oct-01-2025_allEras_BToPhiAnalysis_RooOnly/"; // This is for BToPhi with extension
  //TString inDir = "/ceph/cms/store/user/fernance/Run3ScoutingOutput/outputHistograms_Oct-08-2025_allEras_BToPhiAnalysis_RooOnly_Reweighting/"; // This is for BToPhi reg
  //TString inDir = "/ceph/cms/store/user/fernance/Run3ScoutingOutput/outputHistograms_Oct-17-2025_allEras_forNEWDQCDAnalysis/"; // This is for DQCD

  // Names of the search regions: They should be picked accordingly to the model/fit it has to be done
  vector<TString> dNames = { };
  if (bdtRegions) {
  dNames.push_back("d_Dimuon_lxy0p0to0p2_bdt");
  dNames.push_back("d_Dimuon_lxy0p2to1p0_bdt");
  dNames.push_back("d_Dimuon_lxy1p0to2p4_bdt");
  dNames.push_back("d_Dimuon_lxy2p4to3p1_bdt");
  dNames.push_back("d_Dimuon_lxy3p1to7p0_bdt");
  dNames.push_back("d_Dimuon_lxy7p0to11p0_bdt");
  dNames.push_back("d_Dimuon_lxy11p0to16p0_bdt");
  dNames.push_back("d_Dimuon_lxy16p0to70p0_bdt");
  } else {
  //dNames.push_back("d_FourMu_sep");
  //dNames.push_back("d_Dimuon_lxy0p0to0p2_inclusive");
  //dNames.push_back("d_Dimuon_lxy0p2to1p0_inclusive");
  //dNames.push_back("d_Dimuon_lxy1p0to2p4_inclusive");
  //dNames.push_back("d_Dimuon_lxy2p4to3p1_inclusive");
  //dNames.push_back("d_Dimuon_lxy3p1to7p0_inclusive");
  //dNames.push_back("d_Dimuon_lxy7p0to11p0_inclusive");
  //dNames.push_back("d_Dimuon_lxy11p0to16p0_inclusive");
  //dNames.push_back("d_Dimuon_lxy16p0to70p0_inclusive");
  //
  dNames.push_back("d_FourMu_sep");
  dNames.push_back("d_FourMu_osv");
  dNames.push_back("d_Dimuon_lxy0p0to0p2_iso0_ptlow");
  dNames.push_back("d_Dimuon_lxy0p0to0p2_iso0_pthigh");
  dNames.push_back("d_Dimuon_lxy0p0to0p2_iso1_ptlow");
  dNames.push_back("d_Dimuon_lxy0p0to0p2_iso1_pthigh");
  dNames.push_back("d_Dimuon_lxy0p2to1p0_iso0_ptlow");
  dNames.push_back("d_Dimuon_lxy0p2to1p0_iso0_pthigh");
  dNames.push_back("d_Dimuon_lxy0p2to1p0_iso1_ptlow");
  dNames.push_back("d_Dimuon_lxy0p2to1p0_iso1_pthigh");
  dNames.push_back("d_Dimuon_lxy1p0to2p4_iso0_ptlow"); 
  dNames.push_back("d_Dimuon_lxy1p0to2p4_iso0_pthigh");
  dNames.push_back("d_Dimuon_lxy1p0to2p4_iso1_ptlow");
  dNames.push_back("d_Dimuon_lxy1p0to2p4_iso1_pthigh");
  dNames.push_back("d_Dimuon_lxy2p4to3p1_iso0_ptlow");
  dNames.push_back("d_Dimuon_lxy2p4to3p1_iso0_pthigh");
  dNames.push_back("d_Dimuon_lxy2p4to3p1_iso1_ptlow");
  dNames.push_back("d_Dimuon_lxy2p4to3p1_iso1_pthigh");
  dNames.push_back("d_Dimuon_lxy3p1to7p0_iso0_ptlow");
  dNames.push_back("d_Dimuon_lxy3p1to7p0_iso0_pthigh");
  dNames.push_back("d_Dimuon_lxy3p1to7p0_iso1_ptlow");
  dNames.push_back("d_Dimuon_lxy3p1to7p0_iso1_pthigh");
  dNames.push_back("d_Dimuon_lxy7p0to11p0_iso0_ptlow");
  dNames.push_back("d_Dimuon_lxy7p0to11p0_iso0_pthigh");
  dNames.push_back("d_Dimuon_lxy7p0to11p0_iso1_ptlow");
  dNames.push_back("d_Dimuon_lxy7p0to11p0_iso1_pthigh");
  dNames.push_back("d_Dimuon_lxy11p0to16p0_iso0_ptlow");
  dNames.push_back("d_Dimuon_lxy11p0to16p0_iso0_pthigh");
  dNames.push_back("d_Dimuon_lxy11p0to16p0_iso1_ptlow");
  dNames.push_back("d_Dimuon_lxy11p0to16p0_iso1_pthigh");
  dNames.push_back("d_Dimuon_lxy16p0to70p0_iso0_ptlow");
  dNames.push_back("d_Dimuon_lxy16p0to70p0_iso0_pthigh");
  dNames.push_back("d_Dimuon_lxy16p0to70p0_iso1_ptlow");
  dNames.push_back("d_Dimuon_lxy16p0to70p0_iso1_pthigh");
  dNames.push_back("d_Dimuon_lxy0p0to0p2_non-pointing"); 
  dNames.push_back("d_Dimuon_lxy0p2to1p0_non-pointing");
  dNames.push_back("d_Dimuon_lxy1p0to2p4_non-pointing");
  dNames.push_back("d_Dimuon_lxy2p4to3p1_non-pointing");
  dNames.push_back("d_Dimuon_lxy3p1to7p0_non-pointing");
  dNames.push_back("d_Dimuon_lxy7p0to11p0_non-pointing");
  dNames.push_back("d_Dimuon_lxy11p0to16p0_non-pointing");
  dNames.push_back("d_Dimuon_lxy16p0to70p0_non-pointing");
  //dNames.push_back("");
  }

  // Eras (to be uncommented when adding 2023 and splitting in eras)
  std::vector<std::pair<TString, TString>> periods; // first: year, second: era
  vector<TString> eras;
  if (period=="2022") {
    eras.push_back("2022");
    eras.push_back("2022postEE");
    periods.emplace_back("2022", "2022");
    periods.emplace_back("2022", "2022postEE");
  } else if (period=="2023") {
    eras.push_back("2022postEE");
    eras.push_back("2023");
    eras.push_back("2023BPix");
    periods.emplace_back("2023", "2022postEE");
    periods.emplace_back("2023", "2023");
    periods.emplace_back("2023", "2023BPix");
  } else if (period=="2024") {
    eras.push_back("2024");
    periods.emplace_back("2024", "2024");
  } else if (period=="allEras") {
    eras.push_back("2022");
    eras.push_back("2022postEE");
    eras.push_back("2023");
    eras.push_back("2023BPix");
    periods.emplace_back("2022", "2022");
    periods.emplace_back("2022", "2022postEE");
    periods.emplace_back("2023", "2022postEE");
    periods.emplace_back("2023", "2023");
    periods.emplace_back("2023", "2023BPix");
  }

  vector<TString> samples = { };
  vector<TString> sigmodels = { };
  vector<TString> sigsamples = { };
  vector<TString> sigsamples_outdir = { };
  vector<float> sigmasses_2mu = { };
  vector<float> sigmasses_4mu = { };
  vector<float> sigmasses_ctau = { };

  if ( useData ) {
    samples.push_back("Data");
  }

  // Signals (Should include sigMass, sigCtau and a proper definition for sigmasses_4mu)
  //vector<float> sigMass = {0.5, 0.7, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 14.0, 16.0, 20.0, 22.0, 24.0, 30.0, 34.0, 40.0, 44.0, 50.0};
  if ( model=="HTo2ZdTo2mu2x" ) {
    if ( useSignalMC ) {
        if (!reweighting) {
          vector<float> sigMass = {1.5, 2.0, 2.5, 5.0, 7.0, 8.0, 14.0, 16.0, 20.0, 24.0, 30.0, 34.0, 40.0, 50.0};
          vector<float> sigCtau = {1, 10, 100, 1000};
          for ( unsigned int m=0; m<sigMass.size(); m++ ) {
            TString massString = Form("%.1f",sigMass[m]); 
            massString.ReplaceAll(".", "p");
            for ( unsigned int t=0; t<sigCtau.size(); t++ ) {
              if ( (sigMass[m] < 1.0 && sigCtau[t] > 10) || (sigMass[m] < 1.9 && sigCtau[t] > 100) )
                continue;
              TString ctauString = Form("%.2f",sigCtau[t]); 
              sigsamples.push_back(Form("Signal_HTo2ZdTo2mu2x_MZd-%s_ctau-%smm",massString.Data(),ctauString.Data()));
              sigmasses_2mu.push_back(sigMass[m]);
              sigmasses_4mu.push_back(125.); // Mass of the higgs
              sigmasses_ctau.push_back(sigCtau[t]); // Lifetime
              std::cout << Form("Reading signal sample: Signal_HTo2ZdTo2mu2x_MZd-%s_ctau-%smm",massString.Data(),ctauString.Data()) << std::endl;
            }
          }
        } else { 
          vector<float> sigMass = {1.5, 2.0, 2.5, 5.0, 7.0, 8.0, 14.0, 16.0, 20.0, 22.0, 24.0, 30.0, 40.0, 50.0};
          vector<float> sigCtau = {0.1, 0.16, 0.25, 0.40, 0.63, 1.00, 1.60, 2.50, 4.00, 6.30, 10.00, 16.00, 25.00, 40.00, 63.00, 100.00, 160.00, 250.00, 400.00, 630.00, 1000.00};
          for ( unsigned int m=0; m<sigMass.size(); m++ ) {
            TString massString = Form("%.1f",sigMass[m]);
            massString.ReplaceAll(".", "p");
            for ( unsigned int t=0; t<sigCtau.size(); t++ ) {
              if ( (sigMass[m] < 1.0 && sigCtau[t] > 10) || (sigMass[m] < 1.9 && sigCtau[t] > 100) )
                continue;
                sigsamples.push_back(Form("Signal_HTo2ZdTo2mu2x_MZd-%s_ctau-%.2fmm",massString.Data(),sigCtau[t]));
                sigmasses_2mu.push_back(sigMass[m]);
                sigmasses_4mu.push_back(125.); // Mass of the higgs
                sigmasses_ctau.push_back(sigCtau[t]); // Lifetime
                std::cout << Form("Reading signal sample: Signal_HTo2ZdTo2mu2x_MZd-%s_ctau-%.2fmm",massString.Data(),sigCtau[t]) << std::endl;
            }
          }
          //// Only for validation (leave commented unless you want to use it!!):
          //sigsamples.push_back("Signal_HTo2ZdTo2mu2x_MZd-1p5_rectau-1.00mm");
          //sigsamples.push_back("Signal_HTo2ZdTo2mu2x_MZd-1p5_ctau-1.00mm");
          //sigsamples.push_back("Signal_HTo2ZdTo2mu2x_MZd-5p0_rectau-1.00mm");
          //sigsamples.push_back("Signal_HTo2ZdTo2mu2x_MZd-5p0_ctau-1.00mm");
          //sigsamples.push_back("Signal_HTo2ZdTo2mu2x_MZd-8p0_rectau-1.00mm");
          //sigsamples.push_back("Signal_HTo2ZdTo2mu2x_MZd-8p0_ctau-1.00mm");
          //sigsamples.push_back("Signal_HTo2ZdTo2mu2x_MZd-14p0_rectau-1.00mm");
          //sigsamples.push_back("Signal_HTo2ZdTo2mu2x_MZd-14p0_ctau-1.00mm");
          //sigsamples.push_back("Signal_HTo2ZdTo2mu2x_MZd-22p0_rectau-1.00mm");
          //sigsamples.push_back("Signal_HTo2ZdTo2mu2x_MZd-22p0_ctau-1.00mm");
          //sigsamples.push_back("Signal_HTo2ZdTo2mu2x_MZd-40p0_rectau-1.00mm");
          //sigsamples.push_back("Signal_HTo2ZdTo2mu2x_MZd-40p0_ctau-1.00mm");
          //sigmasses_2mu = {1.5, 1.5, 5.0, 5.0, 8.0, 8.0, 14.0, 14.0, 22.0, 22.0, 40.0, 40.0};
          //sigmasses_4mu = {125.0, 125.0, 125.0, 125.0, 125.0, 125.0, 125.0, 125.0, 125.0, 125.0, 125.0, 125.0};
          //sigmasses_ctau = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
        }
    } else {
        vector<float> sigCtau = {1, 10, 100, 1000};
        vector<float> sigMass;
        std::ifstream infile("data/sigmasses_HTo2ZdTo2mu2x_fine.txt"); // "data/sigmasses_HTo2ZdTo2mu2x_fine.txt"
        std::string line;
        while (std::getline(infile, line)) {
          try {
              float mass = std::stof(line); // Convertir string a float
              if (mass > 0.0) // cut for testing
                sigMass.push_back(mass);
          } catch (const std::invalid_argument& e) {
              std::cerr << "No valid line" << line << std::endl;
          }
        }
        for ( unsigned int m=0; m<sigMass.size(); m++ ) {
            TString massString = Form("%.3f",sigMass[m]); 
            for ( unsigned int t=0; t<sigCtau.size(); t++ ) {
              if ( (sigMass[m] < 1.0 && sigCtau[t] > 10) || (sigMass[m] < 2.00 && sigCtau[t] > 100) )
                continue;
              TString ctauString = Form("%.2f",sigCtau[t]); 
              sigsamples.push_back(Form("Signal_HTo2ZdTo2mu2x_MZd-%s_ctau-%smm",massString.Data(),ctauString.Data()));
              sigmasses_2mu.push_back(sigMass[m]);
              sigmasses_4mu.push_back(125.); // Mass of the higgs
              sigmasses_ctau.push_back(sigCtau[t]); // Lifetime
              std::cout << Form("Reading signal sample: Signal_HTo2ZdTo2mu2x_MZd-%s_ctau-%smm",massString.Data(),ctauString.Data()) << std::endl;
            }
        }

    }
  } else if (model=="ScenarioA" && period=="2024") {
    if ( useSignalMC ) {
        std::ifstream infile("data/ScenarioA_2024_points.txt");
        std::string line;
        while (std::getline(infile, line)) {
          float mpi, mA, ctau;
          std::istringstream iss(line);
          if (!(iss >> mpi >> mA >> ctau)) {
            std::cerr << "No valid line: " << line << std::endl;
            continue;
          }
          TString mpiString = Form("%g", mpi);  mpiString.ReplaceAll(".", "p");
          TString mAString = Form("%g", mA);    mAString.ReplaceAll(".", "p");
          TString ctauString = Form("%g", ctau); ctauString.ReplaceAll(".", "p");
          sigsamples.push_back(Form("Signal_ScenarioA_Mpi-%s_MA-%s_ctau-%smm",mpiString.Data(),mAString.Data(),ctauString.Data()));
          sigmasses_2mu.push_back(mA);
          sigmasses_4mu.push_back(mpi);
          sigmasses_ctau.push_back(ctau);
          std::cout << "Reading signal sample: " << sigsamples.back() << std::endl;
        }
    }
  } else if (model=="ScenarioA") {
    if ( useSignalMC ) {
        std::vector<std::pair<float, float>> sigMass;  // Pairs are std::make_pair(mPi, mA)
        sigMass.push_back(std::make_pair(2.0, 0.67));
        sigMass.push_back(std::make_pair(5.0, 1.67));
        sigMass.push_back(std::make_pair(7.5, 2.50));
        vector<float> sigCtau = {0.1, 0.25, 0.60, 1.00, 2.50, 6.00, 10.00, 25.00, 60.00, 100.00, 1000.00};
        for ( unsigned int m=0; m<sigMass.size(); m++ ) {
            auto massPair = sigMass.at(m); 
            TString massStringPi = Form("%.0f",massPair.first); 
            TString massStringA = Form("%.2f",massPair.second); 
            massStringA.ReplaceAll(".", "p");
            if (massStringPi=="3") massStringPi = "3p33";
            if (massStringPi=="8") massStringPi = "7p50";
            if (massStringA=="2p00") massStringA = "2";
            for ( unsigned int t=0; t<sigCtau.size(); t++ ) {
                TString ctauString = Form("%.2f",sigCtau[t]);
                sigsamples.push_back(Form("Signal_ScenarioA_Mpi-%s_MA-%s_ctau-%smm",massStringPi.Data(),massStringA.Data(),ctauString.Data()));
                sigmasses_2mu.push_back(massPair.second);
                sigmasses_4mu.push_back(massPair.first); // Mass of the mother particle
                sigmasses_ctau.push_back(sigCtau[t]);
                std::cout << Form("Reading signal sample: Signal_ScenarioA_Mpi-%s_MA-%s_ctau-%smm",massStringPi.Data(),massStringA.Data(),ctauString.Data()) << std::endl;
            }
        }
    }
  } else if (model=="ScenarioB1") {
    if ( useSignalMC ) {
        std::vector<std::pair<float, float>> sigMass;
        sigMass.push_back(std::make_pair(2.0, 0.67));
        sigMass.push_back(std::make_pair(5.0, 1.67));
        sigMass.push_back(std::make_pair(7.5, 2.50));
        vector<float> sigCtau = {0.1, 0.25, 0.60, 1.00, 2.50, 6.00, 10.00, 25.00, 60.00, 100.00};
        for ( unsigned int m=0; m<sigMass.size(); m++ ) {
            auto massPair = sigMass.at(m); 
            TString massStringPi = Form("%.0f",massPair.first); 
            TString massStringA = Form("%.2f",massPair.second); 
            massStringA.ReplaceAll(".", "p");
            if (massStringPi=="3") massStringPi = "3p33";
            if (massStringPi=="8") massStringPi = "7p50";
            if (massStringA=="2p00") massStringA = "2";
            for ( unsigned int t=0; t<sigCtau.size(); t++ ) {
                TString ctauString = Form("%.2f",sigCtau[t]);
                sigsamples.push_back(Form("Signal_ScenarioB1_Mpi-%s_MA-%s_ctau-%smm",massStringPi.Data(),massStringA.Data(),ctauString.Data()));
                sigmasses_2mu.push_back(massPair.second);
                sigmasses_4mu.push_back(massPair.first); // Mass of the mother particle
                sigmasses_ctau.push_back(sigCtau[t]);
                std::cout << Form("Reading signal sample: Signal_ScenarioB1_Mpi-%s_MA-%s_ctau-%smm",massStringPi.Data(),massStringA.Data(),ctauString.Data()) << std::endl;
            }
        }
    }
  } 
  else if ( model == "BToPhi") {
    if ( useSignalMC ) {
        //vector<float> sigMass = {0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.9, 1.25, 1.5, 2.0, 2.85, 3.35, 4.0, 5.0};
        //vector<float> sigMass = {0.25, 0.3, 0.4, 0.6, 0.7, 0.9, 1.25, 1.5, 2.85, 3.35};
        //vector<float> sigMass = {5.00};
        vector<float> sigMass = {1.25, 1.5, 2.0};
        vector<float> sigCtau = {0.1, 0.16, 0.25, 0.40, 0.63, 1.00, 1.60, 2.50, 4.00, 6.30, 10.00, 16.00, 25.00, 40.00, 63.00, 100.00};
        //vector<float> sigCtau = {0.1, 1.00, 10.00, 100.00};
        for (unsigned int m = 0; m < sigMass.size(); m++) {
            TString massString = Form("%.2f", sigMass[m]);
            massString.ReplaceAll(".", "p");     
            for (unsigned int t = 0; t < sigCtau.size(); t++) {
                TString ctauString = Form("%.2f",sigCtau[t]);
                sigsamples.push_back(Form("Signal_BToPhi_MPhi-%s_ctau-%smm", massString.Data(), ctauString.Data()));
	              sigmasses_2mu.push_back(sigMass[m]);
                sigmasses_4mu.push_back(100.); // For B-hadron model this doesn't make much sense
                sigmasses_ctau.push_back(sigCtau[t]);
                std::cout << Form("Reading signal sample: Signal_BToPhi_MPhi-%s_ctau-%smm", massString.Data(), ctauString.Data()) << std::endl;
            }
        }
    } else { // BToPhi: useSignalMC=False
        vector<float> sigMass;
        //vector<float> sigCtau = {0.1, 1, 10, 100};
        vector<float> sigCtau = {1, 10, 100};
        std::ifstream infile("data/BToPhi_limitgrid.txt");
        std::string line;
        while (std::getline(infile, line)) {
          try {
              float mass = std::stof(line); // Convertir string a float
              sigMass.push_back(mass);
          } catch (const std::invalid_argument& e) {
              std::cerr << "No valid line" << line << std::endl;
          }
        }
        for ( unsigned int m=0; m<sigMass.size(); m++ ) {
            TString massString = Form("%.3f",sigMass[m]); 
            for ( unsigned int t=0; t<sigCtau.size(); t++ ) {
                TString ctauString = Form("%.2f",sigCtau[t]);
                sigsamples.push_back(Form("Signal_BToPhi_MPhi-%s_ctau-%smm", massString.Data(), ctauString.Data()));
	              sigmasses_2mu.push_back(sigMass[m]);
                sigmasses_4mu.push_back(100.); // For B-hadron model this doesn't make much sense
                sigmasses_ctau.push_back(sigCtau[t]);
                std::cout << Form("Reading signal sample: Signal_BToPhi_MPhi-%s_ctau-%smm", massString.Data(), ctauString.Data()) << std::endl;
            }
        }
    }
  }

 // Loop over datasets
 TSystemDirectory dir(inDir, inDir);
 TList* files = dir.GetListOfFiles(); 
 if (!files) {
   std::cerr << "Error: List of files can't be listed" << inDir << std::endl;
   return;
 }
 vector<RooDataSet> mmumu_bkgs = {};
 vector<vector<RooDataSet>> mmumu_sigs {{}}; 
 vector<vector<RooDataSet>> mmumu_sigs_trg_up {{}}; 
 vector<vector<RooDataSet>> mmumu_sigs_trg_down {{}}; 
 vector<vector<RooDataSet>> mmumu_sigs_sel_up {{}}; 
 vector<vector<RooDataSet>> mmumu_sigs_sel_down {{}}; 
 cout << "Preparing to read datasets..." << endl;
 for ( unsigned int d=0; d<dNames.size(); d++ ) {
   // Loop over datasets
   for ( unsigned int iperiod=0; iperiod<periods.size(); iperiod++ ) {
     //std::cout << "Loading dataset: " << dNames[d] << std::endl;
     vector<TString> dataEras = {};
     TString era = periods[iperiod].second;
     TString year = periods[iperiod].first;
     //if (era.Contains("2022"))
     //  year = "2022";
     //else
     //  year = "2023";
     if (era=="2022" && year=="2022") {
       dataEras.push_back("_DataC_"); dataEras.push_back("_DataD_"); dataEras.push_back("_DataE_");
     } else if (era=="2022postEE" && year=="2022") {
       dataEras.push_back("_DataF_"); dataEras.push_back("_DataG_");
     } else if (era=="2022postEE" && year=="2023") {
       dataEras.push_back("_DataC-triggerV10_");
     } else if (era=="2023" && year=="2023") {
       dataEras.push_back("_DataC_");
     } else if (era=="2023BPix" && year=="2023") {
       dataEras.push_back("_DataD_");
     } else if (era=="2024" && year=="2024") {
       dataEras.push_back("_Data2024_");
     }
     // Loop over data files
     int idata = 0;
     for (const auto& file : *files) {
       TString filename = file->GetName();
       if (!filename.BeginsWith("histograms_Data") || filename.EndsWith("all.root"))
         continue;
       bool toRead = false;
       for ( unsigned int jera=0; jera<dataEras.size(); jera++ ) {
         if (filename.Contains(year) && filename.Contains(dataEras[jera])) {
           toRead = true; break;
         }
       }
       if (!toRead)
         continue;
       TString inFile = Form("%s/%s",inDir.Data(),filename.Data());
       TFile fin(inFile, "READ");
       std::cout << inFile << std::endl;
       if (idata == 0) {
         RooDataSet *tds = (RooDataSet*) fin.Get(dNames[d])->Clone();
         mmumu_bkgs.push_back( *tds );
         //std::cout << "Appended: " << filename.Data() << std::endl;
       } else {
         RooDataSet *tds_other = (RooDataSet*) fin.Get(dNames[d])->Clone();
         mmumu_bkgs[iperiod].append( *tds_other );
         //std::cout << "Appended: " << filename.Data() << std::endl;
       }
       fin.Close();
       idata++;
     }
     mmumu_bkgs[iperiod].SetName(dNames[d]+"_data_"+year+"_"+era);

     // Loop over signals 
     for ( int isample=0; isample<sigsamples.size(); isample++ ) {
       TString sample = sigsamples[isample];
       cout<<"Sample: "<< sample << endl;
       if ( useSignalMC ) {
         TString inFile = Form("%s/histograms_%s_%s_%s_0.root",inDir.Data(),sample.Data(),era.Data(),year.Data());
         TFile fin(inFile);
         RooDataSet *tds = (RooDataSet*) fin.Get(dNames[d])->Clone();
         tds->SetName(dNames[d]+"_"+sample+"_"+year+"_"+era);
         std::cout << "Reading signal file: " <<  inFile << ", with dataset with entries: " << tds->sumEntries() << std::endl;
         if (iperiod == 0) {
           vector<RooDataSet> tds_aux{};
           mmumu_sigs.push_back( tds_aux );
         }
         mmumu_sigs[isample].push_back( *tds );
         if (doUpAndDownVariations) {
           RooDataSet *tds_trg_up = (RooDataSet*) fin.Get(dNames[d]+"_trg_up")->Clone();
           RooDataSet *tds_trg_down = (RooDataSet*) fin.Get(dNames[d]+"_trg_down")->Clone();
           RooDataSet *tds_sel_up = (RooDataSet*) fin.Get(dNames[d]+"_sel_up")->Clone();
           RooDataSet *tds_sel_down = (RooDataSet*) fin.Get(dNames[d]+"_sel_down")->Clone();
           tds_trg_up->SetName(dNames[d]+"_"+sample+"_"+year+"_"+era+"_trg_up");
           tds_trg_down->SetName(dNames[d]+"_"+sample+"_"+year+"_"+era+"_trg_down");
           tds_sel_up->SetName(dNames[d]+"_"+sample+"_"+year+"_"+era+"_sel_up");
           tds_sel_down->SetName(dNames[d]+"_"+sample+"_"+year+"_"+era+"_sel_down");
           if (iperiod == 0) {
             vector<RooDataSet> tds_aux_trg_up{};
             vector<RooDataSet> tds_aux_trg_down{};
             vector<RooDataSet> tds_aux_sel_up{};
             vector<RooDataSet> tds_aux_sel_down{};
             mmumu_sigs_trg_up.push_back( tds_aux_trg_up );
             mmumu_sigs_trg_down.push_back( tds_aux_trg_down );
             mmumu_sigs_sel_up.push_back( tds_aux_sel_up );
             mmumu_sigs_sel_down.push_back( tds_aux_sel_down );
           }
           mmumu_sigs_trg_up[isample].push_back( *tds_trg_up );
           mmumu_sigs_trg_down[isample].push_back( *tds_trg_down );
           mmumu_sigs_sel_up[isample].push_back( *tds_sel_up );
           mmumu_sigs_sel_down[isample].push_back( *tds_sel_down );
         }
         fin.Close();
         cout << "Number of entries for the mmumu_sigs RooDataSet: " << mmumu_sigs[isample][iperiod].numEntries() << endl;
       } //else {
       //  RooDataSet *tds = (RooDataSet*) mmumu_bkgs[iera].emptyClone(dNames[d]+"_"+sample+"_"+year, "");
       //  mmumu_sigs.push_back( *tds );
       //}
     }
   }
   
   if (mergeEras) {
     TString outDir = "fitResults_"+period+"_"+model;
     RooDataSet mmumu_bkg_merged = mmumu_bkgs[0];
     for ( int iperiod=1; iperiod<periods.size(); iperiod++ ) {
         mmumu_bkg_merged.append(mmumu_bkgs[iperiod]);
     }
     mmumu_bkg_merged.SetName(dNames[d]+"_Data_"+period);
     cout << "Merged dataset for background has entries: " << mmumu_bkg_merged.sumEntries() << endl;
     vector<RooDataSet> mmumu_sig_merged = {};
     vector<RooDataSet> mmumu_sig_trg_up_merged = {};
     vector<RooDataSet> mmumu_sig_trg_down_merged = {};
     vector<RooDataSet> mmumu_sig_sel_up_merged = {};
     vector<RooDataSet> mmumu_sig_sel_down_merged = {};
     for (unsigned int isample=0; isample<sigsamples.size(); isample++ ) {
       // Create worksapce, import data and model
       std::cout << "Creating merged workspace" << std::endl;
       RooWorkspace wfit("wfit","workspace"); 
       
       RooWorkspace wfit_trg_up("wfit_trg_up","workspace_trg_up");
       RooWorkspace wfit_trg_down("wfit_trg_down","workspace_trg_down");
       RooWorkspace wfit_sel_up("wfit_sel_up","workspace_sel_up");
       RooWorkspace wfit_sel_down("wfit_sel_down","workspace_sel_down");
       
       //
       if (useSignalMC) {
         for (unsigned int iperiod=0; iperiod<periods.size(); iperiod++ ) {
           if (iperiod==0){
             mmumu_sig_merged.push_back(mmumu_sigs[isample][iperiod]);
           }else{
             mmumu_sig_merged[isample].append(mmumu_sigs[isample][iperiod]);
           }
         }
       } else {
         RooDataSet *tds = (RooDataSet*) mmumu_bkg_merged.emptyClone(dNames[d]+"_"+sigsamples[isample]+"_"+period, "");
         mmumu_sig_merged.push_back( *tds );
       }
       mmumu_sig_merged[isample].SetName(dNames[d]+"_"+sigsamples[isample]+"_"+period);
       if (doUpAndDownVariations) {
         for (unsigned int iperiod=0; iperiod<periods.size(); iperiod++ ) {
           if (iperiod==0) {
             mmumu_sig_trg_up_merged.push_back(mmumu_sigs_trg_up[isample][iperiod]);
             mmumu_sig_trg_down_merged.push_back(mmumu_sigs_trg_down[isample][iperiod]);
             mmumu_sig_sel_up_merged.push_back(mmumu_sigs_sel_up[isample][iperiod]);
             mmumu_sig_sel_down_merged.push_back(mmumu_sigs_sel_down[isample][iperiod]);
           } else {
             mmumu_sig_trg_up_merged[isample].append(mmumu_sigs_trg_up[isample][iperiod]);
             mmumu_sig_trg_down_merged[isample].append(mmumu_sigs_trg_down[isample][iperiod]);
             mmumu_sig_sel_up_merged[isample].append(mmumu_sigs_sel_up[isample][iperiod]);
             mmumu_sig_sel_down_merged[isample].append(mmumu_sigs_sel_down[isample][iperiod]);
           }
         }  
         mmumu_sig_trg_up_merged[isample].SetName(dNames[d]+"_"+sigsamples[isample]+"_"+period+"_trg_up");
         mmumu_sig_trg_down_merged[isample].SetName(dNames[d]+"_"+sigsamples[isample]+"_"+period+"_trg_down");
         mmumu_sig_sel_up_merged[isample].SetName(dNames[d]+"_"+sigsamples[isample]+"_"+period+"_sel_up");
         mmumu_sig_sel_down_merged[isample].SetName(dNames[d]+"_"+sigsamples[isample]+"_"+period+"_sel_down");
       }
       cout << "Merged dataset for signal " << sigsamples[isample] << " with entries " << mmumu_sig_merged[isample].sumEntries() << endl;
       // Fit invariant mass
       std::cout << "Prepare to fit..." << std::endl;
       if (dNames[d].BeginsWith("d_FourMu_")) {
          if (useSignalMC)   
            fitmass(mmumu_sig_merged[isample], "Signal", false, true, true, sigsamples[isample], sigmasses_2mu[isample], sigmasses_4mu[isample], sigmasses_ctau[isample], wfit, true, period, "dcbfastg", outDir);
          else
            fitmass(mmumu_sig_merged[isample], "Signal", false, true, false, sigsamples[isample], sigmasses_2mu[isample], sigmasses_4mu[isample], sigmasses_ctau[isample], wfit, true, period, "dcbfastg", outDir);
          fitmass(mmumu_bkg_merged, "Background", true, false, false, sigsamples[isample], sigmasses_2mu[isample], sigmasses_4mu[isample], sigmasses_ctau[isample], wfit, true, period, "", outDir); 
       } else {
        if (useSignalMC)
          fitmass(mmumu_sig_merged[isample], "Signal", false, true, true, sigsamples[isample], sigmasses_2mu[isample], sigmasses_2mu[isample], sigmasses_ctau[isample], wfit, false, period, "dcbfastg", outDir);
        else
          fitmass(mmumu_sig_merged[isample], "Signal", false, true, false, sigsamples[isample], sigmasses_2mu[isample], sigmasses_2mu[isample], sigmasses_ctau[isample], wfit, false, period, "dcbfastg", outDir);
        fitmass(mmumu_bkg_merged, "Background", true, false, false, sigsamples[isample], sigmasses_2mu[isample], sigmasses_2mu[isample], sigmasses_ctau[isample], wfit, false, period, "", outDir); 
       }
       if (doUpAndDownVariations) { 
         if (dNames[d].BeginsWith("d_FourMu_")) {
           fitmass(mmumu_sig_trg_up_merged[isample], "Signal", false, true, true, sigsamples[isample], sigmasses_2mu[isample], sigmasses_4mu[isample], sigmasses_ctau[isample], wfit_trg_up, true, period, "dcbfastg", outDir);
           fitmass(mmumu_sig_trg_down_merged[isample], "Signal", false, true, true, sigsamples[isample], sigmasses_2mu[isample], sigmasses_4mu[isample], sigmasses_ctau[isample], wfit_trg_down, true, period, "dcbfastg", outDir);
           fitmass(mmumu_sig_sel_up_merged[isample], "Signal", false, true, true, sigsamples[isample], sigmasses_2mu[isample], sigmasses_4mu[isample], sigmasses_ctau[isample], wfit_sel_up, true, period, "dcbfastg", outDir);
           fitmass(mmumu_sig_sel_down_merged[isample], "Signal", false, true, true, sigsamples[isample], sigmasses_2mu[isample], sigmasses_4mu[isample], sigmasses_ctau[isample], wfit_sel_down, true, period, "dcbfastg", outDir);
         } else {
           fitmass(mmumu_sig_trg_up_merged[isample], "Signal", false, true, true, sigsamples[isample], sigmasses_2mu[isample], sigmasses_2mu[isample], sigmasses_ctau[isample], wfit_trg_up, false, period, "dcbfastg", outDir);
           fitmass(mmumu_sig_trg_down_merged[isample], "Signal", false, true, true, sigsamples[isample], sigmasses_2mu[isample], sigmasses_2mu[isample], sigmasses_ctau[isample], wfit_trg_down, false, period, "dcbfastg", outDir);
           fitmass(mmumu_sig_sel_up_merged[isample], "Signal", false, true, true, sigsamples[isample], sigmasses_2mu[isample], sigmasses_2mu[isample], sigmasses_ctau[isample], wfit_sel_up, false, period, "dcbfastg", outDir);
           fitmass(mmumu_sig_sel_down_merged[isample], "Signal", false, true, true, sigsamples[isample], sigmasses_2mu[isample], sigmasses_2mu[isample], sigmasses_ctau[isample], wfit_sel_down, false, period, "dcbfastg", outDir);
         }
       }
       
    
       // Print workspace contents
       std::cout << "Workspace contents: " << std::endl;
       wfit.Print();

       // Save the workspace into a ROOT file
       TString fwsname;
       if ( writeWS ) {
	 fwsname = Form("%s/%s_workspace.root",outDir.Data(),mmumu_sig_merged[isample].GetName());
         TFile *fws = new TFile(fwsname, "RECREATE");
         fws->cd();
         cout << "Writing workspace..." << endl;
         wfit.Write();
         if (doUpAndDownVariations) { 
           wfit_trg_up.Write();
           wfit_trg_down.Write();
           wfit_sel_up.Write();
           wfit_sel_down.Write();
         }
         fws->Close();
       }
       cout<<endl;
     }
   }
   std::cout << "Cleaning the workspace containers..." << std::endl;
   for ( int isample=0; isample<sigsamples.size(); isample++ ) {
     if (useSignalMC)
       mmumu_sigs[isample].clear();
     if (doUpAndDownVariations) { 
       mmumu_sigs_trg_up[isample].clear();
       mmumu_sigs_trg_down[isample].clear();
       mmumu_sigs_sel_up[isample].clear();
       mmumu_sigs_sel_down[isample].clear();
     }
   }
   if (useSignalMC)
     mmumu_sigs.clear();
   if (doUpAndDownVariations) { 
     mmumu_sigs_trg_up.clear();
     mmumu_sigs_trg_down.clear();
     mmumu_sigs_sel_up.clear();
     mmumu_sigs_sel_down.clear();
   }
   mmumu_bkgs.clear();
   
 }
}
