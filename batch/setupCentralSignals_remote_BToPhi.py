import os

## This script put the crab output in the corresponding ceph area and creates the files to run the looper and plotter

# To be modified by user:
model_template = "/ceph/cms/store/group/Run3Scouting/RAWScouting_BToPhi_2022_v5p0/BToPhi-{MASS}_ctau-{CTAU}mm_TuneCP5_13p6TeV_pythia8/crab_centralSkim__BToPhi_2022_m-{MASS}_ctau-{CTAU}mm_5p0/"
#model_template = "/ceph/cms/store/group/Run3Scouting/RAWScouting_BToPhi_2022_v5p0/BToPhi-{MASS}_ctau-{CTAU}mm_TuneCP5_13p6TeV_pythia8/crab_centralSkim__BToPhi_2022postEE_m-{MASS}_ctau-{CTAU}mm_5p0/"
signal_template = "Signal_BToPhi-{MASS}_ctau-{CTAU}mm"
model='BToPhi'

mass_points = []
if model=='BToPhi':
    
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
### File with central datasets fr looper input
fout_ = open("centralDatasets_copy.txt", "w")

### File with condor launcher for looper
f2out_ = open("runScoutingLooper_onCondor_BToPhi.sub", "w") 
f2out_.write("executable      = $ENV(STARTDIR)/condor/condorLooper_executable.sh\n")
f2out_.write("output          = $ENV(STARTDIR)/condor/plotting_logs/job.$(ClusterId).$(ProcId).out\n")
f2out_.write("error           = $ENV(STARTDIR)/condor/plotting_logs/job.$(ClusterId).$(ProcId).err\n")
f2out_.write("log             = $ENV(STARTDIR)/condor/plotting_logs/job.$(ClusterId).$(ProcId).log\n")
f2out_.write("\n")
f2out_.write("getenv = True\n")
f2out_.write("+JobFlavour = 'workday'\n")
f2out_.write("+DESIRED_Sites = 'T2_US_UCSD'\n")
f2out_.write("transfer_input_files = $ENV(STARTDIR)/package.tar.gz\n")
f2out_.write("should_transfer_files = YES\n")
f2out_.write("when_to_transfer_output = ON_EXIT\n")
f2out_.write("x509userproxy=$ENV(X509_USER_PROXY)\n")
f2out_.write("use_x509userproxy = True\n")
f2out_.write("\n")
f2out_.write("queue arguments from (\n")

### File with condor launcher for histos
f3out_ = open("runScoutingHistos_onCondor_BToPhi.sub", "w") 
f3out_.write("executable      = $ENV(STARTDIR)/condor/condorHistos_executable.sh\n")
f3out_.write("output          = $ENV(STARTDIR)/condor/plotting_logs/job.$(ClusterId).$(ProcId).out\n")
f3out_.write("error           = $ENV(STARTDIR)/condor/plotting_logs/job.$(ClusterId).$(ProcId).err\n")
f3out_.write("log             = $ENV(STARTDIR)/condor/plotting_logs/job.$(ClusterId).$(ProcId).log\n")
f3out_.write("\n")
f3out_.write("RequestCpus   = 2\n")
f3out_.write("RequestMemory = 8000\n")
f3out_.write("\n")
f3out_.write("getenv = True\n")
f3out_.write("+JobFlavour = 'workday'\n")
f3out_.write("+DESIRED_Sites = 'T2_US_UCSD'\n")
f3out_.write("transfer_input_files = $ENV(STARTDIR)/package.tar.gz\n")
f3out_.write("should_transfer_files = YES\n")
f3out_.write("when_to_transfer_output = ON_EXIT\n")
f3out_.write("x509userproxy=$ENV(X509_USER_PROXY)\n")
f3out_.write("use_x509userproxy = True\n")
f3out_.write("\n")
f3out_.write("queue arguments from (\n")



for [mass, time] in mass_points:
    if "2022postEE" in model_template:
        era = "2022postEE"
        year = "2022"
    elif "2022" in model_template:
        era = "2022"
        year = "2022"
    elif "2023BPix" in model_template:
        era = "2023BPix"
        year = "2023"
    elif "2023" in model_template:
        era = "2023"
        year = "2023"
    if mass == '0p5' and era == '2022' and time == '100':
        continue
    if mass == '4p0' and era == '2022' and time == '10':
        continue
    #os.system("mv {}*.root {}".format(model_template, destination))
    parent_dir = model_template.format(MASS=mass, CTAU=time)
    first_level_subdirectories = [d for d in os.listdir(parent_dir) if os.path.isdir(os.path.join(parent_dir, d))]
    first_level_path = os.path.join(parent_dir, first_level_subdirectories[0])
    second_level_subdirectories = [d for d in os.listdir(first_level_path) if os.path.isdir(os.path.join(first_level_path, d))]
    full_path = os.path.join(first_level_path, second_level_subdirectories[0])
    dataset = full_path
    #dataset = model_template.format(MASS = mass, CTAU = time)
    signalid = signal_template.format(MASS = mass, CTAU = time) + '_' + era
    line = "{},{}\n".format(signalid, dataset)
    fout_.write(line)
    f2out_.write("$ENV(SCOUTINGOUTPUTDIR) {} {} 0 100 1 1\n".format(year, signalid))
    f3out_.write("$ENV(SCOUTINGINPUTDIR) $ENV(SCOUTINGOUTPUTDIR) --condor --signal --inSample {} --splitIndex 0 --splitPace 1000000 $ENV(SCOUTINGARGS)\n".format(signalid))
f2out_.write(")\n")
f3out_.write(")\n")
fout_.close()
f2out_.close()
