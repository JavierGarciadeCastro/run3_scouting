import os

## This script put the crab output in the corresponding ceph area and creates the files to run the looper and plotter

# To be modified by user:
model_template = "/HTo2ZdTo2mu2x_MZd-{MASS}_ctau-{CTAU}mm_TuneCP5_13p6TeV_madgraph-pythia8/Run3Scouting-crab_centralSkim__HTo2ZdTo2mu2x_2023BPix_m-{MASS}_ctau-{CTAU}mm_5p0-1f253bfad14a992d8887e2452ef8b1e3/USER"
signal_template = "Signal_HTo2ZdTo2mu2x_MZd-{MASS}_ctau-{CTAU}mm"
model='HTo2ZdTo2mu2x'

mass_points = []
if model=='HTo2ZdTo2mu2x':
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



### File with central datasets fr looper input
fout_ = open("centralDatasets.txt", "w")

### File with condor launcher for looper
f2out_ = open("runScoutingLooper_onCondor_CentralSignal.sub", "w") 
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
f3out_ = open("runScoutingHistos_onCondor_CentralSignal.sub", "w") 
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
    #os.system("mv {}*.root {}".format(model_template, destination))
    dataset = model_template.format(MASS = mass, CTAU = time)
    signalid = signal_template.format(MASS = mass, CTAU = time) + '_' + era
    line = "{},{}\n".format(signalid, dataset)
    fout_.write(line)
    f2out_.write("$ENV(SCOUTINGOUTPUTDIR) {} {} 0 100 1 1\n".format(year, signalid))
    f3out_.write("$ENV(SCOUTINGINPUTDIR) $ENV(SCOUTINGOUTPUTDIR) --condor --signal --inSample {} --splitIndex 0 --splitPace 1000000 $ENV(SCOUTINGARGS)\n".format(signalid))
f2out_.write(")\n")
f3out_.write(")\n")
fout_.close()
f2out_.close()
