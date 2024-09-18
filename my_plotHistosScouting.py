import ROOT
import os,sys,json
import argparse
from datetime import date    
import numpy as np
import copy
import math
sys.path.append('utils')
import plotUtils

ROOT.gROOT.SetBatch(1)

user = os.environ.get("USER")
today= date.today().strftime("%b-%d-%Y")

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("--inDir", default=os.environ.get("PWD")+"/outputHistograms_"+today, help="Choose output directory. Default: '"+os.environ.get("PWD")+"/outputHistograms_"+today+"'")
parser.add_argument("--inSamples", default=[], nargs="+", help="Choose sample(s); for all data samples in input directory, choose 'Data'")
parser.add_argument("--inMultiDir", default=[], nargs="+", help="Choose directories for one sample")
parser.add_argument("--inMultiLeg", default=[], nargs="+", help="Choose legends for different flavors of one sample, if --inMultiDir is used")
parser.add_argument("--inMultiWeights", default=[], nargs="+", help="Choose weights for different flavors of one sample, if --inMultiDir is used")
parser.add_argument("--outDir", default=os.environ.get("PWD")+"/plots_"+today, help="Choose output directory. Default: '"+os.environ.get("PWD")+"/plots_"+today+"'")
parser.add_argument("--outSuffix", default="", help="Choose output directory. Default: ''")
parser.add_argument("--data", default=False, action="store_true", help="Plot data")
parser.add_argument("--signal", default=False, action="store_true", help="Plot signal")
parser.add_argument("--generator", default=False, action="store_true", help="Plot GEN-level histograms")
parser.add_argument("--year", default="2022", help="Year to be processes. Default: 2022")
parser.add_argument("--lumi", default=0.0, help="Luminosity in case we want to set one manually")
parser.add_argument("--inYears", default=[], nargs="+", help="Choose years if data/simulation from different years is used (using year if not set)")
parser.add_argument("--partialUnblinding", default=False, action="store_true", help="Have processed x% (default: x=50) of available data")
parser.add_argument("--partialUnblindingFraction", default="0.5", help="Fraction of available data that have been processed")
parser.add_argument("--dimuonMassSel", default=[], nargs="+", help="Selection on dimuon mass: first (or only) value is lower cut, second (optional) value is upper cut")
parser.add_argument("--dimuonPtSel", default=[], nargs="+", help="Selection on dimuon pT: first (or only) value is lower cut, second (optional) value is upper cut")
parser.add_argument("--fourmuonMassSel", default=[], nargs="+", help="Selection on four-muon mass: first (or only) value is lower cut, second (optional) value is upper cut")
parser.add_argument("--fourmuonPtSel", default=[], nargs="+", help="Selection on four-muon pT: first (or only) value is lower cut, second (optional) value is upper cut")
parser.add_argument("--dimuonMassSelForFourMuon", default=[], nargs="+", help="Selection on dimuon mass in four-muon system: first (or only) value is lower cut, second (optional) value is upper cut")
parser.add_argument("--dimuonMassDiffSelForFourMuon", default=["0.05"], nargs="+", help="Selection on dimuon mass difference / mean in four-muon system: first (or only) value is *upper* cut, second (optional) value is *lower* cut")
parser.add_argument("--dimuonPtSelForFourMuon", default=[], nargs="+", help="Selection on dimuon pT in four-muon system: first (or only) value is lower cut, second (optional) value is upper cut")
parser.add_argument("--dimuonMassSelForFourMuonOSV", default=[], nargs="+", help="Selection on dimuon mass in four-muon system from overlapping SV: first (or only) value is lower cut, second (optional) value is upper cut")
parser.add_argument("--dimuonMassDiffSelForFourMuonOSV", default=["0.05"], nargs="+", help="Selection on dimuon mass difference / mean in four-muon system from overlapping SV: first (or only) value is *upper* cut, second (optional) value is *lower* cut")
parser.add_argument("--dimuonPtSelForFourMuonOSV", default=[], nargs="+", help="Selection on dimuon pT in four-muon system from overlapping SV: first (or only) value is lower cut, second (optional) value is upper cut")
parser.add_argument("--lxySel", default=[], nargs="+", help="Selection on lxy: first (or only) value is lower cut, second (optional) value is upper cut")
parser.add_argument("--relaxedSVSel", default=False, action="store_true", help="Extend range for 1D histograms of SV selection variables")
parser.add_argument("--doRatio", default=False, action="store_true", help="Plot ratios for 1D histograms")
parser.add_argument("--shape", default=False, action="store_true", help="Shape normalization")
parser.add_argument("--scaleSignal", default=1.0, help="Signal normalization")
parser.add_argument("--rebinWindow", default=0, help="Rebin mass window, set to 1 for automattic rebinning")
parser.add_argument("--logY", default=False, action="store_true", help="Log-scale for Y axis")
parser.add_argument("--logX", default=False, action="store_true", help="Log-scale for X axis")
parser.add_argument("--zoomMass", default=[], nargs="+", help="Zoom mass [1] < mass < [0] GeV")
parser.add_argument("--zoomLxy", default=[], nargs="+", help="Zoom lxy [1] < lxy < [0] cm")
parser.add_argument("--zoomPixel2D", default=False, action="store_true", help="Zoom -20 <(x,y)< 20 cm")
parser.add_argument("--removeBeampipe2D", default=False, action="store_true", help="Remove points within beampipe")
parser.add_argument("--noPreSel", default=False, action="store_true", help="Do not plot pre-selection/association histograms")
parser.add_argument("--noDiMuon", default=False, action="store_true", help="Do not plot dimuon histograms")
parser.add_argument("--noFourMuon", default=False, action="store_true", help="Do not plot four-muon histograms for four-muon systems")
parser.add_argument("--noFourMuonOSV", default=False, action="store_true", help="Do not plot four-muon histograms for four-muon systems from overlapping SVs")
parser.add_argument("--plotOSV", default=False, action="store_true", help="Plot histograms for (di-/four-) muons from overlapping SVs")
parser.add_argument("--extraLabel", default="", help="Label to put on top of the plot")
parser.add_argument("--extraLabelBold", default="", help="Label to put on top of the plot")
parser.add_argument("--pdf", default=False, action="store_true", help="Output format: .pdf. Default: .png")
parser.add_argument("--dxyScaledAxis", default=False, action="store_true", help="X axis from 0 to 1 in lifetime scaled dxy")
args = parser.parse_args()

isData   = args.data
isSignal = args.signal
samples = args.inSamples
if len(samples)<1:
    print("Please, specify list of samples to be plotted")
    exit()

for s in samples:
    if "Data" in s:
        isData = True
        break

skimFraction = float(args.partialUnblindingFraction)
skimEvents = args.partialUnblinding
if not isData:
    skimEvents = False
luminosity2022 = 35.144834218
luminosity2022B = 0.085456340
luminosity2022C = 5.070714851
luminosity2022D = 3.006332096
luminosity2022E = 5.866801170
luminosity2022F = 18.006671456
luminosity2022G = 3.108858306
luminosity2023 = 27.208114203999997
#luminosity2023 = 17.0604 # Excluding B and D
luminosity2023B = 0.622430830
luminosity2023C = 5.557004785
luminosity2023C_triggerV10 = 11.503479528
luminosity2023D = 9.525199061
if float(args.lumi) > 0.00001:
    luminosity = float(args.lumi)
elif samples[0]=="Data" and args.year=="allYears":
    luminosity = luminosity2022 + luminosity2023
elif samples[0]!="Data" and len(samples)==1 and args.year=="2022":
    ts = samples[0]
    if ts=="DataB":
        luminosity = luminosity2022B
    elif ts=="DataC":
        luminosity = luminosity2022C
    elif ts=="DataD":
        luminosity = luminosity2022D
    elif ts=="DataE":
        luminosity = luminosity2022E
    elif ts=="DataF":
        luminosity = luminosity2022F
    elif ts=="DataG":
        luminosity = luminosity2022G
    else:
        luminosity = luminosity2022
elif samples[0]=="Data" and args.year=="2022":
    luminosity = luminosity2022
elif samples[0]!="Data" and len(samples)==1 and args.year=="2023":
    ts = samples[0]
    if ts=="DataB":
        luminosity = luminosity2023B
    elif ts=="DataC":
        luminosity = luminosity2023C
    elif ts=="DataC-triggerV10":
        luminosity = luminosity2023C_triggerV10
    elif ts=="DataD":
        luminosity = luminosity2023D
elif samples[0]=="Data" and args.year=="2023":
    luminosity = luminosity2023
else:
    luminosity = 1.0
luminosity      = 0.1*luminosity
luminosity2022  = 0.1*luminosity2022
luminosity2022B = 0.1*luminosity2022B
luminosity2022C = 0.1*luminosity2022C
luminosity2022D = 0.1*luminosity2022D
luminosity2022E = 0.1*luminosity2022E
luminosity2022F = 0.1*luminosity2022F
luminosity2022G = 0.1*luminosity2022G
luminosity2023  = 0.1*luminosity2023
luminosity2023B = 0.1*luminosity2023B
luminosity2023C = 0.1*luminosity2023C
luminosity2023C_triggerV10 = 0.1*luminosity2023C_triggerV10
luminosity2023D = 0.1*luminosity2023D
if skimEvents and skimFraction>0.0:
    luminosity      = skimFraction*luminosity
    luminosity2022B = skimFraction*luminosity2022B
    luminosity2022C = skimFraction*luminosity2022C
    luminosity2022D = skimFraction*luminosity2022D
    luminosity2022E = skimFraction*luminosity2022E
    luminosity2022F = skimFraction*luminosity2022F
    luminosity2022G = skimFraction*luminosity2022G



legnames = dict()
legnames["Data"]  = "Data"
legnames["DataB"] = "Run2022B (%.2f/fb)"%luminosity2022B
legnames["DataC"] = "Run2022C (%.2f/fb)"%luminosity2022C
legnames["DataD"] = "Run2022D (%.2f/fb)"%luminosity2022D
legnames["DataE"] = "Run2022E (%.2f/fb)"%luminosity2022E
legnames["DataF"] = "Run2022F (%.2f/fb)"%luminosity2022F
legnames["DataG"] = "Run2022G (%.2f/fb)"%luminosity2022G
legnames["DileptonMinBias"] = "MinBias simulation"
legnames["Signal_HTo2ZdTo2mu2x_MZd10_Epsilon1e-06"] = "HAHM: 10, #epsilon = 1x10^{-6}"
legnames["Signal_HTo2ZdTo2mu2x_MZd10_Epsilon5e-07"] = "HAHM: 10, #epsilon = 5x10^{-7}"
legnames["Signal_HTo2ZdTo2mu2x_MZd10_Epsilon1e-07"] = "HAHM: 10, #epsilon = 1x10^{-7}"
legnames["Signal_HTo2ZdTo2mu2x_MZd10_Epsilon3e-08"] = "HAHM: 10, #epsilon = 3x10^{-8}"
legnames["Signal_ScenB1_30_9p9_4p8_ctau_1mm"] = "#pi_{3}#rightarrowA'A': (m_{#pi_{3}},m_{A'}) = (9.9, 4.8) GeV, c#tau_{#pi_{3}} = 1 mm"
legnames["Signal_ScenB1_30_9p9_4p8_ctau_10mm"] = "#pi_{3}#rightarrowA'A': (m_{#pi_{3}},m_{A'}) = (9.9, 4.8) GeV, c#tau_{#pi_{3}} = 10 mm"
legnames["Signal_ScenB1_30_9p9_4p8_ctau_100mm"] = "#pi_{3}#rightarrowA'A': (m_{#pi_{3}},m_{A'}) = (9.9, 4.8) GeV, c#tau_{#pi_{3}} = 100 mm"

for s in samples:
    if "Signal" in s and s not in legnames.keys():
        legnames[s] = s.replace("_"," ")
    if 'Signal_HTo2ZdTo2mu2x_MZd' in s and 'ctau' in s:
        legtxt = "h#rightarrowZ_{{D}}Z_{{D}}: m = {} GeV, c#tau_{{Z_{{D}}}} = {}".format(s.split('MZd-')[1].split('_')[0].replace('p', '.'), s.split( 'ctau-')[1].split('mm')[0])
        if "integrated" not in s:
            legtxt = legtxt + " mm"
        legnames[s] = legtxt
    if 'Signal_BToPhi-' in s and 'ctau' in s:
        mass = s.split('Phi-')[1].split('_')[0].replace('p', '.')
        ctau_part = s.split('ctau-')[1].split('mm')[0]
        ctau = ctau_part.replace('p', '.') if 'p' in ctau_part else ctau_part

        #ctau = s.split('ctau-')[1].split('mm')[0]
        legtxt = "B#rightarrow#PhiX m = {} GeV, c#tau = {}".format(mass, ctau)
        #legtxt = "B#rightarrow#Phi m = {} GeV, c#tau = {}".format(s.split('Phi-')[1].split('_')[0].replace('p', '.'), s.split( 'ctau-')[1].split('mm')[0])
        if "integrated" not in s:
            legtxt = legtxt + " mm"
        legnames[s] = legtxt

samplecol   = []
for s in samples:
    if "Signal" in s or "signal" in s:
        samplecol.append("signal")
    elif "Data" in s:
        samplecol.append(s)
    else:
        samplecol.append("mc")


indir  = args.inDir
outdir = args.outDir
if args.relaxedSVSel:
    outdir = outdir+"_relaxedSVselection"
if args.outSuffix!="":
    outdir = outdir+"_"+args.outSuffix
if not os.path.exists(outdir):
    os.makedirs(outdir)
os.system('cp '+os.environ.get("PWD")+'/utils/index.php '+outdir)

isMultiDir = (len(samples)==1 and len(args.inMultiDir)>1 or len(samples)==len(args.inMultiDir))
inmultidirs = []
inmultilegs = []
inmultiweights = []
if isMultiDir:
    inmultidirs = args.inMultiDir
    inmultilegs = args.inMultiLeg
if len(inmultilegs)<len(inmultidirs):
    for i in range(len(inmultilegs),len(inmultidirs)):
        inmultilegs.append("Unknown")
if len(args.inMultiWeights)>0:
    if len(args.inMultiWeights)!=len(args.inMultiDir):
        print("Weight structure is wrong! Exiting...")
        exit
    else:
        inmultiweights = [float(x) for x in args.inMultiWeights]
else:
    inmultiweights = [1.0 for n in range(0, len(inmultidirs))]

doRatio = args.doRatio
if len(samples)<=1 and not isMultiDir:
    doRatio = False

roff = 0.0
if not doRatio:
    roff=0.01

infiles = []
inyears = []
if len(args.inYears) > 1:
    inyears = args.inYears
else:
    if not isMultiDir:
        inyears = len(samples)*[args.year]
    else:
        inyears = len(inmultidirs)*[args.year]
hname = "histograms"
if args.generator:
    hname = "histograms_GEN"
if not isMultiDir:
    for i,s in enumerate(samples):
        if not os.path.isfile("%s/%s_%s_%s_all.root"%(indir,hname,s,inyears[i])):
            os.system('hadd '+indir+'/'+hname+'_'+s+'_'+inyears[i]+'_all.root $(find '+indir+' -name "'+hname+'_'+s+'*_'+inyears[i]+'*.root")')
        infiles.append(indir+'/'+hname+'_'+s+'_'+inyears[i]+'_all.root')
else:
    for i,d in enumerate(inmultidirs):
        if len(samples) == len(inmultidirs):
            if not os.path.isfile("%s/%s_%s_%s_all.root"%(d,hname,samples[i],inyears[i])):
                os.system('hadd '+d+'/'+hname+'_'+samples[i]+'_'+inyears[i]+'_all.root $(find '+d+' -name "'+hname+'_'+samples[i]+'*_'+inyears[i]+'_*.root")')
            infiles.append(d+'/'+hname+'_'+samples[i]+'_'+inyears[i]+'_all.root')        
        else:
            if not os.path.isfile("%s/%s_%s_%s_all.root"%(d,hname,samples[0],inyears[i])):
                if inyears[i]=='allYears':
                    os.system('hadd '+d+'/'+hname+'_'+samples[0]+'_'+inyears[i]+'_all.root $(find '+d+' -name "'+hname+'_'+samples[0]+'*_*_*.root")')
                else:
                    os.system('hadd '+d+'/'+hname+'_'+samples[0]+'_'+inyears[i]+'_all.root $(find '+d+' -name "'+hname+'_'+samples[0]+'*_'+inyears[i]+'_*.root")')
            infiles.append(d+'/'+hname+'_'+samples[0]+'_'+inyears[i]+'_all.root')        

print(infiles)
if len(infiles)<1:
    print("No matching input file was found in %s."%indir)
    exit()

colors = dict()
colors["Data"]  = [ROOT.kBlack, ROOT.kBlue]
colors["DataB"] = ROOT.kYellow+1
colors["DataC"] = ROOT.kOrange+1
colors["DataD"] = ROOT.kRed+1
colors["DataE"] = ROOT.kPink+1
colors["DataF"] = ROOT.kMagenta+1
colors["DataG"] = ROOT.kViolet+1
colors["signal"] = [ROOT.kAzure+2, ROOT.kOrange+7, ROOT.kGreen+2, ROOT.kRed+1, ROOT.kMagenta-2]
colors["mc"] = [ROOT.kRed-3]

"""
colorsMultiDir = [ROOT.kBlack,
                  ROOT.TColor.GetColor('#017B93'),
                  ROOT.TColor.GetColor('#10D9DC'),
                  ROOT.TColor.GetColor('#F1B950'),
                  ROOT.TColor.GetColor('#CA6702'),
                  ROOT.TColor.GetColor('#AE2012')]
"""

if isMultiDir > 0:
    """
    if len(inmultidirs) < 8:    
        colorsMultiDir = [ROOT.kBlack,
                          ROOT.TColor.GetColor('#448aff'),
                          ROOT.TColor.GetColor('#009688'),
                          ROOT.TColor.GetColor('#8bc34a'),
                          ROOT.TColor.GetColor('#ffc107'),
                          ROOT.TColor.GetColor('#f44336'),
                          ROOT.TColor.GetColor('#f15bb5')]
    """
    if len(inmultidirs) < 7:    
        colorsMultiDir = [ROOT.kBlack,
                          ROOT.TColor.GetColor('#448aff'),
                          ROOT.TColor.GetColor('#8bc34a'),
                          ROOT.TColor.GetColor('#ffc107'),
                          ROOT.TColor.GetColor('#f44336'),
                          ROOT.TColor.GetColor('#f15bb5')]
    elif len(inmultidirs) < 9:    
        colorsMultiDir = [ROOT.kBlack,
                          ROOT.TColor.GetColor('#448aff'),
                          ROOT.TColor.GetColor('#1565c0'),
                          ROOT.TColor.GetColor('#009688'),
                          ROOT.TColor.GetColor('#8bc34a'),
                          ROOT.TColor.GetColor('#ffc107'),
                          ROOT.TColor.GetColor('#f44336'),
                          ROOT.TColor.GetColor('#f15bb5')]
    else:
        colorsMultiDir = [ROOT.kBlack,
                          ROOT.TColor.GetColor('#448aff'),
                          ROOT.TColor.GetColor('#1565c0'),
                          ROOT.TColor.GetColor('#009688'),
                          ROOT.TColor.GetColor('#8bc34a'),
                          ROOT.TColor.GetColor('#ffc107'),
                          ROOT.TColor.GetColor('#ff9800'),
                          ROOT.TColor.GetColor('#f44336'),
                          ROOT.TColor.GetColor('#ad1457'),
                          ROOT.TColor.GetColor('#9d4edd')]
files_not_working = ['hgenmu_pt','hselmuon_mudphimm','hselmuon_uphi_mudphimm','hselmuon_mutheta','hselmuon_uphi_mutheta','hselmuon_muthetamm','hselmuon_uphi_muthetamm','hdimuon_gen_deltalxy','hselmuon_osv_mudphimm','hselmuon_osv_uphi_mudphimm','hselmuon_osv_mutheta','hselmuon_osv_uphi_mutheta','hselmuon_osv_muthetamm','hselmuon_osv_uphi_muthetamm']
files_not_working_2d = ['hsvsel_xerrvslxy','hsvsel_yerrvslxy','hsvsel_zerrvslxy','hsvsel_xerrvsz','hsvsel_yerrvsz','hsvsel_zerrvsz','hsvselass_xerrvslxy','hsvselass_yerrvslxy','hsvselass_zerrvslxy','hsvselass_xerrvsz','hsvselass_yerrvsz','hsvselass_zerrvsz']
fin = ROOT.TFile.Open(infiles[0],"r")
fin.ls()
listkeys = fin.GetListOfKeys()
size = listkeys.GetSize()
h1dn = []
h2dn = []
for i in range(0,size):
    if "TH1" in listkeys.At(i).GetClassName():
        if listkeys.At(i).GetClassName() in files_not_working:
            continue
        h1dn.append(listkeys.At(i).GetName())
    elif "TH2" in listkeys.At(i).GetClassName():
        if listkeys.At(i).GetClassName() in files_not_working:
            continue
        h2dn.append(listkeys.At(i).GetName())
h1dn = [file for file in h1dn if file not in files_not_working]
h2dn = [file for file in h2dn if file not in files_not_working_2d]
h1d = []
h2d = []
inf = []

# Legend
ncl = 1
xol = 0.30
if len(samples)>5:
    ncl=2
    xol=0.5
if len(inmultidirs)>5:
    longLabel = False
    for label in inmultilegs:
        if len(label) > 30:
            longLabel = True
            break
    if not longLabel:
        ncl=2
        xol=0.5
    else:
        xol=0.35

leg = ROOT.TLegend(0.62-xol, 0.89-0.045*max(len(samples),len(inmultidirs)+1)/ncl, 0.89, 0.89)
leg.SetNColumns(ncl)
leg.SetFillColor(0)
leg.SetFillStyle(0)
leg.SetLineWidth(0)
leg.SetTextSize(0.033)
leg.SetMargin(0.1)

nDataSamples = 0
nSigSamples = 0
nMCSamples = 0
i = -1
#print(h1dn)
#print(len(h1dn))
for fn,f in enumerate(infiles):
    print(f)
    inf.append(ROOT.TFile(f))
    if len(h1dn)>0:
        h1d.append([])
    if len(h2dn)>0:
        h2d.append([])
    for hn in h1dn:
        i += 1
        if not args.plotOSV:
            if "osv" in hn:
                continue
        if args.noPreSel:
            if "hsvsel_" in hn:
                continue
            if "h_nsvsel" in hn and "ass" not in hn:
                continue
            if "hmuon" in hn or "h_nmuons" in hn:
                continue
        if args.noDiMuon:
            if "dimuon" in hn:
                continue
            if "ass" in hn and "fourmu" not in hn:
                continue
            if "selmuon" in hn and "fourmu" not in hn:
                continue
        if args.noFourMuon:
            if "fourmu" in hn and "osv" not in hn:
                continue
        if args.noFourMuonOSV:
            if "fourmu" in hn and "osv" in hn:
                continue
        if hn in files_not_working:
            continue
        
        ht = inf[fn].Get(hn).Clone("%s_%d"%(hn,fn))
        #ht.SetDirectory(0)
        h1d[fn].append(copy.deepcopy(ht))
        if not isMultiDir:
            if ("signal" not in samplecol[fn]) and ("mc" not in samplecol[fn]):
                h1d[fn][len(h1d[fn])-1].SetLineColor  (colors[samplecol[fn]][nDataSamples])
                h1d[fn][len(h1d[fn])-1].SetMarkerColor(colors[samplecol[fn]][nDataSamples])
                h1d[fn][len(h1d[fn])-1].SetMarkerStyle(20)
                h1d[fn][len(h1d[fn])-1].SetMarkerSize(0.5)
            elif "mc" not in samplecol[fn]:
                h1d[fn][len(h1d[fn])-1].SetLineColor  (colors[samplecol[fn]][nSigSamples])
                h1d[fn][len(h1d[fn])-1].SetMarkerColor(colors[samplecol[fn]][nSigSamples])
            else:
                h1d[fn][len(h1d[fn])-1].SetLineColor  (colors[samplecol[fn]][nMCSamples])
                h1d[fn][len(h1d[fn])-1].SetMarkerColor(colors[samplecol[fn]][nMCSamples])
                h1d[fn][len(h1d[fn])-1].SetFillColorAlpha(colors[samplecol[fn]][nMCSamples], 0.5)
        else:
            h1d[fn][len(h1d[fn])-1].SetLineColor  (colorsMultiDir[fn])
            h1d[fn][len(h1d[fn])-1].SetMarkerColor(colorsMultiDir[fn])
    if not isMultiDir:
        if "Data" in samplecol[fn]:
            nDataSamples = nDataSamples+1
        if "signal" in samplecol[fn]:
            nSigSamples = nSigSamples+1
        if "mc" in samplecol[fn]:
            nMCSamples = nMCSamples+1
    for hn in h2dn:
        #print(hn)
        if hn in files_not_working_2d:
            continue
        if not args.plotOSV:
            if "osv" in hn:
                continue
        if args.noPreSel:
            if "hsvsel_" in hn:
                continue
        if args.noDiMuon:
            if "dimuon" in hn:
                continue
            if "hsvselass_" in hn:
                continue
        if args.noFourMuon:
            if ("fourmu" in hn and "osv" not in hn):
                continue
        if args.noFourMuonOSV:
            if ("fourmu" in h and "osv" in h):
                continue
        if hn == 'd_Dimuon_excluded_rawmass':
            continue
        ht = inf[fn].Get(hn).Clone("%s_%d"%(hn,fn))
        #ht.SetDirectory(0)
        h2d[fn].append(copy.deepcopy(ht))
    
    if len(h1d[fn])>0:
        if not isMultiDir:
            if "Data" in samples[fn] and '2022' in inyears and '2023' in inyears:
                leg.AddEntry(h1d[fn][0], legnames[samples[fn]] + " ({})".format(inyears[fn]), "PEL")
            elif "mc" in samplecol[fn]:
                leg.AddEntry(h1d[fn][0], legnames[samples[fn]], "F")
            else:
                leg.AddEntry(h1d[fn][0], legnames[samples[fn]], "L")
        else:
            #leg.SetHeader(legnames[samples[0].split(" ")[0]])
            if "Data" in samples[0] and '2022' in inyears and '2023' in inyears:
                leg.AddEntry(h1d[fn][0], inmultilegs[fn] + " ({})".format(inyears[fn]), "PEL")
            else:
                leg.AddEntry(h1d[fn][0], inmultilegs[fn], "L")
# Draw histograms
unityArea = args.shape
scaleSignal = (not args.shape) and (not args.scaleSignal==1.0)
doLogy = args.logY
doLogx = args.logX
rebinWindow = int(args.rebinWindow)
dxyScaledAxis = args.dxyScaledAxis
# Weights (deactivated as signal should be weighted in filler for better handling of different eras conditions)
# Default xsec of 1 pb for every signal
#nevents = {}
#with open("data/Info_SignalMC.txt", "r") as file_:
#    info = file_.readlines()
#    for s in info:
#        if s[0]=="#": continue
#        s.replace('\n','')
#        nevents[s.split(',')[0]] = float(s.split(',')[1])
#weights = []
#for s_,s in enumerate(samples):
#    if "Signal" in s:
#        sampleTag = s.replace('Signal_', '').split('_202')[0]
#        if s not in nevents.keys():
#            print("ATTENTION: Not weighting due to lack of information in Info_SignalMC.txt")
#            weights.append(1.0) # Default for data (MC not considered)
#        elif "2022" in inyears[s_]:
#            weights.append(1000.0*luminosity2022/nevents[s]*100)
#        elif "2023" in inyears[s_]:
#            weights.append(1000.0*luminosity2023/nevents[s])
#    else:
#        weights.append(1.0) # Default for data (MC not considered)
weights = []
for s_,s in enumerate(samples):
    if "Signal" in s:
        weights.append(float(args.scaleSignal))
    else:
        weights.append(1.0)

# Labels
if args.year!='allYears':
    yearenergy = "%.2f fb^{-1} (%s, 13.6 TeV)"%(luminosity,args.year)
else:
    yearenergy = "%.2f fb^{-1} (2022+2023, 13.6 TeV)"%(luminosity)
cmsExtra = "Preliminary"
if not isData:
    yearenergy = "(13.6 TeV)"
    cmsExtra = "Simulation"

drawCMSOnTop = True
#
latex = ROOT.TLatex()
latex.SetTextFont(42)
latex.SetTextAlign(31)
latex.SetTextSize(0.04)
latex.SetNDC(True)
#
latexCMS = ROOT.TLatex()
latexCMS.SetTextFont(61)
latexCMS.SetTextSize(0.055)
latexCMS.SetNDC(True)
#
latexCMSExtra = ROOT.TLatex()
latexCMSExtra.SetTextFont(52)
latexCMSExtra.SetTextSize(0.04)
latexCMSExtra.SetNDC(True)
#
latexExtra = ROOT.TLatex()
latexExtra.SetTextFont(42)
latexExtra.SetTextSize(0.034)
latexExtra.SetNDC(True)
#
latexExtraBold = ROOT.TLatex()
latexExtraBold.SetTextFont(62)
latexExtraBold.SetTextSize(0.034)
latexExtraBold.SetNDC(True)

h1dr     = []
h1dr_den = []
h_axis       = []
h_axis_ratio = []
minX         = []
maxX         = []
minY         = []
maxY         = []
line         = []
for fn in range(len(infiles)):
    h1dr.append([])

if not args.plotOSV:
    h1dn = [h for h in h1dn if not ("osv" in h)]
    h2dn = [h for h in h2dn if not ("osv" in h)]
if args.noPreSel:
    h1dn = [h for h in h1dn if not ("hsvsel_" in h)]
    h1dn = [h for h in h1dn if not ("h_nsvsel" in h and "ass" not in h)]
    h1dn = [h for h in h1dn if not ("hmuon" in h)]
    h1dn = [h for h in h1dn if not ("h_nmuons" in h)]
    h2dn = [h for h in h2dn if not ("hsvsel_" in h)]
if args.noDiMuon:
    h1dn = [h for h in h1dn if not ("dimuon" in h)]
    h1dn = [h for h in h1dn if not ("ass" in h and "fourmu" not in h)]
    h1dn = [h for h in h1dn if not ("semuon" in h and "fourmu" not in h)]
    h2dn = [h for h in h2dn if not ("hsvselass_" in h)]
if args.noFourMuon:
    h1dn = [h for h in h1dn if not ("fourmu" in h and "osv" not in h)]
if args.noFourMuonOSV:
    h1dn = [h for h in h1dn if not ("fourmu" in h and "osv" in h)]
for hn,hnn in enumerate(h1dn):
    isZoom = False
    tminY = 1e100
    tmaxY = 0.0
    for fn in range(len(infiles)):
        sfSVrange = 1.0
        if args.relaxedSVSel:
            sfSVrange = 5.0
        if "_type" not in hnn:
            if 'selmuon_dxyscaled' in hnn:
                if dxyScaledAxis:
                    h1d[fn][hn].Rebin(3)
            xmin=None
            xmax=None
            if "reldmass" in hnn:
                xmax=2.0
            if "lxy" in hnn:
                if len(args.zoomLxy)>0:
                    xmax=float(args.zoomLxy[0])
                    isZoom = True
                    if len(args.zoomLxy)>1:
                        xmin=float(args.zoomLxy[1])
                if ("selass" in hnn or "dimuon" in hnn) and "sig" not in hnn:
                    if len(args.lxySel)>0:
                        xmin=float(args.lxySel[0])
                        isZoom = True
                        if len(args.lxySel)>1:
                            xmax=float(args.lxySel[1])
            if "mass" in hnn and "reld" not in hnn:
                if len(args.zoomMass)>0:
                    xmax=float(args.zoomMass[0])
                    isZoom = True
                    if len(args.zoomMass)>1:
                        xmin=float(args.zoomMass[1])
                if "dimuon" in hnn or "Dimuon" in hnn:
                    if len(args.dimuonMassSel)>0:
                        xmin=float(args.dimuonMassSel[0])
                        isZoom = True
                        if len(args.dimuonMassSel)>1:
                            xmax=float(args.dimuonMassSel[1])
                if "fourmuon_mass" in hnn or "FourMu" in hnn:
                    if len(args.fourmuonMassSel)>0:
                        xmin=float(args.fourmuonMassSel[0])
                        isZoom = True
                        if len(args.fourmuonMassSel)>1:
                            xmax=float(args.fourmuonMassSel[1])
            if "sv" in hnn and ("xerr" in hnn or "yerr" in hnn):
                xmax=0.05*sfSVrange
            if "sv" in hnn and "zerr" in hnn:
                xmax=0.10*sfSVrange
            if "sv" in hnn and "chi2ndof" in hnn and not "muon" in hnn:
                xmax=5.0*sfSVrange
            if "sv" in hnn and "d3derr" in hnn:
                xmax=0.15
            if not isZoom:
                plotUtils.PutUnderflowInFirstBin(h1d[fn][hn],xmin)
                plotUtils.PutOverflowInLastBin(h1d[fn][hn],xmax)
                if "mass" in hnn and "reld" not in hnn and rebinWindow > 0:
                    h1d[fn][hn].Rebin(rebinWindow)
                elif "mass" in hnn and "reld" not in hnn and rebinWindow < 1:
                    h1d[fn][hn].Rebin(10)
                    ytitle = h1d[fn][hn].GetYaxis().GetTitle()
                    ytitle = ytitle.replace("0.01","0.1")
                    h1d[fn][hn].GetYaxis().SetTitle(ytitle)
                if "sv" in hnn and "sig" in hnn:
                    h1d[fn][hn].Rebin(5)
                    ytitle = h1d[fn][hn].GetYaxis().GetTitle()
                    ytitle = ytitle.replace("0.1","0.5")
                    h1d[fn][hn].GetYaxis().SetTitle(ytitle)
                if "dphisv" in hnn or "3danglesv" in hnn:
                    h1d[fn][hn].Rebin(10)
                    ytitle = h1d[fn][hn].GetYaxis().GetTitle()
                    ytitle = ytitle.replace("0.01","0.1")
                    h1d[fn][hn].GetYaxis().SetTitle(ytitle)
            else:
                if "mass" in hnn and "reld" not in hnn and rebinWindow > 0:
                    if rebinWindow == 1:
                        width = xmax - xmin
                        nbins = int(width/0.01)
                        for n in range(30, 51):
                            if nbins%n==0:
                                h1d[fn][hn].Rebin(int(nbins/n))
                                break
                        ytitle = h1d[fn][hn].GetYaxis().GetTitle()
                        ytitle = ytitle.replace("0.01","{:.2f}".format(int(nbins/n)*0.01))
                        h1d[fn][hn].GetYaxis().SetTitle(ytitle)
                    else: 
                        h1d[fn][hn].Rebin(rebinWindow)
                        ytitle = h1d[fn][hn].GetYaxis().GetTitle()
                        ytitle = ytitle.replace("0.01","{:.2f}".format(rebinWindow*0.01))
                        h1d[fn][hn].GetYaxis().SetTitle(ytitle)
            if "reld" in hnn:
                h1d[fn][hn].Rebin(5)
        if scaleSignal:
            h1d[fn][hn].Scale(weights[fn])
        if isMultiDir:
            h1d[fn][hn].Scale(inmultiweights[fn])
        h1dr[fn].append(h1d[fn][hn].Clone("%s_ratio"%hnn))
        tbm = 1
        tbM = h1d[fn][hn].GetNbinsX()
        if fn==0:
            h1dr_den.append(h1d[fn][hn].Clone("%s_denominator"%hnn))
            if not isMultiDir and samples[fn]=="Data":
                for b in range(1,h1dr[fn][hn].GetNbinsX()+1):
                    h1dr_den[hn].SetBinError(b,0.0)
            if isMultiDir and "lxy" not in inmultidirs[fn]:
                for b in range(1,h1dr[fn][hn].GetNbinsX()+1):
                    h1dr_den[hn].SetBinError(b,0.0)
            if xmin!=None:
                tbm = h1d[fn][hn].GetXaxis().FindBin(xmin)
                h1d [fn][hn].GetXaxis().SetRangeUser(h1d[fn][hn].GetXaxis().GetBinLowEdge(tbm),h1d[fn][hn].GetXaxis().GetBinUpEdge(h1d[fn][hn].GetNbinsX()))
                h1dr[fn][hn].GetXaxis().SetRangeUser(h1d[fn][hn].GetXaxis().GetBinLowEdge(tbm),h1d[fn][hn].GetXaxis().GetBinUpEdge(h1d[fn][hn].GetNbinsX()))
            if xmax!=None:
                tbM = h1d[fn][hn].GetXaxis().FindBin(xmax)
                if (h1d[fn][hn].Integral(tbM, -1)<=0.0 and h1d[fn][hn].Integral(1, tbM-1)>0) or (isZoom and not xmax>h1d[fn][hn].GetXaxis().GetBinLowEdge(tbM)):
                    tbM = tbM-1
                h1d [fn][hn].GetXaxis().SetRangeUser(h1d[fn][hn].GetXaxis().GetBinLowEdge(1),h1d[fn][hn].GetXaxis().GetBinUpEdge(tbM))
                h1dr[fn][hn].GetXaxis().SetRangeUser(h1d[fn][hn].GetXaxis().GetBinLowEdge(1),h1d[fn][hn].GetXaxis().GetBinUpEdge(tbM))
            if xmin!=None and xmax!=None:
                tbm = h1d[fn][hn].GetXaxis().FindBin(xmin)
                tbM = h1d[fn][hn].GetXaxis().FindBin(xmax)
                if (h1d[fn][hn].Integral(tbM, -1)<=0.0 and h1d[fn][hn].Integral(tbm, tbM-1)>0) or (isZoom and not xmax>h1d[fn][hn].GetXaxis().GetBinLowEdge(tbM)):
                    tbM = tbM-1
                h1d [fn][hn].GetXaxis().SetRangeUser(h1d[fn][hn].GetXaxis().GetBinLowEdge(tbm),h1d[fn][hn].GetXaxis().GetBinUpEdge(tbM))
                h1dr[fn][hn].GetXaxis().SetRangeUser(h1d[fn][hn].GetXaxis().GetBinLowEdge(tbm),h1d[fn][hn].GetXaxis().GetBinUpEdge(tbM))
            if not isMultiDir and "Data" in samples[fn]:
                h1d [fn][hn].SetBinErrorOption(ROOT.TH1.kPoisson)
                h1dr[fn][hn].SetBinErrorOption(ROOT.TH1.kPoisson)
            if isMultiDir and "Data" in samples[0]:
                h1d [fn][hn].SetBinErrorOption(ROOT.TH1.kPoisson)
                h1dr[fn][hn].SetBinErrorOption(ROOT.TH1.kPoisson)
            h1d[fn][hn].GetYaxis().SetLabelSize(0.025)
            h1d[fn][hn].GetYaxis().SetMaxDigits(3)
            h1d[fn][hn].SetLineWidth(2)
            ytitle = h1d[fn][hn].GetYaxis().GetTitle()
            """ 
            if not unityArea:
                tmaxY = max(tmaxY, 1.4*h1d[fn][hn].GetMaximum())
                tminY = 0.0
                if doLogy:
                    tmaxY = 100*tmaxY
                    tminY = 0.1
                    if weightSignal:
                        tminY = 0.001
            """
        h1d[fn][hn].SetLineWidth(2)
        if unityArea:
            ytitle = ytitle.replace("Events", "Fraction of events")
            integral = h1d[fn][hn].Integral(0,-1)
            if integral>0.0:
                h1d[fn][hn].Scale(1.0/integral)
                h1dr[fn][hn].Scale(1.0/integral)
                if fn==0:
                    h1dr_den[hn].Scale(1.0/integral)
            tmaxY = max(tmaxY, h1d[fn][hn].GetMaximum())
            if doLogy:
                minb = 1e100
                for b in range(1,h1d[fn][hn].GetNbinsX()+1):
                    if b > tbM or b < tbm:
                        continue
                    if h1d[fn][hn].GetBinContent(b)<minb and h1d[fn][hn].GetBinContent(b)>0.0:
                        minb = h1d[fn][hn].GetBinContent(b)
                tminY = min(0.1*minb,0.1)
                if h1d[fn][hn].GetMaximum() > 0.0:
                    tmaxY = 10**((math.log10(h1d[fn][hn].GetMaximum()) - math.log10(tminY)))
                    #print(h1d[fn][hn].GetName(), tmaxY, h1d[fn][hn].GetMaximum(), tminY)
                else:
                    tmaxY = 1.0
            else:
                tminY = 0.0
                tmaxY = 1.5*h1d[fn][hn].GetMaximum()
        else:
            tminY = 0.0
            if doLogy:
                tmaxY = max(tmaxY, 1000*h1d[fn][hn].GetMaximum())
                tminY = 0.1 # change to 0.01
            else:
                tmaxY = max(tmaxY, 1.4*h1d[fn][hn].GetMaximum())
        h1dr[fn][hn].Divide(h1dr_den[hn])
        if "hsv" in hnn and not "mind" in hnn and not "maxd" in hnn:
            ytitle = ytitle.replace("Events", "Number of SVs")
            ytitle = ytitle.replace("events", "SVs")
        if "hmuon" in hnn and hnn!="hmuon_mindr" and hnn!="hmuon_maxdr":
            ytitle = ytitle.replace("Events", "Number of muons")
            ytitle = ytitle.replace("events", "muons")
        h1d[fn][hn].GetYaxis().SetTitle(ytitle)
        #print(hn)
        if fn==0:
            h_axis.append(h1d[fn][hn].Clone("axis_%s"%hnn))
            h_axis[hn].Reset("ICE")
            h_axis[hn].SetTitle("")
            h_axis[hn].GetYaxis().SetTitleSize(0.04)
            h_axis[hn].GetYaxis().SetTitleOffset(1.35)
            h_axis[hn].GetXaxis().SetTitleSize(0.04)
            h_axis[hn].GetXaxis().SetTitleOffset(1.25-roff*10.0)
            h_axis[hn].GetXaxis().SetLabelSize(0.035)
            h_axis[hn].GetYaxis().SetLabelSize(0.035)
            h_axis_ratio.append(h_axis[hn].Clone("axis_ratio_%s"%hnn))
            if xmax==None:
                maxX.append(h_axis[hn].GetXaxis().GetBinUpEdge(plotUtils.GetLastBin(h_axis[hn])))
            else:
                maxX.append(min(h_axis[hn].GetXaxis().GetBinUpEdge(tbM),
                                h_axis[hn].GetXaxis().GetBinUpEdge(plotUtils.GetLastBin(h_axis[hn]))))
            if xmin==None:
                minX.append(h_axis[hn].GetXaxis().GetBinLowEdge(plotUtils.GetFirstBin(h_axis[hn])))
            else:
                minX.append(max(h_axis[hn].GetXaxis().GetBinLowEdge(tbm),
                                h_axis[hn].GetXaxis().GetBinLowEdge(plotUtils.GetFirstBin(h_axis[hn
]))))
            line.append(ROOT.TLine(minX[hn], 1.0, maxX[hn], 1.0))
            minY.append(tminY)
            maxY.append(tmaxY)
        else:
            if tminY<minY[hn]:
                minY[hn]=tminY
            if tmaxY>maxY[hn]:
                maxY[hn]=tmaxY
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetFrameLineWidth(2);
if doRatio: 
    can = ROOT.TCanvas("can","",600, 600)
else:
    can = ROOT.TCanvas("can","",700, 600)
minR = []
maxR = []
for hn,hnn in enumerate(h1dn):
    #if 'selmuon_dxyscaled' not in hnn and 'muon_dxysig' not in hnn:
    #    continue
    pads = []
    if doRatio:
        pads.append(ROOT.TPad("1","1",0.0,0.3,1.0,1.0))
        pads.append(ROOT.TPad("2","2",0.0,0.01,1.0,0.3))
        pads[0].SetTopMargin(0.08)
        pads[0].SetBottomMargin(0.015)
        pads[0].SetLeftMargin(0.10)
        pads[0].SetRightMargin(0.05)
        pads[1].SetTopMargin(0.05)
        pads[1].SetLeftMargin(0.10)
        pads[1].SetRightMargin(0.05)
        pads[1].SetBottomMargin(0.4)
    else:
        pads.append(ROOT.TPad("1","1",0,0,1,1))
        pads[0].SetTopMargin(0.08)
        pads[0].SetLeftMargin(0.10)
        pads[0].SetRightMargin(0.05)
    if doRatio:
        h_axis[hn].GetXaxis().SetTitleSize(0.0)
        h_axis[hn].GetXaxis().SetLabelSize(0.0)
        h_axis[hn].GetYaxis().SetTitleSize(0.05)
        h_axis[hn].GetYaxis().SetLabelSize(0.045)
        h_axis[hn].GetYaxis().SetTitleOffset(1.05)
        minR.append(0.0)
        maxR.append(2.0)
        ty = np.array([])
        tmax=maxR[hn]
        for fn in range(len(infiles)):
            tg = ROOT.TGraph(h1dr[fn][hn])
            ty = tg.GetY()
            if len(ty)>0:
                tmax = np.amax(ty)
            if tmax>maxR[hn]:
                maxR[hn]=tmax*1.05
            if maxR[hn]>5.0:
                minR[hn]=0.1
        h_axis_ratio[hn].GetYaxis().SetRangeUser(minR[hn],maxR[hn])
        h_axis_ratio[hn].SetMinimum(minR[hn])
        h_axis_ratio[hn].SetMaximum(maxR[hn])
        h_axis_ratio[hn].SetTitle(";;Ratio")
        h_axis_ratio[hn].GetYaxis().SetTitleSize(0.12)
        h_axis_ratio[hn].GetYaxis().SetTitleOffset(0.25)
        if doLogy:
            h_axis_ratio[hn].GetYaxis().SetTitleOffset(0.4)
        h_axis_ratio[hn].GetYaxis().SetNdivisions(505)
        h_axis_ratio[hn].GetYaxis().SetLabelSize(0.12)
        h_axis_ratio[hn].GetYaxis().CenterTitle()
        h_axis_ratio[hn].GetYaxis().SetTickLength(0.02)
        h_axis_ratio[hn].GetXaxis().SetLabelSize(0.12)
        h_axis_ratio[hn].GetXaxis().SetTitleSize(0.12)
        h_axis_ratio[hn].GetXaxis().SetTitle(h1d[0][hn].GetXaxis().GetTitle())
        h_axis_ratio[hn].GetXaxis().SetTickSize(0.06)

        pads[0].Draw()
        pads[1].Draw()
        pads[1].cd()
        if maxR[hn]>5.0:
            pads[1].SetLogy()
        pads[1].SetTickx()
        pads[1].SetTicky()
        h_axis_ratio[hn].Draw("AXIS")
        for fn in range(1,len(infiles)):
            h1dr[fn][hn].Draw("SAME,P,E")
        line[hn].SetLineStyle(2)
        if not isMultiDir:
            line[hn].SetLineColor(h1dr[0][hn].GetLineColor())
        else:
            line[hn].SetLineColor(colorsMultiDir[0])
        line[hn].SetLineWidth(1)
        line[hn].Draw("SAME")
        pads[1].Modified()
        pads[1].Update()

    else:
        h_axis[hn].GetXaxis().SetTitleSize(0.04)
        h_axis[hn].GetXaxis().SetLabelSize(0.04)
        h_axis[hn].GetYaxis().SetTitleSize(0.04)
        h_axis[hn].GetYaxis().SetLabelSize(0.04)
        if 'selmuon_dxyscaled' in hnn:
            if dxyScaledAxis:
                h_axis[hn].GetXaxis().SetRangeUser(0,1)
        #h_axis[hn].GetYaxis().SetTitleOffset(1.05)
        pads[0].Draw()
    
    pads[0].cd()
    #pads[0].SetTickx()
    #pads[0].SetTicky()
    if doLogy:
        pads[0].SetLogy()
    if doLogx:
        pads[0].SetLogx()
        if doRatio:
            pads[1].SetLogx()
    #h_axis[hn].GetYaxis().SetMoreLogLabels()
            
    if doLogy:
        maxY[hn] = maxY[hn]
    else:
        maxY[hn] = maxY[hn]
    h_axis[hn].GetYaxis().SetRangeUser(minY[hn], maxY[hn])
    h_axis[hn].SetMinimum(minY[hn])
    h_axis[hn].SetMaximum(maxY[hn])
    h_axis[hn].Draw("")
    fnd = []
    for fn in range(len(infiles)):
        if (isMultiDir and "Data" in samples[0]):
            #h1d[fn][hn].Draw("SAME,P,E")
            h1d[fn][hn].Draw("SAME,HIST")
        elif (not isMultiDir and "Data" in samples[fn]):
            fnd.append(fn)
        else:
            h1d[fn][hn].Draw("SAME,HIST")
    if (len(fnd) > 0):
        for fn_ in fnd:
            h1d[fnd[fn_]][hn].Draw("SAME,P,E")
    leg.Draw("same")
    pads[0].RedrawAxis()
    pads[0].Update()
    can.cd()
    latex.DrawLatex(0.95, 0.945-roff, yearenergy);
    if drawCMSOnTop:
        latexCMS.DrawLatex(0.11,0.945-roff,"CMS");
        if doRatio:
            latexCMSExtra.DrawLatex(0.23,0.945-roff, cmsExtra);
        else:
            latexCMSExtra.DrawLatex(0.21,0.945-roff, cmsExtra);
    else:
        latexCMS.DrawLatex(0.14,0.875,"CMS");
        latexCMSExtra.DrawLatex(0.14,0.825, cmsExtra);
    latexExtra.DrawLatex(0.14,0.8,args.extraLabel);
    latexExtraBold.DrawLatex(0.14,0.85,args.extraLabelBold);
    outname = "%s"%hnn
    if "d_Dimuon" in outname:
        if "excluded" in outname:
            latexExtra.DrawLatex(0.14,0.86,"Dimuon") 
            latexExtra.DrawLatex(0.14,0.76,"Excluded")
        elif 'inclusive' in outname:
            latexExtra.DrawLatex(0.14,0.86,"Dimuon") 
            latexExtra.DrawLatex(0.14,0.76,"Inclusive")
        elif "non-pointing" not in outname:
            latexExtra.DrawLatex(0.14,0.86,"Dimuon") 
            catnames = outname.split("_")
            lxybin = catnames[2]
            lxybin = (lxybin[3:]).split("to")
            latexExtra.DrawLatex(0.14,0.81,"l_{{xy}} #in [{},{}]".format(lxybin[0].replace("p", "."), lxybin[1].replace("p", ".")))
            isobin = "Isolated" if catnames[3]=="iso1" else "Non isolated"
            latexExtra.DrawLatex(0.14,0.76,isobin)
            ptbin = "High p_{T}^{#mu#mu}" if catnames[4]=="pthigh" else "Low p_{T}^{#mu#mu}"
            latexExtra.DrawLatex(0.14,0.71,ptbin)
        else:
            latexExtra.DrawLatex(0.14,0.86,"Dimuon") 
            catnames = outname.split("_")
            lxybin = catnames[2]
            lxybin = (lxybin[3:]).split("to")
            latexExtra.DrawLatex(0.14,0.81,"l_{{xy}} #in [{},{}]".format(lxybin[0].replace("p", "."), lxybin[1].replace("p", ".")))
            latexExtra.DrawLatex(0.14,0.76,"Non-pointing")
    elif "d_FourMu" in outname:
        latexExtraBold.DrawLatex(0.14,0.86,"Four muon") 
        fourmucat = "Resolved" if "sep" in outname else "Overlapping"
        latexExtra.DrawLatex(0.14,0.81,fourmucat) 
    if "h_" in outname and not "d_" in outname:
        outname = outname.replace("h_","")
    elif "hsv" in outname or "hmuon" in outname or "hsel" in outname:
        outname = outname.replace("hsv","sv")
        outname = outname.replace("hmuon","muon")
        outname = outname.replace("hsel","sel")
    elif "hdi" in outname or "hfour" in outname:
        outname = outname.replace("hdi","di")
        outname = outname.replace("hfour","four")
    if unityArea:
        outname = outname+"_unitArea"
    if doLogy:
        outname = outname+"_logY"
    can.SaveAs("%s/%s.png"%(outdir,outname))
    if args.pdf:
        can.SaveAs("%s/%s.pdf"%(outdir,outname))
    can.Clear()
    del pads

for hn,hnn in enumerate(h2dn):
    continue
    for fn in range(len(infiles)):
        h = h2d[fn][hn].Clone()
        if "lxycomp" in hnn:
            if len(args.lxySel)>0:
                 xmin=float(args.lxySel[0])
                 ymin=float(args.lxySel[0])
                 isZoom = True
                 if len(args.lxySel)>1:
                     xmax=float(args.lxySel[1])
                     ymax=float(args.lxySel[1])
                 h.GetXaxis().SetRangeUser(xmin, xmax)
                 h.GetYaxis().SetRangeUser(ymin, ymax)
        if "yvsx" in hnn:
            maxc      = 100.0
            maxcPixel = 17.5
            zoomPixel2D = args.zoomPixel2D
            if len(args.lxySel)>1 and float(args.lxySel[1])<maxcPixel and "selass" in hnn:
                zoomPixel2D = True
            if len(args.zoomLxy)>0 and float(args.zoomLxy[0])<maxcPixel and "selass" in hnn:
                zoomPixel2D = True
            if zoomPixel2D:
                maxc = maxcPixel
            minc = -1.0*maxc
            h.GetXaxis().SetRangeUser(minc, maxc)
            h.GetYaxis().SetRangeUser(minc, maxc)
        if args.removeBeampipe2D:
            maxbp = 1.0
            for b in range(1,h.GetNbinsX()+1):
                if h.GetXaxis().GetBinLowEdge(b)<=0.0 and h.GetXaxis().GetBinLowEdge(b)>-1.0*maxbp:
                    for bb in range(1,h.GetNbinsY()+1):
                        if h.GetYaxis().GetBinLowEdge(bb)<=0.0 and h.GetYaxis().GetBinLowEdge(bb)>-1.0*maxbp:
                            h.SetBinContent(b,bb,0.0)
                            h.SetBinError  (b,bb,0.0)
                        elif h.GetYaxis().GetBinUpEdge(bb)>0.0 and h.GetYaxis().GetBinUpEdge(bb)<maxbp:
                            h.SetBinContent(b,bb,0.0)
                            h.SetBinError  (b,bb,0.0)
                elif h.GetXaxis().GetBinLowEdge(b)>=0.0 and h.GetXaxis().GetBinUpEdge(b)<maxbp:
                    for bb in range(1,h.GetNbinsY()+1):
                        if h.GetYaxis().GetBinLowEdge(bb)<=0.0 and h.GetYaxis().GetBinLowEdge(bb)>-1.0*maxbp:
                            h.SetBinContent(b,bb,0.0)
                            h.SetBinError  (b,bb,0.0)                
                        elif h.GetYaxis().GetBinUpEdge(bb)>0.0 and h.GetYaxis().GetBinUpEdge(bb)<maxbp:
                            h.SetBinContent(b,bb,0.0)
                            h.SetBinError  (b,bb,0.0)
        can.cd()
        ROOT.gPad.SetLogy(0)
        if doLogy:
            ROOT.gPad.SetLogz()
        h.GetXaxis().SetLabelSize(0.025)
        h.GetYaxis().SetLabelSize(0.025)
        ztitle = h.GetZaxis().GetTitle()
        integral = h.Integral(0,-1,0,-1)
        if unityArea:
            ztitle = ztitle.replace("Events", "Fraction of events")
            ztitle = ztitle.replace("Number", "Fraction")
            if integral>0:
                h.Scale(1.0/integral)
            if not doLogy:
                h.GetZaxis().SetRangeUser(0.0,1.1*h.GetMaximum())
            else:
                if integral>0:
                    h.GetZaxis().SetRangeUser(1.0/integral,1.1*h.GetMaximum())
                else:
                    h.GetZaxis().SetRangeUser(1e-1,1.1*h.GetMaximum())
        if "nmuons" in h.GetName():
            ztitle = ztitle.replace("Events", "Number of muons")
            ztitle = ztitle.replace("events", "muons")
        xtitle = h.GetXaxis().GetTitle()
        ytitle = h.GetYaxis().GetTitle()
        if "yvsx" in hnn:
            xtitle = xtitle.replace("(from PV) ","")
            ytitle = ytitle.replace("(from PV) ","")
        h.GetXaxis().SetTitle(xtitle)
        h.GetYaxis().SetTitle(ytitle)
        h.GetZaxis().SetTitle(ztitle)
        if isMultiDir:
            h.SetTitle(inmultilegs[fn])
        else:
            h.SetTitle("")
        can.cd()
        if unityArea:
            ROOT.gStyle.SetPaintTextFormat(".2f");
        h.Draw("colz")
        can.Update()
        palette = h.GetListOfFunctions().FindObject("palette")
        if palette:
            palette.SetX2NDC(0.925)
        can.Update()
        h.GetZaxis().SetTitle(ztitle)
        h.GetZaxis().SetLabelSize(0.025)
        h.GetZaxis().SetMaxDigits(2)
        if not doLogy:
            h.GetZaxis().SetRangeUser(0.0,1.1*h.GetMaximum())
        else:
            if integral>0.0:
                h.GetZaxis().SetRangeUser(1.0/integral,2.0*h.GetMaximum())
            else:
                h.GetZaxis().SetRangeUser(1e-1,1.1*h.GetMaximum())
        latex.DrawLatex(0.90, 0.91, yearenergy);
        if drawCMSOnTop:
            latexCMS.DrawLatex(0.1,0.91,"CMS");
            latexCMSExtra.DrawLatex(0.19,0.91, cmsExtra);
        else:
            latexCMS.DrawLatex(0.13,0.86,"CMS");
            latexCMSExtra.DrawLatex(0.13,0.815, cmsExtra);
        ROOT.gPad.RedrawAxis()
        outname = ""
        if not isMultiDir:
            outname = "%s_%s"%(hnn,samples[fn])
        else:
            outname = "%s_dir%d"%(hnn,fn)
        if "h_" in outname:
            outname = outname.replace("h_","")
        elif "hsv" in outname:
            outname = outname.replace("hsv","sv")
        if args.zoomPixel2D:
            outname = outname+"_zoomPixel"
        if args.removeBeampipe2D:
            outname = outname+"_noBeampipe"
        if unityArea:
            outname = outname+"_unitArea"
        if doLogy:
            outname = outname+"_logZ"
        can.SaveAs("%s/%s.png"%(outdir,outname))
        if args.pdf:
            can.SaveAs("%s/%s.pdf"%(outdir,outname))
        can.Clear()

del can
