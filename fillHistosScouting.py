import os,sys,json
import argparse
from datetime import date    
import ROOT
import numpy as np
from DataFormats.FWLite import Events, Handle
sys.path.append('utils')
import histDefinition
import math
import csv

ROOT.EnableImplicitMT(2)

user = os.environ.get("USER")
today= date.today().strftime("%b-%d-%Y")

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("--inDir", default="/ceph/cms/store/user/"+user+"/Run3ScoutingOutput/looperOutput_"+today, help="Choose input directory. Default: '/ceph/cms/store/user/"+user+"/Run3ScoutingOutput/looperOutput_"+today+"'")
parser.add_argument("--inSample", default="*", help="Choose sample; for all samples in input directory, choose '*'")
parser.add_argument("--inFile", default="*", help="Choose input file by index (for debug); for all files in input directory, choose '*'")
parser.add_argument("--outDir", default=os.environ.get("PWD")+"/outputHistograms_"+today, help="Choose output directory. Default: '"+os.environ.get("PWD")+"/outputHistograms_"+today+"'")
parser.add_argument("--outSuffix", default="", help="Choose output directory. Default: ''")
parser.add_argument("--condor", default=False, action="store_true", help="Run on condor")
parser.add_argument("--data", default=False, action="store_true", help="Process data")
parser.add_argument("--signal", default=False, action="store_true", help="Process signal")
parser.add_argument("--year", default="2022", help="Year to be processed. Default: 2022")
parser.add_argument("--weightMC", default=True, help="Indicate if MC is weighted")
parser.add_argument("--rooWeight", default="1.00", help="Weight to be used for RooDatasets and Signal Regions (It doesn't weight other histograms)")
parser.add_argument("--partialUnblinding", default=False, action="store_true", help="Process x% (default: x=50) of available data")
parser.add_argument("--partialUnblindingFraction", default="0.5", help="Fraction of available data to be processed")
parser.add_argument("--removeDuplicates", default=False, action="store_true", help="Check for and remove duplicates")
parser.add_argument("--splitIndex", default="-1", help="Split index")
parser.add_argument("--splitPace", default="250000", help="Split pace")
parser.add_argument("--dimuonMassSel", default=[], nargs="+", help="Selection on dimuon mass: first (or only) value is lower cut, second (optional) value is upper cut")
parser.add_argument("--dimuonMassSidebandSel", default=[], nargs="+", help="Selection on dimuon mass sidebands: first pair of values is left sideband, second (optional) is right sideband")
parser.add_argument("--dimuonPtSel", default=[], nargs="+", help="Selection on dimuon pT: first (or only) value is lower cut, second (optional) value is upper cut")
parser.add_argument("--dimuonIsoCatSel", default=-99, help="Selection on dimuon isolation category: 0 (non-iso), 1 (part-iso) or 2 (iso)")
parser.add_argument("--fourmuonMassSel", default=[], nargs="+", help="Selection on four-muon mass: first (or only) value is lower cut, second (optional) value is upper cut")
parser.add_argument("--fourmuonPtSel", default=[], nargs="+", help="Selection on four-muon pT: first (or only) value is lower cut, second (optional) value is upper cut")
parser.add_argument("--dimuonMassSelForFourMuon", default=[], nargs="+", help="Selection on dimuon mass in four-muon system: first (or only) value is lower cut, second (optional) value is upper cut")
parser.add_argument("--dimuonMassDiffSelForFourMuon", default=[], nargs="+", help="Selection on dimuon mass difference / mean in four-muon system: first (or only) value is *upper* cut, second (optional) value is *lower* cut")
parser.add_argument("--dimuonPtSelForFourMuon", default=[], nargs="+", help="Selection on dimuon pT in four-muon system: first (or only) value is lower cut, second (optional) value is upper cut")
parser.add_argument("--dimuonMassSelForFourMuonOSV", default=[], nargs="+", help="Selection on dimuon mass in four-muon system from overlapping SV: first (or only) value is lower cut, second (optional) value is upper cut")
parser.add_argument("--dimuonMassDiffSelForFourMuonOSV", default=[], nargs="+", help="Selection on dimuon mass difference / mean in four-muon system from overlapping SV: first (or only) value is *upper* cut, second (optional) value is *lower* cut")
parser.add_argument("--dimuonPtSelForFourMuonOSV", default=[], nargs="+", help="Selection on dimuon pT in four-muon system from overlapping SV: first (or only) value is lower cut, second (optional) value is upper cut")
parser.add_argument("--lxySel", default=[], nargs="+", help="Selection on lxy: first (or only) value is lower cut, second (optional) value is upper cut")
parser.add_argument("--lzSel", default=[], nargs="+", help="Selection on lz: first (or only) value is lower cut, second (optional) value is upper cut")
parser.add_argument("--lxySelForFourMuon", default=[], nargs="+", help="Selection on lxy: first (or only) value is lower cut, second (optional) value is upper cut")
parser.add_argument("--noMaterialVeto", default=False, action="store_true", help="Do not apply material vertex veto")
parser.add_argument("--noMuonIPSel", default=False, action="store_true", help="Do not apply selection on muon IP (Not applied at four-muon level)")
parser.add_argument("--noMuonHitSel", default=False, action="store_true", help="Do not apply selection on muon hits (Not applied at four-muon level)")
parser.add_argument("--noDiMuonAngularSel", default=False, action="store_true", help="Do not apply selection on dimuon angular variables")
parser.add_argument("--noFourMuonAngularSel", default=False, action="store_true", help="Do not apply selection on fourmuon angular variables")
parser.add_argument("--noFourMuonMassDiffSel", default=False, action="store_true", help="Do not apply selection on fourmuon invariant mass difference")
parser.add_argument("--noPreSel", default=False, action="store_true", help="Do not fill pre-selection/association histograms")
parser.add_argument("--noDiMuon", default=False, action="store_true", help="Do not fill dimuon histograms")
parser.add_argument("--noFourMuon", default=False, action="store_true", help="Do not fill four-muon histograms for four-muon systems")
parser.add_argument("--noFourMuonOSV", default=False, action="store_true", help="Do not fill four-muon histograms for four-muon systems from overlapping SVs")
parser.add_argument("--noSeed", default=[], nargs="+", help="Exclude L1 seeds from the acceptance")
parser.add_argument("--noHistos", default=False, action="store_true", help="Skip histogram filling")
parser.add_argument("--doGen", default=False, action="store_true", help="Fill generation information histograms")
args = parser.parse_args()

# Functions for selection
def applyDiMuonSelection(vec):
    selected = True
    if len(args.dimuonMassSel)>0:
        selected = selected and (vec.M() > float(args.dimuonMassSel[0]))
    if len(args.dimuonMassSel)>1:
        selected = selected and (vec.M() < float(args.dimuonMassSel[1]))
    if len(args.dimuonPtSel)>0:
        selected = selected and (vec.Pt() > float(args.dimuonPtSel[0]))
    if len(args.dimuonPtSel)>1:
        selected = selected and (vec.Pt() < float(args.dimuonPtSel[1]))
    if len(args.dimuonMassSidebandSel)>3:
        selected = selected and ((vec.M() > float(args.dimuonMassSidebandSel[0]) and vec.M() < float(args.dimuonMassSidebandSel[1])) or
                                 (vec.M() > float(args.dimuonMassSidebandSel[2]) and vec.M() < float(args.dimuonMassSidebandSel[3])))
    return selected

def applyFourMuonSelection(vec):
    selected = True
    if len(args.fourmuonMassSel)>0:
        selected = selected and (vec.M() > float(args.fourmuonMassSel[0]))
    if len(args.fourmuonMassSel)>1:
        selected = selected and (vec.M() < float(args.fourmuonMassSel[1]))
    if len(args.fourmuonPtSel)>0:
        selected = selected and (vec.Pt() > float(args.fourmuonPtSel[0]))
    if len(args.fourmuonPtSel)>1:
        selected = selected and (vec.Pt() < float(args.fourmuonPtSel[1]))
    return selected

def applyDiMuonSelectionForFourMuon(vecf,vecs):
    selected = True
    if len(args.dimuonMassSelForFourMuon)>0:
        selected = selected and (vecf.M() > float(args.dimuonMassSelForFourMuon[0])) and (vecs.M() > float(args.dimuonMassSelForFourMuon[0]))
    if len(args.dimuonMassSelForFourMuon)>1:
        selected = selected and (vecf.M() < float(args.dimuonMassSelForFourMuon[1])) and (vecs.M() < float(args.dimuonMassSelForFourMuon[1]))
    if len(args.dimuonPtSelForFourMuon)>0:
        selected = selected and (vecf.Pt() > float(args.dimuonPtSelForFourMuon[0])) and (vecs.Pt() > float(args.dimuonPtSelForFourMuon[0]))
    if len(args.dimuonPtSelForFourMuon)>1:
        selected = selected and (vecf.Pt() < float(args.dimuonPtSelForFourMuon[1])) and (vecs.Pt() < float(args.dimuonPtSelForFourMuon[1]))
    if len(args.dimuonMassDiffSelForFourMuon)>0:
        selected = selected and (abs(vecf.M()-vecs.M())*2.0/(vecf.M()+vecs.M()) < float(args.dimuonMassDiffSelForFourMuon[0]))
    if len(args.dimuonMassDiffSelForFourMuon)>1:
        selected = selected and (abs(vecf.M()-vecs.M())*2.0/(vecf.M()+vecs.M()) > float(args.dimuonMassDiffSelForFourMuon[1]))
    return selected

def applyDiMuonSelectionForFourMuonOSV(vecf,vecs):
    selected = True
    if len(args.dimuonMassSelForFourMuonOSV)>0:
        selected = selected and (vecf.M() > float(args.dimuonMassSelForFourMuonOSV[0])) and (vecs.M() > float(args.dimuonMassSelForFourMuonOSV[0]))
    if len(args.dimuonMassSelForFourMuonOSV)>1:
        selected = selected and (vecf.M() < float(args.dimuonMassSelForFourMuonOSV[1])) and (vecs.M() < float(args.dimuonMassSelForFourMuonOSV[1]))
    if len(args.dimuonPtSelForFourMuonOSV)>0:
        selected = selected and (vecf.Pt() > float(args.dimuonPtSelForFourMuonOSV[0])) and (vecs.Pt() > float(args.dimuonPtSelForFourMuonOSV[0]))
    if len(args.dimuonPtSelForFourMuonOSV)>1:
        selected = selected and (vecf.Pt() < float(args.dimuonPtSelForFourMuonOSV[1])) and (vecs.Pt() < float(args.dimuonPtSelForFourMuonOSV[1]))
    if len(args.dimuonMassDiffSelForFourMuonOSV)>0:
        selected = selected and (abs(vecf.M()-vecs.M())*2.0/(vecf.M()+vecs.M()) < float(args.dimuonMassDiffSelForFourMuonOSV[0]))
    if len(args.dimuonMassDiffSelForFourMuonOSV)>1:
        selected = selected and (abs(vecf.M()-vecs.M())*2.0/(vecf.M()+vecs.M()) > float(args.dimuonMassDiffSelForFourMuonOSV[1]))
    return selected

def applyLxySelection(lxy):
    selected = True
    if len(args.lxySel)>0:
        selected = selected and (lxy > float(args.lxySel[0]))
    if len(args.lxySel)>1:
        selected = selected and (lxy < float(args.lxySel[1]))
    return selected

def applyLzSelection(lz):
    selected = True
    if len(args.lzSel)>0:
        selected = selected and (lxy > float(args.lzSel[0]))
    if len(args.lzSel)>1:
        selected = selected and (lxy < float(args.lzSel[1]))
    return selected

def applyFourMuonLxySelection(lxymin,lxymax):
    selected = True
    if len(args.lxySelForFourMuon)>0:
        selected = selected and (lxymin > float(args.lxySel[0]))
    if len(args.lxySelForFourMuon)>1:
        selected = selected and (lxymax < float(args.lxySel[1]))
    return selected

# Muon type
def muonType(isGlobal, isTracker, isStandAlone):
    if isGlobal and isTracker:
        return 0.5
    elif isGlobal and not isTracker:
        return 1.5
    elif not isGlobal and isTracker:
        return 2.5
    elif not isGlobal and not isTracker and isStandAlone:
        return 3.5
    elif not isGlobal and not isTracker and not isStandAlone:
        return 4.5

# Isolation category
def dimuonIsoCategory(iso1, pt1, iso2, pt2):
    isocat = -99
    #if (iso1>8.0 and iso2>8.0 ):
    if ((iso1>8.0 and iso1/pt1 > 0.2) and (iso2>8.0 and iso2/pt2 > 0.2) ):
        isocat = 0
    elif ((iso1<8.0 or iso1/pt1 < 0.2) and (iso2<8.0 or iso2/pt2 < 0.2) ):
        isocat = 2
    else:
        isocat = 1 
    return isocat

# Evaluate L1 with the possibility of excluding one
def evaluateSeeds(tree, seedList):
    for seed in seedList: 
        if eval('tree.'+seed):
            return True
    return False

# Muon IP sign
def getIPSign(mphi, dmphi):
    dxydir = ROOT.TVector3(-ROOT.TMath.Sin(mphi), ROOT.TMath.Cos(mphi), 0.0)
    dmudir = ROOT.TVector3(ROOT.TMath.Cos(dmphi), ROOT.TMath.Sin(dmphi), 0.0)
    if (dxydir*dmudir > 0):
        return 1.0
    else:
        return -1.0

# Get Signal normalization weight
def getweight(era, ngen, frac=1.0, xsec=1000):
    if era=="2022":
        return frac*8.077046947*xsec/ngen
    if era=="2022postEE":
        return frac*26.982330931*xsec/ngen
    if era=="2023":
        return frac*17.060484313*xsec/ngen
    if era=="2023BPix":
        return frac*9.525199061*xsec/ngen


indir  = args.inDir.replace("/ceph/cms","")
outdir = args.outDir
if args.outSuffix!="":
    outdir = outdir+"_"+args.outSuffix
if not os.path.exists(outdir):
    os.makedirs(outdir)

applyMaterialVeto = not args.noMaterialVeto
applyMuonIPSel = not args.noMuonIPSel
applyDiMuonAngularSel = not args.noDiMuonAngularSel
applyMuonHitSel = not args.noMuonHitSel
applyFourMuonAngularSel = not args.noFourMuonAngularSel
applyFourMuonMassDiffSel = not args.noFourMuonMassDiffSel

isData = args.data
if "Data" in args.inSample:
    isData = True
removeDuplicates = args.removeDuplicates
if not isData:
    removeDuplicates = False
MUON_MASS = 0.10566

if args.signal:
    isData = False

skimFraction = float(args.partialUnblindingFraction)
skimEvents = args.partialUnblinding
if not isData:
    skimEvents = False
rndm_partialUnblinding = ROOT.TRandom3(42)

files = []
prependtodir = ""
if not args.condor:
    prependtodir = "/ceph/cms"
else:
    prependtodir = "davs://redirector.t2.ucsd.edu:1095"
if not args.condor:
    if args.inFile!="*" and args.inSample!="*":
        thisfile="output_%s_%s_%s.root"%(args.inSample,args.year,args.inFile)
        if os.path.isfile("/ceph/cms%s/%s"%(indir,thisfile)):
            files.append("%s%s/%s"%(prependtodir,indir,thisfile))
    elif args.inSample!="*":
        for f in os.listdir("/ceph/cms%s"%indir):
            if ("output_%s_%s_"%(args.inSample,args.year) in f) and os.path.isfile("/ceph/cms%s/%s"%(indir,f)):
                files.append("%s%s/%s"%(prependtodir,indir,f))
    else:
        for f in os.listdir("/ceph/cms%s"%indir):
            if (args.year in f) and (".root" in f) and os.path.isfile("/ceph/cms%s/%s"%(indir,f)):
                files.append("%s%s/%s"%(prependtodir,indir,f))
else:
    os.system('xrdfs redirector.t2.ucsd.edu:1095 ls %s > filein.txt'%indir)
    fin = open("filein.txt","r")
    for f in fin.readlines():
        f = f.strip("\n")
        if args.inFile!="*" and args.inSample!="*":
            thisfile="output_%s_%s_%s.root"%(args.inSample,args.year,args.inFile)
            if thisfile in f:
                files.append("%s/%s"%(prependtodir,thisfile))
        elif args.inSample!="*":
            if "output_%s_%s_"%(args.inSample,args.year) in f:
                files.append("%s/%s"%(prependtodir,f))
        else:
            if (args.year in f) and (".root" in f):
                files.append("%s/%s"%(prependtodir,f))
    fin.close()
    os.system('rm -f filein.txt')
print("Found {} files matching criteria".format(len(files)))
print(files)

index = int(args.splitIndex)
pace  = int(args.splitPace)

# Inputs:
t = ROOT.TChain("tout")
for f in files:
    t.Add(f)

# MC normalization (Need to integrate everything for all signals that we may produce)
ncounts = 1
efilter = 1.0
lumiweight = 1.0
sampleTag = args.inSample.replace('Signal_', '').split('_202')[0]
if not isData and args.weightMC:
    counts = ROOT.TH1F("totals", "", 1, 0, 1)
    print("Simulations: Getting counts")
    if "HTo2ZdTo2mu2x" in sampleTag:
        for _,f in enumerate(files):
            if not args.condor:
                f_ = ROOT.TFile.Open(f.replace('davs://redirector.t2.ucsd.edu:1095//', '/ceph/cms/'))
            else:
                f_ = ROOT.TFile.Open(f)
            h_ = f_.Get("counts").Clone("Clone_{}".format(_))
            counts.Add(h_)
            f_.Close()
        ncounts = counts.GetBinContent(1) 
        with open('data/hahm-request.csv') as mcinfo:
            reader = csv.reader(mcinfo, delimiter=',')
            for row in reader:
                if sampleTag in row[0]:
                    efilter = float(row[-1])
                    break
    if 'BToPhi' in sampleTag:
        for _,f in enumerate(files):
            if not args.condor:
                f_ = ROOT.TFile.Open(f.replace('davs://redirector.t2.ucsd.edu:1095//', '/ceph/cms/'))
            else:
                f_ = ROOT.TFile.Open(f)
            h_ = f_.Get("counts").Clone("Clone_{}".format(_))
            counts.Add(h_)
            f_.Close()
        ncounts = counts.GetBinContent(1)
        with open('data/BToPhi-request.csv') as mcinfo:
            reader = csv.reader(mcinfo, delimiter=',')
            for row in reader:
                if sampleTag in row[0]:
                    efilter = float(row[-1])
                    break
    if "ScenB2" in sampleTag:
        for _,f in enumerate(files):
            if not args.condor:
                f_ = ROOT.TFile.Open(f.replace('davs://redirector.t2.ucsd.edu:1095//', '/ceph/cms/'))
            else:
                f_ = ROOT.TFile.Open(f)
            h_ = f_.Get("counts").Clone("Clone_{}".format(_))
            counts.Add(h_)
            f_.Close()
        ncounts = counts.GetBinContent(1)
        efilter = 1.0
    if "2022postEE" in f:
        lumiweight = getweight("2022postEE", ncounts/efilter, 0.1)
    elif "2022" in f:
        lumiweight = getweight("2022", ncounts/efilter, 0.1)
    elif "2023" in f:
        lumiweight = getweight("2023", ncounts/efilter, 0.1)
    elif "2023BPix" in f:
        lumiweight = getweight("2023BPix", ncounts/efilter, 0.1)
    if "ScenB1" in sampleTag:
        ncounts = 300000.0
        efilter = 1.0
        lumiweight = 0.1*(8.077046947 + 26.982330931)*1000.0/(ncounts/efilter)
    print("Total number of counts: {}".format(ncounts))
    print("Filter efficiency (generation): {}".format(efilter))
    print("Lumiweight: {}".format(lumiweight))

# Histograms:
h1d = dict()
variable1d = dict()
#
h2d = dict()
variable2d = dict()
#
h1d,variable1d,h2d,variable2d = histDefinition.histInitialization(not(args.noPreSel),not(args.noDiMuon),not(args.noFourMuon),not(args.noFourMuonOSV))
if args.noHistos:
    for key in h1d.keys():
        h1d[key] = []
    for key in h2d.keys():
        h2d[key] = []
###
for cat in h1d.keys():
    for h in h1d[cat]:
        if isData:
            h.Sumw2(ROOT.kFALSE)
        else:
            h.Sumw2()
for cat in h2d.keys():
    for h in h2d[cat]:
        if isData:
            h.Sumw2(ROOT.kFALSE)
        else:
            h.Sumw2()
###
L1seeds = []
branch_list = t.GetListOfBranches()
for branch in branch_list:
    if branch.GetName().startswith('L1'):
        L1seeds.append(branch.GetName())
effL1seeds = [s for s in L1seeds if s not in args.noSeed]
print('List of L1 seeds: ', L1seeds)
if args.noSeed:
    print('but excluding: ', args.noSeed)
###

elist = [] # for duplicate removal
print("Starting loop over %d events"%t.GetEntries())
firste = 0
laste  = t.GetEntries()
print(args.inSample, laste)
if index>=0:
    firste = index*pace
    laste  = min((index+1)*pace,t.GetEntries())
if firste >= t.GetEntries():
    exit()

## Init RooDataSets
# Dimuon binning
lxybins = [0.0, 0.2, 1.0, 2.4, 3.1, 7.0, 11.0, 16.0, 70.0]
lxystrs = [str(l).replace('.', 'p') for l in lxybins]
lxybinlabel = ["lxy{}to{}".format(lxystrs[l], lxystrs[l+1]) for l in range(0, len(lxystrs)-1)]
ptcut = 25. # Tried 25, 50 and 100
dphisvcut = 0.02 # Last Run 2 value
# Variables
mfit = ROOT.RooRealVar("mfit", "mfit", 0.4, 140.0)
m4fit = ROOT.RooRealVar("m4fit", "m4fit", 0.4, 140.0)
roow = ROOT.RooRealVar("roow", "roow", -10000.0, 10000.0)
roow4 = ROOT.RooRealVar("roow", "roow", -10000.0, 10000.0)
roods = {}
catmass = {}
dbins = []
rooweight = float(args.rooWeight) # Weight for RooDataset
# Categories
dbins.append("FourMu_sep") # 4mu, multivertex
dbins.append("FourMu_osv") # 4mu, 4mu-vertex
dbins.append("Dimuon_full_inclusive") # Dimuons excluded from categorization
for label in lxybinlabel:
    dbins.append("Dimuon_"+label+"_inclusive")
    dbins.append("Dimuon_"+label+"_iso0_ptlow")
    dbins.append("Dimuon_"+label+"_iso0_pthigh")
    dbins.append("Dimuon_"+label+"_iso1_ptlow")
    dbins.append("Dimuon_"+label+"_iso1_pthigh")
    dbins.append("Dimuon_"+label+"_non-pointing") # non-pointing
dbins.append("Dimuon_excluded") # Dimuons excluded from categorization
#
for dbin in dbins:
    dname = "d_" + dbin
    if 'Dimuon' in dname:
        catmass[dbin] = ROOT.TH1F(dname + "_rawmass","; m_{#mu#mu} [GeV]; Events / 0.01 GeV",15000, 0., 150.)
        roods[dbin] = ROOT.RooDataSet(dname,dname,ROOT.RooArgSet(mfit,roow),"roow")
    else:
        catmass[dbin] = ROOT.TH1F(dname + "_rawmass","; m_{4#mu} [GeV]; Events / 0.01 GeV",15000, 0., 150.)
        roods[dbin] = ROOT.RooDataSet(dname,dname,ROOT.RooArgSet(m4fit,roow4),"roow")


print("From event %d to event %d"%(firste,laste))
for e in range(firste,laste):
    t.GetEntry(e)
    #if e%1000==0:
    #    print("At entry %d"%e)
    if len(args.noSeed) > 0:
        passL1 = evaluateSeeds(t, effL1seeds)
        if not passL1:
            continue
    if removeDuplicates:
        # Run, event number & lumisection
        eid  = t.evtn
        run  = t.run
        lumi = t.lumi
        if (run,lumi,eid) in elist:
            continue
        else:
            elist.append((run,lumi,eid))
    if skimEvents and skimFraction>0.0:
        if rndm_partialUnblinding.Rndm() > skimFraction:
            continue
    ### As advised in LUM POG TWiki 
    ### (https://twiki.cern.ch/twiki/bin/view/CMS/LumiRecommendationsRun3),
    ### exclude runs 359571 + 359661
    if isData and t.run==359571 or t.run==359661:
            continue

    # Event info
    for h in h1d["event"]:
        tn = h.GetName()
        h.Fill(eval(variable1d[h.GetName()]), lumiweight)

    # Gen info
    dmugen = []
    dmumot = []
    if not isData:
        nGEN = len(t.GenPart_pdgId)
        for i in range(nGEN):
            if abs(t.GenPart_pdgId[i]) != 13:            
                continue
            isResonance = False
            for j in range(i+1, nGEN):
                if i==j:
                    continue
                if t.GenPart_motherIndex[i] != t.GenPart_motherIndex[j] or t.GenPart_pdgId[i]*t.GenPart_pdgId[j] > 0:
                    continue
                dmugen.append(i)
                dmugen.append(j)
                isResonance = True
                for k in range(0, nGEN):
                    if t.GenPart_index[k] == t.GenPart_motherIndex[i]:
                        gvec = ROOT.TLorentzVector()
                        gvec.SetPtEtaPhiM(t.GenPart_pt[k], t.GenPart_eta[k], t.GenPart_phi[k], t.GenPart_m[k])
                        dmumot.append(gvec)
                        break 
            '''
            if isResonance:
                for h in h1d["genmu"]:
                    tn = h.GetName()
                    h.Fill(eval(variable1d[h.GetName()]), lumiweight)
            '''
        for g,gp in enumerate(dmumot):
            lxygen = t.GenPart_lxy[dmugen[2*g]] 
            if t.GenPart_motherPdgId[dmugen[2*g]]==443:
                for h in h1d["jpsi"]:
                    tn = h.GetName()
                    h.Fill(eval(variable1d[h.GetName()]), lumiweight)
            for h in h1d["llp"]:
                tn = h.GetName()
                h.Fill(eval(variable1d[h.GetName()]), lumiweight)

    # Loop over SVs
    nSV = len(t.SV_index)
    if nSV<1:
        continue
    nSVsel = 0
    for v in range(nSV):
        if args.noPreSel:
            break
        if not t.SV_selected[v]:
            continue
        if applyMaterialVeto and (t.SV_onModuleWithinUnc[v] or (abs(t.SV_minDistanceFromDet_x[v]) < 0.81 and abs(t.SV_minDistanceFromDet_y[v]) < 3.24 and abs(t.SV_minDistanceFromDet_z[v]) < 0.0145)):
            continue
        nSVsel = nSVsel+1
        lxy = t.SV_lxy[v]
        for h in h1d["svsel"]:
            tn = h.GetName()
            h.Fill(eval(variable1d[h.GetName()]), lumiweight)
        for h in h2d["svsel"]:
            tn = h.GetName()
            h.Fill(eval(variable2d[h.GetName()][0]),eval(variable2d[h.GetName()][1]), lumiweight)
    nSVs = nSVsel
    for h in h1d["nsvsel"]:
        if args.noPreSel:
            break
        tn = h.GetName()
        h.Fill(eval(variable1d[h.GetName()]), lumiweight)

    # Loop over muons
    nMu = len(t.Muon_selected)
    nMuSel = 0
    nMuAss = 0
    nMuAssOverlap = 0
    muselidxs = []
    for m in range(nMu):
        if not t.Muon_selected[m]:
            continue
        nMuSel = nMuSel+1
        if t.Muon_bestAssocSVOverlapIdx[m]>-1:
            nMuAss = nMuAss+1
            nMuAssOverlap = nMuAssOverlap+1
        elif t.Muon_bestAssocSVIdx[m]>-1:
            nMuAss = nMuAss+1
        else:
            continue
        muselidxs.append(m)
        if args.noPreSel:
            continue
        for h in h1d["muon"]:
            if args.noPreSel:
                break
            tn = h.GetName()
            h.Fill(eval(variable1d[h.GetName()]), lumiweight)
        for h in h2d["muon"]:
            if args.noPreSel:
                break
            tn = h.GetName()
            h.Fill(eval(variable2d[h.GetName()][0]),eval(variable2d[h.GetName()][1]), lumiweight)
    for h in h1d["nmuon"]:
        if args.noPreSel:
            break
        tn = h.GetName()
        h.Fill(eval(variable1d[h.GetName()]), lumiweight)
    # Select events witb at least two muons associated to a SV
    if nMuAss<2:
        continue

    # Muon pairing
    dmuvec = []
    dmu_muvecdp = []
    dmuidxs = []
    svvec = []
    svidx = []
    dmuvec_osv = []
    dmu_muvecdp_osv = []
    dmuidxs_osv = []
    osvvec = []
    osvidx = []
    qmuvec_osv = []
    qmuidxs_osv = []
    qmuidxs_osv_sel = []
    qmu_dmuvec_osv = []
    qmu_muvecdp_osv = []
    osvvec_qmu = []
    osvidx_qmu = []
    qmuvec = []
    qmuidxs = []
    qmuidxs_sel = []
    qmuidxsminlxy = []
    qmuidxsmaxlxy = []
    qmu_muvecdp = []
    qmu_muvecdpminlxy = []
    qmu_muvecdpmaxlxy = []
    qmu_dmuvecminlxy = []
    qmu_dmuvecmaxlxy = []
    qmu_dmuvecdpminlxy = []
    qmu_dmuvecdpmaxlxy = []
    svvecminlxy_qmu = []
    svidxminlxy_qmu = []
    svvecmaxlxy_qmu = []
    svidxmaxlxy_qmu = []
    for m in muselidxs:
        chg   = t.Muon_ch[m]
        ovidx = -1
        vidx  = -1
        ovpos = -1
        vpos = -1
        # First, identify muons from overlapping SVs
        if t.Muon_bestAssocSVOverlapIdx[m]>-1:
            if applyMaterialVeto and (t.SV_onModuleWithinUnc[t.SVOverlap_vtxIdxs[t.Muon_bestAssocSVOverlapIdx[m]][0]] or (abs(t.SV_minDistanceFromDet_x[t.SVOverlap_vtxIdxs[t.Muon_bestAssocSVOverlapIdx[m]][0]]) < 0.81 and abs(t.SV_minDistanceFromDet_y[t.SVOverlap_vtxIdxs[t.Muon_bestAssocSVOverlapIdx[m]][0]]) < 3.24 and abs(t.SV_minDistanceFromDet_z[t.SVOverlap_vtxIdxs[t.Muon_bestAssocSVOverlapIdx[m]][0]]) < 0.0145)):
                continue
            ovidx = t.Muon_bestAssocSVOverlapIdx[m]
            ovpos = ovidx
        # Then, identify muons from non-overlapping SVs
        elif t.Muon_bestAssocSVIdx[m]>-1:
            vidx = t.Muon_bestAssocSVIdx[m]
            for v in range(len(t.SV_index)):
                if applyMaterialVeto and (t.SV_onModuleWithinUnc[v] or (abs(t.SV_minDistanceFromDet_x[v]) < 0.81 and abs(t.SV_minDistanceFromDet_y[v]) < 3.24 and abs(t.SV_minDistanceFromDet_z[v]) < 0.0145)):
                    continue
                if t.SV_index[v]==vidx:
                    vpos = v
                    break
        # Loop over muons, and do pairing
        for mm in muselidxs:
            if mm==m:
                continue
            if abs(chg+t.Muon_ch[mm])>0:
                continue
            # First, identify muon pairs from overlapping SVs
            if ovidx>-1 and t.Muon_bestAssocSVOverlapIdx[mm]==ovidx and ovpos>-1:
                if not (m in dmuidxs_osv or mm in dmuidxs_osv):
                    dmuvec_osv.append(t.Muon_vec[m])
                    dmuvec_osv[len(dmuvec_osv)-1] = dmuvec_osv[len(dmuvec_osv)-1] + t.Muon_vec[mm]
                    dmuidxs_osv.append(m)
                    dmuidxs_osv.append(mm)
                    dmu_muvecdp_osv.append(ROOT.TLorentzVector())
                    dmu_muvecdp_osv[-1].SetPtEtaPhiM(t.Muon_pt[m],t.Muon_eta[m],t.Muon_phi[m],MUON_MASS)
                    dmu_muvecdp_osv.append(ROOT.TLorentzVector())
                    dmu_muvecdp_osv[-1].SetPtEtaPhiM(t.Muon_pt[mm],t.Muon_eta[mm],t.Muon_phi[mm],MUON_MASS)
                    osvvec.append(ROOT.TVector3())
                    osvvec[len(osvvec)-1].SetXYZ(t.SVOverlap_x[ovpos]-t.PV_x, t.SVOverlap_y[ovpos]-t.PV_y, t.SVOverlap_z[ovpos]-t.PV_z)
                    osvidx.append(ovpos)
            # Then, identify muon pairs from non-overlapping SVs
            elif vidx>-1 and t.Muon_bestAssocSVIdx[mm]==vidx and vpos>-1:
                if not (m in dmuidxs or mm in dmuidxs):
                    dmuvec.append(t.Muon_vec[m])
                    dmuvec[len(dmuvec)-1] = dmuvec[len(dmuvec)-1] + t.Muon_vec[mm]
                    dmuidxs.append(m)
                    dmuidxs.append(mm)
                    dmu_muvecdp.append(ROOT.TLorentzVector())
                    dmu_muvecdp[-1].SetPtEtaPhiM(t.Muon_pt[m],t.Muon_eta[m],t.Muon_phi[m],MUON_MASS)
                    dmu_muvecdp.append(ROOT.TLorentzVector())
                    dmu_muvecdp[-1].SetPtEtaPhiM(t.Muon_pt[mm],t.Muon_eta[mm],t.Muon_phi[mm],MUON_MASS)
                    svvec.append(ROOT.TVector3())
                    svvec[len(svvec)-1].SetXYZ(t.SV_x[vpos]-t.PV_x, t.SV_y[vpos]-t.PV_y, t.SV_z[vpos]-t.PV_z)
                    svidx.append(vpos)

    # If multiple muon pairs from overlapping SVs are found, create four-muon system
    if len(dmuidxs_osv)>2 and float(len(dmuidxs_osv))/float(len(set(osvidx)))>2:
        for m in range(len(dmuidxs_osv)):
            if m%2>0:
                continue
            if len(qmuidxs_osv)>3:
                break
            for mm in range(m+2,len(dmuidxs_osv)):
                if mm%2>0:
                    continue
                if len(qmuidxs_osv)>3:
                    break
                if osvidx[int(m/2)]==osvidx[int(mm/2)]:
                    if (dmuidxs_osv[m] in qmuidxs_osv) or (dmuidxs_osv[mm] in qmuidxs_osv):
                        continue
                    else:
                        qmuidxs_osv.append(dmuidxs_osv[m])
                        qmuidxs_osv.append(dmuidxs_osv[m+1])
                        qmuidxs_osv.append(dmuidxs_osv[mm])
                        qmuidxs_osv.append(dmuidxs_osv[mm+1])
                        qmu_muvecdp_osv.append(dmu_muvecdp_osv[m])
                        qmu_muvecdp_osv.append(dmu_muvecdp_osv[m+1])
                        qmu_muvecdp_osv.append(dmu_muvecdp_osv[mm])
                        qmu_muvecdp_osv.append(dmu_muvecdp_osv[mm+1])
                        qmuvec_osv.append(dmuvec_osv[int(m/2)])
                        qmuvec_osv[len(qmuvec_osv)-1] = qmuvec_osv[len(qmuvec_osv)-1]+dmuvec_osv[int(mm/2)]
                        qmu_dmuvec_osv.append(dmuvec_osv[int(m/2)])
                        qmu_dmuvec_osv.append(dmuvec_osv[int(mm/2)])
                        osvvec_qmu.append(ROOT.TVector3())
                        osvvec_qmu[len(osvvec_qmu)-1].SetXYZ(t.SVOverlap_x[osvidx[int(m/2)]]-t.PV_x, t.SVOverlap_y[osvidx[int(m/2)]]-t.PV_y, t.SVOverlap_z[osvidx[int(m/2)]]-t.PV_z)
                        osvidx_qmu.append(osvidx[int(m/2)])

    dmuidxs_all = dmuidxs_osv+dmuidxs
    dmuvec_all = dmuvec_osv+dmuvec
    dmu_muvecdp_all = dmu_muvecdp_osv+dmu_muvecdp
    svidx_all = osvidx+svidx
    svvec_all = osvvec+svvec
    # If multiple muon pairs are found not from overlapping SVs, create four-muon system from non-overlapping SVs
    if len(dmuidxs_all)>=4:
        for m in range(len(dmuidxs_all)):
            if m%2>0:
                continue
            if len(qmuidxs)>3:
                break
            if m in dmuidxs_osv and m in qmuidxs_osv and len(qmuidxs_osv)>3:
                continue
            for mm in range(m+2,len(dmuidxs_all)):
                if mm%2>0:
                    continue
                if len(qmuidxs)>3:
                    break
                if mm in dmuidxs_osv and mm in qmuidxs_osv and len(qmuidxs_osv)>3:
                    continue
                if svidx_all[int(m/2)] in osvidx_qmu or svidx_all[int(mm/2)] in osvidx_qmu:
                    continue
                else:
                    qmuidxs.append(dmuidxs_all[m])
                    qmuidxs.append(dmuidxs_all[m+1])
                    qmuidxs.append(dmuidxs_all[mm])
                    qmuidxs.append(dmuidxs_all[mm+1])
                    qmu_muvecdp.append(dmu_muvecdp_all[m])
                    qmu_muvecdp.append(dmu_muvecdp_all[m+1])
                    qmu_muvecdp.append(dmu_muvecdp_all[mm])
                    qmu_muvecdp.append(dmu_muvecdp_all[mm+1])
                    qmuvec.append(dmuvec_all[int(m/2)])
                    qmuvec[len(qmuvec)-1] = qmuvec[len(qmuvec)-1]+dmuvec_all[int(mm/2)]
                    if svvec_all[int(m/2)].Perp() < svvec_all[int(mm/2)].Perp():
                        svvecminlxy_qmu.append(svvec_all[int(m/2)])
                        svidxminlxy_qmu.append(svidx_all[int(m/2)])
                        qmu_dmuvecminlxy.append(dmuvec_all[int(m/2)])
                        svvecmaxlxy_qmu.append(svvec_all[int(mm/2)])
                        svidxmaxlxy_qmu.append(svidx_all[int(mm/2)])
                        qmu_dmuvecmaxlxy.append(dmuvec_all[int(mm/2)])
                        qmu_dmuvecdpminlxy.append(dmu_muvecdp_all[m]+dmu_muvecdp_all[m+1])
                        qmu_dmuvecdpmaxlxy.append(dmu_muvecdp_all[mm]+dmu_muvecdp_all[mm+1])
                        qmuidxsminlxy.append(dmuidxs_all[m])
                        qmuidxsminlxy.append(dmuidxs_all[m+1])
                        qmuidxsmaxlxy.append(dmuidxs_all[mm])
                        qmuidxsmaxlxy.append(dmuidxs_all[mm+1])
                        qmu_muvecdpminlxy.append(dmu_muvecdp_all[m])
                        qmu_muvecdpminlxy.append(dmu_muvecdp_all[m+1])
                        qmu_muvecdpmaxlxy.append(dmu_muvecdp_all[mm])
                        qmu_muvecdpmaxlxy.append(dmu_muvecdp_all[mm+1])
                    else:
                        svvecminlxy_qmu.append(svvec_all[int(mm/2)])
                        svidxminlxy_qmu.append(svidx_all[int(mm/2)])
                        qmu_dmuvecminlxy.append(dmuvec_all[int(mm/2)])
                        svvecmaxlxy_qmu.append(svvec_all[int(m/2)])
                        svidxmaxlxy_qmu.append(svidx_all[int(m/2)])
                        qmu_dmuvecmaxlxy.append(dmuvec_all[int(m/2)])
                        qmu_dmuvecdpminlxy.append(dmu_muvecdp_all[mm]+dmu_muvecdp_all[mm+1])
                        qmu_dmuvecdpmaxlxy.append(dmu_muvecdp_all[m]+dmu_muvecdp_all[m+1])
                        qmuidxsminlxy.append(dmuidxs_all[mm])
                        qmuidxsminlxy.append(dmuidxs_all[mm+1])
                        qmuidxsmaxlxy.append(dmuidxs_all[m])
                        qmuidxsmaxlxy.append(dmuidxs_all[m+1])
                        qmu_muvecdpminlxy.append(dmu_muvecdp_all[mm])
                        qmu_muvecdpminlxy.append(dmu_muvecdp_all[mm+1])
                        qmu_muvecdpmaxlxy.append(dmu_muvecdp_all[m])
                        qmu_muvecdpmaxlxy.append(dmu_muvecdp_all[m+1])


    ### Scan analysis initialization 
    # Cat selection:
    filledcat4musep = False
    filledcat4muosv = False
    filledcat2mu = False


    # Apply selections and fill histograms for four-muon systems from non-overlapping SVs
    selqmusvidxs = []
    selqmuvecs = []
    mindrmm, mindpmm, mindemm, mindedpmm, mina3dmm = 1e6, 1e6, 1e6, 1e6, 1e6
    maxdrmm, maxdpmm, maxdemm, maxdedpmm, maxa3dmm = -1., -1., -1., -1., -1.
    mindrmmu, mindpmmu, mindemmu, mindedpmmu, mina3dmmu = 1e6, 1e6, 1e6, 1e6, 1e6
    maxdrmmu, maxdpmmu, maxdemmu, maxdedpmmu, maxa3dmmu = -1., -1., -1., -1., -1.
    for vn,v in enumerate(qmuvec):
        if args.noFourMuon:
            break
        if not applyDiMuonSelectionForFourMuon(qmu_dmuvecminlxy[vn], qmu_dmuvecmaxlxy[vn]):
            continue
        if not applyFourMuonSelection(v):
            continue
        minlxy  = svvecminlxy_qmu[vn].Perp()
        maxlxy  = svvecmaxlxy_qmu[vn].Perp()
        minl3d  = svvecminlxy_qmu[vn].Mag()
        maxl3d  = svvecmaxlxy_qmu[vn].Mag()
        if not applyFourMuonLxySelection(minlxy,maxlxy):
            continue
        qmuidxs_sel.append(qmuidxs[vn*4])
        qmuidxs_sel.append(qmuidxs[vn*4+1])
        qmuidxs_sel.append(qmuidxs[vn*4+2])
        qmuidxs_sel.append(qmuidxs[vn*4+3])
        selqmusvidxs.append(svidxminlxy_qmu[vn])
        selqmusvidxs.append(svidxmaxlxy_qmu[vn])
        selqmuvecs.append(vn)
        mass    = v.M()
        minmass = min(qmu_dmuvecminlxy[vn].M(), qmu_dmuvecmaxlxy[vn].M())
        maxmass = max(qmu_dmuvecminlxy[vn].M(), qmu_dmuvecmaxlxy[vn].M())
        avgmass = 0.5*(minmass+maxmass)
        reldmass= (maxmass-minmass)/avgmass
        pt    = v.Pt()
        minpt = min(qmu_dmuvecminlxy[vn].Pt(), qmu_dmuvecmaxlxy[vn].Pt())
        maxpt = max(qmu_dmuvecminlxy[vn].Pt(), qmu_dmuvecmaxlxy[vn].Pt())
        nhitsbeforesvtotal = t.Muon_nhitsbeforesv[qmuidxs[vn*4]] + t.Muon_nhitsbeforesv[qmuidxs[vn*4+1]] + t.Muon_nhitsbeforesv[qmuidxs[vn*4+2]] + t.Muon_nhitsbeforesv[qmuidxs[vn*4+3]]
        for m in range(vn*4,vn*4+4):
            if not m%2==0:
                continue
            drmm = t.Muon_vec[qmuidxs[m]].DeltaR(t.Muon_vec[qmuidxs[m+1]])
            dpmm = abs(t.Muon_vec[qmuidxs[m]].DeltaPhi(t.Muon_vec[qmuidxs[m+1]]))
            demm = abs(t.Muon_vec[qmuidxs[m]].Eta()-t.Muon_vec[qmuidxs[m+1]].Eta())
            dedpmm = 1e6
            if dpmm>0.0:
                dedpmm = demm/dpmm
            a3dmm = abs(t.Muon_vec[qmuidxs[m]].Angle(t.Muon_vec[qmuidxs[m+1]].Vect()))
            if drmm<mindrmm:
                mindrmm = drmm
            if drmm>maxdrmm:
                maxdrmm = drmm
            if dpmm<mindpmm:
                mindpmm = dpmm
            if dpmm>maxdpmm:
                maxdpmm = dpmm
            if demm<mindemm:
                mindemm = demm
            if demm>maxdemm:
                maxdemm = demm
            if dedpmm<mindedpmm:
                mindedpmm = dedpmm
            if dedpmm>maxdedpmm:
                maxdedpmm = dedpmm
            if a3dmm<mina3dmm:
                mina3dmm = a3dmm
            if a3dmm>maxa3dmm:
                maxa3dmm = a3dmm
            #
            drmmu = qmu_muvecdp[m].DeltaR(qmu_muvecdp[m+1])
            dpmmu = abs(qmu_muvecdp[m].DeltaPhi(qmu_muvecdp[m+1]))
            demmu = abs(qmu_muvecdp[m].Eta()-qmu_muvecdp[m+1].Eta())
            dedpmmu = 1e6
            if dpmmu>0.0:
                dedpmmu = demmu/dpmmu
            a3dmmu = abs(qmu_muvecdp[m].Angle(qmu_muvecdp[m+1].Vect()))
            if drmmu<mindrmmu:
                mindrmmu = drmmu
            if drmmu>maxdrmmu:
                maxdrmmu = drmmu
            if dpmmu<mindpmmu:
                mindpmmu = dpmmu
            if dpmmu>maxdpmmu:
                maxdpmmu = dpmmu
            if demmu<mindemmu:
                mindemmu = demmu
            if demmu>maxdemmu:
                maxdemmu = demmu
            if dedpmmu<mindedpmmu:
                mindedpmmu = dedpmmu
            if dedpmmu>maxdedpmmu:
                maxdedpmmu = dedpmmu
            if a3dmmu<mina3dmmu:
                mina3dmmu = a3dmmu
            if a3dmmu>maxa3dmmu:
                maxa3dmmu = a3dmmu
        mindphisv = min(abs(qmu_dmuvecminlxy[vn].Vect().DeltaPhi(svvecminlxy_qmu[vn])),abs(qmu_dmuvecmaxlxy[vn].Vect().DeltaPhi(svvecmaxlxy_qmu[vn])))
        maxdphisv = max(abs(qmu_dmuvecminlxy[vn].Vect().DeltaPhi(svvecminlxy_qmu[vn])),abs(qmu_dmuvecmaxlxy[vn].Vect().DeltaPhi(svvecmaxlxy_qmu[vn])))
        mindetasv = min(abs(qmu_dmuvecminlxy[vn].Vect().Eta()-svvecminlxy_qmu[vn].Eta()),abs(qmu_dmuvecmaxlxy[vn].Vect().Eta()-svvecmaxlxy_qmu[vn].Eta()))
        maxdetasv = max(abs(qmu_dmuvecminlxy[vn].Vect().Eta()-svvecminlxy_qmu[vn].Eta()),abs(qmu_dmuvecmaxlxy[vn].Vect().Eta()-svvecmaxlxy_qmu[vn].Eta()))
        mindetadphisv = min(abs(qmu_dmuvecminlxy[vn].Vect().Eta()-svvecminlxy_qmu[vn].Eta())/abs(qmu_dmuvecminlxy[vn].Vect().DeltaPhi(svvecminlxy_qmu[vn])),
                            abs(qmu_dmuvecmaxlxy[vn].Vect().Eta()-svvecmaxlxy_qmu[vn].Eta())/abs(qmu_dmuvecmaxlxy[vn].Vect().DeltaPhi(svvecmaxlxy_qmu[vn])))
        maxdetadphisv = max(abs(qmu_dmuvecminlxy[vn].Vect().Eta()-svvecminlxy_qmu[vn].Eta())/abs(qmu_dmuvecminlxy[vn].Vect().DeltaPhi(svvecminlxy_qmu[vn])),
                            abs(qmu_dmuvecmaxlxy[vn].Vect().Eta()-svvecmaxlxy_qmu[vn].Eta())/abs(qmu_dmuvecmaxlxy[vn].Vect().DeltaPhi(svvecmaxlxy_qmu[vn])))
        mina3dsv = min(abs(qmu_dmuvecminlxy[vn].Vect().Angle(svvecminlxy_qmu[vn])),abs(qmu_dmuvecmaxlxy[vn].Vect().Angle(svvecmaxlxy_qmu[vn])))
        maxa3dsv = max(abs(qmu_dmuvecminlxy[vn].Vect().Angle(svvecminlxy_qmu[vn])),abs(qmu_dmuvecmaxlxy[vn].Vect().Angle(svvecmaxlxy_qmu[vn])))
        #
        mindphisvu = min(abs(qmu_dmuvecdpminlxy[vn].Vect().DeltaPhi(svvecminlxy_qmu[vn])),abs(qmu_dmuvecdpmaxlxy[vn].Vect().DeltaPhi(svvecmaxlxy_qmu[vn])))
        maxdphisvu = max(abs(qmu_dmuvecdpminlxy[vn].Vect().DeltaPhi(svvecminlxy_qmu[vn])),abs(qmu_dmuvecdpmaxlxy[vn].Vect().DeltaPhi(svvecmaxlxy_qmu[vn])))
        mindetasvu = min(abs(qmu_dmuvecdpminlxy[vn].Vect().Eta()-svvecminlxy_qmu[vn].Eta()),abs(qmu_dmuvecdpmaxlxy[vn].Vect().Eta()-svvecmaxlxy_qmu[vn].Eta()))
        maxdetasvu = max(abs(qmu_dmuvecdpminlxy[vn].Vect().Eta()-svvecminlxy_qmu[vn].Eta()),abs(qmu_dmuvecdpmaxlxy[vn].Vect().Eta()-svvecmaxlxy_qmu[vn].Eta()))
        mindetadphisvu = min(abs(qmu_dmuvecdpminlxy[vn].Vect().Eta()-svvecminlxy_qmu[vn].Eta())/abs(qmu_dmuvecdpminlxy[vn].Vect().DeltaPhi(svvecminlxy_qmu[vn])),
                            abs(qmu_dmuvecdpmaxlxy[vn].Vect().Eta()-svvecmaxlxy_qmu[vn].Eta())/abs(qmu_dmuvecdpmaxlxy[vn].Vect().DeltaPhi(svvecmaxlxy_qmu[vn])))
        maxdetadphisvu = max(abs(qmu_dmuvecdpminlxy[vn].Vect().Eta()-svvecminlxy_qmu[vn].Eta())/abs(qmu_dmuvecdpminlxy[vn].Vect().DeltaPhi(svvecminlxy_qmu[vn])),
                            abs(qmu_dmuvecdpmaxlxy[vn].Vect().Eta()-svvecmaxlxy_qmu[vn].Eta())/abs(qmu_dmuvecdpmaxlxy[vn].Vect().DeltaPhi(svvecmaxlxy_qmu[vn])))
        mina3dsvu = min(abs(qmu_dmuvecdpminlxy[vn].Vect().Angle(svvecminlxy_qmu[vn])),abs(qmu_dmuvecdpmaxlxy[vn].Vect().Angle(svvecmaxlxy_qmu[vn])))
        maxa3dsvu = max(abs(qmu_dmuvecdpminlxy[vn].Vect().Angle(svvecminlxy_qmu[vn])),abs(qmu_dmuvecdpmaxlxy[vn].Vect().Angle(svvecmaxlxy_qmu[vn])))
        #
        minlxy_desvdpsvmin = 1e6
        minlxy_demmdpsvmin = 1e6
        minlxy_dpsvmin = 1e6
        minlxy_desvdpsvminu = 1e6
        minlxy_demmdpsvminu = 1e6
        minlxy_dpsvminu = 1e6
        tdpsv = abs(t.Muon_vec[qmuidxsminlxy[vn*2]].Vect().DeltaPhi(svvecminlxy_qmu[vn]))
        if tdpsv < minlxy_dpsvmin:
            minlxy_dpsvmin = tdpsv
        tdpsv = abs(t.Muon_vec[qmuidxsminlxy[vn*2+1]].Vect().DeltaPhi(svvecminlxy_qmu[vn]))
        if tdpsv < minlxy_dpsvmin:
            minlxy_dpsvmin = tdpsv
        if minlxy_dpsvmin>0.0:
            tdemm = abs(t.Muon_vec[qmuidxsminlxy[vn*2]].Eta()-t.Muon_vec[qmuidxsminlxy[vn*2+1]].Eta())
            tdetasv = abs(qmu_dmuvecminlxy[vn].Vect().Eta()-svvecminlxy_qmu[vn].Eta())
            minlxy_demmdpsvmin = tdemm/minlxy_dpsvmin
            minlxy_desvdpsvmin = tdetasv/minlxy_dpsvmin
        #
        tdpsvu = abs(qmu_muvecdpminlxy[vn*2].Vect().DeltaPhi(svvecminlxy_qmu[vn]))
        if tdpsvu < minlxy_dpsvminu:
            minlxy_dpsvminu = tdpsvu
        tdpsvu = abs(qmu_muvecdpminlxy[vn*2+1].Vect().DeltaPhi(svvecminlxy_qmu[vn]))
        if tdpsvu < minlxy_dpsvminu:
            minlxy_dpsvminu = tdpsvu
        if minlxy_dpsvminu>0.0:
            tdemmu = abs(qmu_muvecdpminlxy[vn*2].Eta()-qmu_muvecdpminlxy[vn*2+1].Eta())
            tdetasvu = abs(qmu_dmuvecdpminlxy[vn].Vect().Eta()-svvecminlxy_qmu[vn].Eta())
            minlxy_demmdpsvminu = tdemmu/minlxy_dpsvminu
            minlxy_desvdpsvminu = tdetasvu/minlxy_dpsvminu
        #
        maxlxy_desvdpsvmin = 1e6
        maxlxy_demmdpsvmin = 1e6
        maxlxy_dpsvmin = 1e6
        tdpsv = abs(t.Muon_vec[qmuidxsmaxlxy[vn*2]].Vect().DeltaPhi(svvecmaxlxy_qmu[vn]))
        if tdpsv < maxlxy_dpsvmin:
            maxlxy_dpsvmin = tdpsv
        tdpsv = abs(t.Muon_vec[qmuidxsmaxlxy[vn*2+1]].Vect().DeltaPhi(svvecmaxlxy_qmu[vn]))
        if tdpsv < maxlxy_dpsvmin:
            maxlxy_dpsvmin = tdpsv
        if maxlxy_dpsvmin>0.0:
            tdemm = abs(t.Muon_vec[qmuidxsmaxlxy[vn*2]].Eta()-t.Muon_vec[qmuidxsmaxlxy[vn*2+1]].Eta())
            tdetasv = abs(qmu_dmuvecmaxlxy[vn].Vect().Eta()-svvecmaxlxy_qmu[vn].Eta())
            maxlxy_demmdpsvmin = tdemm/maxlxy_dpsvmin
            maxlxy_desvdpsvmin = tdetasv/maxlxy_dpsvmin
        #
        maxlxy_desvdpsvminu = 1e6
        maxlxy_demmdpsvminu = 1e6
        maxlxy_dpsvminu = 1e6
        tdpsvu = abs(qmu_muvecdpmaxlxy[vn*2].Vect().DeltaPhi(svvecmaxlxy_qmu[vn]))
        if tdpsvu < maxlxy_dpsvminu:
            maxlxy_dpsvminu = tdpsvu
        tdpsvu = abs(qmu_muvecdpmaxlxy[vn*2+1].Vect().DeltaPhi(svvecmaxlxy_qmu[vn]))
        if tdpsvu < maxlxy_dpsvminu:
            maxlxy_dpsvminu = tdpsvu
        if maxlxy_dpsvminu>0.0:
            tdemmu = abs(qmu_muvecdpmaxlxy[vn*2].Eta()-qmu_muvecdpmaxlxy[vn*2+1].Eta())
            tdetasvu = abs(qmu_dmuvecdpmaxlxy[vn].Vect().Eta()-svvecmaxlxy_qmu[vn].Eta())
            maxlxy_demmdpsvminu = tdemmu/maxlxy_dpsvminu
            maxlxy_desvdpsvminu = tdetasvu/maxlxy_dpsvminu
        #
        mindesvdpsvmin = min(minlxy_desvdpsvmin, maxlxy_desvdpsvmin)
        maxdesvdpsvmin = max(minlxy_desvdpsvmin, maxlxy_desvdpsvmin)
        mindesvdpsvminu = min(minlxy_desvdpsvminu, maxlxy_desvdpsvminu)
        maxdesvdpsvminu = max(minlxy_desvdpsvminu, maxlxy_desvdpsvminu)
        mindemmdpsvmin = min(minlxy_demmdpsvmin, maxlxy_demmdpsvmin)
        maxdemmdpsvmin = max(minlxy_demmdpsvmin, maxlxy_demmdpsvmin)
        mindemmdpsvminu = min(minlxy_demmdpsvminu, maxlxy_demmdpsvminu)
        maxdemmdpsvminu = max(minlxy_demmdpsvminu, maxlxy_demmdpsvminu)
        #
        minlxy_sindpsvlxy = abs(ROOT.TMath.Sin(qmu_dmuvecminlxy[vn].Vect().DeltaPhi(svvecminlxy_qmu[vn])))*minlxy
        minlxy_sindpsvlxyu = abs(ROOT.TMath.Sin(qmu_dmuvecdpminlxy[vn].Vect().DeltaPhi(svvecminlxy_qmu[vn])))*minlxy
        minlxy_sina3dsvl3d = abs(ROOT.TMath.Sin(qmu_dmuvecminlxy[vn].Vect().Angle(svvecminlxy_qmu[vn])))*minl3d
        minlxy_sina3dsvl3du = abs(ROOT.TMath.Sin(qmu_dmuvecdpminlxy[vn].Vect().Angle(svvecminlxy_qmu[vn])))*minl3d
        #
        maxlxy_sindpsvlxy = abs(ROOT.TMath.Sin(qmu_dmuvecmaxlxy[vn].Vect().DeltaPhi(svvecmaxlxy_qmu[vn])))*maxlxy
        maxlxy_sindpsvlxyu = abs(ROOT.TMath.Sin(qmu_dmuvecdpmaxlxy[vn].Vect().DeltaPhi(svvecmaxlxy_qmu[vn])))*maxlxy
        maxlxy_sina3dsvl3d = abs(ROOT.TMath.Sin(qmu_dmuvecmaxlxy[vn].Vect().Angle(svvecmaxlxy_qmu[vn])))*maxl3d
        maxlxy_sina3dsvl3du = abs(ROOT.TMath.Sin(qmu_dmuvecdpmaxlxy[vn].Vect().Angle(svvecmaxlxy_qmu[vn])))*maxl3d
        #
        minsindpsvlxy = min(minlxy_sindpsvlxy, maxlxy_sindpsvlxy)
        maxsindpsvlxy = max(minlxy_sindpsvlxy, maxlxy_sindpsvlxy)        
        minsindpsvlxyu = min(minlxy_sindpsvlxyu, maxlxy_sindpsvlxyu)
        maxsindpsvlxyu = max(minlxy_sindpsvlxyu, maxlxy_sindpsvlxyu)
        minsina3dsvl3d = min(minlxy_sina3dsvl3d, maxlxy_sina3dsvl3d)
        maxsina3dsvl3d = max(minlxy_sina3dsvl3d, maxlxy_sina3dsvl3d)
        minsina3dsvl3du = min(minlxy_sina3dsvl3du, maxlxy_sina3dsvl3du)
        maxsina3dsvl3du = max(minlxy_sina3dsvl3du, maxlxy_sina3dsvl3du)
        #
        if applyFourMuonAngularSel:
            if maxa3dmmu>0.9*(ROOT.TMath.Pi()) or maxa3dmmu>0.9*(ROOT.TMath.Pi()):
                continue
            if maxdpmmu>0.9*(ROOT.TMath.Pi()) or mindpmmu>0.9*(ROOT.TMath.Pi()):
                continue
            if maxdphisvu>ROOT.TMath.PiOver2() or mindphisvu>ROOT.TMath.PiOver2():
                continue
            if maxa3dsvu>ROOT.TMath.PiOver2() or mina3dsvu>ROOT.TMath.PiOver2():
                continue
        #
        if applyFourMuonMassDiffSel:
            if reldmass > 0.1:
                continue
        #
        for h in h1d["fourmuon"]:
            tn = h.GetName()
            h.Fill(eval(variable1d[h.GetName()]), lumiweight)
        for h in h2d["fourmuon"]:
            tn = h.GetName()
            h.Fill(eval(variable2d[h.GetName()][0]),eval(variable2d[h.GetName()][1]), lumiweight)
        # Scan:
        if ((not filledcat4musep) and (not filledcat4muosv) and (not filledcat2mu)): 
            m4fit.setVal(mass)
            roow4.setVal(lumiweight*rooweight);
            roods["FourMu_sep"].add(ROOT.RooArgSet(m4fit,roow4),roow4.getVal());
            catmass["FourMu_sep"].Fill(mass, lumiweight*rooweight);
            #mfit.setVal(avgmass)
            #roow.setVal(lumiweight);
            #roods["FourMu_sep"].add(ROOT.RooArgSet(mfit,roow),roow.getVal());
            #catmass["FourMu_sep"].Fill(avgmass, lumiweight);
            filledcat4musep = True
        else:
            filledcat4musep = False

    # Fill histograms for selected non-overlapping SVs (with a selected four muon system)
    nSVselass_qmu = 0
    for v in set(selqmusvidxs):
        nSVselass_qmu = nSVselass_qmu+1
        lxy = t.SV_lxy[v]
        for h in h1d["svselass_fourmu"]:
            tn = h.GetName()
            h.Fill(eval(variable1d[h.GetName()]), lumiweight)
        for h in h2d["svselass_fourmu"]:
            tn = h.GetName()
            h.Fill(eval(variable2d[h.GetName()][0]),eval(variable2d[h.GetName()][1]), lumiweight)
    nSVs = nSVselass_qmu
    for h in h1d["nsvselass_fourmu"]:
        if args.noFourMuon:
            break
        tn = h.GetName()
        h.Fill(eval(variable1d[h.GetName()]), lumiweight)

    # Apply selections and fill histograms for four-muon systems from overlapping SVs
    selqmusvidxs_osv = []
    selqmuvecs_osv = []
    mindrmm, mindpmm, mindemm, mindedpmm, mina3dmm = 1e6, 1e6, 1e6, 1e6, 1e6
    maxdrmm, maxdpmm, maxdemm, maxdedpmm, maxa3dmm = -1., -1., -1., -1., -1.
    mindrmmu, mindpmmu, mindemmu, mindedpmmu, mina3dmmu = 1e6, 1e6, 1e6, 1e6, 1e6
    maxdrmmu, maxdpmmu, maxdemmu, maxdedpmmu, maxa3dmmu = -1., -1., -1., -1., -1.
    for vn,v in enumerate(qmuvec_osv):
        if args.noFourMuonOSV:
            break
        if not applyDiMuonSelectionForFourMuonOSV(qmu_dmuvec_osv[vn*2], qmu_dmuvec_osv[vn*2+1]):
            continue
        if not applyFourMuonSelection(v):
            continue
        lxy  = t.SVOverlap_lxy[osvidx_qmu[vn]]
        if not applyLxySelection(lxy):
            continue
        qmuidxs_osv_sel.append(qmuidxs_osv[vn*4])
        qmuidxs_osv_sel.append(qmuidxs_osv[vn*4+1])
        qmuidxs_osv_sel.append(qmuidxs_osv[vn*4+2])
        qmuidxs_osv_sel.append(qmuidxs_osv[vn*4+3])
        selqmusvidxs_osv.append(t.SVOverlap_vtxIdxs[osvidx_qmu[vn]][0])
        selqmuvecs_osv.append(vn)
        mass = v.M()
        pt   = v.Pt()
        nhitsbeforesvtotal = t.Muon_nhitsbeforesv[qmuidxs_osv[vn*4]] + t.Muon_nhitsbeforesv[qmuidxs_osv[vn*4+1]] + t.Muon_nhitsbeforesv[qmuidxs_osv[vn*4+2]] + t.Muon_nhitsbeforesv[qmuidxs_osv[vn*4+3]]
        for m in range(vn*4,vn*4+4):
            for mm in range(m+1,vn*4+4):
                drmm = t.Muon_vec[qmuidxs_osv[m]].DeltaR(t.Muon_vec[qmuidxs_osv[mm]])
                dpmm = abs(t.Muon_vec[qmuidxs_osv[m]].DeltaPhi(t.Muon_vec[qmuidxs_osv[mm]]))
                demm = abs(t.Muon_vec[qmuidxs_osv[m]].Eta()-t.Muon_vec[qmuidxs_osv[mm]].Eta())
                dedpmm = 1e6
                if dpmm>0.0:
                    dedpmm = demm/dpmm
                a3dmm = abs(t.Muon_vec[qmuidxs_osv[m]].Angle(t.Muon_vec[qmuidxs_osv[mm]].Vect()))
                if drmm<mindrmm:
                    mindrmm = drmm
                if drmm>maxdrmm:
                    maxdrmm = drmm
                if dpmm<mindpmm:
                    mindpmm = dpmm
                if dpmm>maxdpmm:
                    maxdpmm = dpmm
                if demm<mindemm:
                    mindemm = demm
                if demm>maxdemm:
                    maxdemm = demm
                if dedpmm<mindedpmm:
                    mindedpmm = dedpmm
                if dedpmm>maxdedpmm:
                    maxdedpmm = dedpmm
                if a3dmm<mina3dmm:
                    mina3dmm = a3dmm
                if a3dmm>maxa3dmm:
                    maxa3dmm = a3dmm
                #
                drmmu = qmu_muvecdp_osv[m].DeltaR(qmu_muvecdp_osv[mm])
                dpmmu = abs(qmu_muvecdp_osv[m].DeltaPhi(qmu_muvecdp_osv[mm]))
                demmu = abs(qmu_muvecdp_osv[m].Eta()-qmu_muvecdp_osv[mm].Eta())
                dedpmmu = 1e6
                if dpmmu>0.0:
                    dedpmmu = demmu/dpmmu
                a3dmmu = abs(qmu_muvecdp_osv[m].Angle(qmu_muvecdp_osv[mm].Vect()))
                if drmmu<mindrmmu:
                    mindrmmu = drmmu
                if drmmu>maxdrmmu:
                    maxdrmmu = drmmu
                if dpmmu<mindpmmu:
                    mindpmmu = dpmmu
                if dpmmu>maxdpmmu:
                    maxdpmmu = dpmmu
                if demmu<mindemmu:
                    mindemmu = demmu
                if demmu>maxdemmu:
                    maxdemmu = demmu
                if dedpmmu<mindedpmmu:
                    mindedpmmu = dedpmmu
                if dedpmmu>maxdedpmmu:
                    maxdedpmmu = dedpmmu
                if a3dmmu<mina3dmmu:
                    mina3dmmu = a3dmmu
                if a3dmmu>maxa3dmmu:
                    maxa3dmmu = a3dmmu
        #
        dphisv = abs(v.Vect().DeltaPhi(osvvec_qmu[vn]))
        detasv = abs(v.Vect().Eta()-osvvec_qmu[vn].Eta())
        detadphisv = 1e6
        if dphisv>0.0:
            detadphisv = detasv/dphisv
        a3dsv  = abs(v.Vect().Angle(osvvec_qmu[vn]))
        #
        desvdpsvmin = 1e6
        demmdpsvmin = 1e6
        dpsvmin = 1e6
        for m in range(vn*4,vn*4+4):
            tdpsv = abs(t.Muon_vec[qmuidxs_osv[m]].Vect().DeltaPhi(osvvec_qmu[vn]))
            if tdpsv < dpsvmin:
                dpsvmin = tdpsv
        if dpsvmin>0.0:
            demmdpsvmin = maxdemm/dpsvmin
            desvdpsvmin = detasv/dpsvmin
        #
        vu = qmu_muvecdp_osv[vn*4]
        for mm in range(vn*4+1,vn*4+4):
            vu = vu + qmu_muvecdp_osv[mm]
        dphisvu = abs(vu.Vect().DeltaPhi(osvvec_qmu[vn]))
        detasvu = abs(vu.Vect().Eta()-osvvec_qmu[vn].Eta())
        detadphisvu = 1e6
        if dphisvu>0.0:
            detadphisvu = detasvu/dphisvu
        a3dsvu  = abs(vu.Vect().Angle(osvvec_qmu[vn]))
        #
        desvdpsvminu = 1e6
        demmdpsvminu = 1e6
        dpsvminu = 1e6
        for m in range(vn*4,vn*4+4):
            tdpsvu = abs(qmu_muvecdp_osv[mm].Vect().DeltaPhi(osvvec_qmu[vn]))
            if tdpsvu < dpsvminu:
                dpsvminu = tdpsvu
        if dpsvminu>0.0:
            demmdpsvminu = demmu/dpsvminu
            desvdpsvminu = detasvu/dpsvminu
        #
        sindpsvlxy = abs(ROOT.TMath.Sin(v.Vect().DeltaPhi(osvvec_qmu[vn])))*lxy
        sindpsvlxyu = abs(ROOT.TMath.Sin(vu.Vect().DeltaPhi(osvvec_qmu[vn])))*lxy
        sina3dsvl3d = abs(ROOT.TMath.Sin(v.Vect().Angle(osvvec_qmu[vn])))*(osvvec_qmu[vn].Mag())
        sina3dsvl3du = abs(ROOT.TMath.Sin(vu.Vect().Angle(osvvec_qmu[vn])))*(osvvec_qmu[vn].Mag())
        #
        if applyFourMuonAngularSel:
            if maxa3dmmu>0.9*(ROOT.TMath.Pi()) or maxa3dmmu>0.9*(ROOT.TMath.Pi()):
                continue
            if maxdpmmu>0.9*(ROOT.TMath.Pi()) or mindpmmu>0.9*(ROOT.TMath.Pi()):
                continue
            if a3dsvu>ROOT.TMath.PiOver2():
                continue
            if dphisvu>ROOT.TMath.PiOver2():
                continue
        #
        for h in h1d["fourmuon_osv"]:
            tn = h.GetName()
            h.Fill(eval(variable1d[h.GetName()]), lumiweight)
        for h in h2d["fourmuon_osv"]:
            tn = h.GetName()
            h.Fill(eval(variable2d[h.GetName()][0]),eval(variable2d[h.GetName()][1]), lumiweight)
        # Scan:
        if ( (not filledcat4musep) and (not filledcat4muosv) and (not filledcat2mu) ): 
            m4fit.setVal(mass)
            roow4.setVal(lumiweight*rooweight);
            roods["FourMu_osv"].add(ROOT.RooArgSet(m4fit,roow4),roow4.getVal());
            catmass["FourMu_osv"].Fill(mass, lumiweight*rooweight);
            filledcat4muosv = True
        else:
            filledcat4muosv = False

    # Fill histograms for selected overlapping SVs (with a selected four muon system)
    nSVselass_qmu_osv = 0
    for v in set(selqmusvidxs_osv):
        nSVselass_qmu_osv = nSVselass_qmu_osv+1
        lxy = t.SV_lxy[v]
        for h in h1d["svselass_fourmu_osv"]:
            tn = h.GetName()
            h.Fill(eval(variable1d[h.GetName()]), lumiweight)
        for h in h2d["svselass_fourmu_osv"]:
            tn = h.GetName()
            h.Fill(eval(variable2d[h.GetName()][0]),eval(variable2d[h.GetName()][1]), lumiweight)
    nSVs = nSVselass_qmu_osv
    for h in h1d["nsvselass_fourmu_osv"]:
        if args.noFourMuonOSV:
            break
        tn = h.GetName()
        h.Fill(eval(variable1d[h.GetName()]), lumiweight)

    # Apply selections and fill histograms for muon pairs from non-overlapping SVs
    seldmuidxs = []
    seldmusvidxs = []
    seldmuvecs = []
    for vn,v in enumerate(dmuvec):
        if args.noDiMuon:
            break
        if not applyDiMuonSelection(v):
            continue
        lxy  = t.SV_lxy[svidx[vn]]
        if not applyLxySelection(lxy):
            continue
        lz  = abs(t.SV_z[svidx[vn]])
        if not applyLzSelection(lz):
            continue
        mass = v.M()
        pt   = v.Pt()
        nhitsbeforesvtotal = t.Muon_nhitsbeforesv[dmuidxs[int(vn*2)]] + t.Muon_nhitsbeforesv[dmuidxs[int(vn*2)+1]]
        #
        # Check if the dimuon is gen-matched
        isgen = False
        drgen = 9999.
        lxygen = -999.
        idgen = 0
        for gn,g in enumerate(dmumot):
            tmp_drgen = g.DeltaR(v)
            deltapt = abs(g.Pt() - v.Pt())/g.Pt()
            if tmp_drgen < drgen:
                drgen = tmp_drgen 
        if drgen < 0.1:
            isgen = True
            lxygen = t.GenPart_lxy[dmugen[2*gn]]
            idgen = t.GenPart_motherPdgId[dmugen[2*gn]]
        deltalxy = (lxy - lxygen)/lxygen
        #if not isgen or abs(deltalxy) < 0.1 or lxy < 2:
        #    continue
        #
        #  Apply selection on muon lifetime-scaled dxy
        if applyMuonIPSel:
            if ( abs(t.Muon_dxyCorr[dmuidxs[int(vn*2)]]  )/(lxy*mass/pt)<0.1 or
                 abs(t.Muon_dxyCorr[dmuidxs[int(vn*2)+1]])/(lxy*mass/pt)<0.1 ):
                continue
            if ( abs(t.Muon_dxyCorr[dmuidxs[int(vn*2)]]/t.Muon_dxye[dmuidxs[int(vn*2)]])<2.0 or abs(t.Muon_dxyCorr[dmuidxs[int(vn*2)+1]]/t.Muon_dxye[dmuidxs[int(vn*2)+1]])<2.0 ):
                continue
        if applyMuonHitSel:
            if lxy < 11.0:
                if ( nhitsbeforesvtotal > 0):
                    continue
            elif lxy > 11.0 and lxy < 16.0:
                if ( nhitsbeforesvtotal > 1):
                    continue
            elif lxy > 16.0:
                if ( nhitsbeforesvtotal > 2):
                    continue
        drmm = t.Muon_vec[dmuidxs[int(vn*2)]].DeltaR(t.Muon_vec[dmuidxs[int(vn*2)+1]])
        dpmm = abs(t.Muon_vec[dmuidxs[int(vn*2)]].DeltaPhi(t.Muon_vec[dmuidxs[int(vn*2)+1]]))
        demm = abs(t.Muon_vec[dmuidxs[int(vn*2)]].Eta()-t.Muon_vec[dmuidxs[int(vn*2)+1]].Eta())
        dedpmm = 1e6
        if dpmm>0.0:
            dedpmm = demm/dpmm
        a3dmm = abs(t.Muon_vec[dmuidxs[int(vn*2)]].Angle(t.Muon_vec[dmuidxs[int(vn*2)+1]].Vect()))
        #
        drmmu = dmu_muvecdp[int(vn*2)].DeltaR(dmu_muvecdp[int(vn*2)+1])
        dpmmu = abs(dmu_muvecdp[int(vn*2)].DeltaPhi(dmu_muvecdp[int(vn*2)+1]))
        demmu = abs(dmu_muvecdp[int(vn*2)].Eta()-dmu_muvecdp[int(vn*2)+1].Eta())
        dedpmmu = 1e6
        if dpmmu>0.0:
            dedpmmu = demmu/dpmmu
        a3dmmu = abs(dmu_muvecdp[int(vn*2)].Angle(dmu_muvecdp[int(vn*2)+1].Vect()))
        #
        dphisv = abs(v.Vect().DeltaPhi(svvec[vn]))
        dphisv1 = abs(t.Muon_vec[dmuidxs[int(vn*2)]].Vect().DeltaPhi(svvec[vn]))
        dphisv2 = abs(t.Muon_vec[dmuidxs[int(vn*2)+1]].Vect().DeltaPhi(svvec[vn]))
        detasv = abs(v.Vect().Eta()-svvec[vn].Eta())
        detadphisv = 1e6
        if dphisv>0.0:
            detadphisv = detasv/dphisv
        a3dsv  = abs(v.Vect().Angle(svvec[vn]))
        #
        desvdpsvmin = 1e6
        demmdpsvmin = 1e6
        dpsvmin = 1e6
        tdpsv = abs(t.Muon_vec[dmuidxs[int(vn*2)]].Vect().DeltaPhi(svvec[vn]))
        if tdpsv < dpsvmin:
            dpsvmin = tdpsv
        tdpsv = abs(t.Muon_vec[dmuidxs[int(vn*2)+1]].Vect().DeltaPhi(svvec[vn]))
        if tdpsv < dpsvmin:
            dpsvmin = tdpsv
        if dpsvmin>0.0:
            demmdpsvmin = demm/dpsvmin
            desvdpsvmin = detasv/dpsvmin
        #
        vu = (dmu_muvecdp[int(vn*2)]+dmu_muvecdp[int(vn*2)+1])
        dphisvu = abs(vu.Vect().DeltaPhi(svvec[vn]))
        dphisv1u = abs(dmu_muvecdp[int(vn*2)].Vect().DeltaPhi(svvec[vn]))
        dphisv2u = abs(dmu_muvecdp[int(vn*2)+1].Vect().DeltaPhi(svvec[vn]))
        detasvu = abs(vu.Vect().Eta()-svvec[vn].Eta())
        detadphisvu = 1e6
        if dphisvu>0.0:
            detadphisvu = detasvu/dphisvu
        a3dsvu  = abs(vu.Vect().Angle(svvec[vn]))
        #
        desvdpsvminu = 1e6
        demmdpsvminu = 1e6
        dpsvminu = 1e6
        tdpsvu = abs(dmu_muvecdp[int(vn*2)].Vect().DeltaPhi(svvec[vn]))
        if tdpsvu < dpsvminu:
            dpsvminu = tdpsvu
        tdpsvu = abs(dmu_muvecdp[int(vn*2)+1].Vect().DeltaPhi(svvec[vn]))
        if tdpsvu < dpsvminu:
            dpsvminu = tdpsvu
        if dpsvminu>0.0:
            demmdpsvminu = demmu/dpsvminu
            desvdpsvminu = detasvu/dpsvminu
        #
        sindpsvlxy = abs(ROOT.TMath.Sin(v.Vect().DeltaPhi(svvec[vn])))*lxy
        sindpsvlxyu = abs(ROOT.TMath.Sin(vu.Vect().DeltaPhi(svvec[vn])))*lxy
        sina3dsvl3d = abs(ROOT.TMath.Sin(v.Vect().Angle(svvec[vn])))*(svvec[vn].Mag())
        sina3dsvl3du = abs(ROOT.TMath.Sin(vu.Vect().Angle(svvec[vn])))*(svvec[vn].Mag())
        #
        # Apply selection on angular distances
        if applyDiMuonAngularSel:
            if ROOT.TMath.Log10(dedpmmu)>1.25:
                continue
            if dpmmu>0.9*(ROOT.TMath.Pi()) or a3dmmu>0.9*(ROOT.TMath.Pi()):
                continue
            #if dphisv1u>ROOT.TMath.PiOver2() or dphisv2u>ROOT.TMath.PiOver2() or a3dsvu>ROOT.TMath.PiOver2(): # uncomment if we want to apply this cut per muon
            if dphisvu>ROOT.TMath.PiOver2() or a3dsvu>ROOT.TMath.PiOver2():
            #if dphisvu>0.02 or a3dsvu>ROOT.TMath.PiOver2():
                continue
        #
        #  Apply categorization and selection on muon isolation
        isocat = dimuonIsoCategory(t.Muon_PFIsoAll0p4[dmuidxs[int(vn*2)]], t.Muon_pt[dmuidxs[int(vn*2)]], t.Muon_PFIsoAll0p4[dmuidxs[int(vn*2)+1]], t.Muon_pt[dmuidxs[int(vn*2)+1]])
        if float(args.dimuonIsoCatSel) > -1 and isocat!=float(args.dimuonIsoCatSel):
            continue
        #
        seldmuidxs.append(dmuidxs[int(vn*2)])
        seldmuidxs.append(dmuidxs[int(vn*2)+1])
        seldmusvidxs.append(svidx[vn])        
        seldmuvecs.append(vn)        
        for h in h1d["dimuon"]:
            tn = h.GetName()
            if "dimuon_gen_" in tn and not isgen:
                continue
            if "dimuon_genjpsi_" in tn and not isgen and not idgen==443:
                continue
            h.Fill(eval(variable1d[h.GetName()]), lumiweight)
        for h in h2d["dimuon"]:
            tn = h.GetName()
            if "dimuon_gen_" in tn and not isgen:
                continue
            if "dimuon_genjpsi_" in tn and not isgen and not idgen==443:
                continue
            h.Fill(eval(variable2d[h.GetName()][0]),eval(variable2d[h.GetName()][1]), lumiweight)
        # Scan:
        if ( (not filledcat4musep) and (not filledcat4muosv) and (not filledcat2mu) ): 
            slice = ""
            if dphisvu < dphisvcut:
                for l in range(len(lxybins)-1):
                    label = lxybinlabel[l]  
                    if lxy > lxybins[l] and lxy < lxybins[l+1]:
                        break
                if isocat==2:
                    if pt < ptcut:
                        slice = "Dimuon_"+label+"_iso1_ptlow"
                    else:
                        slice = "Dimuon_"+label+"_iso1_pthigh"
                else:
                    if pt < ptcut:
                        slice = "Dimuon_"+label+"_iso0_ptlow"
                    else:
                        slice = "Dimuon_"+label+"_iso0_pthigh"
            else:
                for l in range(len(lxybins)-1):
                    label = lxybinlabel[l]  
                    if lxy > lxybins[l] and lxy < lxybins[l+1]:
                        slice = "Dimuon_"+label+"_non-pointing"
                        break
            if slice!="":
                mfit.setVal(mass)
                roow.setVal(lumiweight*rooweight);
                roods[slice].add(ROOT.RooArgSet(mfit,roow),roow.getVal());
                catmass[slice].Fill(mass, lumiweight*rooweight);
                roods["Dimuon_"+label+"_inclusive"].add(ROOT.RooArgSet(mfit,roow),roow.getVal());
                roods["Dimuon_full_inclusive"].add(ROOT.RooArgSet(mfit,roow),roow.getVal());
                filledcat2mu = True

    # Apply selections and fill histograms for muon pairs from overlapping SVs
    seldmuidxs_osv = []
    seldmusvidxs_osv = []
    seldmuvecs_osv = []
    for vn,v in enumerate(dmuvec_osv):
        if args.noDiMuon:
            break
        if not applyDiMuonSelection(v):
            continue
        lxy  = t.SVOverlap_lxy[osvidx[vn]]
        if not applyLxySelection(lxy):
            continue
        lz  = abs(t.SVOverlap_z[osvidx[vn]])
        if not applyLzSelection(lz):
            continue
        mass = v.M()
        pt   = v.Pt()
        nhitsbeforesvtotal = t.Muon_nhitsbeforesv[dmuidxs_osv[int(vn*2)]] + t.Muon_nhitsbeforesv[dmuidxs_osv[int(vn*2)+1]]
        #  Apply selection on muon lifetime-scaled dxy
        if applyMuonIPSel:
            if ( abs(t.Muon_dxyCorr[dmuidxs_osv[int(vn*2)]]  )/(lxy*mass/pt)<0.1 or
                 abs(t.Muon_dxyCorr[dmuidxs_osv[int(vn*2)+1]])/(lxy*mass/pt)<0.1 ):
                continue
            if ( abs(t.Muon_dxyCorr[dmuidxs_osv[int(vn*2)]]/t.Muon_dxye[dmuidxs_osv[int(vn*2)]])<2.0 or abs(t.Muon_dxyCorr[dmuidxs_osv[int(vn*2)+1]]/t.Muon_dxye[dmuidxs_osv[int(vn*2)+1]])<2.0 ):
                continue
        if applyMuonHitSel:
            if lxy < 11.0:
                if ( nhitsbeforesvtotal > 0):
                    continue
            elif lxy > 11.0 and lxy < 16.0:
                if ( nhitsbeforesvtotal > 1):
                    continue
            elif lxy > 16.0:
                if ( nhitsbeforesvtotal > 2):
                    continue
        drmm = t.Muon_vec[dmuidxs_osv[int(vn*2)]].DeltaR(t.Muon_vec[dmuidxs_osv[int(vn*2)+1]])
        dpmm = abs(t.Muon_vec[dmuidxs_osv[int(vn*2)]].DeltaPhi(t.Muon_vec[dmuidxs_osv[int(vn*2)+1]]))
        demm = abs(t.Muon_vec[dmuidxs_osv[int(vn*2)]].Eta()-t.Muon_vec[dmuidxs_osv[int(vn*2)+1]].Eta())
        dedpmm = 1e6
        if dpmm>0.0:
            dedpmm = demm/dpmm
        a3dmm = abs(t.Muon_vec[dmuidxs_osv[int(vn*2)]].Angle(t.Muon_vec[dmuidxs_osv[int(vn*2)+1]].Vect()))
        #
        drmmu = dmu_muvecdp_osv[int(vn*2)].DeltaR(dmu_muvecdp_osv[int(vn*2)+1])
        dpmmu = abs(dmu_muvecdp_osv[int(vn*2)].DeltaPhi(dmu_muvecdp_osv[int(vn*2)+1]))
        demmu = abs(dmu_muvecdp_osv[int(vn*2)].Eta()-dmu_muvecdp_osv[int(vn*2)+1].Eta())
        dedpmmu = 1e6
        if dpmmu>0.0:
            dedpmmu = demmu/dpmmu
        a3dmmu = abs(dmu_muvecdp_osv[int(vn*2)].Angle(dmu_muvecdp_osv[int(vn*2)+1].Vect()))
        #
        dphisv = abs(v.Vect().DeltaPhi(osvvec[vn]))
        dphisv1 = abs(t.Muon_vec[dmuidxs_osv[int(vn*2)]].Vect().DeltaPhi(osvvec[vn]))
        dphisv2 = abs(t.Muon_vec[dmuidxs_osv[int(vn*2)+1]].Vect().DeltaPhi(osvvec[vn]))
        detasv = abs(v.Vect().Eta()-osvvec[vn].Eta())
        detadphisv = 1e6
        if dphisv>0.0:
            detadphisv = detasv/dphisv
        a3dsv  = abs(v.Vect().Angle(osvvec[vn]))
        #
        desvdpsvmin = 1e6
        demmdpsvmin = 1e6
        dpsvmin = 1e6
        tdpsv = abs(t.Muon_vec[dmuidxs_osv[int(vn*2)]].Vect().DeltaPhi(osvvec[vn]))
        if tdpsv < dpsvmin:
            dpsvmin = tdpsv
        tdpsv = abs(t.Muon_vec[dmuidxs_osv[int(vn*2)+1]].Vect().DeltaPhi(osvvec[vn]))
        if tdpsv < dpsvmin:
            dpsvmin = tdpsv
        if dpsvmin>0.0:
            demmdpsvmin = demm/dpsvmin
            desvdpsvmin = detasv/dpsvmin
        #
        vu = (dmu_muvecdp_osv[int(vn*2)]+dmu_muvecdp_osv[int(vn*2)+1])
        dphisvu = abs(vu.Vect().DeltaPhi(osvvec[vn]))
        dphisv1u = abs(dmu_muvecdp_osv[int(vn*2)].Vect().DeltaPhi(osvvec[vn]))
        dphisv2u = abs(dmu_muvecdp_osv[int(vn*2)+1].Vect().DeltaPhi(osvvec[vn]))
        detasvu = abs(vu.Vect().Eta()-osvvec[vn].Eta())
        detadphisvu = 1e6
        if dphisvu>0.0:
            detadphisvu = detasvu/dphisvu
        a3dsvu  = abs(vu.Vect().Angle(osvvec[vn]))
        #
        desvdpsvminu = 1e6
        demmdpsvminu = 1e6
        dpsvminu = 1e6
        tdpsvu = abs(dmu_muvecdp_osv[int(vn*2)].Vect().DeltaPhi(osvvec[vn]))
        if tdpsvu < dpsvminu:
            dpsvminu = tdpsvu
        tdpsvu = abs(dmu_muvecdp_osv[int(vn*2)+1].Vect().DeltaPhi(osvvec[vn]))
        if tdpsvu < dpsvminu:
            dpsvminu = tdpsvu
        if dpsvminu>0.0:
            demmdpsvminu = demmu/dpsvminu
            desvdpsvminu = detasvu/dpsvminu
        #
        sindpsvlxy = abs(ROOT.TMath.Sin(v.Vect().DeltaPhi(osvvec[vn])))*lxy
        sindpsvlxyu = abs(ROOT.TMath.Sin(vu.Vect().DeltaPhi(osvvec[vn])))*lxy
        sina3dsvl3d = abs(ROOT.TMath.Sin(v.Vect().Angle(osvvec[vn])))*(osvvec[vn].Mag())
        sina3dsvl3du = abs(ROOT.TMath.Sin(vu.Vect().Angle(osvvec[vn])))*(osvvec[vn].Mag())
        #
        # Apply selection on angular distances
        if applyDiMuonAngularSel:
            if ROOT.TMath.Log10(dedpmmu)>1.25:
                continue
            if dpmmu>0.9*(ROOT.TMath.Pi()) or a3dmmu>0.9*(ROOT.TMath.Pi()):
                continue
            #if dphisv1u>ROOT.TMath.PiOver2() or dphisv2u>ROOT.TMath.PiOver2() or a3dsvu>ROOT.TMath.PiOver2():
            if dphisvu>ROOT.TMath.PiOver2() or a3dsvu>ROOT.TMath.PiOver2():
            #if dphisvu>0.02 or a3dsvu>ROOT.TMath.PiOver2():
                continue
         #
        #  Apply categorization and selection on muon isolation
        isocat = dimuonIsoCategory(t.Muon_PFIsoAll0p4[dmuidxs_osv[int(vn*2)]], t.Muon_pt[dmuidxs_osv[int(vn*2)]], t.Muon_PFIsoAll0p4[dmuidxs_osv[int(vn*2)+1]], t.Muon_pt[dmuidxs_osv[int(vn*2)+1]])
        if float(args.dimuonIsoCatSel) > -1 and isocat!=float(args.dimuonIsoCatSel):
            continue
        #
        # Check if the dimuon is gen-matched
        isgen = False
        drgen = 9999.
        lxygen = -999.
        idgen = 0
        for gn,g in enumerate(dmumot):
            tmp_drgen = g.DeltaR(v)
            if tmp_drgen < drgen:
                drgen = tmp_drgen 
        if drgen < 0.1:
            isgen = True
            lxygen = t.GenPart_lxy[dmugen[2*gn]]
            idgen = t.GenPart_motherPdgId[dmugen[2*gn]]
        #
        seldmuidxs_osv.append(dmuidxs_osv[int(vn*2)])
        seldmuidxs_osv.append(dmuidxs_osv[int(vn*2)+1])
        seldmusvidxs_osv.append(t.SVOverlap_vtxIdxs[osvidx[vn]][0])
        seldmuvecs_osv.append(vn)
        for h in h1d["dimuon"] + h1d["dimuon_osv"]:
            tn = h.GetName()
            if "dimuon_gen_" in tn and not isgen:
                continue
            if "dimuon_genjpsi_" in tn and not isgen and not idgen==443:
                continue
            h.Fill(eval(variable1d[h.GetName()]), lumiweight)
        for h in h2d["dimuon"] + h2d["dimuon_osv"]:
            tn = h.GetName()
            if "dimuon_gen_" in tn and not isgen:
                continue
            if "dimuon_genjpsi_" in tn and not isgen and not idgen==443:
                continue
            h.Fill(eval(variable2d[h.GetName()][0]),eval(variable2d[h.GetName()][1]), lumiweight)
        # Scan:
        if ( (not filledcat4musep) and (not filledcat4muosv) and (not filledcat2mu) ): 
            slice = "Dimuon_excluded"
            if dphisvu < dphisvcut:
                for l in range(len(lxybins)-1):
                    label = lxybinlabel[l]  
                    if lxy > lxybins[l] and lxy < lxybins[l+1]:
                        break
                if isocat==2:
                    if pt < ptcut:
                        slice = "Dimuon_"+label+"_iso1_ptlow"
                    else:
                        slice = "Dimuon_"+label+"_iso1_pthigh"
                else:
                    if pt < ptcut:
                        slice = "Dimuon_"+label+"_iso0_ptlow"
                    else:
                        slice = "Dimuon_"+label+"_iso0_pthigh"
            else:
                for l in range(len(lxybins)-1):
                    label = lxybinlabel[l]  
                    if lxy > lxybins[l] and lxy < lxybins[l+1]:
                        slice = "Dimuon_"+label+"_non-pointing"
                        break
            mfit.setVal(mass)
            roow.setVal(lumiweight*rooweight);
            roods[slice].add(ROOT.RooArgSet(mfit,roow),roow.getVal());
            catmass[slice].Fill(mass, lumiweight*rooweight);
            roods["Dimuon_"+label+"_inclusive"].add(ROOT.RooArgSet(mfit,roow),roow.getVal());
            roods["Dimuon_full_inclusive"].add(ROOT.RooArgSet(mfit,roow),roow.getVal());
            filledcat2mu = True

    # Fill histograms for selected SVs (with a selected muon pair)
    seldmusvidxs_all = seldmusvidxs_osv+seldmusvidxs
    nSVselass = 0
    nSVselass_osv = 0
    for v in set(seldmusvidxs_all):
        nSVselass = nSVselass+1
        lxy = t.SV_lxy[v]
        for h in h1d["svselass"]:
            tn = h.GetName()
            h.Fill(eval(variable1d[h.GetName()]), lumiweight)
        for h in h2d["svselass"]:
            tn = h.GetName()
            h.Fill(eval(variable2d[h.GetName()][0]),eval(variable2d[h.GetName()][1]), lumiweight)
        if v in seldmusvidxs_osv:
            nSVselass_osv = nSVselass_osv+1
            for h in h1d["svselass_osv"]:
                tn = h.GetName()
                h.Fill(eval(variable1d[h.GetName()]), lumiweight)
            for h in h2d["svselass_osv"]:
                tn = h.GetName()
                h.Fill(eval(variable2d[h.GetName()][0]),eval(variable2d[h.GetName()][1]), lumiweight)
    nSVs = nSVselass
    for h in h1d["nsvselass"]:
        if args.noDiMuon:
            break
        tn = h.GetName()
        h.Fill(eval(variable1d[h.GetName()]), lumiweight)
    nSVs = nSVselass_osv
    for h in h1d["nsvselass_osv"]:
        if args.noDiMuon:
            break
        tn = h.GetName()
        h.Fill(eval(variable1d[h.GetName()]), lumiweight)

    # Fill histograms for muons from selected dimuon and four-muon systems
    selmuidxs_dmu = seldmuidxs+seldmuidxs_osv
    selvecidxs_dmu = seldmuvecs+seldmuvecs_osv
    for m_,m in enumerate(selmuidxs_dmu):
        vn = selvecidxs_dmu[m_//2]
        if m in seldmuidxs:
            lxy  = t.SV_lxy[svidx[vn]]
            mass = (dmuvec[vn]).M()
            pt =  (dmuvec[vn]).Pt()
            phi =  (dmuvec[vn]).Phi()
            ssvphi = t.Muon_vec[m].Vect().DeltaPhi(svvec[vn])
            ssvphiu = dmu_muvecdp[dmuidxs.index(m)].Vect().DeltaPhi(svvec[vn])
            svphi = abs(ssvphi)
            svphiu = abs(ssvphiu)
            mmuphi = t.Muon_vec[m].Vect().DeltaPhi(dmuvec[vn].Vect())
            mmuphiu = dmu_muvecdp[dmuidxs.index(m)].Vect().DeltaPhi(dmuvec[vn].Vect())
            mmtheta = abs(mmuphi - math.copysign(ROOT.TMath.PiOver2(), ssvphi)) 
            mmthetau = abs(mmuphiu - math.copysign(ROOT.TMath.PiOver2(), ssvphiu)) 
            mmtheta_mmsign = abs(mmuphi - math.copysign(ROOT.TMath.PiOver2(), mmuphi)) 
            mmthetau_mmsign = abs(mmuphiu - math.copysign(ROOT.TMath.PiOver2(), mmuphiu)) 
        else:
            lxy = t.SVOverlap_lxy[osvidx[vn]]
            mass = (dmuvec_osv[vn]).M()
            pt =  (dmuvec_osv[vn]).Pt()
            phi =  (dmuvec_osv[vn]).Phi()
            ssvphi = t.Muon_vec[m].Vect().DeltaPhi(osvvec[vn])
            ssvphiu = dmu_muvecdp_osv[dmuidxs_osv.index(m)].Vect().DeltaPhi(osvvec[vn])
            svphi = abs(ssvphi)
            svphiu = abs(ssvphiu)
            mmuphi = t.Muon_vec[m].Vect().DeltaPhi(dmuvec_osv[vn].Vect())
            mmuphiu = dmu_muvecdp_osv[dmuidxs_osv.index(m)].Vect().DeltaPhi(dmuvec_osv[vn].Vect())
            mmtheta = abs(mmuphi - math.copysign(ROOT.TMath.PiOver2(), ssvphi))
            mmthetau = abs(mmuphiu - math.copysign(ROOT.TMath.PiOver2(), ssvphiu))
            mmtheta_mmsign = abs(mmuphi - math.copysign(ROOT.TMath.PiOver2(), mmuphi)) 
            mmthetau_mmsign = abs(mmuphiu - math.copysign(ROOT.TMath.PiOver2(), mmuphiu)) 
        dxysign = getIPSign(t.Muon_phiCorr[m], phi)
        for h in h1d["selmuon"]:
            tn = h.GetName()
            h.Fill(eval(variable1d[h.GetName()]), lumiweight)
        for h in h2d["selmuon"]:
            tn = h.GetName()
            h.Fill(eval(variable2d[h.GetName()][0]),eval(variable2d[h.GetName()][1]), lumiweight)
        if m in seldmuidxs_osv:
            for h in h1d["selmuon_osv"]:
                tn = h.GetName()
                h.Fill(eval(variable1d[h.GetName()]), lumiweight)
            for h in h2d["selmuon_osv"]:
                tn = h.GetName()
                h.Fill(eval(variable2d[h.GetName()][0]),eval(variable2d[h.GetName()][1]), lumiweight)

    if len(qmuidxs_sel)>3:
        for m_,m in enumerate(qmuidxs_sel):
            vn = selqmuvecs[m_//4]
            if m in qmuidxsminlxy:
                lxy = svvecminlxy_qmu[vn].Perp()
                pt = (qmu_dmuvecminlxy[vn]).Pt()
                mass = (qmu_dmuvecminlxy[vn]).M()
            else:
                lxy = svvecmaxlxy_qmu[vn].Perp()
                pt = (qmu_dmuvecmaxlxy[vn]).Pt()
                mass = (qmu_dmuvecmaxlxy[vn]).M()
            for h in h1d["selmuon_fourmu"]:
                tn = h.GetName()
                h.Fill(eval(variable1d[h.GetName()]), lumiweight)
            for h in h2d["selmuon_fourmu"]:
                tn = h.GetName()
                h.Fill(eval(variable2d[h.GetName()][0]),eval(variable2d[h.GetName()][1]), lumiweight)

    if len(qmuidxs_osv_sel)>3:
        for m_,m in enumerate(qmuidxs_osv_sel):
            vn = selqmuvecs_osv[m_//4]
            lxy  = t.SVOverlap_lxy[osvidx_qmu[vn]]
            mass = (qmuvec_osv[vn]).M()
            pt = (qmuvec_osv[vn]).Pt()
            for h in h1d["selmuon_fourmu_osv"]:
                tn = h.GetName()
                h.Fill(eval(variable1d[h.GetName()]), lumiweight)
            for h in h2d["selmuon_fourmu_osv"]:
                tn = h.GetName()
                h.Fill(eval(variable2d[h.GetName()][0]),eval(variable2d[h.GetName()][1]), lumiweight)


### Write histograms
foname = "%s/histograms_%s_all.root"%(outdir,args.year)
if args.inSample!="*":
    if args.inFile!="*":
        foname = "%s/histograms_file%s_%s_%s"%(outdir,args.inFile,args.inSample,args.year)
    else:
        foname = "%s/histograms_%s_%s"%(outdir,args.inSample,args.year)
if index>=0:
    foname = foname+("_%d"%index)
fout = ROOT.TFile(foname+".root","RECREATE")
fout.cd()
for cat in h1d.keys():
    for h in h1d[cat]:
        h.Write()
for cat in h2d.keys():
    for h in h2d[cat]:
        h.Write()
for dbin in dbins:
    print("RooDataSet {}  with {} entries".format(dbin, roods[dbin].sumEntries()))
    roods[dbin].Write()
    catmass[dbin].Write()
fout.Close()
