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
import correctionlib
import awkward as ak


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
parser.add_argument("--unblind", default=False, action="store_true", help="Unblind data")
parser.add_argument("--year", default="2022", help="Year to be processed. Default: 2022")
parser.add_argument("--weightMC", default=True, help="Indicate if MC is weighted")
parser.add_argument("--weightB", default="True", help="Indicate if BToPhi MC is reweighted for b-hadron pt")
parser.add_argument("--doPUreweighting", default=False, action="store_true", help="Apply PU reweighting")
parser.add_argument("--reweightFrom", default=-1, help="Indicate ctau of the sample")
parser.add_argument("--reweightTo", default=-1, help="Indicate ctau to reweight to")
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
parser.add_argument("--noDiMuonResonanceMasking", default=False, action="store_true", help="Do not apply dimuon resonance masking in dimuon events")
parser.add_argument("--noFourMuonResonanceMasking", default=False, action="store_true", help="Do not apply dimuon resonance masking in dimuon events")
parser.add_argument("--noFourMuonAngularSel", default=False, action="store_true", help="Do not apply selection on fourmuon angular variables")
parser.add_argument("--noFourMuonIPSel", default=False, action="store_true", help="Do not apply selection on muon IP (Not applied at four-muon level)")
parser.add_argument("--noFourMuonMassDiffSel", default=False, action="store_true", help="Do not apply selection on fourmuon invariant mass difference")
parser.add_argument("--noPreSel", default=False, action="store_true", help="Do not fill pre-selection/association histograms")
parser.add_argument("--noMuonChi2Sel", default=False, action="store_true", help="Drop the muon normChi2<3.0 cut (keep only pt>3.0 and |eta|<2.4) for both Muon and Muon_vtx collections")
parser.add_argument("--noDiMuon", default=False, action="store_true", help="Do not fill dimuon histograms")
parser.add_argument("--noFourMuon", default=False, action="store_true", help="Do not fill four-muon histograms for four-muon systems")
parser.add_argument("--noFourMuonOSV", default=False, action="store_true", help="Do not fill four-muon histograms for four-muon systems from overlapping SVs")
parser.add_argument("--noOverlappingSV", default=False, action="store_true", help="Do not use overlapping vertices for anything and just use simple SVs")
parser.add_argument("--noSeed", default=[], nargs="+", help="Exclude L1 seeds from the acceptance")
parser.add_argument("--requireL1", default=False, action="store_true", help="Require the event to pass L1 (looper 'passL1' branch) before filling. Apply consistently to data and MC.")
parser.add_argument("--tupleDir", default="", help="Output subdirectory name for the tuples (under STARTDIR); default tuples_parking[_nochi2] per --noMuonChi2Sel.")
parser.add_argument("--noHistos", default=False, action="store_true", help="Skip histogram filling")
parser.add_argument("--noSF", default=False, action="store_true", help="Dont apply SF")
parser.add_argument("--doGen", default=False, action="store_true", help="Fill generation information histograms")
parser.add_argument("--collection", default="OR", choices=["noVtx", "Vtx", "both", "OR"], help="Vertex collection to use: 'noVtx', 'Vtx', 'both' (union), or 'OR' (Vtx-preferred, unique muons via dR matching)")
args = parser.parse_args()
use_novtx = args.collection in ('noVtx', 'both', 'OR')
use_vtx   = args.collection in ('Vtx', 'both', 'OR')
use_OR    = (args.collection == 'OR')
OR_SV_DLXY_MAX = 0.1

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
        print(f"Getting {era}: Normalize to {frac*8.064921449}") # Leonardo: 8.077046947
        return frac*8.077046947*xsec/ngen
    if era=="2022postEE":
        print(f"Getting {era}: Normalize to {frac*26.658464830}") # Leonardo: 26.982330931
        return frac*26.982330931*xsec/ngen
    if era=="2023postEE":
        print(f"Getting {era}: Normalize to {frac*5.613751282}") # Leonardo: 5.557004785
        return frac*5.557004785*xsec/ngen # Total C: 17.060484313 
    if era=="2023":
        print(f"Getting {era}: Normalize to {frac*11.503479528}") # Leonardo: 11.503479528
        return frac*11.503479528*xsec/ngen # Total C: 17.060484313
    if era=="2023BPix":
        print(f"Getting {era}: Normalize to {frac*9.525199061}")
        return frac*9.525199061*xsec/ngen
    if era=="2024":
        print(f"Getting {era}: Normalize to {frac*109.95}")
        return frac*109.95*xsec/ngen

# Get closest from list
def getClosest(val, values):
    iv = 0
    for v,vv in enumerate(values):
        if abs(val-vv) < abs(val-values[iv]):
            iv = v
    return iv, abs(val-values[iv])
#
def getClosestAngular(vals = [0,0,0], values_lxy = [], values_eta = [], values_phi = []):
    iv = 0
    ivt = ROOT.TVector3()
    ivt.SetPtEtaPhi(values_lxy[0], values_eta[0], values_phi[0])
    v1 = ROOT.TVector3()
    v1.SetPtEtaPhi(vals[0], vals[1], vals[2])
    for v in range(0, len(values_lxy)):
        v2 = ROOT.TVector3()
        v2.SetPtEtaPhi(values_lxy[v], values_eta[v], values_phi[v])
        #if abs(val-vv) < abs(val-values[iv]):
        if v1.DeltaR(v2) < v1.DeltaR(ivt):
            iv = v
            ivt = v2
    return iv


# Get trigger SF (2022)
def getTriggerSF(subpt, lxy):
    sf = 1.0
    sfup = 1.0
    sfdown = 1.0
    if lxy > 0.0 and lxy < 0.2:
        if subpt < 5.0:
            sf = 1.12
            sfup = 1.12 + 0.01
            sfdown = 1.12 - 0.01
        elif subpt > 5. and subpt < 10.:
            sf = 0.996
            sfup = 0.996 + 0.004
            sfdown = 0.996 - 0.003
        elif subpt > 10. and subpt < 15.:
            sf = 0.974
            sfup = 0.974 + 0.007
            sfdown = 0.974 - 0.008
        elif subpt > 15.0:
            sf = 0.915
            sfup = 0.915 + 0.008
            sfdown = 0.915 - 0.008
    elif lxy > 0.2 and lxy < 1.0:
        if subpt < 5.0:
            sf = 1.19
            sfup = 1.19 + 0.04
            sfdown = 1.19 - 0.05
        elif subpt > 5. and subpt < 10.:
            sf = 1.008
            sfup = 1.008 + 0.007
            sfdown = 1.008 - 0.008
        elif subpt > 10. and subpt < 15.:
            sf = 0.969
            sfup = 0.969 + 0.006
            sfdown = 0.969 - 0.008
        elif subpt > 15.0:
            sf = 0.96
            sfup = 0.96 + 0.01
            sfdown = 0.96 - 0.01
    elif lxy > 1.0 and lxy < 2.4:
        if subpt < 5.0:
            sf = 1.32
            sfup = 1.32 + 0.1
            sfdown = 1.32 - 0.2
        elif subpt > 5. and subpt < 10.:
            sf = 1.03
            sfup = 1.03 + 0.02
            sfdown = 1.03 - 0.02
        elif subpt > 10. and subpt < 15.:
            sf = 0.986
            sfup = 0.986 + 0.008
            sfdown = 0.986 - 0.017
        elif subpt > 15.0:
            sf = 0.98
            sfup = 0.98 + 0.01
            sfdown = 0.98 - 0.04
    elif lxy > 2.4 and subpt < 5.0:
        sf = 1.00
        sfup = 1.00 + 0.3
        sfdown = 1.00 - 0.3
    elif lxy > 2.4 and lxy < 3.1:
        sf = 1.04
        sfup = 1.04 + 0.01
        sfdown = 1.04 - 0.02
    elif lxy > 3.1 and lxy < 7.0:
        sf = 1.02
        sfup = 1.02 + 0.01
        sfdown = 1.02 - 0.03
    elif lxy > 7.0:
        sf = 1.00
        sfup = 1.00 + 0.05
        sfdown = 1.00 - 0.05
    return sf, sfup, sfdown

# Get Selection SF
def getSelectionSF(lxy):
    sf = 1.0
    sfup = 1.0
    sfdown = 1.0
    if lxy > 0.0 and lxy < 0.2:
        sf = 1.10
        sfup = 1.10 + 0.01
        sfdown = 1.10 - 0.01
    elif lxy > 0.2 and lxy < 1.0:
        sf = 1.05
        sfup = 1.05 + 0.02
        sfdown = 1.05 - 0.02
    elif lxy > 1.0 and lxy < 2.4:
        sf = 1.03
        sfup = 1.03 + 0.03
        sfdown = 1.03 - 0.03
    elif lxy > 2.4 and lxy < 3.1:
        sf = 0.97
        sfup = 0.97 + 0.06
        sfdown = 0.97 - 0.06
    elif lxy > 3.1 and lxy < 7.0:
        sf = 0.82
        sfup = 0.82 + 0.1
        sfdown = 0.82 - 0.1
    elif lxy > 7.0: # This weight is applied later with datacards are filled.
        sf = 1.0
        sfup = 1.0 + 0.3
        sfdown = 1.0 - 0.3
    return sf, sfup, sfdown

## Running settings
indir  = args.inDir.replace("/ceph/cms","")
is_local = os.path.isdir(args.inDir)
outdir = args.outDir
if args.outSuffix!="":
    outdir = outdir+"_"+args.outSuffix
if not os.path.exists(outdir):
    os.makedirs(outdir)

applyMaterialVeto = not args.noMaterialVeto
applyMuonIPSel = False            # moved to cut-and-count
applyDiMuonAngularSel = False     # moved to cut-and-count
applyMuonHitSel = False           # moved to cut-and-count
applyDiMuonResonanceMasking = False  # moved to cut-and-count
applyFourMuonResonanceMasking = not args.noFourMuonResonanceMasking
applyFourMuonIPSel = not args.noFourMuonIPSel
applyFourMuonAngularSel = not args.noFourMuonAngularSel
applyFourMuonMassDiffSel = not args.noFourMuonMassDiffSel
applyPUreweighting = args.doPUreweighting
applyBweights = eval(args.weightB)
doOverlappingSV = not args.noOverlappingSV

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
    prependtodir = "" if is_local else "/ceph/cms"
else:
    prependtodir = "root://redir.t2.ucsd.edu:1095"
if not args.condor:
    cephbase = args.inDir if is_local else "/ceph/cms%s"%indir
    if args.inFile!="*" and args.inSample!="*":
        thisfile="output_%s_%s_%s.root"%(args.inSample,args.year,args.inFile)
        if os.path.isfile("%s/%s"%(cephbase,thisfile)):
            files.append("%s%s/%s"%(prependtodir,indir,thisfile))
    elif args.inSample!="*":
        for f in os.listdir(cephbase):
            print("output_%s_%s_"%(args.inSample,args.year))
            if ("output_%s_%s_"%(args.inSample,args.year) in f) and os.path.isfile("%s/%s"%(cephbase,f)):
                files.append("%s%s/%s"%(prependtodir,indir,f))
    else:
        for f in os.listdir(cephbase):
            if (args.year in f) and (".root" in f) and os.path.isfile("%s/%s"%(cephbase,f)):
                files.append("%s%s/%s"%(prependtodir,indir,f))
else:
    os.system('xrdfs redir.t2.ucsd.edu:1095 ls %s > filein.txt'%indir)
    fin = open("filein.txt","r")
    for f in fin.readlines():
        f = f.strip("\n")
        if args.inFile!="*" and args.inSample!="*":
            thisfile="output_%s_%s_%s.root"%(args.inSample,args.year,args.inFile)
            if thisfile in f:
                files.append("%s/%s"%(prependtodir,f))
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
unblind_frac = 1.0 if args.unblind else 0.1

if not isData and args.weightMC and 'MinBias' not in args.inSample and 'QCD' not in args.inSample:
    counts = ROOT.TH1F("totals", "", 1, 0, 1)
    print("Simulations: Getting counts")
    if "HTo2ZdTo2mu2x" in sampleTag:
        for _,f in enumerate(files):
            if not args.condor:
                f_ = ROOT.TFile.Open(f.replace('davs://redir.t2.ucsd.edu:1095//', '/ceph/cms/'))
            else:
                f_ = ROOT.TFile.Open(f)
            h_ = f_.Get("counts").Clone("Clone_{}".format(_))
            counts.Add(h_)
            f_.Close()
        ncounts = counts.GetBinContent(1)
        hahm_csv = 'data/hahm-request_2024.csv' if args.year == "2024" else 'data/hahm-request.csv'
        efilter = None
        with open(hahm_csv) as mcinfo:
            reader = csv.reader(mcinfo, delimiter=',')
            for row in reader:
                if sampleTag in row[0]:
                    efilter = float(row[-1])
                    break
        if efilter is None:
            sys.exit(f"ERROR: no filter efficiency for {sampleTag} in {hahm_csv}")
    if 'BToPhi' in sampleTag:
        for _,f in enumerate(files):
            if not args.condor:
                f_ = ROOT.TFile.Open(f.replace('davs://redir.t2.ucsd.edu:1095//', '/ceph/cms/'))
            else:
                f_ = ROOT.TFile.Open(f)
            h_ = f_.Get("counts").Clone("Clone_{}".format(_))
            counts.Add(h_)
            f_.Close()
        ncounts = counts.GetBinContent(1)
        efilter = 1.
        #with open('data/BToPhi-request.csv') as mcinfo:
        #    reader = csv.reader(mcinfo, delimiter=',')
        #    for row in reader:
        #        if sampleTag.replace('_MPhi','') in row[0]:
        #            efilter = float(row[-1])
        #            break
    if "ScenarioB1" in sampleTag:
        for _,f in enumerate(files):
            if not args.condor:
                f_ = ROOT.TFile.Open(f.replace('davs://redir.t2.ucsd.edu:1095//', '/ceph/cms/'))
            else:
                f_ = ROOT.TFile.Open(f)
            h_ = f_.Get("counts").Clone("Clone_{}".format(_))
            counts.Add(h_)
            f_.Close()
        ncounts = counts.GetBinContent(1)
        with open('data/DQCD-request.csv') as mcinfo:
            reader = csv.reader(mcinfo, delimiter=',')
            for row in reader:
                if sampleTag in row[0]:
                    efilter = float(row[-1])
                    break
    if "ScenarioA" in sampleTag:
        for _,f in enumerate(files):
            if not args.condor:
                f_ = ROOT.TFile.Open(f.replace('davs://redir.t2.ucsd.edu:1095//', '/ceph/cms/'))
            else:
                f_ = ROOT.TFile.Open(f)
            h_ = f_.Get("counts").Clone("Clone_{}".format(_))
            counts.Add(h_)
            f_.Close()
        ncounts = counts.GetBinContent(1)
        with open('data/DQCD-request.csv') as mcinfo:
            reader = csv.reader(mcinfo, delimiter=',')
            for row in reader:
                if sampleTag in row[0]:
                    efilter = float(row[-1])
                    break
    if "2022postEE" in files[0] and args.year=="2022":
        lumiweight = getweight("2022postEE", ncounts/efilter, unblind_frac)
    elif "2022" in files[0] and args.year=="2022":
        lumiweight = getweight("2022", ncounts/efilter, unblind_frac)
    if "2022postEE" in files[0] and args.year=="2023":
        lumiweight = getweight("2023postEE", ncounts/efilter, unblind_frac)
    elif "2023BPix" in files[0] and args.year=="2023":
        lumiweight = getweight("2023BPix", ncounts/efilter, unblind_frac)
    elif "2023" in files[0] and args.year=="2023":
        lumiweight = getweight("2023", ncounts/efilter, unblind_frac)
    elif "2024" in files[0] and args.year=="2024":
        lumiweight = getweight("2024", ncounts/efilter, unblind_frac)

    print("Total number of counts: {}".format(ncounts))
    print("Filter efficiency (generation): {}".format(efilter))
    print("Lumiweight: {}".format(lumiweight))

    if applyBweights and 'BToPhi' in sampleTag:
        print("B-hadron reweighting is being applied")

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

# Lifetime-reweighting
reweightFrom = float(args.reweightFrom)
reweightTo = float(args.reweightTo)

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
mfit = ROOT.RooRealVar("mfit", "mfit", 0.1, 140.0)
m4fit = ROOT.RooRealVar("m4fit", "m4fit", 0.1, 140.0)
roow = ROOT.RooRealVar("roow", "roow", -10000.0, 10000.0)
roow_trg_up = ROOT.RooRealVar("roow_trg_up", "roow_trg_up", -10000.0, 10000.0)
roow_trg_down = ROOT.RooRealVar("roow_trg_down", "roow_trg_down", -10000.0, 10000.0)
roow_sel_up = ROOT.RooRealVar("roow_sel_up", "roow_sel_up", -10000.0, 10000.0)
roow_sel_down = ROOT.RooRealVar("roow_sel_down", "roow_sel_down", -10000.0, 10000.0)
roow4 = ROOT.RooRealVar("roow", "roow", -10000.0, 10000.0)
roow4_trg_up = ROOT.RooRealVar("roow_trg_up", "roow_trg_up", -10000.0, 10000.0)
roow4_trg_down = ROOT.RooRealVar("roow_trg_down", "roow_trg_down", -10000.0, 10000.0)
roow4_sel_up = ROOT.RooRealVar("roow_sel_up", "roow_sel_up", -10000.0, 10000.0)
roow4_sel_down = ROOT.RooRealVar("roow_sel_down", "roow_sel_down", -10000.0, 10000.0)
roods = {}
roods_trg_up = {}
roods_trg_down = {}
roods_sel_up = {}
roods_sel_down = {}
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
        roods_trg_up[dbin] = ROOT.RooDataSet(dname + "_trg_up",dname,ROOT.RooArgSet(mfit,roow_trg_up),"roow_trg_up")
        roods_trg_down[dbin] = ROOT.RooDataSet(dname + "_trg_down",dname,ROOT.RooArgSet(mfit,roow_trg_down),"roow_trg_down")
        roods_sel_up[dbin] = ROOT.RooDataSet(dname + "_sel_up",dname,ROOT.RooArgSet(mfit,roow_sel_up),"roow_sel_up")
        roods_sel_down[dbin] = ROOT.RooDataSet(dname + "_sel_down",dname,ROOT.RooArgSet(mfit,roow_sel_down),"roow_sel_down")
    else:
        catmass[dbin] = ROOT.TH1F(dname + "_rawmass","; m_{4#mu} [GeV]; Events / 0.01 GeV",2500, 0., 250.)
        roods[dbin] = ROOT.RooDataSet(dname,dname,ROOT.RooArgSet(m4fit,roow4),"roow")
        roods_trg_up[dbin] = ROOT.RooDataSet(dname + "_trg_up",dname,ROOT.RooArgSet(m4fit,roow4_trg_up),"roow_trg_up")
        roods_trg_down[dbin] = ROOT.RooDataSet(dname + "_trg_down",dname,ROOT.RooArgSet(m4fit,roow4_trg_down),"roow_trg_down")
        roods_sel_up[dbin] = ROOT.RooDataSet(dname + "_sel_up",dname,ROOT.RooArgSet(m4fit,roow4_sel_up),"roow_sel_up")
        roods_sel_down[dbin] = ROOT.RooDataSet(dname + "_sel_down",dname,ROOT.RooArgSet(m4fit,roow4_sel_down),"roow_sel_down")
#
catmass["FourMu_osv_dimuonmass"] = ROOT.TH1F("d_FourMu_osv_dimuonmass_rawmass","; m_{#mu#mu} [GeV]; Events / 0.01 GeV",1500, 0., 150.)
catmass["FourMu_sep_dimuonmass"] = ROOT.TH1F("d_FourMu_sep_dimuonmass_rawmass","; m_{#mu#mu} [GeV]; Events / 0.01 GeV",1500, 0., 150.)


#
#
# Parquet branches
branches = {}
#
SV_BRANCHES = [
    "ptmm", "chi2", "prob", "x", "y", "z", "lxy",
    "xErr", "yErr", "zErr", "dphi_mumu_SV", "d3d_mumu_SV", "a3d_mumu",
    "dr_mumu", "dphi_mumu", "deta_mumu",
    "dphisv_mu1", "dphisv_mu2", "deta_mumu_SV",
    "sina3d_l3d", "sindphi_lxy",
    "mu1_pt", "mu1_eta", "mu1_phi", "mu1_isvtx", "mu1_normChi2", "mu1_dxy", "mu1_dxyErr",
    "mu1_dxysig", "mu1_dz", "mu1_dze", "mu1_dzsig", "mu1_nhitsbeforesv",
    "mu1_isGlobal", "mu1_isTracker", "mu1_isStandAlone",
    "mu1_pixHits", "mu1_stripHits", "mu1_pixLayers", "mu1_trkLayers",
    "mu1_saHits", "mu1_saMatchedStats",
    "mu1_muHits", "mu1_muChambs", "mu1_muCSCDT", "mu1_muMatch", "mu1_muMatchedStats", "mu1_muExpMatchedStats", "mu1_muMatchedRPC",
    "mu1_ecalIso", "mu1_hcalIso", "mu1_trackIso", "mu1_ecalRelIso", "mu1_hcalRelIso", "mu1_trackRelIso",
    "mu1_PFIsoChg0p3", "mu1_PFIsoAll0p3", "mu1_PFRelIsoChg0p3", "mu1_PFRelIsoAll0p3", "mu1_mindrPF0p3",
    "mu1_PFIsoChg0p4", "mu1_PFIsoAll0p4", "mu1_PFRelIsoChg0p4", "mu1_PFRelIsoAll0p4", "mu1_mindrPF0p4",
    "mu1_mindr", "mu1_maxdr", "mu1_mindrJet", "mu1_mindphiJet", "mu1_mindetaJet",
    "mu1_ncompatible", "mu1_ncompatibletotal",
    "mu1_nexpectedhits", "mu1_nexpectedhitsmultiple", "mu1_nexpectedhitsmultipletotal", "mu1_nexpectedhitstotal",
    "mu1_phiCorr",
    "mu2_pt", "mu2_eta", "mu2_phi", "mu2_isvtx", "mu2_normChi2", "mu2_dxy", "mu2_dxyErr",
    "mu2_dxysig", "mu2_dz", "mu2_dze", "mu2_dzsig", "mu2_nhitsbeforesv",
    "mu2_isGlobal", "mu2_isTracker", "mu2_isStandAlone",
    "mu2_pixHits", "mu2_stripHits", "mu2_pixLayers", "mu2_trkLayers",
    "mu2_saHits", "mu2_saMatchedStats",
    "mu2_muHits", "mu2_muChambs", "mu2_muCSCDT", "mu2_muMatch", "mu2_muMatchedStats", "mu2_muExpMatchedStats", "mu2_muMatchedRPC",
    "mu2_ecalIso", "mu2_hcalIso", "mu2_trackIso", "mu2_ecalRelIso", "mu2_hcalRelIso", "mu2_trackRelIso",
    "mu2_PFIsoChg0p3", "mu2_PFIsoAll0p3", "mu2_PFRelIsoChg0p3", "mu2_PFRelIsoAll0p3", "mu2_mindrPF0p3",
    "mu2_PFIsoChg0p4", "mu2_PFIsoAll0p4", "mu2_PFRelIsoChg0p4", "mu2_PFRelIsoAll0p4", "mu2_mindrPF0p4",
    "mu2_mindr", "mu2_maxdr", "mu2_mindrJet", "mu2_mindphiJet", "mu2_mindetaJet",
    "mu2_ncompatible", "mu2_ncompatibletotal",
    "mu2_nexpectedhits", "mu2_nexpectedhitsmultiple", "mu2_nexpectedhitsmultipletotal", "mu2_nexpectedhitstotal",
    "mu2_phiCorr",
    "mass",
    "ndof", "chi2Ndof", "l3d",
    "mindx", "mindy", "mindz", "mindxy", "mind3d",
    "maxdx", "maxdy", "maxdz", "maxdxy", "maxd3d",
    "onModule", "onModuleWithinUnc",
    "minDistanceFromDet", "minDistanceFromDet_x", "minDistanceFromDet_y", "minDistanceFromDet_z",
    "closestDet_x", "closestDet_y", "closestDet_z",
]
#
for sv in ["SV1", "SV2"]:
    for var in SV_BRANCHES:
        branches[f"{sv}_{var}"] = []

branches["passL1"] = []
EVENT_COUNT_BRANCHES = ["nMuons_NoVtx", "nMuons_Vtx", "nSVs_NoVtx", "nSVs_Vtx"]
for _cb in EVENT_COUNT_BRANCHES:
    branches[_cb] = []
EVENT_INFO_BRANCHES = ["run", "lumi", "evtn", "nPUTrue", "wPU"]
for _eb in EVENT_INFO_BRANCHES:
    branches[_eb] = []
GEN_B_BRANCHES = ["nGenBplus"]
for _gb in GEN_B_BRANCHES:
    branches[_gb] = []
HAS_GEN_B = bool(t.GetBranch("nGenBplus"))



def muattr(t, attr, idx, is_vtx):
    """Get muon attribute from the appropriate vertex collection."""
    if is_vtx:
        return getattr(t, 'Muon_vtx_' + attr)[idx]
    return getattr(t, 'Muon_' + attr)[idx]

def mu_nhits(t, idx, is_vtx):
    """Return nhitsbeforesv for the noVtx or Vtx collection."""
    if is_vtx:
        return t.Muon_vtx_nhitsbeforesv[idx]
    return t.Muon_nhitsbeforesv[idx]

def muon_ip_pv(t, idx, is_vtx):
    phi = muattr(t, 'phi', idx, is_vtx)
    eta = muattr(t, 'eta', idx, is_vtx)
    dxy_pv = muattr(t, 'dxy', idx, is_vtx) + t.PV_x*np.sin(phi) - t.PV_y*np.cos(phi)
    dz_pv = muattr(t, 'dz', idx, is_vtx) - t.PV_z + (t.PV_x*np.cos(phi) + t.PV_y*np.sin(phi))*np.sinh(eta)
    return dxy_pv, dz_pv

def vtx_muvec_corr(t, idx):
    vec = ROOT.TLorentzVector()
    vec.SetPtEtaPhiM(t.Muon_vtx_pt[idx], t.Muon_vtx_eta[idx], t.Muon_vtx_phiCorr[idx], MUON_MASS)
    return vec

def svattr(t, attr, idx, is_vtx):
    """Get SV attribute from the appropriate vertex collection."""
    if is_vtx:
        return getattr(t, 'SV_vtx_' + attr)[idx]
    return getattr(t, 'SV_' + attr)[idx]

# Event loop
print("From event %d to event %d"%(firste,laste))
for e in range(firste,laste):
    # Init parquet branches
    for sv in ["SV1", "SV2"]:
        for var in SV_BRANCHES:
            branches[f"{sv}_{var}"].append(-1.)
    branches["passL1"].append(-1.)   # event-level; overwritten with the real value below
    for _cb in EVENT_COUNT_BRANCHES:
        branches[_cb].append(-1.)    # event-level; overwritten with the real value below
    for _eb in EVENT_INFO_BRANCHES:
        branches[_eb].append(-1.)    # event-level; overwritten with the real value below
    for _gb in GEN_B_BRANCHES:
        branches[_gb].append(-1.)

    # Access event
    t.GetEntry(e)
    #if e%1000==0:
    #    print("At entry %d"%e)
    # Store the looper L1 decision per event so it can be cut on downstream.
    branches["passL1"][-1] = float(t.passL1)
    for _cb in EVENT_COUNT_BRANCHES:
        branches[_cb][-1] = float(getattr(t, _cb))
    for _eb in EVENT_INFO_BRANCHES:
        branches[_eb][-1] = float(getattr(t, _eb))
    if HAS_GEN_B:
        for _gb in GEN_B_BRANCHES:
            branches[_gb][-1] = float(getattr(t, _gb))
    if args.requireL1 and not t.passL1:
        continue
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

    # PU reweighting
    puweight = 1.
    if applyPUreweighting and not isData:
        puweight = t.wPU

    # B-hadron reweighting
    bweight = 1.0
    if 'BToPhi' in sampleTag and applyBweights and not isData:
        cset = correctionlib.CorrectionSet.from_file("data/b-hadron_weights.json")
        bweight = cset['PT_weight'].evaluate(t.GenB_pt, 'nominal')

    # Event info
    event_weight = lumiweight*bweight*puweight

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
            if isResonance:
                for h in h1d["genmu"]:
                    tn = h.GetName()
                    h.Fill(eval(variable1d[h.GetName()]), lumiweight*bweight*puweight)
        for g,gp in enumerate(dmumot):
            lxygen = t.GenPart_lxy[dmugen[2*g]] 
            if t.GenPart_motherPdgId[dmugen[2*g]]==443:
                for h in h1d["jpsi"]:
                    tn = h.GetName()
                    h.Fill(eval(variable1d[h.GetName()]), lumiweight*bweight*puweight)
            for h in h1d["llp"]:
                tn = h.GetName()
                h.Fill(eval(variable1d[h.GetName()]), lumiweight*bweight*puweight)
        # Identify the LLPs
        LLPs = []
        LLPs_lxy = []
        LLPs_ct = []
        LLPs_eta = []
        LLPs_phi = []
        if (reweightTo > 0 and reweightFrom > 0):
            for i in range(0, len(t.GenPart_pdgId)):
                if t.GenPart_pdgId[i] not in [1023, 9900015, 6000211]:
                    continue
                for j in range(0, len(t.GenPart_pdgId)):
                    if (t.GenPart_pdgId[j]==13 and t.GenPart_motherIndex[j]==t.GenPart_index[i]):
                        LLPs.append(i)
                        LLPs_lxy.append(t.GenPart_lxy[j]) # from the muon
                        LLPs_eta.append(t.GenPart_eta[i])
                        LLPs_phi.append(t.GenPart_phi[i])
                        LLPs_ct.append(t.GenPart_ct[i])
                        break

    # Loop over SVs (noVtx)
    nSV_noVtx = len(t.SV_index)
    if use_novtx and not use_vtx and nSV_noVtx<1:
        continue
    nSVsel_noVtx = 0
    for v in range(nSV_noVtx):
        if args.noPreSel:
            break
        if not t.SV_selected[v]:
            continue
        if applyMaterialVeto and (t.SV_onModuleWithinUnc[v] or (abs(t.SV_minDistanceFromDet_x[v]) < 0.81 and abs(t.SV_minDistanceFromDet_y[v]) < 3.24 and abs(t.SV_minDistanceFromDet_z[v]) < 0.0145)):
            continue
        nSVsel_noVtx = nSVsel_noVtx + 1
        lxy = t.SV_lxy[v]
    nSVs_noVtx = nSVsel_noVtx

    # Loop over muons (noVtx)
    nMu_noVtx = len(t.Muon_selected)
    nMuSel_noVtx = 0
    nMuAss_noVtx = 0
    nMuAssOverlap_noVtx = 0
    muselidxs_noVtx = []
    for m in range(nMu_noVtx):
        if args.noMuonChi2Sel:
            if not (t.Muon_pt[m] > 3.0 and abs(t.Muon_eta[m]) < 2.4):
                continue
        elif not t.Muon_selected[m]:
            continue
        nMuSel_noVtx = nMuSel_noVtx+1
        if t.Muon_bestAssocSVOverlapIdx[m]>-1:
            if applyMaterialVeto:
                _ovsv = t.SVOverlap_vtxIdxs[t.Muon_bestAssocSVOverlapIdx[m]][0]
                if t.SV_onModuleWithinUnc[_ovsv] or (abs(t.SV_minDistanceFromDet_x[_ovsv]) < 0.81 and abs(t.SV_minDistanceFromDet_y[_ovsv]) < 3.24 and abs(t.SV_minDistanceFromDet_z[_ovsv]) < 0.0145):
                    continue
            nMuAss_noVtx = nMuAss_noVtx+1
            nMuAssOverlap_noVtx = nMuAssOverlap_noVtx+1
        elif t.Muon_bestAssocSVIdx[m]>-1:
            if applyMaterialVeto:
                _vidx = t.Muon_bestAssocSVIdx[m]
                _sv_on_module = False
                for _v in range(len(t.SV_index)):
                    if t.SV_index[_v] == _vidx:
                        _sv_on_module = bool(t.SV_onModuleWithinUnc[_v] or (abs(t.SV_minDistanceFromDet_x[_v]) < 0.81 and abs(t.SV_minDistanceFromDet_y[_v]) < 3.24 and abs(t.SV_minDistanceFromDet_z[_v]) < 0.0145))
                        break
                if _sv_on_module:
                    continue
            nMuAss_noVtx = nMuAss_noVtx+1
        else:
            continue
        muselidxs_noVtx.append(m)
        if args.noPreSel:
            continue

    # Select events with at least two muons associated to a SV (noVtx)
    if use_novtx and not use_vtx and nMuAss_noVtx < 2:
        continue

    # Muon pairing
    dmuvec_noVtx = []
    dmu_muvecdp_noVtx = []
    dmuidxs_noVtx = []
    svvec_noVtx = []
    svidx_noVtx = []
    dmuvec_osv_noVtx = []
    dmu_muvecdp_osv_noVtx = []
    dmuidxs_osv_noVtx = []
    osvvec_noVtx = []
    osvidx_noVtx = []
    qmuvec_osv_noVtx = []
    qmuidxs_osv_noVtx = []
    qmuidxs_osv_sel_noVtx = []
    qmu_dmuvec_osv_noVtx = []
    qmu_muvecdp_osv_noVtx = []
    osvvec_qmu_noVtx = []
    osvidx_qmu_noVtx = []
    qmuvec_noVtx = []
    qmuidxs_noVtx = []
    qmuidxs_sel_noVtx = []
    qmuidxsminlxy_noVtx = []
    qmuidxsmaxlxy_noVtx = []
    qmu_muvecdp_noVtx = []
    qmu_muvecdpminlxy_noVtx = []
    qmu_muvecdpmaxlxy_noVtx = []
    qmu_dmuvecminlxy_noVtx = []
    qmu_dmuvecmaxlxy_noVtx = []
    qmu_dmuvecdpminlxy_noVtx = []
    qmu_dmuvecdpmaxlxy_noVtx = []
    svvecminlxy_qmu_noVtx = []
    svidxminlxy_qmu_noVtx = []
    svvecmaxlxy_qmu_noVtx = []
    svidxmaxlxy_qmu_noVtx = []
    for m in muselidxs_noVtx:
        chg   = t.Muon_ch[m]
        ovidx = -1
        vidx  = -1
        ovpos = -1
        vpos = -1
        # First, identify muons from overlapping SVs
        if t.Muon_bestAssocSVOverlapIdx[m]>-1 and doOverlappingSV:
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
        for mm in muselidxs_noVtx:
            if mm==m:
                continue
            if abs(chg+t.Muon_ch[mm])>0:
                continue
            # First, identify muon pairs from overlapping SVs
            if ovidx>-1 and t.Muon_bestAssocSVOverlapIdx[mm]==ovidx and ovpos>-1:
                if not (m in dmuidxs_osv_noVtx or mm in dmuidxs_osv_noVtx):
                    dmuvec_osv_noVtx.append(t.Muon_vec[m])
                    dmuvec_osv_noVtx[len(dmuvec_osv_noVtx)-1] = dmuvec_osv_noVtx[len(dmuvec_osv_noVtx)-1] + t.Muon_vec[mm]
                    dmuidxs_osv_noVtx.append(m)
                    dmuidxs_osv_noVtx.append(mm)
                    dmu_muvecdp_osv_noVtx.append(ROOT.TLorentzVector())
                    dmu_muvecdp_osv_noVtx[-1].SetPtEtaPhiM(t.Muon_pt[m],t.Muon_eta[m],t.Muon_phi[m],MUON_MASS)
                    dmu_muvecdp_osv_noVtx.append(ROOT.TLorentzVector())
                    dmu_muvecdp_osv_noVtx[-1].SetPtEtaPhiM(t.Muon_pt[mm],t.Muon_eta[mm],t.Muon_phi[mm],MUON_MASS)
                    osvvec_noVtx.append(ROOT.TVector3())
                    osvvec_noVtx[len(osvvec_noVtx)-1].SetXYZ(t.SVOverlap_x[ovpos]-t.PV_x, t.SVOverlap_y[ovpos]-t.PV_y, t.SVOverlap_z[ovpos]-t.PV_z)
                    osvidx_noVtx.append(ovpos)
            # Then, identify muon pairs from non-overlapping SVs
            elif vidx>-1 and t.Muon_bestAssocSVIdx[mm]==vidx and vpos>-1:
                if not (m in dmuidxs_noVtx or mm in dmuidxs_noVtx):
                    dmuvec_noVtx.append(t.Muon_vec[m])
                    dmuvec_noVtx[len(dmuvec_noVtx)-1] = dmuvec_noVtx[len(dmuvec_noVtx)-1] + t.Muon_vec[mm]
                    dmuidxs_noVtx.append(m)
                    dmuidxs_noVtx.append(mm)
                    dmu_muvecdp_noVtx.append(ROOT.TLorentzVector())
                    dmu_muvecdp_noVtx[-1].SetPtEtaPhiM(t.Muon_pt[m],t.Muon_eta[m],t.Muon_phi[m],MUON_MASS)
                    dmu_muvecdp_noVtx.append(ROOT.TLorentzVector())
                    dmu_muvecdp_noVtx[-1].SetPtEtaPhiM(t.Muon_pt[mm],t.Muon_eta[mm],t.Muon_phi[mm],MUON_MASS)
                    svvec_noVtx.append(ROOT.TVector3())
                    svvec_noVtx[len(svvec_noVtx)-1].SetXYZ(t.SV_x[vpos]-t.PV_x, t.SV_y[vpos]-t.PV_y, t.SV_z[vpos]-t.PV_z)
                    svidx_noVtx.append(vpos)


    # -----------------------------------------------------------------------
    # Vtx SV loop
    # -----------------------------------------------------------------------
    nSVsel_vtx = 0
    _vtx_sv_selected_idxs = []
    if use_vtx:
        nSV_vtx = len(t.SV_vtx_selected)
        for v in range(nSV_vtx):
            if args.noPreSel:
                break
            if not t.SV_vtx_selected[v]:
                continue
            if applyMaterialVeto and (t.SV_vtx_onModuleWithinUnc[v] or (abs(t.SV_vtx_minDistanceFromDet_x[v]) < 0.81 and abs(t.SV_vtx_minDistanceFromDet_y[v]) < 3.24 and abs(t.SV_vtx_minDistanceFromDet_z[v]) < 0.0145)):
                continue
            nSVsel_vtx += 1
            _vtx_sv_selected_idxs.append(v)
    nSVs_vtx = nSVsel_vtx

    # -----------------------------------------------------------------------
    # Vtx muon selection
    # -----------------------------------------------------------------------
    nMu_vtx = 0
    nMuSel_vtx = 0
    nMuAss_vtx = 0
    muselidxs_vtx = []
    if use_vtx:
        nMu_vtx = len(t.Muon_vtx_selected)
        for m in range(nMu_vtx):
            if args.noMuonChi2Sel:
                if not (t.Muon_vtx_pt[m] > 3.0 and abs(t.Muon_vtx_eta[m]) < 2.4):
                    continue
            elif not t.Muon_vtx_selected[m]:
                continue
            nMuSel_vtx += 1
            if t.Muon_vtx_bestAssocSVOverlapVtxIdx[m] > -1:
                if applyMaterialVeto:
                    _ovsv = t.SVOverlap_vtx_vtxIdxs[t.Muon_vtx_bestAssocSVOverlapVtxIdx[m]][0]
                    if t.SV_vtx_onModuleWithinUnc[_ovsv] or (abs(t.SV_vtx_minDistanceFromDet_x[_ovsv]) < 0.81 and abs(t.SV_vtx_minDistanceFromDet_y[_ovsv]) < 3.24 and abs(t.SV_vtx_minDistanceFromDet_z[_ovsv]) < 0.0145):
                        continue
                nMuAss_vtx += 1
            elif t.Muon_vtx_bestAssocSVVtxIdx[m] > -1:
                if applyMaterialVeto:
                    _vpos = t.Muon_vtx_bestAssocSVVtxIdx[m]
                    if t.SV_vtx_onModuleWithinUnc[_vpos] or (abs(t.SV_vtx_minDistanceFromDet_x[_vpos]) < 0.81 and abs(t.SV_vtx_minDistanceFromDet_y[_vpos]) < 3.24 and abs(t.SV_vtx_minDistanceFromDet_z[_vpos]) < 0.0145):
                        continue
                nMuAss_vtx += 1
            else:
                continue
            muselidxs_vtx.append(m)

    # -----------------------------------------------------------------------
    # Vtx muon pairing (mirror of noVtx block, including overlapping SVs)
    # -----------------------------------------------------------------------
    dmuvec_vtx = []
    dmu_muvecdp_vtx = []
    dmuidxs_vtx = []
    svvec_vtx = []
    svidx_vtx = []
    dmuvec_osv_vtx = []
    dmu_muvecdp_osv_vtx = []
    dmuidxs_osv_vtx = []
    osvvec_vtx = []
    osvidx_vtx = []
    if use_vtx and nMuAss_vtx >= 2:
        for m in muselidxs_vtx:
            chg   = t.Muon_vtx_ch[m]
            ovidx = -1
            vidx  = -1
            ovpos = -1
            vpos  = -1
            # Overlapping SVs
            if t.Muon_vtx_bestAssocSVOverlapVtxIdx[m] > -1 and doOverlappingSV:
                if applyMaterialVeto and (t.SV_vtx_onModuleWithinUnc[t.SVOverlap_vtx_vtxIdxs[t.Muon_vtx_bestAssocSVOverlapVtxIdx[m]][0]] or (abs(t.SV_vtx_minDistanceFromDet_x[t.SVOverlap_vtx_vtxIdxs[t.Muon_vtx_bestAssocSVOverlapVtxIdx[m]][0]]) < 0.81 and abs(t.SV_vtx_minDistanceFromDet_y[t.SVOverlap_vtx_vtxIdxs[t.Muon_vtx_bestAssocSVOverlapVtxIdx[m]][0]]) < 3.24 and abs(t.SV_vtx_minDistanceFromDet_z[t.SVOverlap_vtx_vtxIdxs[t.Muon_vtx_bestAssocSVOverlapVtxIdx[m]][0]]) < 0.0145)):
                    continue
                ovidx = t.Muon_vtx_bestAssocSVOverlapVtxIdx[m]
                ovpos = ovidx
            # Non-overlapping SVs
            if t.Muon_vtx_bestAssocSVVtxIdx[m] > -1:
                vidx = t.Muon_vtx_bestAssocSVVtxIdx[m]
                vpos = vidx  # direct index for Vtx
                if applyMaterialVeto and (t.SV_vtx_onModuleWithinUnc[vpos] or (abs(t.SV_vtx_minDistanceFromDet_x[vpos]) < 0.81 and abs(t.SV_vtx_minDistanceFromDet_y[vpos]) < 3.24 and abs(t.SV_vtx_minDistanceFromDet_z[vpos]) < 0.0145)):
                    continue
            for mm in muselidxs_vtx:
                if mm == m:
                    continue
                if abs(chg + t.Muon_vtx_ch[mm]) > 0:
                    continue
                # Overlapping SV pair
                if ovidx > -1 and t.Muon_vtx_bestAssocSVOverlapVtxIdx[mm] == ovidx and ovpos > -1:
                    if not (m in dmuidxs_osv_vtx or mm in dmuidxs_osv_vtx):
                        dmuvec_osv_vtx.append(vtx_muvec_corr(t, m))
                        dmuvec_osv_vtx[-1] = dmuvec_osv_vtx[-1] + vtx_muvec_corr(t, mm)
                        dmuidxs_osv_vtx.append(m)
                        dmuidxs_osv_vtx.append(mm)
                        dmu_muvecdp_osv_vtx.append(ROOT.TLorentzVector())
                        dmu_muvecdp_osv_vtx[-1].SetPtEtaPhiM(t.Muon_vtx_pt[m], t.Muon_vtx_eta[m], t.Muon_vtx_phi[m], MUON_MASS)
                        dmu_muvecdp_osv_vtx.append(ROOT.TLorentzVector())
                        dmu_muvecdp_osv_vtx[-1].SetPtEtaPhiM(t.Muon_vtx_pt[mm], t.Muon_vtx_eta[mm], t.Muon_vtx_phi[mm], MUON_MASS)
                        osvvec_vtx.append(ROOT.TVector3())
                        osvvec_vtx[-1].SetXYZ(t.SVOverlap_vtx_x[ovpos]-t.PV_x, t.SVOverlap_vtx_y[ovpos]-t.PV_y, t.SVOverlap_vtx_z[ovpos]-t.PV_z)
                        osvidx_vtx.append(ovpos)
                # Non-overlapping SV pair
                elif vidx > -1 and t.Muon_vtx_bestAssocSVVtxIdx[mm] == vidx and vpos > -1:
                    if not (m in dmuidxs_vtx or mm in dmuidxs_vtx):
                        dmuvec_vtx.append(vtx_muvec_corr(t, m))
                        dmuvec_vtx[-1] = dmuvec_vtx[-1] + vtx_muvec_corr(t, mm)
                        dmuidxs_vtx.append(m)
                        dmuidxs_vtx.append(mm)
                        dmu_muvecdp_vtx.append(ROOT.TLorentzVector())
                        dmu_muvecdp_vtx[-1].SetPtEtaPhiM(t.Muon_vtx_pt[m], t.Muon_vtx_eta[m], t.Muon_vtx_phi[m], MUON_MASS)
                        dmu_muvecdp_vtx.append(ROOT.TLorentzVector())
                        dmu_muvecdp_vtx[-1].SetPtEtaPhiM(t.Muon_vtx_pt[mm], t.Muon_vtx_eta[mm], t.Muon_vtx_phi[mm], MUON_MASS)
                        svvec_vtx.append(ROOT.TVector3())
                        svvec_vtx[-1].SetXYZ(t.SV_vtx_x[vpos]-t.PV_x, t.SV_vtx_y[vpos]-t.PV_y, t.SV_vtx_z[vpos]-t.PV_z)
                        svidx_vtx.append(vpos)

    # -----------------------------------------------------------------------
    # Combined event guard: skip if no usable dimuon candidate from any active collection
    # -----------------------------------------------------------------------
    n_dimu_novtx = len(dmuvec_noVtx) + len(dmuvec_osv_noVtx)
    n_dimu_vtx   = len(dmuvec_vtx)   + len(dmuvec_osv_vtx)
    if use_novtx and not use_vtx and n_dimu_novtx == 0:
        continue
    if use_vtx and not use_novtx and n_dimu_vtx == 0:
        continue
    if use_novtx and use_vtx and n_dimu_novtx == 0 and n_dimu_vtx == 0:
        continue

    # -----------------------------------------------------------------------
    # OR SV deduplication: build Vtx SV positions so NoVtx SVs that are
    # geometrically matched to a Vtx SV can be dropped in OR mode.
    # -----------------------------------------------------------------------
    _or_vtx_sv_vecs = []
    _or_vtx_sv_lxy  = []
    if use_OR:
        for vi in _vtx_sv_selected_idxs:
            _or_vtx_sv_vecs.append(ROOT.TVector3(
                t.SV_vtx_x[vi] - t.PV_x,
                t.SV_vtx_y[vi] - t.PV_y,
                t.SV_vtx_z[vi] - t.PV_z
            ))
            _or_vtx_sv_lxy.append(t.SV_vtx_lxy[vi])

        def _matched_sv_to_vtx(sv_vec, sv_lxy):
            for vsv, vsv_lxy in zip(_or_vtx_sv_vecs, _or_vtx_sv_lxy):
                deta = sv_vec.PseudoRapidity() - vsv.PseudoRapidity()
                dphi_val = (sv_vec.Phi() - vsv.Phi() + math.pi) % (2 * math.pi) - math.pi
                if math.sqrt(deta**2 + dphi_val**2) < 0.1 and abs(sv_lxy - vsv_lxy) < OR_SV_DLXY_MAX:
                    return True
            return False

    # -----------------------------------------------------------------------
    # Build combined collection (dmuvec_all, etc.) based on --collection
    # Vtx candidates are always added first so OR deduplication can compare
    # NoVtx muons against the already-selected Vtx muon set.
    # -----------------------------------------------------------------------
    dmuvec_all      = []
    dmu_muvecdp_all = []
    dmuidxs_all     = []
    svvec_all       = []
    svidx_all       = []
    dimu_type_all   = []  # 'noVtx', 'noVtx_osv', 'vtx', 'vtx_osv'

    if use_vtx:
        for i in range(len(dmuvec_vtx)):
            dmuvec_all.append(dmuvec_vtx[i])
            dmu_muvecdp_all.append(dmu_muvecdp_vtx[2*i])
            dmu_muvecdp_all.append(dmu_muvecdp_vtx[2*i+1])
            dmuidxs_all.append(dmuidxs_vtx[2*i])
            dmuidxs_all.append(dmuidxs_vtx[2*i+1])
            svvec_all.append(svvec_vtx[i])
            svidx_all.append(svidx_vtx[i])
            dimu_type_all.append('vtx')
        for i in range(len(dmuvec_osv_vtx)):
            dmuvec_all.append(dmuvec_osv_vtx[i])
            dmu_muvecdp_all.append(dmu_muvecdp_osv_vtx[2*i])
            dmu_muvecdp_all.append(dmu_muvecdp_osv_vtx[2*i+1])
            dmuidxs_all.append(dmuidxs_osv_vtx[2*i])
            dmuidxs_all.append(dmuidxs_osv_vtx[2*i+1])
            svvec_all.append(osvvec_vtx[i])
            svidx_all.append(osvidx_vtx[i])
            dimu_type_all.append('vtx_osv')

    if use_novtx:
        # In OR mode: build (eta,phi) for every Vtx muon used in any Vtx dimuon,
        # then skip NoVtx candidates whose muons dR-match a Vtx muon (dR < 0.1).
        if use_OR:
            vtx_mu_coords = [(t.Muon_vtx_eta[mi], t.Muon_vtx_phi[mi])
                             for mi in dmuidxs_vtx + dmuidxs_osv_vtx]
            def _matched_to_vtx(eta, phi):
                for veta, vphi in vtx_mu_coords:
                    dphi_val = (phi - vphi + math.pi) % (2*math.pi) - math.pi
                    if math.sqrt((eta - veta)**2 + dphi_val**2) < 0.1:
                        return True
                return False

        for i in range(len(dmuvec_noVtx)):
            mi  = dmuidxs_noVtx[2*i]
            mii = dmuidxs_noVtx[2*i+1]
            if use_OR and (_matched_to_vtx(t.Muon_eta[mi],  t.Muon_phi[mi]) or
                           _matched_to_vtx(t.Muon_eta[mii], t.Muon_phi[mii])):
                continue
            if use_OR and _matched_sv_to_vtx(svvec_noVtx[i], t.SV_lxy[svidx_noVtx[i]]):
                continue
            dmuvec_all.append(dmuvec_noVtx[i])
            dmu_muvecdp_all.append(dmu_muvecdp_noVtx[2*i])
            dmu_muvecdp_all.append(dmu_muvecdp_noVtx[2*i+1])
            dmuidxs_all.append(dmuidxs_noVtx[2*i])
            dmuidxs_all.append(dmuidxs_noVtx[2*i+1])
            svvec_all.append(svvec_noVtx[i])
            svidx_all.append(svidx_noVtx[i])
            dimu_type_all.append('noVtx')
        for i in range(len(dmuvec_osv_noVtx)):
            mi  = dmuidxs_osv_noVtx[2*i]
            mii = dmuidxs_osv_noVtx[2*i+1]
            if use_OR and (_matched_to_vtx(t.Muon_eta[mi],  t.Muon_phi[mi]) or
                           _matched_to_vtx(t.Muon_eta[mii], t.Muon_phi[mii])):
                continue
            if use_OR and _matched_sv_to_vtx(osvvec_noVtx[i], t.SVOverlap_lxy[osvidx_noVtx[i]]):
                continue
            dmuvec_all.append(dmuvec_osv_noVtx[i])
            dmu_muvecdp_all.append(dmu_muvecdp_osv_noVtx[2*i])
            dmu_muvecdp_all.append(dmu_muvecdp_osv_noVtx[2*i+1])
            dmuidxs_all.append(dmuidxs_osv_noVtx[2*i])
            dmuidxs_all.append(dmuidxs_osv_noVtx[2*i+1])
            svvec_all.append(osvvec_noVtx[i])
            svidx_all.append(osvidx_noVtx[i])
            dimu_type_all.append('noVtx_osv')

    ### Scan analysis initialization
    ## If you put one of the to True, it won't fill.
    # Cat selection:
    filledcat4musep = False # False
    filledcat4muosv = False # False
    filledcat2mu = False # False

    # Apply selections and fill histograms for all dimuon candidates (noVtx, noVtx_osv, Vtx)
    for vn,v in enumerate(dmuvec_all):
        if args.noDiMuon and dimu_type_all[vn] in ('noVtx_osv', 'vtx_osv'):
            continue
        if not applyDiMuonSelection(v):
            continue
        # --- collection-aware dispatch ---
        dimu_type  = dimu_type_all[vn]
        is_vtx     = dimu_type in ('vtx', 'vtx_osv')
        mu1idx     = dmuidxs_all[int(vn*2)]
        mu2idx     = dmuidxs_all[int(vn*2)+1]
        _svvec     = svvec_all[vn]
        _muvecdp1  = dmu_muvecdp_all[int(vn*2)]
        _muvecdp2  = dmu_muvecdp_all[int(vn*2)+1]
        if dimu_type == 'vtx':
            lxy    = t.SV_vtx_lxy[svidx_all[vn]]
            lz     = abs(t.SV_vtx_z[svidx_all[vn]])
            _sv_x  = t.SV_vtx_x[svidx_all[vn]]
            _sv_y  = t.SV_vtx_y[svidx_all[vn]]
            _sv_z  = t.SV_vtx_z[svidx_all[vn]]
        elif dimu_type == 'vtx_osv':
            lxy    = t.SVOverlap_vtx_lxy[svidx_all[vn]]
            lz     = abs(t.SVOverlap_vtx_z[svidx_all[vn]])
            _sv_x  = t.SVOverlap_vtx_x[svidx_all[vn]]
            _sv_y  = t.SVOverlap_vtx_y[svidx_all[vn]]
            _sv_z  = t.SVOverlap_vtx_z[svidx_all[vn]]
        elif dimu_type == 'noVtx_osv':
            lxy    = t.SVOverlap_lxy[svidx_all[vn]]
            lz     = abs(t.SVOverlap_z[svidx_all[vn]])
            _sv_x  = t.SVOverlap_x[svidx_all[vn]]
            _sv_y  = t.SVOverlap_y[svidx_all[vn]]
            _sv_z  = t.SVOverlap_z[svidx_all[vn]]
        else:
            lxy    = t.SV_lxy[svidx_all[vn]]
            lz     = abs(t.SV_z[svidx_all[vn]])
            _sv_x  = t.SV_x[svidx_all[vn]]
            _sv_y  = t.SV_y[svidx_all[vn]]
            _sv_z  = t.SV_z[svidx_all[vn]]
        _mu1_pt    = muattr(t, 'pt',           mu1idx, is_vtx)
        _mu2_pt    = muattr(t, 'pt',           mu2idx, is_vtx)
        _mu1_dxyCr = muattr(t, 'dxyCorr',     mu1idx, is_vtx)
        _mu2_dxyCr = muattr(t, 'dxyCorr',     mu2idx, is_vtx)
        _mu1_dxye  = muattr(t, 'dxye',        mu1idx, is_vtx)
        _mu2_dxye  = muattr(t, 'dxye',        mu2idx, is_vtx)
        _mu1_vec   = muattr(t, 'vec',          mu1idx, is_vtx)
        _mu2_vec   = muattr(t, 'vec',          mu2idx, is_vtx)
        _mu1_iso   = muattr(t, 'PFIsoAll0p4',  mu1idx, is_vtx)
        _mu2_iso   = muattr(t, 'PFIsoAll0p4',  mu2idx, is_vtx)
        _nhits1    = mu_nhits(t, mu1idx, is_vtx)
        _nhits2    = mu_nhits(t, mu2idx, is_vtx)
        # ---
        if not applyLxySelection(lxy):
            continue
        if not applyLzSelection(lz):
            continue
        mass = v.M()
        pt   = v.Pt()
        nhitsbeforesvtotal = _nhits1 + _nhits2
        subpt = _mu1_pt if _mu1_pt < _mu2_pt else _mu2_pt
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
        #
        #  Apply selection on muon lifetime-scaled dxy
        if applyMuonIPSel:
            if ( abs(_mu1_dxyCr)/(lxy*mass/pt)<0.1 or
                 abs(_mu2_dxyCr)/(lxy*mass/pt)<0.1 ):
                continue
            if ( abs(_mu1_dxyCr/_mu1_dxye)<2.0 or abs(_mu2_dxyCr/_mu2_dxye)<2.0 ):
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
        if applyDiMuonResonanceMasking:
            if (mass > 0.41 and mass < 0.50):
                continue
            if (mass > 0.51 and mass < 0.59):
                continue
            if (mass > 0.731 and mass < 0.83):
                continue
            if (mass > 0.96 and mass < 1.08):
                continue
            if (mass > 2.91 and mass < 3.27):
                continue
            if (mass > 3.47 and mass < 3.89):
                continue
            if (mass > 8.99 and mass < 9.91):
                continue
            if (mass > 9.64 and mass < 10.56):
                continue
            if (mass > 9.90 and mass < 10.78):
                continue
            if (mass > 81 and mass < 101):
                continue
        #

        drmm = _mu1_vec.DeltaR(_mu2_vec)
        dpmm = abs(_mu1_vec.DeltaPhi(_mu2_vec))
        demm = abs(_mu1_vec.Eta()-_mu2_vec.Eta())
        dedpmm = 1e6
        if dpmm>0.0:
            dedpmm = demm/dpmm
        a3dmm = abs(_mu1_vec.Angle(_mu2_vec.Vect()))
        #
        drmmu = _muvecdp1.DeltaR(_muvecdp2)
        dpmmu = abs(_muvecdp1.DeltaPhi(_muvecdp2))
        demmu = abs(_muvecdp1.Eta()-_muvecdp2.Eta())
        dedpmmu = 1e6
        if dpmmu>0.0:
            dedpmmu = demmu/dpmmu
        a3dmmu = abs(_muvecdp1.Angle(_muvecdp2.Vect()))
        #
        dphisv = abs(v.Vect().DeltaPhi(_svvec))
        dphisv1 = abs(_mu1_vec.Vect().DeltaPhi(_svvec))
        dphisv2 = abs(_mu2_vec.Vect().DeltaPhi(_svvec))
        detasv = abs(v.Vect().Eta()-_svvec.Eta())
        detadphisv = 1e6
        if dphisv>0.0:
            detadphisv = detasv/dphisv
        a3dsv  = abs(v.Vect().Angle(_svvec))
        #
        desvdpsvmin = 1e6
        demmdpsvmin = 1e6
        dpsvmin = 1e6
        tdpsv = abs(_mu1_vec.Vect().DeltaPhi(_svvec))
        if tdpsv < dpsvmin:
            dpsvmin = tdpsv
        tdpsv = abs(_mu2_vec.Vect().DeltaPhi(_svvec))
        if tdpsv < dpsvmin:
            dpsvmin = tdpsv
        if dpsvmin>0.0:
            demmdpsvmin = demm/dpsvmin
            desvdpsvmin = detasv/dpsvmin
        #
        vu = (_muvecdp1+_muvecdp2)
        dphisvu = abs(vu.Vect().DeltaPhi(_svvec))
        dphisv1u = abs(_muvecdp1.Vect().DeltaPhi(_svvec))
        dphisv2u = abs(_muvecdp2.Vect().DeltaPhi(_svvec))
        detasvu = abs(vu.Vect().Eta()-_svvec.Eta())
        detadphisvu = 1e6
        if dphisvu>0.0:
            detadphisvu = detasvu/dphisvu
        a3dsvu  = abs(vu.Vect().Angle(_svvec))
        #
        desvdpsvminu = 1e6
        demmdpsvminu = 1e6
        dpsvminu = 1e6
        tdpsvu = abs(_muvecdp1.Vect().DeltaPhi(_svvec))
        if tdpsvu < dpsvminu:
            dpsvminu = tdpsvu
        tdpsvu = abs(_muvecdp2.Vect().DeltaPhi(_svvec))
        if tdpsvu < dpsvminu:
            dpsvminu = tdpsvu
        if dpsvminu>0.0:
            demmdpsvminu = demmu/dpsvminu
            desvdpsvminu = detasvu/dpsvminu
        #
        sindpsvlxy = abs(ROOT.TMath.Sin(v.Vect().DeltaPhi(_svvec)))*lxy
        sindpsvlxyu = abs(ROOT.TMath.Sin(vu.Vect().DeltaPhi(_svvec)))*lxy
        sina3dsvl3d = abs(ROOT.TMath.Sin(v.Vect().Angle(_svvec)))*(_svvec.Mag())
        sina3dsvl3du = abs(ROOT.TMath.Sin(vu.Vect().Angle(_svvec)))*(_svvec.Mag())
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
        isocat = dimuonIsoCategory(_mu1_iso, _mu1_pt, _mu2_iso, _mu2_pt)
        if float(args.dimuonIsoCatSel) > -1 and isocat!=float(args.dimuonIsoCatSel):
            continue

        #
        # Branch filling
        if branches["SV1_x"][-1] < 0.:
            _sv_slot = "SV1"
        elif branches["SV2_x"][-1] < 0.:
            _sv_slot = "SV2"
        else:
            break
        _mu1_dxyPV, _mu1_dzPV = muon_ip_pv(t, mu1idx, is_vtx)
        _mu2_dxyPV, _mu2_dzPV = muon_ip_pv(t, mu2idx, is_vtx)
        _mu1_eta_ntup  = muattr(t, 'eta',      mu1idx, is_vtx)
        _mu2_eta_ntup  = muattr(t, 'eta',      mu2idx, is_vtx)
        _mu1_phi_ntup  = muattr(t, 'phi',      mu1idx, is_vtx)
        _mu2_phi_ntup  = muattr(t, 'phi',      mu2idx, is_vtx)
        _mu1_nchi2     = muattr(t, 'chi2Ndof', mu1idx, is_vtx)
        _mu2_nchi2     = muattr(t, 'chi2Ndof', mu2idx, is_vtx)
        _sv_direct     = dimu_type in ('vtx', 'noVtx')
        _svi           = svidx_all[vn]
        _sv_chi2   = svattr(t, 'chi2',  _svi, is_vtx) if _sv_direct else -1.
        _sv_prob   = svattr(t, 'prob',  _svi, is_vtx) if _sv_direct else -1.
        _sv_xErr   = svattr(t, 'xe', _svi, is_vtx) if _sv_direct else -1.
        _sv_yErr   = svattr(t, 'ye', _svi, is_vtx) if _sv_direct else -1.
        _sv_zErr   = svattr(t, 'ze', _svi, is_vtx) if _sv_direct else -1.
        branches[f"{_sv_slot}_x"][-1]            = _sv_x
        branches[f"{_sv_slot}_y"][-1]            = _sv_y
        branches[f"{_sv_slot}_z"][-1]            = _sv_z
        branches[f"{_sv_slot}_lxy"][-1]          = lxy
        branches[f"{_sv_slot}_ptmm"][-1]         = pt
        branches[f"{_sv_slot}_chi2"][-1]         = _sv_chi2
        branches[f"{_sv_slot}_prob"][-1]         = _sv_prob
        branches[f"{_sv_slot}_xErr"][-1]         = _sv_xErr
        branches[f"{_sv_slot}_yErr"][-1]         = _sv_yErr
        branches[f"{_sv_slot}_zErr"][-1]         = _sv_zErr
        branches[f"{_sv_slot}_dphi_mumu_SV"][-1]  = dphisvu
        branches[f"{_sv_slot}_d3d_mumu_SV"][-1]  = a3dsvu
        branches[f"{_sv_slot}_a3d_mumu"][-1]      = a3dmmu
        branches[f"{_sv_slot}_dr_mumu"][-1]       = drmmu
        branches[f"{_sv_slot}_dphi_mumu"][-1]     = dpmmu
        branches[f"{_sv_slot}_deta_mumu"][-1]     = demmu
        branches[f"{_sv_slot}_dphisv_mu1"][-1]    = dphisv1u
        branches[f"{_sv_slot}_dphisv_mu2"][-1]    = dphisv2u
        branches[f"{_sv_slot}_deta_mumu_SV"][-1]  = detasvu
        branches[f"{_sv_slot}_sina3d_l3d"][-1]    = sina3dsvl3du
        branches[f"{_sv_slot}_sindphi_lxy"][-1]   = sindpsvlxyu
        branches[f"{_sv_slot}_mu1_pt"][-1]                         = _mu1_pt
        branches[f"{_sv_slot}_mu1_eta"][-1]                        = _mu1_eta_ntup
        branches[f"{_sv_slot}_mu1_phi"][-1]                        = _mu1_phi_ntup
        branches[f"{_sv_slot}_mu1_isvtx"][-1]                      = float(is_vtx)
        branches[f"{_sv_slot}_mu1_normChi2"][-1]                   = _mu1_nchi2
        branches[f"{_sv_slot}_mu1_dxy"][-1]                        = _mu1_dxyCr
        branches[f"{_sv_slot}_mu1_dxyErr"][-1]                     = _mu1_dxye
        branches[f"{_sv_slot}_mu1_dxysig"][-1]                     = _mu1_dxyPV/muattr(t, 'dxye', mu1idx, is_vtx)
        branches[f"{_sv_slot}_mu1_dz"][-1]                         = _mu1_dzPV
        branches[f"{_sv_slot}_mu1_dze"][-1]                        = muattr(t, 'dze',                     mu1idx, is_vtx)
        branches[f"{_sv_slot}_mu1_dzsig"][-1]                      = _mu1_dzPV/muattr(t, 'dze', mu1idx, is_vtx)
        branches[f"{_sv_slot}_mu1_nhitsbeforesv"][-1]              = mu_nhits(t, mu1idx, is_vtx)
        branches[f"{_sv_slot}_mu1_isGlobal"][-1]                   = muattr(t, 'isGlobal',                mu1idx, is_vtx)
        branches[f"{_sv_slot}_mu1_isTracker"][-1]                  = muattr(t, 'isTracker',               mu1idx, is_vtx)
        branches[f"{_sv_slot}_mu1_isStandAlone"][-1]               = muattr(t, 'isStandAlone',            mu1idx, is_vtx)
        branches[f"{_sv_slot}_mu1_pixHits"][-1]                    = muattr(t, 'pixHits',                 mu1idx, is_vtx)
        branches[f"{_sv_slot}_mu1_stripHits"][-1]                  = muattr(t, 'stripHits',               mu1idx, is_vtx)
        branches[f"{_sv_slot}_mu1_pixLayers"][-1]                  = muattr(t, 'pixLayers',               mu1idx, is_vtx)
        branches[f"{_sv_slot}_mu1_trkLayers"][-1]                  = muattr(t, 'trkLayers',               mu1idx, is_vtx)
        branches[f"{_sv_slot}_mu1_saHits"][-1]                     = muattr(t, 'saHits',                  mu1idx, is_vtx)
        branches[f"{_sv_slot}_mu1_saMatchedStats"][-1]             = muattr(t, 'saMatchedStats',          mu1idx, is_vtx)
        branches[f"{_sv_slot}_mu1_muHits"][-1]                     = muattr(t, 'muHits',                  mu1idx, is_vtx)
        branches[f"{_sv_slot}_mu1_muChambs"][-1]                   = muattr(t, 'muChambs',                mu1idx, is_vtx)
        branches[f"{_sv_slot}_mu1_muCSCDT"][-1]                    = muattr(t, 'muCSCDT',                 mu1idx, is_vtx)
        branches[f"{_sv_slot}_mu1_muMatch"][-1]                    = muattr(t, 'muMatch',                 mu1idx, is_vtx)
        branches[f"{_sv_slot}_mu1_muMatchedStats"][-1]             = muattr(t, 'muMatchedStats',          mu1idx, is_vtx)
        branches[f"{_sv_slot}_mu1_muExpMatchedStats"][-1]          = muattr(t, 'muExpMatchedStats',       mu1idx, is_vtx)
        branches[f"{_sv_slot}_mu1_muMatchedRPC"][-1]               = muattr(t, 'muMatchedRPC',            mu1idx, is_vtx)
        branches[f"{_sv_slot}_mu1_ecalIso"][-1]                    = muattr(t, 'ecalIso',                 mu1idx, is_vtx)
        branches[f"{_sv_slot}_mu1_hcalIso"][-1]                    = muattr(t, 'hcalIso',                 mu1idx, is_vtx)
        branches[f"{_sv_slot}_mu1_trackIso"][-1]                   = muattr(t, 'trackIso',                mu1idx, is_vtx)
        branches[f"{_sv_slot}_mu1_ecalRelIso"][-1]                 = muattr(t, 'ecalRelIso',              mu1idx, is_vtx)
        branches[f"{_sv_slot}_mu1_hcalRelIso"][-1]                 = muattr(t, 'hcalRelIso',              mu1idx, is_vtx)
        branches[f"{_sv_slot}_mu1_trackRelIso"][-1]                = muattr(t, 'trackRelIso',             mu1idx, is_vtx)
        branches[f"{_sv_slot}_mu1_PFIsoChg0p3"][-1]               = muattr(t, 'PFIsoChg0p3',             mu1idx, is_vtx)
        branches[f"{_sv_slot}_mu1_PFIsoAll0p3"][-1]               = muattr(t, 'PFIsoAll0p3',             mu1idx, is_vtx)
        branches[f"{_sv_slot}_mu1_PFRelIsoChg0p3"][-1]            = muattr(t, 'PFRelIsoChg0p3',          mu1idx, is_vtx)
        branches[f"{_sv_slot}_mu1_PFRelIsoAll0p3"][-1]            = muattr(t, 'PFRelIsoAll0p3',          mu1idx, is_vtx)
        branches[f"{_sv_slot}_mu1_mindrPF0p3"][-1]                = muattr(t, 'mindrPF0p3',              mu1idx, is_vtx)
        branches[f"{_sv_slot}_mu1_PFIsoChg0p4"][-1]               = muattr(t, 'PFIsoChg0p4',             mu1idx, is_vtx)
        branches[f"{_sv_slot}_mu1_PFIsoAll0p4"][-1]               = muattr(t, 'PFIsoAll0p4',             mu1idx, is_vtx)
        branches[f"{_sv_slot}_mu1_PFRelIsoChg0p4"][-1]            = muattr(t, 'PFRelIsoChg0p4',          mu1idx, is_vtx)
        branches[f"{_sv_slot}_mu1_PFRelIsoAll0p4"][-1]            = muattr(t, 'PFRelIsoAll0p4',          mu1idx, is_vtx)
        branches[f"{_sv_slot}_mu1_mindrPF0p4"][-1]                = muattr(t, 'mindrPF0p4',              mu1idx, is_vtx)
        branches[f"{_sv_slot}_mu1_mindr"][-1]                      = muattr(t, 'mindr',                   mu1idx, is_vtx)
        branches[f"{_sv_slot}_mu1_maxdr"][-1]                      = muattr(t, 'maxdr',                   mu1idx, is_vtx)
        branches[f"{_sv_slot}_mu1_mindrJet"][-1]                   = muattr(t, 'mindrJet',                mu1idx, is_vtx)
        branches[f"{_sv_slot}_mu1_mindphiJet"][-1]                 = muattr(t, 'mindphiJet',              mu1idx, is_vtx)
        branches[f"{_sv_slot}_mu1_mindetaJet"][-1]                 = muattr(t, 'mindetaJet',              mu1idx, is_vtx)
        branches[f"{_sv_slot}_mu1_ncompatible"][-1]                = muattr(t, 'ncompatible',             mu1idx, is_vtx)
        branches[f"{_sv_slot}_mu1_ncompatibletotal"][-1]           = muattr(t, 'ncompatibletotal',        mu1idx, is_vtx)
        branches[f"{_sv_slot}_mu1_nexpectedhits"][-1]              = muattr(t, 'nexpectedhits',           mu1idx, is_vtx)
        branches[f"{_sv_slot}_mu1_nexpectedhitsmultiple"][-1]      = muattr(t, 'nexpectedhitsmultiple',   mu1idx, is_vtx)
        branches[f"{_sv_slot}_mu1_nexpectedhitsmultipletotal"][-1] = muattr(t, 'nexpectedhitsmultipletotal', mu1idx, is_vtx)
        branches[f"{_sv_slot}_mu1_nexpectedhitstotal"][-1]         = muattr(t, 'nexpectedhitstotal',      mu1idx, is_vtx)
        branches[f"{_sv_slot}_mu1_phiCorr"][-1]                    = muattr(t, 'phiCorr',                 mu1idx, is_vtx)
        branches[f"{_sv_slot}_mu2_pt"][-1]                         = _mu2_pt
        branches[f"{_sv_slot}_mu2_eta"][-1]                        = _mu2_eta_ntup
        branches[f"{_sv_slot}_mu2_phi"][-1]                        = _mu2_phi_ntup
        branches[f"{_sv_slot}_mu2_isvtx"][-1]                      = float(is_vtx)
        branches[f"{_sv_slot}_mu2_normChi2"][-1]                   = _mu2_nchi2
        branches[f"{_sv_slot}_mu2_dxy"][-1]                        = _mu2_dxyCr
        branches[f"{_sv_slot}_mu2_dxyErr"][-1]                     = _mu2_dxye
        branches[f"{_sv_slot}_mu2_dxysig"][-1]                     = _mu2_dxyPV/muattr(t, 'dxye', mu2idx, is_vtx)
        branches[f"{_sv_slot}_mu2_dz"][-1]                         = _mu2_dzPV
        branches[f"{_sv_slot}_mu2_dze"][-1]                        = muattr(t, 'dze',                     mu2idx, is_vtx)
        branches[f"{_sv_slot}_mu2_dzsig"][-1]                      = _mu2_dzPV/muattr(t, 'dze', mu2idx, is_vtx)
        branches[f"{_sv_slot}_mu2_nhitsbeforesv"][-1]              = mu_nhits(t, mu2idx, is_vtx)
        branches[f"{_sv_slot}_mu2_isGlobal"][-1]                   = muattr(t, 'isGlobal',                mu2idx, is_vtx)
        branches[f"{_sv_slot}_mu2_isTracker"][-1]                  = muattr(t, 'isTracker',               mu2idx, is_vtx)
        branches[f"{_sv_slot}_mu2_isStandAlone"][-1]               = muattr(t, 'isStandAlone',            mu2idx, is_vtx)
        branches[f"{_sv_slot}_mu2_pixHits"][-1]                    = muattr(t, 'pixHits',                 mu2idx, is_vtx)
        branches[f"{_sv_slot}_mu2_stripHits"][-1]                  = muattr(t, 'stripHits',               mu2idx, is_vtx)
        branches[f"{_sv_slot}_mu2_pixLayers"][-1]                  = muattr(t, 'pixLayers',               mu2idx, is_vtx)
        branches[f"{_sv_slot}_mu2_trkLayers"][-1]                  = muattr(t, 'trkLayers',               mu2idx, is_vtx)
        branches[f"{_sv_slot}_mu2_saHits"][-1]                     = muattr(t, 'saHits',                  mu2idx, is_vtx)
        branches[f"{_sv_slot}_mu2_saMatchedStats"][-1]             = muattr(t, 'saMatchedStats',          mu2idx, is_vtx)
        branches[f"{_sv_slot}_mu2_muHits"][-1]                     = muattr(t, 'muHits',                  mu2idx, is_vtx)
        branches[f"{_sv_slot}_mu2_muChambs"][-1]                   = muattr(t, 'muChambs',                mu2idx, is_vtx)
        branches[f"{_sv_slot}_mu2_muCSCDT"][-1]                    = muattr(t, 'muCSCDT',                 mu2idx, is_vtx)
        branches[f"{_sv_slot}_mu2_muMatch"][-1]                    = muattr(t, 'muMatch',                 mu2idx, is_vtx)
        branches[f"{_sv_slot}_mu2_muMatchedStats"][-1]             = muattr(t, 'muMatchedStats',          mu2idx, is_vtx)
        branches[f"{_sv_slot}_mu2_muExpMatchedStats"][-1]          = muattr(t, 'muExpMatchedStats',       mu2idx, is_vtx)
        branches[f"{_sv_slot}_mu2_muMatchedRPC"][-1]               = muattr(t, 'muMatchedRPC',            mu2idx, is_vtx)
        branches[f"{_sv_slot}_mu2_ecalIso"][-1]                    = muattr(t, 'ecalIso',                 mu2idx, is_vtx)
        branches[f"{_sv_slot}_mu2_hcalIso"][-1]                    = muattr(t, 'hcalIso',                 mu2idx, is_vtx)
        branches[f"{_sv_slot}_mu2_trackIso"][-1]                   = muattr(t, 'trackIso',                mu2idx, is_vtx)
        branches[f"{_sv_slot}_mu2_ecalRelIso"][-1]                 = muattr(t, 'ecalRelIso',              mu2idx, is_vtx)
        branches[f"{_sv_slot}_mu2_hcalRelIso"][-1]                 = muattr(t, 'hcalRelIso',              mu2idx, is_vtx)
        branches[f"{_sv_slot}_mu2_trackRelIso"][-1]                = muattr(t, 'trackRelIso',             mu2idx, is_vtx)
        branches[f"{_sv_slot}_mu2_PFIsoChg0p3"][-1]               = muattr(t, 'PFIsoChg0p3',             mu2idx, is_vtx)
        branches[f"{_sv_slot}_mu2_PFIsoAll0p3"][-1]               = muattr(t, 'PFIsoAll0p3',             mu2idx, is_vtx)
        branches[f"{_sv_slot}_mu2_PFRelIsoChg0p3"][-1]            = muattr(t, 'PFRelIsoChg0p3',          mu2idx, is_vtx)
        branches[f"{_sv_slot}_mu2_PFRelIsoAll0p3"][-1]            = muattr(t, 'PFRelIsoAll0p3',          mu2idx, is_vtx)
        branches[f"{_sv_slot}_mu2_mindrPF0p3"][-1]                = muattr(t, 'mindrPF0p3',              mu2idx, is_vtx)
        branches[f"{_sv_slot}_mu2_PFIsoChg0p4"][-1]               = muattr(t, 'PFIsoChg0p4',             mu2idx, is_vtx)
        branches[f"{_sv_slot}_mu2_PFIsoAll0p4"][-1]               = muattr(t, 'PFIsoAll0p4',             mu2idx, is_vtx)
        branches[f"{_sv_slot}_mu2_PFRelIsoChg0p4"][-1]            = muattr(t, 'PFRelIsoChg0p4',          mu2idx, is_vtx)
        branches[f"{_sv_slot}_mu2_PFRelIsoAll0p4"][-1]            = muattr(t, 'PFRelIsoAll0p4',          mu2idx, is_vtx)
        branches[f"{_sv_slot}_mu2_mindrPF0p4"][-1]                = muattr(t, 'mindrPF0p4',              mu2idx, is_vtx)
        branches[f"{_sv_slot}_mu2_mindr"][-1]                      = muattr(t, 'mindr',                   mu2idx, is_vtx)
        branches[f"{_sv_slot}_mu2_maxdr"][-1]                      = muattr(t, 'maxdr',                   mu2idx, is_vtx)
        branches[f"{_sv_slot}_mu2_mindrJet"][-1]                   = muattr(t, 'mindrJet',                mu2idx, is_vtx)
        branches[f"{_sv_slot}_mu2_mindphiJet"][-1]                 = muattr(t, 'mindphiJet',              mu2idx, is_vtx)
        branches[f"{_sv_slot}_mu2_mindetaJet"][-1]                 = muattr(t, 'mindetaJet',              mu2idx, is_vtx)
        branches[f"{_sv_slot}_mu2_ncompatible"][-1]                = muattr(t, 'ncompatible',             mu2idx, is_vtx)
        branches[f"{_sv_slot}_mu2_ncompatibletotal"][-1]           = muattr(t, 'ncompatibletotal',        mu2idx, is_vtx)
        branches[f"{_sv_slot}_mu2_nexpectedhits"][-1]              = muattr(t, 'nexpectedhits',           mu2idx, is_vtx)
        branches[f"{_sv_slot}_mu2_nexpectedhitsmultiple"][-1]      = muattr(t, 'nexpectedhitsmultiple',   mu2idx, is_vtx)
        branches[f"{_sv_slot}_mu2_nexpectedhitsmultipletotal"][-1] = muattr(t, 'nexpectedhitsmultipletotal', mu2idx, is_vtx)
        branches[f"{_sv_slot}_mu2_nexpectedhitstotal"][-1]         = muattr(t, 'nexpectedhitstotal',      mu2idx, is_vtx)
        branches[f"{_sv_slot}_mu2_phiCorr"][-1]                    = muattr(t, 'phiCorr',                 mu2idx, is_vtx)
        branches[f"{_sv_slot}_mass"][-1]         = mass
        branches[f"{_sv_slot}_ndof"][-1]                    = svattr(t, 'ndof',                  _svi, is_vtx) if _sv_direct else -1.
        branches[f"{_sv_slot}_chi2Ndof"][-1]               = svattr(t, 'chi2Ndof',              _svi, is_vtx) if _sv_direct else -1.
        branches[f"{_sv_slot}_l3d"][-1]                     = svattr(t, 'l3d',                   _svi, is_vtx) if _sv_direct else -1.
        branches[f"{_sv_slot}_mindx"][-1]                   = svattr(t, 'mindx',                 _svi, is_vtx) if _sv_direct else -1.
        branches[f"{_sv_slot}_mindy"][-1]                   = svattr(t, 'mindy',                 _svi, is_vtx) if _sv_direct else -1.
        branches[f"{_sv_slot}_mindz"][-1]                   = svattr(t, 'mindz',                 _svi, is_vtx) if _sv_direct else -1.
        branches[f"{_sv_slot}_mindxy"][-1]                  = svattr(t, 'mindxy',                _svi, is_vtx) if _sv_direct else -1.
        branches[f"{_sv_slot}_mind3d"][-1]                  = svattr(t, 'mind3d',                _svi, is_vtx) if _sv_direct else -1.
        branches[f"{_sv_slot}_maxdx"][-1]                   = svattr(t, 'maxdx',                 _svi, is_vtx) if _sv_direct else -1.
        branches[f"{_sv_slot}_maxdy"][-1]                   = svattr(t, 'maxdy',                 _svi, is_vtx) if _sv_direct else -1.
        branches[f"{_sv_slot}_maxdz"][-1]                   = svattr(t, 'maxdz',                 _svi, is_vtx) if _sv_direct else -1.
        branches[f"{_sv_slot}_maxdxy"][-1]                  = svattr(t, 'maxdxy',                _svi, is_vtx) if _sv_direct else -1.
        branches[f"{_sv_slot}_maxd3d"][-1]                  = svattr(t, 'maxd3d',                _svi, is_vtx) if _sv_direct else -1.
        branches[f"{_sv_slot}_onModule"][-1]                = svattr(t, 'onModule',              _svi, is_vtx) if _sv_direct else -1.
        branches[f"{_sv_slot}_onModuleWithinUnc"][-1]       = svattr(t, 'onModuleWithinUnc',     _svi, is_vtx) if _sv_direct else -1.
        branches[f"{_sv_slot}_minDistanceFromDet"][-1]      = svattr(t, 'minDistanceFromDet',    _svi, is_vtx) if _sv_direct else -1.
        branches[f"{_sv_slot}_minDistanceFromDet_x"][-1]    = svattr(t, 'minDistanceFromDet_x',  _svi, is_vtx) if _sv_direct else -1.
        branches[f"{_sv_slot}_minDistanceFromDet_y"][-1]    = svattr(t, 'minDistanceFromDet_y',  _svi, is_vtx) if _sv_direct else -1.
        branches[f"{_sv_slot}_minDistanceFromDet_z"][-1]    = svattr(t, 'minDistanceFromDet_z',  _svi, is_vtx) if _sv_direct else -1.
        branches[f"{_sv_slot}_closestDet_x"][-1]            = svattr(t, 'closestDet_x',          _svi, is_vtx) if _sv_direct else -1.
        branches[f"{_sv_slot}_closestDet_y"][-1]            = svattr(t, 'closestDet_y',          _svi, is_vtx) if _sv_direct else -1.
        branches[f"{_sv_slot}_closestDet_z"][-1]            = svattr(t, 'closestDet_z',          _svi, is_vtx) if _sv_direct else -1.

        # Lifetime reweighting
        tweight = 1.0
        if (reweightTo > 0 and reweightFrom > 0 and not isData):
            #if "Scenario" in sampleTag and len(LLPs_lxy) < 1:
            #    tweight = 1.0
            if "Scenario" in sampleTag:
                for ict in range(len(LLPs_ct)):
                    tweight = tweight * (reweightFrom/reweightTo) * math.exp(LLPs_ct[ict]/reweightFrom - LLPs_ct[ict]/reweightTo)
            else:
                ict = getClosestAngular([lxy, v.Eta(), v.Phi()], LLPs_lxy, LLPs_eta, LLPs_phi)
                tweight = (reweightFrom/reweightTo) * math.exp(LLPs_ct[ict]/reweightFrom - LLPs_ct[ict]/reweightTo)
        # Scan:
        sf_trg, sf_trg_up, sf_trg_down  = 1., 1., 1.
        sf_sel, sf_sel_up, sf_sel_down  = 1., 1., 1.
        if not isData and not args.noSF:
            sf_trg, sf_trg_up, sf_trg_down = getTriggerSF(subpt, lxy)
            sf_sel, sf_sel_up, sf_sel_down = getSelectionSF(lxy)
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
                roow.setVal(lumiweight*rooweight*sf_trg*sf_sel*tweight*bweight*puweight);
                roow_trg_up.setVal(lumiweight*rooweight*sf_sel*sf_trg_up*tweight*bweight*puweight);
                roow_trg_down.setVal(lumiweight*rooweight*sf_sel*sf_trg_down*tweight*bweight*puweight);
                roow_sel_up.setVal(lumiweight*rooweight*sf_trg*sf_sel_up*tweight*bweight*puweight);
                roow_sel_down.setVal(lumiweight*rooweight*sf_trg*sf_sel_down*tweight*bweight*puweight);
                roods[slice].add(ROOT.RooArgSet(mfit,roow),roow.getVal());
                roods_trg_up[slice].add(ROOT.RooArgSet(mfit,roow_trg_up),roow_trg_up.getVal());
                roods_trg_down[slice].add(ROOT.RooArgSet(mfit,roow_trg_down),roow_trg_down.getVal());
                roods_sel_up[slice].add(ROOT.RooArgSet(mfit,roow_sel_up),roow_sel_up.getVal());
                roods_sel_down[slice].add(ROOT.RooArgSet(mfit,roow_sel_down),roow_sel_down.getVal());
                roods["Dimuon_"+label+"_inclusive"].add(ROOT.RooArgSet(mfit,roow),roow.getVal());
                roods_trg_up["Dimuon_"+label+"_inclusive"].add(ROOT.RooArgSet(mfit,roow_trg_up),roow_trg_up.getVal());
                roods_trg_down["Dimuon_"+label+"_inclusive"].add(ROOT.RooArgSet(mfit,roow_trg_down),roow_trg_down.getVal());
                roods_sel_up["Dimuon_"+label+"_inclusive"].add(ROOT.RooArgSet(mfit,roow_sel_up),roow_sel_up.getVal());
                roods_sel_down["Dimuon_"+label+"_inclusive"].add(ROOT.RooArgSet(mfit,roow_sel_down),roow_sel_down.getVal());
                roods["Dimuon_full_inclusive"].add(ROOT.RooArgSet(mfit,roow),roow.getVal());
                roods_trg_up["Dimuon_full_inclusive"].add(ROOT.RooArgSet(mfit,roow_trg_up),roow_trg_up.getVal());
                roods_trg_down["Dimuon_full_inclusive"].add(ROOT.RooArgSet(mfit,roow_trg_down),roow_trg_down.getVal());
                roods_sel_up["Dimuon_full_inclusive"].add(ROOT.RooArgSet(mfit,roow_sel_up),roow_sel_up.getVal());
                roods_sel_down["Dimuon_full_inclusive"].add(ROOT.RooArgSet(mfit,roow_sel_down),roow_sel_down.getVal());
                catmass[slice].Fill(mass, lumiweight*rooweight*sf_trg*sf_sel*tweight*bweight*puweight);
                catmass["Dimuon_"+label+"_inclusive"].Fill(mass, lumiweight*rooweight*sf_trg*sf_sel*tweight*bweight*puweight);
                catmass["Dimuon_full_inclusive"].Fill(mass, lumiweight*rooweight*sf_trg*sf_sel*tweight*bweight*puweight);
                filledcat2mu = True

    # Apply selections and fill histograms for muon pairs from overlapping SVs
    # (now merged into the dmuvec_all loop above; kept here as dead code for reference)
    for vn,v in enumerate([]):
        if args.noDiMuon:
            break
        if not applyDiMuonSelection(v):
            continue
        lxy  = t.SVOverlap_lxy[osvidx_noVtx[vn]]
        if not applyLxySelection(lxy):
            continue
        lz  = abs(t.SVOverlap_z[osvidx_noVtx[vn]])
        if not applyLzSelection(lz):
            continue
        mass = v.M()
        pt   = v.Pt()
        nhitsbeforesvtotal = t.Muon_nhitsbeforesv[dmuidxs_osv_noVtx[int(vn*2)]] + t.Muon_nhitsbeforesv[dmuidxs_osv_noVtx[int(vn*2)+1]]
        subpt = t.Muon_pt[dmuidxs_osv_noVtx[int(vn*2)]] if t.Muon_pt[dmuidxs_osv_noVtx[int(vn*2)]] < t.Muon_pt[dmuidxs_osv_noVtx[int(vn*2)+1]] else t.Muon_pt[dmuidxs_osv_noVtx[int(vn*2)+1]]
        #  Apply selection on muon lifetime-scaled dxy
        if applyMuonIPSel:
            if ( abs(t.Muon_dxyCorr[dmuidxs_osv_noVtx[int(vn*2)]]  )/(lxy*mass/pt)<0.1 or
                 abs(t.Muon_dxyCorr[dmuidxs_osv_noVtx[int(vn*2)+1]])/(lxy*mass/pt)<0.1 ):
                continue
            if ( abs(t.Muon_dxyCorr[dmuidxs_osv_noVtx[int(vn*2)]]/t.Muon_dxye[dmuidxs_osv_noVtx[int(vn*2)]])<2.0 or abs(t.Muon_dxyCorr[dmuidxs_osv_noVtx[int(vn*2)+1]]/t.Muon_dxye[dmuidxs_osv_noVtx[int(vn*2)+1]])<2.0 ):
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
        if applyDiMuonResonanceMasking:
            if (mass > 0.41 and mass < 0.50):
                continue
            if (mass > 0.51 and mass < 0.59):
                continue
            if (mass > 0.731 and mass < 0.83): # 0.73
                continue
            if (mass > 0.96 and mass < 1.08):
                continue
            if (mass > 2.91 and mass < 3.27):
                continue
            if (mass > 3.47 and mass < 3.89):
                continue
            if (mass > 8.99 and mass < 9.91):
                continue
            if (mass > 9.64 and mass < 10.56):
                continue
            if (mass > 9.90 and mass < 10.78):
                continue
            if (mass > 81 and mass < 101):
                continue
        #
        drmm = t.Muon_vec[dmuidxs_osv_noVtx[int(vn*2)]].DeltaR(t.Muon_vec[dmuidxs_osv_noVtx[int(vn*2)+1]])
        dpmm = abs(t.Muon_vec[dmuidxs_osv_noVtx[int(vn*2)]].DeltaPhi(t.Muon_vec[dmuidxs_osv_noVtx[int(vn*2)+1]]))
        demm = abs(t.Muon_vec[dmuidxs_osv_noVtx[int(vn*2)]].Eta()-t.Muon_vec[dmuidxs_osv_noVtx[int(vn*2)+1]].Eta())
        dedpmm = 1e6
        if dpmm>0.0:
            dedpmm = demm/dpmm
        a3dmm = abs(t.Muon_vec[dmuidxs_osv_noVtx[int(vn*2)]].Angle(t.Muon_vec[dmuidxs_osv_noVtx[int(vn*2)+1]].Vect()))
        #
        drmmu = dmu_muvecdp_osv_noVtx[int(vn*2)].DeltaR(dmu_muvecdp_osv_noVtx[int(vn*2)+1])
        dpmmu = abs(dmu_muvecdp_osv_noVtx[int(vn*2)].DeltaPhi(dmu_muvecdp_osv_noVtx[int(vn*2)+1]))
        demmu = abs(dmu_muvecdp_osv_noVtx[int(vn*2)].Eta()-dmu_muvecdp_osv_noVtx[int(vn*2)+1].Eta())
        dedpmmu = 1e6
        if dpmmu>0.0:
            dedpmmu = demmu/dpmmu
        a3dmmu = abs(dmu_muvecdp_osv_noVtx[int(vn*2)].Angle(dmu_muvecdp_osv_noVtx[int(vn*2)+1].Vect()))
        #
        dphisv = abs(v.Vect().DeltaPhi(osvvec_noVtx[vn]))
        dphisv1 = abs(t.Muon_vec[dmuidxs_osv_noVtx[int(vn*2)]].Vect().DeltaPhi(osvvec_noVtx[vn]))
        dphisv2 = abs(t.Muon_vec[dmuidxs_osv_noVtx[int(vn*2)+1]].Vect().DeltaPhi(osvvec_noVtx[vn]))
        detasv = abs(v.Vect().Eta()-osvvec_noVtx[vn].Eta())
        detadphisv = 1e6
        if dphisv>0.0:
            detadphisv = detasv/dphisv
        a3dsv  = abs(v.Vect().Angle(osvvec_noVtx[vn]))
        #
        desvdpsvmin = 1e6
        demmdpsvmin = 1e6
        dpsvmin = 1e6
        tdpsv = abs(t.Muon_vec[dmuidxs_osv_noVtx[int(vn*2)]].Vect().DeltaPhi(osvvec_noVtx[vn]))
        if tdpsv < dpsvmin:
            dpsvmin = tdpsv
        tdpsv = abs(t.Muon_vec[dmuidxs_osv_noVtx[int(vn*2)+1]].Vect().DeltaPhi(osvvec_noVtx[vn]))
        if tdpsv < dpsvmin:
            dpsvmin = tdpsv
        if dpsvmin>0.0:
            demmdpsvmin = demm/dpsvmin
            desvdpsvmin = detasv/dpsvmin
        #
        vu = (dmu_muvecdp_osv_noVtx[int(vn*2)]+dmu_muvecdp_osv_noVtx[int(vn*2)+1])
        dphisvu = abs(vu.Vect().DeltaPhi(osvvec_noVtx[vn]))
        dphisv1u = abs(dmu_muvecdp_osv_noVtx[int(vn*2)].Vect().DeltaPhi(osvvec_noVtx[vn]))
        dphisv2u = abs(dmu_muvecdp_osv_noVtx[int(vn*2)+1].Vect().DeltaPhi(osvvec_noVtx[vn]))
        detasvu = abs(vu.Vect().Eta()-osvvec_noVtx[vn].Eta())
        detadphisvu = 1e6
        if dphisvu>0.0:
            detadphisvu = detasvu/dphisvu
        a3dsvu  = abs(vu.Vect().Angle(osvvec_noVtx[vn]))
        #
        desvdpsvminu = 1e6
        demmdpsvminu = 1e6
        dpsvminu = 1e6
        tdpsvu = abs(dmu_muvecdp_osv_noVtx[int(vn*2)].Vect().DeltaPhi(osvvec_noVtx[vn]))
        if tdpsvu < dpsvminu:
            dpsvminu = tdpsvu
        tdpsvu = abs(dmu_muvecdp_osv_noVtx[int(vn*2)+1].Vect().DeltaPhi(osvvec_noVtx[vn]))
        if tdpsvu < dpsvminu:
            dpsvminu = tdpsvu
        if dpsvminu>0.0:
            demmdpsvminu = demmu/dpsvminu
            desvdpsvminu = detasvu/dpsvminu
        #
        sindpsvlxy = abs(ROOT.TMath.Sin(v.Vect().DeltaPhi(osvvec_noVtx[vn])))*lxy
        sindpsvlxyu = abs(ROOT.TMath.Sin(vu.Vect().DeltaPhi(osvvec_noVtx[vn])))*lxy
        sina3dsvl3d = abs(ROOT.TMath.Sin(v.Vect().Angle(osvvec_noVtx[vn])))*(osvvec_noVtx[vn].Mag())
        sina3dsvl3du = abs(ROOT.TMath.Sin(vu.Vect().Angle(osvvec_noVtx[vn])))*(osvvec_noVtx[vn].Mag())
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
        isocat = dimuonIsoCategory(t.Muon_PFIsoAll0p4[dmuidxs_osv_noVtx[int(vn*2)]], t.Muon_pt[dmuidxs_osv_noVtx[int(vn*2)]], t.Muon_PFIsoAll0p4[dmuidxs_osv_noVtx[int(vn*2)+1]], t.Muon_pt[dmuidxs_osv_noVtx[int(vn*2)+1]])
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

        #
        # Branch filling
        if branches["SV1_x"][-1] < 0.:
            branches["SV1_x"][-1]    =  t.SV_x[osvidx_noVtx[vn]]
            branches["SV1_y"][-1]    =  t.SV_y[osvidx_noVtx[vn]]
            branches["SV1_z"][-1]    =  t.SV_z[osvidx_noVtx[vn]]
        elif branches["SV2_x"][-1] < 0.:
            branches["SV2_x"][-1]    =  t.SV_x[osvidx_noVtx[vn]]
            branches["SV2_y"][-1]    =  t.SV_y[osvidx_noVtx[vn]]
            branches["SV2_z"][-1]    =  t.SV_z[osvidx_noVtx[vn]]
        else:
            break


        # Lifetime reweighting
        tweight = 1.0
        #if (reweightTo > 0 and reweightFrom > 0 and not isData): # provisional, could be adapted, but better go for the lower option...
        #    ict, ival = getClosest(lxy, LLPs_lxy)
        #    tweight = (reweightFrom/reweightTo) * math.exp(LLPs_ct[ict]/reweightFrom - LLPs_ct[ict]/reweightTo)
        if (reweightTo > 0 and reweightFrom > 0 and not isData):
            #if "Scenario" in sampleTag and len(LLPs_lxy) < 1:
            #    tweight = 1.0
            if "Scenario" in sampleTag:
                for ict in range(len(LLPs_ct)):
                    tweight = tweight * (reweightFrom/reweightTo) * math.exp(LLPs_ct[ict]/reweightFrom - LLPs_ct[ict]/reweightTo)
            else:
                ict = getClosestAngular([lxy, v.Eta(), v.Phi()], LLPs_lxy, LLPs_eta, LLPs_phi)
                tweight = (reweightFrom/reweightTo) * math.exp(LLPs_ct[ict]/reweightFrom - LLPs_ct[ict]/reweightTo)
        # Scan:
        sf_trg, sf_trg_up, sf_trg_down  = 1., 1., 1.
        sf_sel, sf_sel_up, sf_sel_down  = 1., 1., 1.
        if not isData and not args.noSF:
            sf_trg, sf_trg_up, sf_trg_down = getTriggerSF(subpt, lxy)
            sf_sel, sf_sel_up, sf_sel_down = getSelectionSF(lxy)
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
            roow.setVal(lumiweight*rooweight*sf_trg*sf_sel*tweight*bweight*puweight);
            roow_trg_up.setVal(lumiweight*rooweight*sf_sel*sf_trg_up*tweight*bweight*puweight);
            roow_trg_down.setVal(lumiweight*rooweight*sf_sel*sf_trg_down*tweight*bweight*puweight);
            roow_sel_up.setVal(lumiweight*rooweight*sf_trg*sf_sel_up*tweight*bweight*puweight);
            roow_sel_down.setVal(lumiweight*rooweight*sf_trg*sf_sel_down*tweight*bweight*puweight);
            roods[slice].add(ROOT.RooArgSet(mfit,roow),roow.getVal());
            roods_trg_up[slice].add(ROOT.RooArgSet(mfit,roow_trg_up),roow_trg_up.getVal());
            roods_trg_down[slice].add(ROOT.RooArgSet(mfit,roow_trg_down),roow_trg_down.getVal());
            roods_sel_up[slice].add(ROOT.RooArgSet(mfit,roow_sel_up),roow_sel_up.getVal());
            roods_sel_down[slice].add(ROOT.RooArgSet(mfit,roow_sel_down),roow_sel_down.getVal());
            roods["Dimuon_"+label+"_inclusive"].add(ROOT.RooArgSet(mfit,roow),roow.getVal());
            roods_trg_up["Dimuon_"+label+"_inclusive"].add(ROOT.RooArgSet(mfit,roow_trg_up),roow_trg_up.getVal());
            roods_trg_down["Dimuon_"+label+"_inclusive"].add(ROOT.RooArgSet(mfit,roow_trg_down),roow_trg_down.getVal());
            roods_sel_up["Dimuon_"+label+"_inclusive"].add(ROOT.RooArgSet(mfit,roow_sel_up),roow_sel_up.getVal());
            roods_sel_down["Dimuon_"+label+"_inclusive"].add(ROOT.RooArgSet(mfit,roow_sel_down),roow_sel_down.getVal());
            roods["Dimuon_full_inclusive"].add(ROOT.RooArgSet(mfit,roow),roow.getVal());
            roods_trg_up["Dimuon_full_inclusive"].add(ROOT.RooArgSet(mfit,roow_trg_up),roow_trg_up.getVal());
            roods_trg_down["Dimuon_full_inclusive"].add(ROOT.RooArgSet(mfit,roow_trg_down),roow_trg_down.getVal());
            roods_sel_up["Dimuon_full_inclusive"].add(ROOT.RooArgSet(mfit,roow_sel_up),roow_sel_up.getVal());
            roods_sel_down["Dimuon_full_inclusive"].add(ROOT.RooArgSet(mfit,roow_sel_down),roow_sel_down.getVal());
            catmass[slice].Fill(mass, lumiweight*rooweight*sf_trg*sf_sel*tweight*bweight*puweight);
            catmass["Dimuon_"+label+"_inclusive"].Fill(mass, lumiweight*rooweight*sf_trg*sf_sel*tweight*bweight*puweight);
            catmass["Dimuon_full_inclusive"].Fill(mass, lumiweight*rooweight*sf_trg*sf_sel*tweight*bweight*puweight);
            filledcat2mu = True
    
    ######### Fill the nutples


### Write histograms
foname = "%s/histograms_%s_all.root"%(outdir,args.year)
if args.inSample!="*":
    if args.inFile!="*":
        foname = "%s/histograms_file%s_%s_%s"%(outdir,args.inFile,args.inSample,args.year)
    else:
        foname = "%s/histograms_%s_%s"%(outdir,args.inSample,args.year)
if index>=0:
    foname = foname+("_%d"%index)

### Write ntuple
import uproot

sv1_mask = np.array(branches["SV1_x"]) != -1.
filtered_branches = {k: np.array(v)[sv1_mask] for k, v in branches.items()}
_tuples_subdir = args.tupleDir if args.tupleDir else ("tuples_parking_nochi2" if args.noMuonChi2Sel else "tuples_parking")
_tuples_dir = os.path.join(
    os.environ.get("STARTDIR", "/home/users/garciaja/fullRun3/CMSSW_15_0_2/src/run3_scouting"),
    _tuples_subdir
)
os.makedirs(_tuples_dir, exist_ok=True)
if args.inSample != "*":
    ntuple_name = f"{_tuples_dir}/tuples_{args.inSample}_{args.year}"
    if args.inFile != "*":
        ntuple_name += f"_{args.inFile}"
    if index >= 0:
        ntuple_name += f"_{index}"
    ntuple_name += ".root"
else:
    ntuple_name = f"{_tuples_dir}/tuples_{args.year}_all.root"
with uproot.recreate(ntuple_name) as f:
    f["tuples"] = filtered_branches


if not isData:
    ## Convention for ctau (uniform for every sample)
    if reweightFrom > 0 and reweightTo > 0:
        foname = foname.replace("ctau-%imm"%reweightFrom, "ctau-%.2fmm"%(reweightTo))
        foname = foname.replace("ctau-%ip0mm"%reweightFrom, "ctau-%.2fmm"%(reweightTo)) # In case we have p0 in the name
    elif 'ctau-' in foname:
        ctau_string = foname.split('ctau-')[1].split('mm')[0]
        if 'p' in ctau_string:
            ctau_string_ = ctau_string.replace('p','.')
        else:
            ctau_string_ = ctau_string
        foname = foname.replace(ctau_string+'mm', "%.2fmm"%(float(ctau_string_)))
    ## Convention for uniform mass for BToPhi
    if 'BToPhi' in foname:
        mass_string = foname.split('MPhi-')[1].split('_ctau')[0]
        mass_value = float(mass_string.replace('p','.'))
        new_mass_string = ('%.2f'%(mass_value)).replace('.','p')
        foname = foname.replace(mass_string, new_mass_string)
        
fout = ROOT.TFile(foname+".root","RECREATE")
fout.cd()
for cat in h1d.keys():
    for h in h1d[cat]:
        h.Write()
for cat in h2d.keys():
    for h in h2d[cat]:
        h.Write()
for dbin in dbins:
    print("RooDataSet {}  with {}({}) entries".format(dbin, roods[dbin].sumEntries(), roods[dbin].numEntries()))
    roods[dbin].Write()
    print("RooDataSet (up) {}  with {} entries".format(dbin, roods_trg_up[dbin].sumEntries()))
    roods_trg_up[dbin].Write()
    print("RooDataSet (down) {}  with {} entries".format(dbin, roods_trg_down[dbin].sumEntries()))
    roods_trg_down[dbin].Write()
    print("RooDataSet (up) {}  with {} entries".format(dbin, roods_sel_up[dbin].sumEntries()))
    roods_sel_up[dbin].Write()
    print("RooDataSet (down) {}  with {} entries".format(dbin, roods_sel_down[dbin].sumEntries()))
    roods_sel_down[dbin].Write()
    catmass[dbin].Write()
catmass["FourMu_osv_dimuonmass"].Write()
catmass["FourMu_sep_dimuonmass"].Write()
#
fout.Close()