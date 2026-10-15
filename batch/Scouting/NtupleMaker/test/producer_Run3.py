import FWCore.ParameterSet.VarParsing as VarParsing

opts     = VarParsing.VarParsing('python')
vpbool   = VarParsing.VarParsing.varType.bool
vpint    = VarParsing.VarParsing.varType.int
vpstring = VarParsing.VarParsing.varType.string

opts.register('data',       True,          mytype = vpbool)
opts.register('monitor',    False,         mytype = vpbool)
opts.register('background', False,         mytype = vpbool) # apply skim filters for background MC (e.g. QCD)
opts.register('era',        "2022D",       mytype = vpstring)
opts.register('output',     "output.root", mytype = vpstring)
opts.register('inputs',     "",            mytype = vpstring) # comma separated list of input files
opts.register('nevents',    -1,            mytype = vpint)
opts.register('testL1',     False,         mytype = vpbool)
opts.register('isMiniAOD',  True,          mytype = vpbool) # set False for private RAW/HLT samples
opts.parseArguments()

def convert_fname(fname):
    if fname.startswith("/store") or fname.startswith("root://"): return fname
    return "file:" + fname.replace("file:","",1)

inputs = []
if opts.inputs:
    inputs = map(convert_fname,opts.inputs.split(","))

lin=list(inputs)

print(""" Running with options:data = {} era = {} output = {} inputs = {} nevents = {} """.format(opts.data,opts.era,opts.output,lin,opts.nevents))

import FWCore.ParameterSet.Config as cms
process = cms.Process('SLIM')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

gtag=""
# from https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions
if opts.data:
    if '2022' in opts.era:
        gtag="124X_dataRun3_Prompt_v10" # latest prompt RECO GT
        #gtag="124X_dataRun3_Prompt_v4"
        #gtag="124X_dataRun3_HLT_v7" # latest HLT GT
    elif '2023' in opts.era:
        #gtag="130X_dataRun3_Prompt_frozen_v3" # latest prompt RECO GT (CMSSW>=13_0_10)
        gtag="130X_dataRun3_Prompt_v4"
        #gtag="130X_dataRun3_HLT_frozen_v3" # latest HLT GT (CMSSW>=13_0_10)
    elif '2024' in opts.era:
        gtag = "140X_dataRun3_Prompt_v4"
    elif '2025' in opts.era:
        gtag = "150X_dataRun3_Prompt_v1"
else:
    if '2022' in opts.era:
        if not 'postEE' in opts.era:
            #gtag="124X_mcRun3_2022_realistic_v12" # latest MC GT
            gtag="130X_mcRun3_2022_realistic_v5" # Found for central production
        else:
            #gtag="124X_mcRun3_2022_realistic_postEE_v1"
            gtag="130X_mcRun3_2022_realistic_postEE_v6" # Found for central production
    elif '2024' in opts.era:
        gtag="150X_mcRun3_2024_realistic_v2" # Found for central production (RunIII2024Summer24DRPremix)
    else:
        #gtag="130X_mcRun3_2023_realistic_v9" # latest MC GT (=phase1_2023_realistic, in CMSSW_13_1_0)
        if not 'BPix' in opts.era:
            #gtag="124X_mcRun3_2022_realistic_v12" # latest MC GT
            gtag="130X_mcRun3_2023_realistic_v14" # Found for central production
        else:
            #gtag="124X_mcRun3_2022_realistic_postEE_v1"
            gtag="130X_mcRun3_2023_realistic_postBPix_v2" # Found for central production
process.GlobalTag.globaltag = gtag

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(opts.nevents))

if '2024' in opts.era or '2025' in opts.era:
    process.options = cms.untracked.PSet(TryToContinue = cms.untracked.vstring('ProductNotFound'))
else:
    process.options = cms.untracked.PSet(SkipEvent = cms.untracked.vstring('ProductNotFound'))


if opts.data:
    process.MessageLogger.cerr.FwkReport.reportEvery = 1000
else:
    process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource",
    dropDescendantsOfDroppedBranches = cms.untracked.bool(True),
    fileNames = cms.untracked.vstring(),
    bypassVersionCheck = cms.untracked.bool(True),
    inputCommands = cms.untracked.vstring(
        'keep *',
        'drop *_hltScoutingTrackPacker_*_*',
    ),
)
process.source.fileNames = lin

if not opts.data:
    process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

# skim data and background MC; keep all events for signal MC acceptance calculations
do_skim = opts.data or opts.background

# Build output keep commands based on era (2024+ uses Vtx/NoVtx split collections)
if '2024' in opts.era or '2025' in opts.era:
    out_keep_muon = ["keep *_hltScoutingMuonPackerVtx_*_*", "keep *_hltScoutingMuonPackerNoVtx_*_*"]
    out_keep_hit  = ["keep *_vertexMakerVtx_*_*", "keep *_hitMakerVtx_*_*", "keep *_hitMakerNoVtx_*_*"]
else:
    out_keep_muon = ["keep *_hltScoutingMuonPacker_*_*"]
    out_keep_hit  = ["keep *_hitMaker_*_*"]

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('file:output.root'),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('skimpath'),
        ),
    outputCommands = cms.untracked.vstring(
        ["drop *"] + out_keep_muon + [
        "keep *_hltScoutingPFPacker_*_*",
        "keep *_hltScoutingPrimaryVertexPacker_*_*",
        "keep *_triggerMaker_*_*",
        ] + out_keep_hit + [
        "keep *_beamSpotMaker_*_*",
        "keep *_genParticles_*_HLT",        # private-HLT inputs (un-pruned)
        "keep *_prunedGenParticles_*_PAT",  # central MiniAOD inputs
        "keep *_addPileupInfo_*_*",         # private-HLT PU info
        "keep *_slimmedAddPileupInfo_*_PAT",
        ]),
     basketSize = cms.untracked.int32(128*1024), # 128kb basket size instead of ~30kb default
)

if len(opts.output):
    process.out.fileName = "file:" + str(opts.output).strip().replace("file:","")

process.outpath = cms.EndPath(process.out)

process.Timing = cms.Service("Timing",
        summaryOnly = cms.untracked.bool(True)
        )

if '2024' in opts.era or '2025' in opts.era:
    process.countmuVtx = cms.EDFilter("ScoutingMuonCountFilter",
        src = cms.InputTag("hltScoutingMuonPackerVtx"),
        minNumber = cms.uint32(2)
    )
    process.countmuNoVtx = cms.EDFilter("ScoutingMuonCountFilter",
        src = cms.InputTag("hltScoutingMuonPackerNoVtx"),
        minNumber = cms.uint32(2)
    )
    process.countvtxNoVtx = cms.EDFilter("ScoutingVertexCountFilter",
        src = cms.InputTag("hltScoutingMuonPackerNoVtx","displacedVtx"),
        minNumber = cms.uint32(1)
    )
else:
    process.countmu = cms.EDFilter("ScoutingMuonCountFilter",
        src = cms.InputTag("hltScoutingMuonPacker"),
        minNumber = cms.uint32(2)
    )
    process.countvtx = cms.EDFilter("ScoutingVertexCountFilter",
        src = cms.InputTag("hltScoutingMuonPacker","displacedVtx"),
        minNumber = cms.uint32(1)
    )

# To recover trigger info: https://hlt-config-editor-confdbv3.app.cern.ch/ + https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmVRun3Analysis + https://cmsoms.cern.ch/cms/triggers/report?cms_run=<run>
L1Info  = []
HLTInfo = []
if '2022' in opts.era or (opts.data and '2023B' in opts.era) or '2023C-triggerV10' in opts.era:
    L1Info = ["L1_DoubleMu_12_5","L1_DoubleMu_15_7","L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7","L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18","L1_DoubleMu4_SQ_OS_dR_Max1p2","L1_DoubleMu4p5_SQ_OS_dR_Max1p2"]
    if opts.testL1:
        L1Info = L1Info + ["L1_DoubleMu0_Upt15_Upt7", "L1_DoubleMu0_Upt6_IP_Min1_Upt4"]
    HLTInfo=[ [ 'Run3_PFScouting', 'DST_Run3_PFScoutingPixelTracking_v*' ] ]
    if opts.monitor:
        HLTInfo = HLTInfo + [
            ['HLT_Ele115_CaloIdVT_GsfTrkIdT', 'HLT_Ele115_CaloIdVT_GsfTrkIdT_v*'],
            ['HLT_Ele35_WPTight_Gsf', 'HLT_Ele35_WPTight_Gsf_v*'],
            ['HLT_IsoMu27', 'HLT_IsoMu27_v*'],
            ['HLT_Mu50', 'HLT_Mu50_v*'],
            ['HLT_PFHT1050', 'HLT_PFHT1050_v*'],
            ['HLT_Photon200', 'HLT_Photon200_v*']
        ]
        L1Info = L1Info + [
            "L1_HTT200er", "L1_HTT255er", "L1_HTT280er", "L1_HTT320er", "L1_HTT360er", "L1_HTT400er", "L1_HTT450er",
            "L1_ETT2000",
            "L1_SingleJet180", "L1_SingleJet200",
            "L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5", "L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5", "L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5",
            "L1_SingleLooseIsoEG28er2p1", "L1_SingleLooseIsoEG28er1p5", "L1_SingleLooseIsoEG30er1p5", "L1_SingleIsoEG28er2p1", "L1_SingleIsoEG30er2p1", "L1_SingleIsoEG32er2p1",
            "L1_DoubleEG_LooseIso16_LooseIso12_er1p5", "L1_DoubleEG_LooseIso18_LooseIso12_er1p5", "L1_DoubleEG_LooseIso20_LooseIso12_er1p5", "L1_DoubleEG_LooseIso22_LooseIso12_er1p5"
        ]
        L1Info = L1Info + [
            "L1_SingleEG34er2p5", "L1_SingleEG36er2p5", "L1_SingleEG38er2p5", "L1_SingleEG40er2p5", "L1_SingleJet160er2p5", "L1_SingleJet180", "L1_SingleJet200", "L1_SingleTau120er2p1", "L1_SingleTau130er2p1", "L1_SingleEG42er2p5", "L1_SingleEG45er2p5", "L1_SingleEG60"
        ]
        L1Info = L1Info + [
            "L1_SingleLooseIsoEG26er2p5", "L1_SingleLooseIsoEG26er1p5", "L1_SingleLooseIsoEG28er2p5", "L1_SingleLooseIsoEG28er2p1", "L1_SingleLooseIsoEG28er1p5", "L1_SingleLooseIsoEG30er2p5", "L1_SingleLooseIsoEG30er1p5", "L1_SingleEG26er2p5", "L1_SingleEG38er2p5", "L1_SingleEG40er2p5", "L1_SingleEG42er2p5", "L1_SingleEG45er2p5", "L1_SingleEG60", "L1_SingleEG34er2p5", "L1_SingleEG36er2p5", "L1_SingleIsoEG24er2p1", "L1_SingleIsoEG26er2p1", "L1_SingleIsoEG28er2p1", "L1_SingleIsoEG30er2p1", "L1_SingleIsoEG32er2p1", "L1_SingleIsoEG26er2p5", "L1_SingleIsoEG28er2p5", "L1_SingleIsoEG30er2p5", "L1_SingleIsoEG32er2p5", "L1_SingleIsoEG34er2p5"
        ]
        L1Info = L1Info + [
            "L1_SingleMu22", "L1_SingleMu25"
        ]
        L1Info = L1Info + [
            "L1_HTT120er", "L1_HTT160er", "L1_HTT200er", "L1_HTT255er", "L1_HTT280er", "L1_HTT280er_QuadJet_70_55_40_35_er2p5", "L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3", "L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3", "L1_HTT320er", "L1_HTT360er", "L1_ETT2000", "L1_HTT400er", "L1_HTT450er"
        ]
        L1Info = L1Info + [
            "L1_SingleEG34er2p5", "L1_SingleEG36er2p5", "L1_SingleEG38er2p5", "L1_SingleEG40er2p5", "L1_SingleJet160er2p5", "L1_SingleJet180", "L1_SingleJet200", "L1_SingleTau120er2p1", "L1_SingleTau130er2p1", "L1_SingleEG42er2p5", "L1_SingleEG45er2p5", "L1_SingleEG60"
        ]
        L1Info = list(set(L1Info))
elif ('2023' in opts.era) and ('2023C-triggerV10' not in opts.era):
    # for run >=367621 (during era Run2023C)
    L1Info = ["L1_DoubleMu_12_5","L1_DoubleMu_15_7","L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7","L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18","L1_DoubleMu4_SQ_OS_dR_Max1p2","L1_DoubleMu4p5_SQ_OS_dR_Max1p2","L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4","L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4","L1_DoubleMu8_SQ"]
    HLTInfo = [ [ 'Run3_DoubleMu3_PFScouting', 'DST_Run3_DoubleMu3_PFScoutingPixelTracking_v*' ] ]
    if opts.monitor:
        HLTInfo = HLTInfo + [
            ['DST_Run3_EG16_EG12_PFScoutingPixelTracking', 'DST_Run3_EG16_EG12_PFScoutingPixelTracking_v*'],
            ['DST_Run3_EG30_PFScoutingPixelTracking', 'DST_Run3_EG30_PFScoutingPixelTracking_v*'],
            ['DST_Run3_JetHT_PFScoutingPixelTracking', 'DST_Run3_JetHT_PFScoutingPixelTracking_v*'],
            ['HLT_Ele115_CaloIdVT_GsfTrkIdT', 'HLT_Ele115_CaloIdVT_GsfTrkIdT_v*'],
            ['HLT_Ele115_CaloIdVT_GsfTrkIdT', 'HLT_Ele115_CaloIdVT_GsfTrkIdT_v*'],
            ['HLT_Ele35_WPTight_Gsf', 'HLT_Ele35_WPTight_Gsf_v*'],
            ['HLT_IsoMu27', 'HLT_IsoMu27_v*'],
            ['HLT_Mu50', 'HLT_Mu50_v*'],
            ['HLT_PFHT1050', 'HLT_PFHT1050_v*'],
            ['HLT_Photon200', 'HLT_Photon200_v*']
        ]
        L1Info = L1Info + [
            "L1_DoubleEG_LooseIso16_LooseIso12_er1p5", "L1_DoubleEG_LooseIso18_LooseIso12_er1p5", "L1_DoubleEG_LooseIso20_LooseIso12_er1p5", "L1_DoubleEG_LooseIso22_LooseIso12_er1p5"
        ]
        L1Info = L1Info + [
            "L1_SingleLooseIsoEG28er2p1", "L1_SingleLooseIsoEG28er1p5", "L1_SingleLooseIsoEG30er1p5", "L1_SingleIsoEG28er2p1", "L1_SingleIsoEG30er2p1", "L1_SingleIsoEG32er2p1", "L1_SingleEG40er2p5", "L1_SingleEG42er2p5", "L1_SingleEG45er2p5"
        ]
        L1Info = L1Info + [
            "L1_HTT200er", "L1_HTT255er", "L1_HTT280er", "L1_HTT320er", "L1_HTT360er", "L1_HTT400er", "L1_HTT450er", "L1_ETT2000", "L1_SingleJet180", "L1_SingleJet200", "L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5", "L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5", "L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5", "L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5"
        ]
        L1Info = L1Info + [
            "L1_HTT200er", "L1_HTT255er", "L1_HTT280er", "L1_HTT320er", "L1_HTT360er", "L1_HTT400er", "L1_HTT450er",
            "L1_ETT2000",
            "L1_SingleJet180", "L1_SingleJet200",
            "L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5", "L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5", "L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5",
            "L1_SingleLooseIsoEG28er2p1", "L1_SingleLooseIsoEG28er1p5", "L1_SingleLooseIsoEG30er1p5", "L1_SingleIsoEG28er2p1", "L1_SingleIsoEG30er2p1", "L1_SingleIsoEG32er2p1",
            "L1_DoubleEG_LooseIso16_LooseIso12_er1p5", "L1_DoubleEG_LooseIso18_LooseIso12_er1p5", "L1_DoubleEG_LooseIso20_LooseIso12_er1p5", "L1_DoubleEG_LooseIso22_LooseIso12_er1p5"
        ]
        L1Info = L1Info + [
            "L1_SingleEG34er2p5", "L1_SingleEG36er2p5", "L1_SingleEG38er2p5", "L1_SingleEG40er2p5", "L1_SingleJet160er2p5", "L1_SingleJet180", "L1_SingleJet200", "L1_SingleTau120er2p1", "L1_SingleTau130er2p1", "L1_SingleEG42er2p5", "L1_SingleEG45er2p5", "L1_SingleEG60"
        ]
        L1Info = L1Info + [
            "L1_SingleLooseIsoEG26er2p5", "L1_SingleLooseIsoEG26er1p5", "L1_SingleLooseIsoEG28er2p5", "L1_SingleLooseIsoEG28er2p1", "L1_SingleLooseIsoEG28er1p5", "L1_SingleLooseIsoEG30er2p5", "L1_SingleLooseIsoEG30er1p5", "L1_SingleEG26er2p5", "L1_SingleEG38er2p5", "L1_SingleEG40er2p5", "L1_SingleEG42er2p5", "L1_SingleEG45er2p5", "L1_SingleEG60", "L1_SingleEG34er2p5", "L1_SingleEG36er2p5", "L1_SingleIsoEG24er2p1", "L1_SingleIsoEG26er2p1", "L1_SingleIsoEG28er2p1", "L1_SingleIsoEG30er2p1", "L1_SingleIsoEG32er2p1", "L1_SingleIsoEG26er2p5", "L1_SingleIsoEG28er2p5", "L1_SingleIsoEG30er2p5", "L1_SingleIsoEG32er2p5", "L1_SingleIsoEG34er2p5"
        ]
        L1Info = L1Info + [
            "L1_SingleMu22", "L1_SingleMu25"
        ]
        L1Info = L1Info + [
            "L1_HTT120er", "L1_HTT160er", "L1_HTT200er", "L1_HTT255er", "L1_HTT280er", "L1_HTT280er_QuadJet_70_55_40_35_er2p5", "L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3", "L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3", "L1_HTT320er", "L1_HTT360er", "L1_ETT2000", "L1_HTT400er", "L1_HTT450er"
        ]
        L1Info = L1Info + [
            "L1_SingleEG34er2p5", "L1_SingleEG36er2p5", "L1_SingleEG38er2p5", "L1_SingleEG40er2p5", "L1_SingleJet160er2p5", "L1_SingleJet180", "L1_SingleJet200", "L1_SingleTau120er2p1", "L1_SingleTau130er2p1", "L1_SingleEG42er2p5", "L1_SingleEG45er2p5", "L1_SingleEG60"
        ]        
        L1Info = list(set(L1Info))
elif '2025' in opts.era:
    L1Info = [
        "L1_DoubleMu_15_7",
        "L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7",
        "L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18",
        "L1_DoubleMu8_SQ",
        "L1_DoubleMu4er2p0_SQ_OS_dR_Max1p6",
        "L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4",
        "L1_DoubleMu4p5_SQ_OS_dR_Max1p2",
        "L1_DoubleMu0_Upt15_Upt7",
        "L1_DoubleMu0_Upt6_IP_Min1_Upt4",
        "L1_DoubleMu0_Upt6",
        "L1_DoubleMu0_Upt7",
        "L1_DoubleMu0_Upt8",
        "L1_SingleMu10_SQ14_BMTF",
        "L1_SingleMu11_SQ14_BMTF",
    ]
    HLTInfo = [
        [ 'Run3_DoubleMu_PFScouting',      'DST_PFScouting_DoubleMuon_v*' ],
        [ 'Run3_DoubleMuNoVtx_PFScouting', 'DST_PFScouting_DoubleMuonNoVtx_v*' ],
        [ 'Run3_DoubleMuVtx_PFScouting',   'DST_PFScouting_DoubleMuonVtx_v*' ],
        [ 'Run3_SingleMu_PFScouting',      'DST_PFScouting_SingleMuon_v*' ],
    ]
elif '2024' in opts.era or '2025' in opts.era:
    print("Running here")
    # 2024+: unprescaled DoubleMuon (no PixelTracking in path name) + SingleMuon L1/HLT seeds
    L1Info = [
        "L1_DoubleMu_15_7",
        "L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7",
        "L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18",
        "L1_DoubleMu8_SQ",
        "L1_DoubleMu4er2p0_SQ_OS_dR_Max1p6",
        "L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4",
        "L1_DoubleMu4p5_SQ_OS_dR_Max1p2",
        "L1_DoubleMu0_Upt15_Upt7",
        "L1_DoubleMu0_Upt6_IP_Min1_Upt4"
    ]
    # unconstrained pt seeds added in 2024F
    if opts.era in ["2024F", "2024G", "2024H", "2024I"]:
        L1Info = L1Info + ["L1_DoubleMu0_Upt6", "L1_DoubleMu0_Upt7", "L1_DoubleMu0_Upt8"]
    # Single muon seeds (only unprescaled)
    L1Info = L1Info + ["L1_SingleMu10_SQ14_BMTF", "L1_SingleMu11_SQ14_BMTF"]
    #
    HLTInfo = [
        [ 'Run3_DoubleMu_PFScouting', 'DST_PFScouting_DoubleMuon_v*' ],
        [ 'Run3_SingleMu_PFScouting', 'DST_PFScouting_SingleMuon_v*' ],
    ]
    #L1Info = list(set(L1Info))

do_trigger_objects = not (opts.data or opts.monitor)

process.triggerMaker = cms.EDProducer("TriggerMaker",
        triggerAlias = cms.vstring(list(zip(*HLTInfo))[0]),
        triggerSelection = cms.vstring(list(zip(*HLTInfo))[1]),
        triggerConfiguration = cms.PSet(
            hltResults            = cms.InputTag('TriggerResults','','HLT'),
            l1tResults            = cms.InputTag('','',''),
            l1tIgnoreMaskAndPrescale = cms.bool(False),
            throw                 = cms.bool(True),
            usePathStatus = cms.bool(False),
            ),
        doL1 = cms.bool(True),
        doTriggerObjects = cms.bool(do_trigger_objects and opts.isMiniAOD),
        isMiniAOD = cms.bool(('2024' in opts.era or '2025' in opts.era) and not opts.data and opts.isMiniAOD),
        #AlgInputTag = cms.InputTag("gtStage2Digis"),
        #AlgInputTag = cms.InputTag("gtStage2Digis","","RECO"),
        #l1tAlgBlkInputTag = cms.InputTag("gtStage2Digis"),
        #l1tAlgBlkInputTag = cms.InputTag("gtStage2Digis","","RECO"),
        #l1tExtBlkInputTag = cms.InputTag("gtStage2Digis"),
        #l1tExtBlkInputTag = cms.InputTag("gtStage2Digis","","RECO"),
        AlgInputTag        = cms.InputTag("gtStage2Digis", "", "RECO"),
        l1tAlgBlkInputTag  = cms.InputTag("gtStage2Digis", "", "RECO"),
        l1tExtBlkInputTag  = cms.InputTag("gtStage2Digis", "", "RECO"),
        ReadPrescalesFromFile = cms.bool(False),
        l1Seeds = cms.vstring(L1Info),
        )

if '2025' in opts.era:
    process.triggerMaker.triggerConfiguration.throw = cms.bool(False)
    if opts.data:
        process.triggerMaker.AlgInputTag       = cms.InputTag("gtStage2Digis")
        process.triggerMaker.l1tAlgBlkInputTag = cms.InputTag("gtStage2Digis")
        process.triggerMaker.l1tExtBlkInputTag = cms.InputTag("gtStage2Digis")

if '2024' in opts.era or '2025' in opts.era:

    from TrackingTools.TransientTrack.TransientTrackBuilder_cfi import *
    process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")

    process.vertexMakerVtx = cms.EDProducer("VertexMaker",
            muonVtxInputTag = cms.InputTag("hltScoutingMuonPackerVtx")
            )
    
    process.hitMakerVtx = cms.EDProducer("HitMaker",
            muonInputTag = cms.InputTag("hltScoutingMuonPackerVtx"),
            dvInputTag = cms.InputTag("vertexMakerVtx:verticesVtx"),
            vtxIndxInputTag = cms.InputTag("vertexMakerVtx", "vtxIndxVtx"),
            measurementTrackerEventInputTag = cms.InputTag("MeasurementTrackerEvent"),
            )
    process.hitMakerNoVtx = cms.EDProducer("HitMaker",
            muonInputTag = cms.InputTag("hltScoutingMuonPackerNoVtx"),
            dvInputTag = cms.InputTag("hltScoutingMuonPackerNoVtx:displacedVtx"),
            vtxIndxInputTag = cms.InputTag(""),
            measurementTrackerEventInputTag = cms.InputTag("MeasurementTrackerEvent"),
            )
else:
    process.hitMaker = cms.EDProducer("HitMaker",
            muonInputTag = cms.InputTag("hltScoutingMuonPacker"),
            dvInputTag = cms.InputTag("hltScoutingMuonPacker:displacedVtx"),
            measurementTrackerEventInputTag = cms.InputTag("MeasurementTrackerEvent"),
            )

process.beamSpotMaker = cms.EDProducer("BeamSpotMaker")

from RecoTracker.MeasurementDet.measurementTrackerEventDefault_cfi import measurementTrackerEventDefault as _measurementTrackerEventDefault
process.MeasurementTrackerEvent = _measurementTrackerEventDefault.clone()

process.load("EventFilter.L1TRawToDigi.gtStage2Digis_cfi")
if opts.monitor:
    process.gtStage2Digis.InputLabel = cms.InputTag( "rawDataCollector", "", "LHC" )
else:
    process.gtStage2Digis.InputLabel = cms.InputTag( "hltFEDSelectorL1" )

process.load("PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cfi")
#process.patTrigger.triggerResults = cms.InputTag("TriggerResults","","HLT")
#process.patTrigger.triggerEvent = cms.InputTag("hltTriggerSummaryAOD","","HLT")
process.patTrigger.stageL1Trigger = cms.uint32(2)

process.offlineBeamSpot = cms.EDProducer("BeamSpotProducer")

if '2024' in opts.era or '2025' in opts.era:
    if do_skim:
        process.skimpath_vtx   = cms.Path(process.countmuVtx+process.gtStage2Digis+process.triggerMaker+process.offlineBeamSpot+process.beamSpotMaker+process.MeasurementTrackerEvent+process.vertexMakerVtx+process.hitMakerVtx+process.hitMakerNoVtx)
        process.skimpath_novtx = cms.Path(process.countmuNoVtx+process.gtStage2Digis+process.countvtxNoVtx+process.triggerMaker+process.offlineBeamSpot+process.beamSpotMaker+process.MeasurementTrackerEvent+process.vertexMakerVtx+process.hitMakerVtx+process.hitMakerNoVtx)
        process.out.SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('skimpath_vtx', 'skimpath_novtx'))
    else: 
        process.skimpath = cms.Path(process.gtStage2Digis+process.triggerMaker+process.offlineBeamSpot+process.beamSpotMaker+process.MeasurementTrackerEvent+process.vertexMakerVtx+process.hitMakerVtx+process.hitMakerNoVtx)
else:
    if do_skim:
        process.skimpath = cms.Path(process.countmu+process.countvtx+process.gtStage2Digis+process.triggerMaker+process.offlineBeamSpot+process.beamSpotMaker+process.MeasurementTrackerEvent+process.hitMaker)
    else:
        process.skimpath = cms.Path(process.gtStage2Digis+process.patTrigger+process.triggerMaker+process.offlineBeamSpot+process.beamSpotMaker+process.MeasurementTrackerEvent+process.hitMaker)
