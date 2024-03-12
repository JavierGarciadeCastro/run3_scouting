#!/usr/bin/env python3
import os
import socket
import time

import dask
import glob
from distributed import Client, LocalCluster
from dask.distributed import progress
import uproot

from coffea.nanoevents.methods import vector
import awkward as ak
ak.behavior.update(vector.behavior)

import hist
import matplotlib.pyplot as plt
import mplhep as hep #matplotlib wrapper for easy plotting in HEP
plt.style.use(hep.style.CMS)

from utils import das_wrapper, redirectors, choose

def get_muons(events, branch="Run3ScoutingMuons_hltScoutingMuonPacker__HLT.obj."):
    from coffea.nanoevents.methods import vector
    import awkward as ak
    ak.behavior.update(vector.behavior)
    return ak.zip({
        'pt': events[f'{branch}pt_'],
        'eta': events[f'{branch}eta_'],
        'phi': events[f'{branch}phi_'],
        'mass': events[f'{branch}m_'],
        'charge': events[f'{branch}charge_'],
        'trackIso': events[f'{branch}trackIso_'],
        'ecalIso': events[f'{branch}ecalIso_'],
        'hcalIso': events[f'{branch}hcalIso_'],
        }, with_name="PtEtaPhiMLorentzVector")

def get_tracks(events, branch="Run3ScoutingTracks_hltScoutingTrackPacker__HLT.obj."):
    from coffea.nanoevents.methods import vector
    import awkward as ak
    ak.behavior.update(vector.behavior)
    return ak.zip({
        'pt': events[f'{branch}pt_'],
        'eta': events[f'{branch}eta_'],
        'phi': events[f'{branch}phi_'],
        'mass': events[f'{branch}m_'],
        'charge': events[f'{branch}charge_'],
        #'trackIso': events[f'{branch}trackIso'],
        }, with_name="PtEtaPhiMLorentzVector")

def get_pfcands(events, branch="Run3ScoutingTracks_hltScoutingTrackPacker__HLT.obj."):
    from coffea.nanoevents.methods import vector
    import awkward as ak
    ak.behavior.update(vector.behavior)
    return ak.zip({
        'pt': events[f'{branch}pt_'],
        'eta': events[f'{branch}eta_'],
        'phi': events[f'{branch}phi_'],
        'mass': events[f'{branch}m_'],
        'charge': events[f'{branch}charge_'],
        #'trackIso': events[f'{branch}trackIso'],
        }, with_name="PtEtaPhiMLorentzVector")

def get_vertices(events, branch="Run3ScoutingVertexs_hltScoutingMuonPacker_displacedVtx_HLT.obj."):
    from coffea.nanoevents.methods import vector
    import awkward as ak
    ak.behavior.update(vector.behavior)
    return ak.zip({
        'x': events[f'{branch}x_'],
        'y': events[f'{branch}y_'],
        'z': events[f'{branch}z_'],
        #'trackIso': events[f'{branch}trackIso'],
        }, with_name="ThreeVector")


test_files = [
    '/store/data/Run2022C/ScoutingPFRun3/RAW/v1/000/357/479/00000/e68468b2-1e5a-4f1c-8f68-1dac20389864.root',
    '/store/data/Run2022C/ScoutingPFRun3/RAW/v1/000/357/479/00000/1d967be3-0308-4000-a6d4-75997faff68c.root',
    '/store/data/Run2022C/ScoutingPFRun3/RAW/v1/000/357/479/00000/632bc65a-f8c7-40d2-a172-0fbe1de8d0e1.root',
    '/store/data/Run2022C/ScoutingPFRun3/RAW/v1/000/357/479/00000/4cc9d927-4d88-4d44-8c51-763ddea5f439.root',
    '/store/data/Run2022C/ScoutingPFRun3/RAW/v1/000/357/479/00000/c009bea1-ff3f-4316-a61d-ef9c69c5bfd9.root',
    '/store/data/Run2022C/ScoutingPFRun3/RAW/v1/000/357/479/00000/80b42ed3-9a22-419b-a36b-f833dc9dbdf1.root',
    '/store/data/Run2022C/ScoutingPFRun3/RAW/v1/000/357/482/00000/16a0f548-52e3-4638-ba20-9777d853bbcf.root',
    '/store/data/Run2022C/ScoutingPFRun3/RAW/v1/000/357/479/00000/6202e446-445d-448e-9b25-9dba5356bb5d.root',
    '/store/data/Run2022C/ScoutingPFRun3/RAW/v1/000/357/479/00000/abf8278d-eaf6-4c11-84b2-7fddf98ef999.root',
    '/store/data/Run2022C/ScoutingPFRun3/RAW/v1/000/357/479/00000/b905b6d6-f9ac-46ef-979e-620e5a153407.root',
    '/store/data/Run2022C/ScoutingPFRun3/RAW/v1/000/357/479/00000/57c46af4-c899-4b98-bd13-6f2dcd11bc1d.root',
    '/store/data/Run2022C/ScoutingPFRun3/RAW/v1/000/357/479/00000/48e8802d-53ce-445d-806c-84821b2fb479.root',
    '/store/data/Run2022C/ScoutingPFRun3/RAW/v1/000/357/479/00000/cd3932b2-bcfc-4dad-9be6-022daf9a6550.root',
    '/store/data/Run2022C/ScoutingPFRun3/RAW/v1/000/357/479/00000/710fff40-67bd-4a41-af0f-9e6e3ee47732.root',
    '/store/data/Run2022C/ScoutingPFRun3/RAW/v1/000/357/479/00000/c30d2157-fe2c-414f-b6a5-51a5925867ea.root',
]

test_files = [redirectors["fnal"] + x for x in test_files]

def get_nevents(f_in):
    with uproot.open(f_in) as f:
        events = f["Events"]
        return len(events["double_hltScoutingPFPacker_pfMetPt_HLT./double_hltScoutingPFPacker_pfMetPt_HLT.obj"].arrays())

def make_simple_hist(f_in):
    #pt_axis = hist.axis.Regular(20, 0.0, 200, name="pt", label=r"$p_{T}^{miss\ (GeV)$")
    results = {}
    mass_axis = hist.axis.Regular(500, 0.0, 50, name="mass", label=r"$M(\mu\mu)\ (GeV)$")
    N_axis = hist.axis.Regular(15, -0.5, 14.5, name="n", label=r"$N$")
    dataset_axis = hist.axis.StrCategory([], name="dataset", label="Dataset", growth=True)
    results['mass'] = hist.Hist(
        dataset_axis,
        #pt_axis,
        mass_axis,
    )
    results['nmuons'] = hist.Hist(
        dataset_axis,
        N_axis,
    )
    results['vertices'] = hist.Hist(
        dataset_axis,
        hist.axis.Regular(1000, -10, 10, name="x", label=r"$x (cm)$"),
        hist.axis.Regular(1000, -10, 10, name="y", label=r"$y (cm)$"),
        hist.axis.Regular(100, -5, 5, name="z", label=r"$z (cm)$"),
    )
    try:
        with uproot.open(f_in, timeout=300) as f:
            events = f["Events"]
            muons = events["Run3ScoutingMuons_hltScoutingMuonPacker__HLT./Run3ScoutingMuons_hltScoutingMuonPacker__HLT.obj"].arrays()
            vertices = events["Run3ScoutingVertexs_hltScoutingMuonPacker_displacedVtx_HLT./Run3ScoutingVertexs_hltScoutingMuonPacker_displacedVtx_HLT.obj"].arrays()
            muons4 = get_muons(muons)
            vert = get_vertices(vertices)
            dimuon = choose(muons4, 2)
            OS_dimuon   = dimuon[((dimuon['0'].charge*dimuon['1'].charge)<0)]
            results['mass'].fill(
                dataset=f_in,
                mass=ak.flatten(OS_dimuon.mass, axis=1),
                #pt=ak.flatten(OS_dimuon.pt, axis=1),
                #pt=events["double_hltScoutingPFPacker_pfMetPt_HLT./double_hltScoutingPFPacker_pfMetPt_HLT.obj"].arrays(),
            )
            results['vertices'].fill(
                dataset=f_in,
                x=ak.flatten(vert.x, axis=1),
                y=ak.flatten(vert.y, axis=1),
                z=ak.flatten(vert.z, axis=1),
            )
            results['nmuons'].fill(
                dataset=f_in,
                n=ak.num(muons4, axis=1),
            )
    except OSError:
        print(f"Could not open file {f_in}, skipping.")
        #raise
    except ValueError:
        print(f"Array missing")
    return results


    
if __name__ == '__main__':

    import argparse

    argParser = argparse.ArgumentParser(description = "Argument parser")
    argParser.add_argument('--cluster', action='store_true', default=False, help="Run on a cluster")
    argParser.add_argument('--small', action='store_true', default=False, help="Run on a small subset")
    argParser.add_argument('--workers', action='store', default=4,  help="Set the number of workers for the DASK cluster")
    argParser.add_argument('--input', action='store', default=None,  help="Provide input file")
    args = argParser.parse_args()

    print("Preparing")

    local = not args.cluster
    workers = int(args.workers)
    host = socket.gethostname()

    print("Starting cluster")
    if local:
        cluster = LocalCluster(
            n_workers=workers,
            threads_per_worker=1,
        )
    else:
        if host.count('ucsd'):
            raise NotImplementedError("Can't yet use a condor cluster on UAF. Please run locally.")
        from lpcjobqueue import LPCCondorCluster
        cluster = LPCCondorCluster(
            transfer_input_files="scouting",
        )
        dask.config.set({'distributed.scheduler.allowed-failures': '20'})

    cluster.adapt(minimum=0, maximum=workers)
    client = Client(cluster)

    if args.input:
        all_files = glob.glob("/ceph/cms/store/user/legianni/testRAWScouting_0/ScoutingPFRun3/crab_skim_2022D_0/230613_184336/0000/*.root")[:100]
        #all_files = [args.input]
        n_events = [350248]  # this is for my example /ceph/cms/store/user/legianni/testRAWScouting_0/ScoutingPFRun3/crab_skim_2022D_0/230613_184336/0000/output_91.root
    else:

        print("Querying files. This should get cached.")
        files = das_wrapper('/ScoutingPFRun3/Run2022B-v1/RAW', query='file', mask=' | grep file.name, file.nevents')  # can't filter only files on Caltech...
        # adding the redirector + weed out useless empty files
        nmax = 10 if args.small else int(1e7)
        all_files = [redirectors["fnal"] + x.split()[0] for x in files if int(x.split()[1])>10][:nmax]
        n_events = [int(x.split()[1]) for x in files if int(x.split()[1])>10][:nmax]

    print("Computing the results")
    tic = time.time()
    futures = client.map(make_simple_hist, all_files)  # .reduction does not work
    # NOTE: accumulation below is potentially memory intensive
    # This is WIP
    # https://docs.dask.org/en/stable/bag.html
    # reduce? https://stackoverflow.com/questions/70563132/distributed-chained-computing-with-dask-on-a-high-failure-rate-cluster
    progress(futures)
    print("Gathering")
    print(client.gather(futures))
    results = client.gather(futures)

    elapsed = time.time() - tic
    print(f"Finished in {elapsed:.1f}s")
    print(f"Total events {round(sum(n_events)/1e6,2)}M")
    print(f"Events/s: {sum(n_events) / elapsed:.0f}")

    print("Accumulating")
    # transpose the results list of dicts
    results_t = {k: [dic[k] for dic in results] for k in results[0]}
    total_hist = sum(results_t['mass'])
    #total_hist[{"dataset":sum, "pt":sum}].show(columns=100)
    total_hist[{"dataset":sum}].show(columns=100)

    plot_dir = os.path.expandvars('./plots/')
    fig, ax = plt.subplots(figsize=(8, 8))

    total_hist[{"dataset":sum}].plot1d(
        histtype="step",
        ax=ax,
    )

    hep.cms.label(
        "Preliminary",
        data=True,
        lumi='0.086',
        com=13.6,
        loc=0,
        ax=ax,
        fontsize=15,
    )
    fig.savefig(f'{plot_dir}/dimuon_mass.png')

    total_hist = sum(results_t['nmuons'])
    vert_hist = sum(results_t['vertices'])
    #total_hist[{"dataset":sum, "pt":sum}].show(columns=100)
    total_hist[{"dataset":sum}].show(columns=50)


    fig, ax = plt.subplots(figsize=(8, 8))
    total_hist[{"dataset":sum}].plot1d(
        histtype="step",
        ax=ax,
    )

    hep.cms.label(
        "Preliminary",
        data=True,
        lumi='0.086',
        com=13.6,
        loc=0,
        ax=ax,
        fontsize=15,
    )
    ax.set_yscale('log')
    fig.savefig(f'{plot_dir}/nmuons.png')

    from matplotlib.colors import LogNorm
    fig, ax = plt.subplots(figsize=(8, 8))
    vert_hist[:, -10.0j:10.0j:2j, -10.0j:10.0j:2j, :][{"dataset":sum, 'z':sum}].plot2d(
        #histtype="step",
        ax=ax,
        norm=LogNorm(),
    )

    hep.cms.label(
        "Preliminary",
        data=True,
        lumi='0.086',
        com=13.6,
        loc=0,
        ax=ax,
        fontsize=15,
    )
    #ax.set_zscale('log')
    fig.savefig(f'{plot_dir}/vertices.png')
