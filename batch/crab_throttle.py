#!/usr/bin/env python3
import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

TERMINAL = {"COMPLETED", "FAILED", "KILLED", "SUBMITFAILED", "RESUBMITFAILED", "KILLFAILED"}


DATASETS = {
    "DQCD_sig": "data/datasets_dqcd_2024_signal.txt",
    "DQCD_sig_private": "data/datasets_dqcd_2024_signal_private.txt",
    "HTo2ZdTo2mu2x": "data/datasets_hahm_2024_signal.txt",
}


def task_dir(ds, era, tag):
    if "HTo2Zd" in ds:
        m = ds.split("MZd-")[1].split("-Ctau")[0]
        ctau = ds.split("-Ctau-")[1].split("mm_")[0]
        return "crab_centralSkim__HTo2ZdTo2mu2x_%s_m-%s_ctau-%smm_%s" % (era, m, ctau, tag)
    model = ds.split("-")[1]
    ctau = ds.split("ctau-")[1].split("-")[0]
    mA = ds.split("mA-")[1].split("-")[0]
    mpi = re.split("[_/]", ds.split("mpi-")[1])[0]
    return "crab_centralSkim__%s_%s_mpi-%s_mA-%s_ctau-%smm%s" % (model, era, mpi, mA, ctau, tag)


def status(d):
    r = subprocess.run(["crab", "status", "-d", str(d)], capture_output=True, text=True)
    out = r.stdout + r.stderr
    m = re.search(r"Status on the scheduler:\s+(\S+)", out)
    if m:
        return m.group(1)
    m = re.search(r"Status on the CRAB server:\s+(\S+)", out)
    return m.group(1) if m else "UNKNOWN"


def read_lines(p):
    return [l.strip() for l in open(p) if l.strip() and not l.startswith("#")]


def main():
    p = argparse.ArgumentParser(description="Submit the 2024 signal CRAB tasks keeping at most --max in flight.")
    p.add_argument("model", choices=list(DATASETS))
    p.add_argument("--datasets", default=None, help="dataset list (default: the one of the model in data/)")
    p.add_argument("--match", default="", help="only submit datasets whose name contains this text, e.g. ScenarioA")
    p.add_argument("--only-dirs", default=None,
                   help="only submit datasets whose task dir name is listed in this file (e.g. B1_kill.txt)")
    p.add_argument("--inflight", nargs="*", default=[],
                   help="files listing task dirs already in flight that count toward --max (e.g. B1_keep.txt)")
    p.add_argument("--also", nargs="*", default=[], help="extra task dirs that count toward --max")
    p.add_argument("--max", type=int, default=10)
    p.add_argument("--sleep", type=int, default=3600, help="seconds between status rounds")
    p.add_argument("--era", default="2024")
    p.add_argument("--tag", default="_dqcd_2024_L1seeds", help="ntuple_version suffix of the request names")
    p.add_argument("--cfg", default="crabcfg_run3_centralmc.py")
    p.add_argument("--state", default=None, help="state file (default: crab_throttle_<model><tag>.json)")
    p.add_argument("--killed-dir", default="killed")
    args = p.parse_args()

    work = Path(os.environ["CMSSW_BASE"]).resolve().parent
    datasets = [d for d in read_lines(args.datasets or DATASETS[args.model]) if args.match in d]
    if args.only_dirs:
        wanted = {Path(x).name for x in read_lines(args.only_dirs)}
        datasets = [d for d in datasets if task_dir(d, args.era, args.tag) in wanted]

    state_path = Path(args.state or "crab_throttle_%s%s.json" % (args.model, args.tag))
    if state_path.is_file():
        state = json.load(open(state_path))
    else:
        inflight = [Path(x).name for f in args.inflight for x in read_lines(f)] + [Path(x).name for x in args.also]
        state = {"inflight": inflight, "submitted": [], "finished": {}}
    queue = [d for d in datasets if d not in state["submitted"]]
    print(f"work area {work}: {len(queue)} datasets queued, {len(state['inflight'])} tasks in flight, max {args.max}")

    try:
        while True:
            still = []
            for d in state["inflight"]:
                s = status(work / d)
                print(f"  {d}: {s}")
                if s in TERMINAL:
                    state["finished"][d] = s
                else:
                    still.append(d)
            state["inflight"] = still
            free = args.max - len(still)
            if free > 0 and queue:
                batch = queue[:free]
                for ds in batch:
                    old = work / task_dir(ds, args.era, args.tag)
                    if old.is_dir():
                        (work / args.killed_dir).mkdir(exist_ok=True)
                        dst = work / args.killed_dir / old.name
                        if dst.exists():
                            dst = dst.with_name(dst.name + "_%d" % int(time.time()))
                        shutil.move(str(old), str(dst))
                with tempfile.NamedTemporaryFile("w", suffix=".txt") as tmp:
                    tmp.write("".join(ds + "\n" for ds in batch))
                    tmp.flush()
                    subprocess.run([sys.executable, args.cfg, args.era, args.model],
                                   env=dict(os.environ, CRAB_DATASETS=tmp.name))
                for ds in batch:
                    d = task_dir(ds, args.era, args.tag)
                    if (work / d / ".requestcache").is_file():
                        state["inflight"].append(d)
                        state["submitted"].append(ds)
                        queue.remove(ds)
                    else:
                        print(f"  submit failed, will retry next round: {d}")
            json.dump(state, open(state_path, "w"), indent=1)
            print(time.strftime("%F %T"), f"in flight {len(state['inflight'])}, queued {len(queue)}, "
                  f"finished {len(state['finished'])}", flush=True)
            if not queue and not state["inflight"]:
                break
            subprocess.run(["kinit", "-R"], capture_output=True)
            left = subprocess.run(["voms-proxy-info", "-timeleft"], capture_output=True, text=True).stdout.strip()
            if left.isdigit() and int(left) < 2 * args.sleep:
                print(f"WARNING: proxy has only {left} s left: voms-proxy-init -voms cms -valid 192:00", flush=True)
            time.sleep(args.sleep)
    except KeyboardInterrupt:
        sys.exit("stopped; run the same command to continue")
    print("all done:", json.dumps(state["finished"], indent=1))


if __name__ == "__main__":
    main()
