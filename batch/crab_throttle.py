#!/usr/bin/env python3
import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

TERMINAL = {"COMPLETED", "FAILED", "KILLED", "SUBMITFAILED", "RESUBMITFAILED", "KILLFAILED"}


def task_dir(ds, era, tag):
    model = ds.split("-")[1]
    ctau = ds.split("ctau-")[1].split("-")[0]
    mA = ds.split("mA-")[1].split("-")[0]
    mpi = ds.split("mpi-")[1].split("/")[0]
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
    p = argparse.ArgumentParser(description="Submit DQCD_sig_private CRAB tasks keeping at most --max in flight.")
    p.add_argument("datasets", help="dataset list (one DAS name per line)")
    p.add_argument("--only-dirs", default=None,
                   help="only submit datasets whose task dir name is listed in this file (e.g. B1_kill.txt)")
    p.add_argument("--inflight", nargs="*", default=[],
                   help="files listing task dirs already in flight that count toward --max (e.g. B1_keep.txt)")
    p.add_argument("--also", nargs="*", default=[], help="extra task dirs that count toward --max")
    p.add_argument("--max", type=int, default=10)
    p.add_argument("--sleep", type=int, default=3600, help="seconds between status rounds")
    p.add_argument("--era", default="2024")
    p.add_argument("--tag", default="_dqcd_2024", help="ntuple_version suffix of the request names")
    p.add_argument("--cfg", default="crabcfg_run3_centralmc.py")
    p.add_argument("--data-file", default="data/datasets_dqcd_2024_signal_private.txt")
    p.add_argument("--state", default="crab_throttle_state.json")
    p.add_argument("--killed-dir", default="killed")
    args = p.parse_args()

    work = Path(os.environ["CMSSW_BASE"]).resolve().parent
    data_file = Path(args.data_file)
    backup = data_file.with_name(data_file.name + ".throttle_bak")
    if data_file.is_file() and not backup.is_file():
        shutil.copy(data_file, backup)

    datasets = read_lines(args.datasets)
    if args.only_dirs:
        wanted = {Path(x).name for x in read_lines(args.only_dirs)}
        datasets = [d for d in datasets if task_dir(d, args.era, args.tag) in wanted]

    state_path = Path(args.state)
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
                data_file.write_text("".join(ds + "\n" for ds in batch))
                subprocess.run([sys.executable, args.cfg, args.era, "DQCD_sig_private"])
                for ds in batch:
                    d = task_dir(ds, args.era, args.tag)
                    if (work / d).is_dir():
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
    finally:
        if backup.is_file():
            shutil.copy(backup, data_file)
    print("all done:", json.dumps(state["finished"], indent=1))


if __name__ == "__main__":
    main()
