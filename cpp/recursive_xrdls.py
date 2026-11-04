#!/usr/bin/env python3
import subprocess
import sys
import os

HOST = "redir.t2.ucsd.edu:1095"
SKIP_DIRS = {'log', 'failed', 'temp'}

def xrdls(path, retries=5, timeout=600):
    for attempt in range(retries):
        print(f"[xrdls] listing: {path}", file=sys.stderr, flush=True)
        try:
            result = subprocess.run(
                ["xrdfs", HOST, "ls", path],
                capture_output=True, text=True, timeout=timeout
            )
        except subprocess.TimeoutExpired:
            print(f"[xrdls] TIMEOUT on: {path} (attempt {attempt+1}/{retries})", file=sys.stderr, flush=True)
            continue
        if result.returncode != 0:
            print(f"[xrdls] ERROR on: {path} — {result.stderr.strip()}", file=sys.stderr, flush=True)
            return []
        entries = [l.strip() for l in result.stdout.strip().split('\n') if l.strip()]
        print(f"[xrdls] got {len(entries)} entries", file=sys.stderr, flush=True)
        return entries
    return []

def find_roots(path, depth=0):
    if depth > 6:
        return
    for entry in xrdls(path):
        if entry.endswith('.root'):
            print(entry)
        elif entry and os.path.basename(entry) not in SKIP_DIRS:
            find_roots(entry, depth + 1)

if __name__ == "__main__":
    find_roots(sys.argv[1])
