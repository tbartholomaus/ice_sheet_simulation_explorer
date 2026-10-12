"""
Profile the Dash app's memory and CPU usage.

Launches app.py (or any command given after --) as a child process, samples
the whole process tree (the app plus any children it spawns, e.g. gunicorn
workers) every --interval seconds, and appends one line per sample to
profile_log.txt. When the app exits -- on its own, or via Ctrl+C / SIGTERM --
a summary with the max, mean and median memory and CPU is written to the log and
printed.

Usage:
    python profile_app.py                        # profiles `python app.py`
    python profile_app.py --interval 0.5
    python profile_app.py -- gunicorn app:server --timeout 120

Memory is RSS (resident set size) summed over the process tree. CPU % is
like `top`: 100% = one fully busy core, so it can exceed 100% on multi-core
machines.

Requires psutil (pip install psutil).
"""

import argparse
import os
import signal
import statistics
import subprocess
import sys
import time
from datetime import datetime

import psutil

HERE = os.path.dirname(os.path.abspath(__file__))


def _tree(root, cache):
    """Return the live processes in root's tree, reusing cached Process
    objects so cpu_percent() measures since the previous sample."""
    try:
        procs = [root] + root.children(recursive=True)
    except psutil.NoSuchProcess:
        return []
    live = []
    for p in procs:
        live.append(cache.setdefault(p.pid, p))
    for pid in set(cache) - {p.pid for p in live}:
        del cache[pid]
    return live


def _sample(procs):
    rss = cpu = 0.0
    for p in procs:
        try:
            with p.oneshot():
                rss += p.memory_info().rss
                cpu += p.cpu_percent(None)
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            pass
    return rss / 2**20, cpu


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--interval", type=float, default=1.0, help="seconds between samples (default 1)")
    ap.add_argument("--log", default=os.path.join(HERE, "profile_log.txt"), help="log file path")
    ap.add_argument("cmd", nargs=argparse.REMAINDER, help="command to profile (after --)")
    args = ap.parse_args()
    cmd = args.cmd[1:] if args.cmd[:1] == ["--"] else args.cmd
    cmd = cmd or [sys.executable, "app.py"]

    child = subprocess.Popen(cmd, cwd=HERE)
    root = psutil.Process(child.pid)
    cache = {}
    _sample(_tree(root, cache))  # prime cpu_percent baselines

    stop = False

    def _on_signal(signum, _frame):
        nonlocal stop
        stop = True
        try:
            child.send_signal(signum)
        except ProcessLookupError:
            pass

    signal.signal(signal.SIGINT, _on_signal)
    signal.signal(signal.SIGTERM, _on_signal)

    start = time.time()
    peak_mem = peak_cpu = 0.0
    peak_mem_t = peak_cpu_t = start
    mems, cpus = [], []

    with open(args.log, "a", buffering=1) as log:
        log.write(f"\n=== profiling started {datetime.now():%Y-%m-%d %H:%M:%S} | "
                  f"cmd: {' '.join(cmd)} | pid {child.pid} | interval {args.interval}s | "
                  f"{psutil.cpu_count()} CPUs ===\n")
        log.write(f"{'timestamp':<19}  {'elapsed_s':>9}  {'rss_mb':>9}  {'cpu_pct':>7}  {'nprocs':>6}\n")

        while child.poll() is None:
            time.sleep(args.interval)
            procs = _tree(root, cache)
            if not procs or child.poll() is not None:
                break
            mem, cpu = _sample(procs)
            now = time.time()
            mems.append(mem)
            cpus.append(cpu)
            if mem > peak_mem:
                peak_mem, peak_mem_t = mem, now
            if cpu > peak_cpu:
                peak_cpu, peak_cpu_t = cpu, now
            log.write(f"{datetime.fromtimestamp(now):%Y-%m-%d %H:%M:%S}  {now - start:9.1f}  "
                      f"{mem:9.1f}  {cpu:7.1f}  {len(procs):6d}\n")

        # After a Ctrl+C, give the app a moment to shut down before reporting.
        try:
            child.wait(timeout=10)
        except subprocess.TimeoutExpired:
            child.kill()
            child.wait()

        fmt = lambda t: f"{datetime.fromtimestamp(t):%H:%M:%S}, {t - start:.1f}s in"
        avg = lambda xs: statistics.mean(xs) if xs else 0.0
        med = lambda xs: statistics.median(xs) if xs else 0.0
        summary = (
            f"--- summary ({'interrupted' if stop else 'app exited'}, exit code {child.returncode}) ---\n"
            f"duration:   {time.time() - start:.1f} s, {len(mems)} samples\n"
            f"memory (MB RSS): max {peak_mem:.1f} ({fmt(peak_mem_t)}), "
            f"mean {avg(mems):.1f}, median {med(mems):.1f}\n"
            f"CPU (%, 100% = one core): max {peak_cpu:.1f} ({fmt(peak_cpu_t)}), "
            f"mean {avg(cpus):.1f}, median {med(cpus):.1f}\n"
        )
        log.write(summary)

    print("\n" + summary + f"log: {args.log}")


if __name__ == "__main__":
    main()
