#!/usr/bin/env python3
"""Compare a benchmark run with the shipped reference (reference/<case>.json and reference/<case>/*.npy).

    python tests/check_reference.py --case case1 --dir OUT_case1_py
    python tests/check_reference.py --case case1 --dir RUN_case1/OUTPUT_SLM     (Fortran example output)
    python tests/check_reference.py --case case3 --dir OUT_case3_py --topo-pres INPUT_FILES/others/topo_pres
    python tests/check_reference.py --case case1 --run --input /path/to/INPUT_FILES --workdir RUNDIR

--run first runs the reference configuration (Python, ducc0, --opt all --inner extrap, --threads,
default 8) in --workdir and then compares. Tolerance --tol (default 1e-6 m) applies to the max |diff|
of every full field, to the area-weighted RMS and max |value| summaries and to the site values;
inner-iteration counts must be identical (use --iter-slack N to allow differences of up to N per step,
e.g. on another machine or BLAS). Prints one line per quantity and exits 1 if any differs.
Reference fields whose source files are not in --dir are reported as SKIP and listed at the end (the
Fortran Case 1 example writes only tgrid550, so dtopo_550 and the inner-iteration counts from the text
file OUTPUT_SLM/numiter are compared and dtopo_1, dtopo_100 and the dsl_* fields are skipped).
"""
import argparse
import json
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
sys.path.insert(0, os.path.join(ROOT, "python"))
from slevolve.tools import make_reference as mr  # noqa: E402

LOVE = {"case1": "prem_512.l60K2C.sum18p6.dum19p2.tz19p4.lm22", "case2": "prem_512.l120C.ump5.lm5",
        "case3": "prem_512.l96C.ump3.lm1"}


def run(case, inp, work, threads):
    py, pd = sys.executable, os.path.join(ROOT, "python")
    earth = os.path.join(ROOT, "lovenum", LOVE[case])
    out = os.path.join(work, f"OUT_{case}_py")
    common = ["--inner", "extrap", "--opt", "all", "--lmax", "511", "--threads", str(threads)]
    if case == "case1":
        cmd = [[py, f"{pd}/run_ais.py", "--steps", "550", "--sht", "ducc", *common, "--icedir",
                f"{inp}/icemodel", "--topo", f"{inp}/others/topo0_archive", "--earth", earth,
                "--save-fields", "tgrid,ocean", "--out", out]]
    elif case == "case2":
        cmd = [[py, f"{pd}/run_outer.py", "--iters", "1", *common, "--inputs", inp, "--icepfx", "iceload",
                "--nfiles", "1201", "--dt", "200", "--topo-pres", f"{inp}/others/etopo2_512_ice6gC",
                "--earth", earth, "--save-steps", "1,600,1200", "--out", out]]
    else:
        t0 = os.path.join(work, "T0_segmented.npy")
        cmd = [[py, f"{pd}/init_topo_segmented.py", "--inputs", inp, "--icepfx", "nadi_", "--nfiles", "51",
                "--dt", "500", "--earth", earth, "--lmax", "511", "--workers", str(threads), "--out", t0],
               [py, f"{pd}/run_outer.py", "--iters", "5", "--reuse", "warm", *common, "--inputs", inp,
                "--earth", earth, "--topo0", t0, "--save-steps", "50", "--out", out]]
    os.makedirs(work, exist_ok=True)
    for c in cmd:
        print("+", " ".join(c), flush=True)
        subprocess.run(c, check=True)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", required=True, choices=("case1", "case2", "case3"))
    ap.add_argument("--dir", help="output directory of an existing run")
    ap.add_argument("--run", action="store_true", help="run the reference configuration first")
    ap.add_argument("--input", help="INPUT_FILES directory (with --run; case3 also default for --topo-pres)")
    ap.add_argument("--workdir", default=".", help="where --run writes its output")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--topo-pres", default=None)
    ap.add_argument("--refdir", default=os.path.join(ROOT, "reference"))
    ap.add_argument("--tol", type=float, default=1e-6, help="m (default 1e-6)")
    ap.add_argument("--iter-slack", type=int, default=0)
    a = ap.parse_args()
    d = run(a.case, a.input, a.workdir, a.threads) if a.run else a.dir
    if d is None:
        ap.error("give --dir or --run --input")
    tp = a.topo_pres or (os.path.join(a.input, "others", "topo_pres") if a.input else None)
    ref = json.load(open(os.path.join(a.refdir, a.case + ".json")))
    F, meta = mr.fields(a.case, d, tp)
    new = mr.summary(a.case, F, meta)
    _, _, w = mr.grid()
    bad = []

    def report(name, diff, tol):
        ok = diff <= tol
        print(f"{'ok  ' if ok else 'FAIL'} {name:45s} {diff:.3e} (tol {tol:g})")
        if not ok:
            bad.append(name)

    skipped = []
    for name in ref["fields"]:
        if name not in F:
            print(f"SKIP {name:45s} not in {d}")
            skipped.append(name)
            continue
        f = F[name][0]
        r = np.load(os.path.join(a.refdir, a.case, name + ".npy"))
        dd = f - r
        report(f"{name} field max|diff|", float(np.abs(dd).max()), a.tol)
        report(f"{name} field area-weighted RMS diff",
               float(np.sqrt((w[:, None] * dd * dd).sum() / (w.sum() * dd.shape[1]))), a.tol)
        for q in ("rms", "maxabs"):
            report(f"{name} {q}", abs(new["fields"][name][q] - ref["fields"][name][q]), a.tol)
        sd = [abs(x - y) for x, y in zip(new["fields"][name]["sites"], ref["fields"][name]["sites"])]
        k = int(np.argmax(sd))
        report(f"{name} sites (worst: {ref['sites'][k]['name']})", sd[k], a.tol)
    if len(skipped) == len(ref["fields"]):
        print(f"FAIL no reference field found in {d}")
        bad.append("fields")
    rn, nn = np.array(ref["numiter"]).ravel(), np.array(new.get("numiter", [])).ravel()
    if nn.size == 0:
        print("SKIP numiter                                       no numiter.npy / numiter file")
        skipped.append("numiter")
    elif rn.shape != nn.shape:
        print(f"FAIL numiter: {nn.size} steps vs {rn.size} in the reference")
        bad.append("numiter")
    else:
        dif = np.abs(rn - nn)
        report(f"numiter per step (total {nn.sum()} vs {rn.sum()}, {int((dif > 0).sum())} steps differ)",
               float(dif.max()), a.iter_slack)
    for key in ("misfit_rms", "misfit_max"):
        if key in ref:
            report(key, float(np.max(np.abs(np.array(new[key]) - np.array(ref[key])))), a.tol)
    if skipped:
        print(f"skipped (not in the run output): {', '.join(skipped)}")
    print(f"check_reference {a.case}: " + ("PASS" if not bad else f"FAIL ({len(bad)}: {', '.join(bad)})")
          + (f", {len(skipped)} item(s) skipped" if skipped else ""))
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
