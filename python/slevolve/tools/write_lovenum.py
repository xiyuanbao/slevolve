#!/usr/bin/env python3
"""write_lovenum.py -- write Love numbers in the text layout read by SLEvolve.

SLEvolve does not compute Love numbers. It reads them from a text file (the reader is
`LoveNumbers` in python/slevolve/solver.py and the block "Read in Love numbers" in
fortran/src/sl_model_mod.f90). This tool writes that file from arrays, so that
normal-mode or collocation results from any Earth-model code can be used.

Layout (every group below starts on a new line; a group of n values may span several lines):

    line 1                 header, free text (skipped by the readers)
    then, for l = 1, 2, ..., L (consecutive, starting at 1; degree 0 is not in the file):
      l  K_l               degree and number of modes K_l (may differ between degrees; K_l <= 500)
      s_1 .. s_K           inverse relaxation times of the modes, 1/kyr (positive)
      h^E  l*l^E  0  0  l*k^E          loading, elastic
      h^F  l*l^F  0  0  l*k^F          loading, fluid limit (t -> infinity); read, not used
      r_1 .. r_K           loading h residues, 1/kyr
      l*r'_1 .. l*r'_K     loading l residues times l, 1/kyr
      l*q_1 .. l*q_K       loading k residues times l, 1/kyr
      and, only for l >= 2, the same five groups for the tidal Love numbers:
      h^TE  l*l^TE  0  0  l*k^TE
      h^TF  l*l^TF  0  0  l*k^TF
      tidal h, l*(tidal l), l*(tidal k) residues
    The two middle columns of the elastic and fluid lines are read and discarded (zero in the
    shipped files).

Time-domain response to a Heaviside load at degree l:
    h_l(t) = h^E_l + sum_k (r_k / s_k) (1 - exp(-s_k t)),  t in kyr,
and the same for l and k. With t in years the code uses exp(-s_k dt / 1000).
k and l are multiplied by the degree l in the file (h is not); the readers divide by l.
The solver reads degrees 1..lmax, so the file must hold at least lmax degrees (more are ignored).
Degree 1 must be in the frame you want the degree-1 outputs in; the shipped files use the CE
frame, k^E_1 = 0 (lovenum/DEGREE1_NOTES.md).

Python use (arrays indexed by degree, entry 0 unused; per-degree mode lists):

    import numpy as np, sys
    from slevolve.tools.write_lovenum import write_lovenum
    L = 2
    s  = [None, np.array([0.5, 2.0]), np.array([0.5, 2.0])]          # 1/kyr
    zr = [None, np.zeros(2), np.zeros(2)]
    write_lovenum("my_earth.txt", "my model, CE frame",
                  hE=np.array([0, -0.29, -0.99]), lE=np.array([0, 0.10, 0.024]),
                  kE=np.array([0, 0.0, -0.31]),
                  hF=np.array([0, -0.69, -2.8]), lF=np.array([0, 0.31, -0.26]),
                  kF=np.array([0, 0.0, -0.99]),
                  s=s, rh=zr, rl=zr, rk=zr,
                  hTE=np.array([0, 0, 0.60]), lTE=np.array([0, 0, 0.084]),
                  kTE=np.array([0, 0, 0.30]), hTF=np.array([0, 0, 1.9]),
                  lTF=np.array([0, 0, 0.45]), kTF=np.array([0, 0, 0.92]),
                  rhT=zr, rlT=zr, rkT=zr)

(The numbers above only illustrate the call; they are not an Earth model.) k, l and residues
are passed unscaled; the tool multiplies by l. Tidal arguments may be omitted (then zeros
are written); they are used only with rotation (TPW) on, at degree 2.

Command line: rewrite an existing file through the reader (round-trip check):

    python python/tools/write_lovenum.py --from lovenum/prem_512.l96C.ump3.lm1 --lmax 512 --out copy.txt
"""
import argparse
import os
import sys

import numpy as np


def _fmt(v, digits):
    return " % .*E" % (digits, v)


def _group(fh, vals, digits, per_line=5):
    vals = list(vals)
    if not vals:
        fh.write("\n")
        return
    for i in range(0, len(vals), per_line):
        fh.write("".join(_fmt(float(x), digits) for x in vals[i:i + per_line]) + "\n")


def write_lovenum(path, header, hE, lE, kE, hF, lF, kF, s, rh, rl, rk,
                  hTE=None, lTE=None, kTE=None, hTF=None, lTF=None, kTF=None,
                  rhT=None, rlT=None, rkT=None, lmax=None, digits=15):
    """Write the Love-number text file.

    hE, lE, kE, hF, lF, kF (and the tidal hTE ... kTF): arrays indexed by degree (0..L, entry 0
    unused), unscaled. s, rh, rl, rk (and rhT, rlT, rkT): sequences indexed by degree whose
    entry l is an array of K_l poles (1/kyr) or residues (1/kyr), unscaled. lmax defaults to
    len(hE) - 1. digits: significant digits after the point (default 15, i.e. 16 in total).
    """
    L = int(lmax if lmax is not None else len(hE) - 1)
    z = np.zeros(L + 1)
    hTE, lTE, kTE, hTF, lTF, kTF = (z if a is None else a for a in (hTE, lTE, kTE, hTF, lTF, kTF))
    with open(path, "w") as fh:
        fh.write(" " + header.replace("\n", " ") + "\n")
        for l in range(1, L + 1):
            sl = np.asarray(s[l], float)
            K = len(sl)
            if K > 500:
                raise ValueError(f"degree {l}: {K} modes, the readers allow at most 500")
            fh.write("%10d%6d\n" % (l, K))
            _group(fh, sl, digits)
            _group(fh, [hE[l], l * lE[l], 0.0, 0.0, l * kE[l]], digits)
            _group(fh, [hF[l], l * lF[l], 0.0, 0.0, l * kF[l]], digits)
            _group(fh, np.asarray(rh[l], float), digits)
            _group(fh, l * np.asarray(rl[l], float), digits)
            _group(fh, l * np.asarray(rk[l], float), digits)
            if l < 2:
                continue
            zK = np.zeros(K)
            _group(fh, [hTE[l], l * lTE[l], 0.0, 0.0, l * kTE[l]], digits)
            _group(fh, [hTF[l], l * lTF[l], 0.0, 0.0, l * kTF[l]], digits)
            _group(fh, zK if rhT is None else np.asarray(rhT[l], float), digits)
            _group(fh, zK if rlT is None else l * np.asarray(rlT[l], float), digits)
            _group(fh, zK if rkT is None else l * np.asarray(rkT[l], float), digits)


def write_from_lovenumbers(path, love, header, digits=15):
    """Write a LoveNumbers object (python/slevolve/solver.py) back to the text layout."""
    L = love.norder
    n = love.nmod
    sl = lambda a: [None] + [a[: n[l], l] for l in range(1, L + 1)]  # noqa: E731
    write_lovenum(path, header, love.he, love.le, love.ke, love.h_fluid, love.l_fluid, love.k_fluid,
                  sl(love.s), sl(love.r), sl(love.rl), sl(love.rprime),
                  love.hTE, love.lTE, love.kTE, love.hT_fluid, love.lT_fluid, love.kT_fluid,
                  sl(love.rT), sl(love.rlT), sl(love.rprimeT), lmax=L, digits=digits)


def main():
    ap = argparse.ArgumentParser(description="rewrite a Love-number file through the SLEvolve reader")
    ap.add_argument("--from", dest="src", required=True, help="existing Love-number file")
    ap.add_argument("--lmax", type=int, required=True, help="number of degrees to read and write")
    ap.add_argument("--out", required=True)
    ap.add_argument("--digits", type=int, default=15)
    a = ap.parse_args()
    from ..solver import LoveNumbers

    with open(a.src) as fh:
        header = fh.readline().strip()
    write_from_lovenumbers(a.out, LoveNumbers(a.src, a.lmax), header, a.digits)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
