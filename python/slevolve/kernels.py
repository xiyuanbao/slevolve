"""
kernels.py -- numba-compiled fused kernels for the hot passes of ``solver.py``.

Every kernel here reproduces the *exact* arithmetic of the NumPy statements it
replaces, operation by operation and (for the reductions) accumulation order by
accumulation order, so the results are **identical** to the NumPy path.
``fastmath`` is deliberately **off** everywhere: with it LLVM would contract
``a*b + c`` into an FMA and results would change at round-off level.

The point of the kernels is not flops -- these passes are all memory bound --
but *fusion*: the NumPy versions write and re-read a 2 MB (grid) or 4 MB
(spectral) temporary between every binary operation, and the hereditary advance
streams 109 MB of state through the cache twelve times per relaxation mode.

Each kernel was checked against the NumPy expression it replaces on random
data.
"""

from __future__ import annotations

import numpy as np
from numba import njit, prange

_JIT = dict(cache=True, fastmath=False, nogil=True)


# ===================================================== hereditary advance ===
@njit(parallel=True, **_JIT)
def advance_and_read(R, decay, omd, a_sl, a_rr, nmod, P, visc, rad):
    """R[l,k,m] <- d R + (1-d) P ; visc = sum_k a_sl R ; rad = sum_k a_rr R.

    Replaces, with identical results, the k-loop of ``ExactMemory.advance_and_read``

        Rk *= decay[k][:, None]
        Rk += P * omd[k][:, None]
        visc += Rk * a_sl[k][:, None]
        rad  += Rk * a_rr[k][:, None]

    The degree axis ``l`` is the parallel axis (each ``l`` writes a disjoint
    slice); the accumulation over the relaxation modes ``k`` keeps its original
    order inside every row, hence identical results.  Modes ``k >= nmod[l]`` do not
    exist for that degree: ``a_sl = a_rr = 0`` there, so they contribute
    exactly 0.0 and are skipped.

    ``R`` is stored as ``[l, k, m]`` (not ``[k, l, m]``) so that the 26 mode
    rows of one degree are 208 kB of *contiguous* memory, and the innermost
    loop runs over ``m`` so that it is unit-stride and vectorisable while the
    two accumulator rows ``visc[l, :]`` and ``rad[l, :]`` (8 kB each) stay in
    L1.  Putting ``k`` innermost instead costs 3x (6.4 ms vs 2.1 ms at 4
    threads) for exactly the same arithmetic.
    """
    nl = R.shape[0]
    for l in prange(nl):
        nk = nmod[l]
        for m in range(nl):
            visc[l, m] = 0.0
            rad[l, m] = 0.0
        for k in range(nk):
            d = decay[k, l]
            o = omd[k, l]
            s = a_sl[k, l]
            q = a_rr[k, l]
            for m in range(l + 1):
                x = R[l, k, m] * d + P[l, m] * o
                R[l, k, m] = x
                visc[l, m] += x * s
                rad[l, m] += x * q


# ============================================================ grid kernels ===
@njit(parallel=True, **_JIT)
def icestar(ice, marine_pos, thr, out):
    """``out = ice * (marine_pos | (ice > thr))`` (checkmarine, eq. 44)."""
    n0, n1 = ice.shape
    for i in prange(n0):
        for j in range(n1):
            v = ice[i, j]
            out[i, j] = v if (marine_pos[i, j] or v > thr[i, j]) else 0.0


@njit(parallel=True, **_JIT)
def beta_cstar_tO(ice, cxy, cstar0, tinit, eps0, beta, cstarxy, tOxy):
    """beta = (ice < eps0), cstarxy = cxy * beta, tOxy = tinit * (cstarxy - cstar0) in one pass
    (``--checkmarine F``; with checkmarine beta is the grounded-ice mask, set in solver_opt)."""
    n0, n1 = ice.shape
    for i in prange(n0):
        for j in range(n1):
            b = 1.0 if ice[i, j] < eps0 else 0.0
            beta[i, j] = b
            c = cxy[i, j] * b
            cstarxy[i, j] = c
            tOxy[i, j] = tinit[i, j] * (c - cstar0[i, j])


@njit(parallel=True, **_JIT)
def add_const(src, c, out):
    """``out = src + c`` for a scalar ``c``."""
    n0, n1 = src.shape
    for i in prange(n0):
        for j in range(n1):
            out[i, j] = src[i, j] + c


@njit(parallel=True, **_JIT)
def rO_only(dslxy, cstarxy, rOxy):
    """``rOxy = dslxy * cstarxy``."""
    n0, n1 = dslxy.shape
    for i in prange(n0):
        for j in range(n1):
            rOxy[i, j] = dslxy[i, j] * cstarxy[i, j]


@njit(parallel=True, **_JIT)
def topo_ocean_cstar_tO(tinit, deltaslxy, beta, cstar0, topoxy, cxy, cstarxy, tOxy):
    """The inner-loop topography / ocean-function refresh, fused.

        topoxy  = tinit - deltaslxy
        cxy     = (topoxy < 0)
        cstarxy = cxy * beta
        tOxy    = tinit * (cstarxy - cstar0)
    """
    n0, n1 = tinit.shape
    for i in prange(n0):
        for j in range(n1):
            t = tinit[i, j] - deltaslxy[i, j]
            topoxy[i, j] = t
            o = 1.0 if t < 0.0 else 0.0
            cxy[i, j] = o
            c = o * beta[i, j]
            cstarxy[i, j] = c
            tOxy[i, j] = tinit[i, j] * (c - cstar0[i, j])


@njit(parallel=True, **_JIT)
def topo_ocean(tinit, deltaslxy, topoxy, cxy):
    """``topoxy = tinit - deltaslxy`` ; ``cxy = (topoxy < 0)``."""
    n0, n1 = tinit.shape
    for i in prange(n0):
        for j in range(n1):
            t = tinit[i, j] - deltaslxy[i, j]
            topoxy[i, j] = t
            cxy[i, j] = 1.0 if t < 0.0 else 0.0


@njit(parallel=True, **_JIT)
def gmsl_rowsums(deltaslxy, cxy, beta, s_slc_ocn, s_c, s_slc_beta, s_cstar):
    """Row sums of the four post-loop GMSL / area grids, in one pass.

        s_slc_ocn  <- rowsum(deltaslxy * cxy)
        s_c        <- rowsum(cxy)
        s_slc_beta <- rowsum(deltaslxy * cxy * beta)
        s_cstar    <- rowsum(cxy * beta)

    Only the degree-0 coefficient of those four fields is ever used, and that
    coefficient is a quadrature average -- see ``DuccSHT.mean00_weights``.
    Summation is left to right within a row, exactly as ``np.add.reduce`` would
    not do (it is pairwise), so these four numbers differ from a full analysis
    at round-off level; they only enter printed diagnostics.
    """
    n0, n1 = cxy.shape
    for i in prange(n0):
        a = 0.0
        b = 0.0
        c = 0.0
        d = 0.0
        for j in range(n1):
            o = cxy[i, j]
            sl = deltaslxy[i, j] * o
            ob = o * beta[i, j]
            a += sl
            b += o
            c += sl * beta[i, j]
            d += ob
        s_slc_ocn[i] = a
        s_c[i] = b
        s_slc_beta[i] = c
        s_cstar[i] = d


# ======================================================== spectral kernels ===
@njit(parallel=True, **_JIT)
def b18(deltaicestar, deltaS_prev, load_inc_lm, viscous, rhoi, rhow, ttlekhl, ttl, out):
    """Eq. B18 assembly, fused.

        out = ttlekhl * ((rhoi*di + rhow*dSprev) + rhow*dS) + ttl * viscous

    which is exactly the operation order of the five NumPy statements it
    replaces (``A = rhoi*di; A += rhow*dSprev; A += rhow*dS;
    out = ttlekhl*A; out += ttl*viscous``).
    """
    nl = out.shape[0]
    for l in prange(nl):
        te = ttlekhl[l]
        t = ttl[l]
        for m in range(nl):
            a = (deltaicestar[l, m] * rhoi + deltaS_prev[l, m] * rhow) + load_inc_lm[l, m] * rhow
            out[l, m] = a * te + viscous[l, m] * t


@njit(parallel=True, **_JIT)
def dS_update(deltaS_prev, rOlm, cstarlm, tOlm, conserv, new, deltaS_cur):
    """The inner-loop ocean-load update, fused.

        new        = (((-dSprev) + rOlm) + conserv*cstarlm) - tOlm
        deltaS_cur = dSprev + new
    """
    nl = new.shape[0]
    for l in prange(nl):
        for m in range(nl):
            p = deltaS_prev[l, m]
            v = ((-p + rOlm[l, m]) + cstarlm[l, m] * conserv) - tOlm[l, m]
            new[l, m] = v
            deltaS_cur[l, m] = p + v


@njit(parallel=True, **_JIT)
def current_load(dicestar, load_inc_lm, rhoi, rhow, out):
    """``out = rhoi*dicestar + rhow*load_inc_lm`` with row ``l = 0`` zeroed."""
    nl = out.shape[0]
    for l in prange(nl):
        if l == 0:
            for m in range(nl):
                out[0, m] = 0.0
        else:
            for m in range(nl):
                out[l, m] = dicestar[l, m] * rhoi + load_inc_lm[l, m] * rhow


@njit(parallel=True, **_JIT)
def commit_add(p_load, current_load):
    """``p_load += current_load``."""
    nl = p_load.shape[0]
    for l in prange(nl):
        for m in range(nl):
            p_load[l, m] += current_load[l, m]


@njit(parallel=True, **_JIT)
def copy3(src_ice, dst_ice, src_ocn, dst_ocn, src_prev, dst_prev):
    """The three post-commit spectral copies, in one parallel pass."""
    nl = dst_ice.shape[0]
    for l in prange(nl):
        for m in range(nl):
            dst_ice[l, m] = src_ice[l, m]
            dst_ocn[l, m] = src_ocn[l, m]
            dst_prev[l, m] = src_prev[l, m]


@njit(parallel=True, **_JIT)
def dice_and_deltaice(icestarlm, oldicestarlm, p_ice, dicestarlm, deltaicestarlm):
    """``dice = icestar - oldicestar`` ; ``deltaice = p_ice + dice``."""
    nl = icestarlm.shape[0]
    for l in prange(nl):
        for m in range(nl):
            d = icestarlm[l, m] - oldicestarlm[l, m]
            dicestarlm[l, m] = d
            deltaicestarlm[l, m] = p_ice[l, m] + d


# ============================================ second-order load placement ===
@njit(parallel=True, **_JIT)
def ramp_commit(R, omphi, nmod, L):
    """``R[l,k,m] += omphi[k,l] * L[l,m]`` (ramp commit, ``ExactMemory.commit`` with ramp on)."""
    nl = R.shape[0]
    for l in prange(nl):
        for k in range(nmod[l]):
            c = omphi[k, l]
            for m in range(l + 1):
                R[l, k, m] += L[l, m] * c


@njit(parallel=True, **_JIT)
def rate_sums(R, w_sl, w_rr, nmod, out_sl, out_rr):
    """``out_sl = sum_k w_sl[k,l] R[l,k,m]``, ``out_rr`` likewise (rate read-out; state unchanged)."""
    nl = R.shape[0]
    for l in prange(nl):
        for m in range(nl):
            out_sl[l, m] = 0.0
            out_rr[l, m] = 0.0
        for k in range(nmod[l]):
            a = w_sl[k, l]
            b = w_rr[k, l]
            for m in range(l + 1):
                x = R[l, k, m]
                out_sl[l, m] += x * a
                out_rr[l, m] += x * b


@njit(parallel=True, **_JIT)
def ramp_add(dicestar, load_inc_lm, rhoi, rhow, w, out):
    """``out += w[l] * (rhoi*dicestar + rhow*load_inc_lm)`` (same-step viscous term of the ramp)."""
    nl = out.shape[0]
    for l in prange(nl):
        c = w[l]
        for m in range(nl):
            out[l, m] += (dicestar[l, m] * rhoi + load_inc_lm[l, m] * rhow) * c


# NOTE: the "load_inc_lm first guess" statement ``(cstarlm / cstarlm[0,0]) * scale`` is
# deliberately NOT compiled here.  NumPy and LLVM use different complex-division
# algorithms (Smith's vs. the naive formula), so a compiled version differs from
# the NumPy one at round-off level (measured: 9.2e-16 relative).  It is a single
# pass over one 4 MB array and is left to NumPy so the path stays exact.


def warmup(nl=8, ng=(8, 16)):
    """Compile every kernel once on tiny arrays (call before timing)."""
    n0, n1 = ng
    g = np.zeros((n0, n1))
    gb = np.zeros((n0, n1), dtype=np.bool_)
    s = np.zeros((nl, nl), dtype=np.complex128)
    r1 = np.zeros(nl)
    R = np.zeros((nl, 2, nl), dtype=np.complex128)
    w = np.zeros((2, nl))
    nm = np.full(nl, 2, dtype=np.int64)
    advance_and_read(R, w, w, w, w, nm, s, s.copy(), s.copy())
    icestar(g, gb, g, g.copy())
    beta_cstar_tO(g, g, g, g, 1e-16, g.copy(), g.copy(), g.copy())
    add_const(g, 1.0, g.copy())
    rO_only(g, g, g.copy())
    topo_ocean_cstar_tO(g, g, g, g, g.copy(), g.copy(), g.copy(), g.copy())
    topo_ocean(g, g, g.copy(), g.copy())
    gmsl_rowsums(g, g, g, np.zeros(n0), np.zeros(n0), np.zeros(n0), np.zeros(n0))
    b18(s, s, s, s, 1.0, 1.0, r1, r1, s.copy())
    dS_update(s, s, s, s, 1.0, s.copy(), s.copy())
    current_load(s, s, 1.0, 1.0, s.copy())
    commit_add(s.copy(), s)
    copy3(s, s.copy(), s, s.copy(), s, s.copy())
    dice_and_deltaice(s, s, s, s.copy(), s.copy())
    ramp_commit(R, w, nm, s)
    rate_sums(R, w, w, nm, s.copy(), s.copy())
    ramp_add(s, s, 1.0, 1.0, r1, s.copy())
