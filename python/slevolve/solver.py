"""
solver.py -- NumPy implementation of the exact-memory (history-free) pseudo-spectral
sea-level solver.

Corresponds to the Fortran sources
  fortran/src/sl_model_mod.f90               (sl_solver_init, sl_solver)
  fortran/src/modal_state_mod.f90    (R_k <- d_k R_k + (1-d_k) P)
  the legacy spharmt transform with corrected fft99b   (via python/slevolve/sht.py)
  fortran/src/support_mod.f90        (constants_mod / planets_mod)

The step length dt (years) is a parameter: run_ais.py passes --dt (default 1 yr; ice file k is
timestep k), and run_outer.py passes its --dt (500 yr for Case 3).

Everything the Fortran round-trips through text files between timesteps
(tgrid0/ocean0/beta0, the previous ocean function, the hereditary checkpoint,
the Love-number file) is kept in memory here; the physics, the iteration order,
the convergence test and the advance/commit semantics are unchanged.
"""

from __future__ import annotations

import math
import time

import numpy as np

from .sht import SHT

# --------------------------------------------------------------------------- constants
# constants_mod (note: pi is the *truncated* literal in the legacy source, not math.pi)
PI = 3.1415926535898
II = 1j
EPS0 = float(np.finfo(np.float64).eps)  # Fortran epsilon(0.0) under -fdefault-real-8
# Floating-ice convention: beta = 1 where the grounded ice I* is zero (Kendall et al. 2005,
# eqs 43-44), flotation on the current topography iterate.  A cell at flotation can make the grounded mask alternate between two states
# from one inner iterate to the next, so the mask is updated at most GR_MAXACC times per step and
# then kept for the rest of the step (the ocean function is still re-guessed every iterate).
GR_MAXACC = 3
NPAM = 500


class Planet:
    """planets_mod::earth_init"""

    def __init__(self):
        self.radius = 6.371e6
        self.mass = 5.976e24
        self.rhoi = 920.0
        self.rhow = 1000.0
        self.gacc = 9.80665
        self.omega = 7.292e-5
        self.moiA = 0.3296145 * self.mass * self.radius**2
        self.moiC = 0.3307007 * self.mass * self.radius**2
        self.kf = 0.9342 + 0.008


# ------------------------------------------------------------------- Love number reader
class _ListDirectedReader:
    """Fortran list-directed READ semantics: every READ statement starts a new record."""

    def __init__(self, path):
        with open(path, "r") as fh:
            self._lines = fh.read().split("\n")
        self._i = 0

    def skip_record(self):
        self._i += 1

    def read(self, n):
        vals = []
        while len(vals) < n:
            vals.extend(self._lines[self._i].split())
            self._i += 1
        return vals[:n]  # surplus items in the last record are discarded

    def read_f(self, n):
        return [float(v) for v in self.read(n)]


class LoveNumbers:
    """Reader for the Love-number text file, mirroring the loop in sl_solver."""

    def __init__(self, path, norder):
        self.norder = norder
        rd = _ListDirectedReader(path)
        rd.skip_record()  # header line

        nmod = np.zeros(norder + 1, dtype=np.int64)
        s = np.zeros((NPAM, norder + 1))
        elast = np.zeros((3, norder + 1))
        asymv = np.zeros((3, norder + 1))
        telast = np.zeros((3, norder + 1))
        tasymv = np.zeros((3, norder + 1))
        resh = np.zeros((NPAM, norder + 1))
        resl = np.zeros((NPAM, norder + 1))
        resk = np.zeros((NPAM, norder + 1))
        tresh = np.zeros((NPAM, norder + 1))
        tresl = np.zeros((NPAM, norder + 1))
        tresk = np.zeros((NPAM, norder + 1))

        for _ in range(norder):
            tok = rd.read(2)
            legord = int(tok[0])
            nm = int(tok[1])
            ll = legord
            nmod[ll] = nm
            s[:nm, ll] = rd.read_f(nm)
            v = rd.read_f(5)
            elast[0, ll], elast[1, ll], elast[2, ll] = v[0], v[1], v[4]
            v = rd.read_f(5)
            asymv[0, ll], asymv[1, ll], asymv[2, ll] = v[0], v[1], v[4]
            resh[:nm, ll] = rd.read_f(nm)
            resl[:nm, ll] = rd.read_f(nm)
            resk[:nm, ll] = rd.read_f(nm)
            if float(legord) < 1.5:
                continue
            v = rd.read_f(5)
            telast[0, ll], telast[1, ll], telast[2, ll] = v[0], v[1], v[4]
            v = rd.read_f(5)
            tasymv[0, ll], tasymv[1, ll], tasymv[2, ll] = v[0], v[1], v[4]
            tresh[:nm, ll] = rd.read_f(nm)
            tresl[:nm, ll] = rd.read_f(nm)
            tresk[:nm, ll] = rd.read_f(nm)

        ell = np.arange(1, norder + 1, dtype=np.float64)
        resl[:, 1:] /= ell
        resk[:, 1:] /= ell
        tresl[:, 1:] /= ell
        tresk[:, 1:] /= ell

        self.nmod = nmod
        self.s = s
        self.he = np.zeros(norder + 1)
        self.ke = np.zeros(norder + 1)
        self.hTE = np.zeros(norder + 1)
        self.kTE = np.zeros(norder + 1)
        self.he[1:] = elast[0, 1:]
        self.ke[1:] = elast[2, 1:] / ell
        self.hTE[1:] = telast[0, 1:]
        self.kTE[1:] = telast[2, 1:] / ell
        self.rprime = resk
        self.r = resh
        self.rprimeT = tresk
        self.rT = tresh
        self.max_modes = int(nmod.max())
        # Tangential-displacement (l) Love numbers, used only by the optional
        # horizontal read-out: the file stores l*l^E and l*l_k, divided above
        # (residues) and here (elastic), exactly as for k.
        self.le = np.zeros(norder + 1)
        self.lTE = np.zeros(norder + 1)
        self.le[1:] = elast[1, 1:] / ell
        self.lTE[1:] = telast[1, 1:] / ell
        self.rl = resl
        self.rlT = tresl
        # Fluid (t -> infinity) limits, read from the second line of each block.  Not used by
        # the solver; kept so that python/tools/write_lovenum.py can write a file back.
        # Same scaling as the elastic line: h as stored, l and k divided by the degree.
        self.h_fluid = np.zeros(norder + 1)
        self.l_fluid = np.zeros(norder + 1)
        self.k_fluid = np.zeros(norder + 1)
        self.hT_fluid = np.zeros(norder + 1)
        self.lT_fluid = np.zeros(norder + 1)
        self.kT_fluid = np.zeros(norder + 1)
        self.h_fluid[1:] = asymv[0, 1:]
        self.l_fluid[1:] = asymv[1, 1:] / ell
        self.k_fluid[1:] = asymv[2, 1:] / ell
        self.hT_fluid[1:] = tasymv[0, 1:]
        self.lT_fluid[1:] = tasymv[1, 1:] / ell
        self.kT_fluid[1:] = tasymv[2, 1:] / ell


# ------------------------------------------------------------- second-order placement
def ramp_factor(x):
    """phi(x) = (1 - e^{-x}) / x and 1 - phi(x), stable as x -> 0 (phi(0) = 1, 1 - phi(0) = 0).

    With the optional second-order load placement (``ramp``) a load increment accrues linearly over
    its step of length dt; mode k then carries phi(s_k dt) of the step response as decaying memory
    and 1 - phi(s_k dt) as the response already reached at the end of the step."""
    x = np.asarray(x, dtype=np.float64)
    small = np.abs(x) < 1.0e-4
    xs = np.where(small, 1.0, x)
    phi = np.where(small, 1.0 - x / 2.0 + x * x / 6.0 - x**3 / 24.0, -np.expm1(-xs) / xs)
    omp = np.where(small, x / 2.0 - x * x / 6.0 + x**3 / 24.0, (xs + np.expm1(-xs)) / xs)
    return phi, omp


# ------------------------------------------------------------- exact hereditary memory
class ExactMemory:
    """History-free hereditary state, port of ``modal_state_mod``.

    Stored as dense ``[k, l, m]`` blocks (upper triangle m>l stays identically
    zero); ``R_k <- d_k R_k + (1-d_k) P`` advanced once per timestep, read while
    frozen through the inner loop, and committed once after convergence.
    """

    def __init__(self, love: LoveNumbers, lmax: int, dt_years: float, tpw: bool, nthreads: int = 1):
        self.lmax = lmax
        self.tpw = tpw
        kmax = love.max_modes
        self.kmax = kmax
        nl = lmax + 1

        s = np.zeros((kmax, nl))
        a_sl = np.zeros((kmax, nl))
        a_rr = np.zeros((kmax, nl))
        a_il = np.zeros(kmax)
        a_m = np.zeros(kmax)
        a_tsl = np.zeros(kmax)
        a_trr = np.zeros(kmax)
        for l in range(1, lmax + 1):
            nm = int(love.nmod[l])
            sl = love.s[:nm, l]
            s[:nm, l] = sl
            a_sl[:nm, l] = (love.rprime[:nm, l] - love.r[:nm, l]) / sl
            a_rr[:nm, l] = love.r[:nm, l] / sl
            if l == 2:
                a_il[:nm] = love.rprime[:nm, l] / sl
                a_m[:nm] = love.rprimeT[:nm, l] / sl
                a_tsl[:nm] = (love.rprimeT[:nm, l] - love.rT[:nm, l]) / sl
                a_trr[:nm] = love.rT[:nm, l] / sl
        self.a_sl, self.a_rr = a_sl, a_rr
        self.a_il, self.a_m, self.a_tsl, self.a_trr = a_il, a_m, a_tsl, a_trr
        self.nmod2 = int(love.nmod[2])

        # decay(k,l) = exp(-s*dt/1000), zero where the mode does not exist
        self._s = s
        self._nmod = [int(love.nmod[l]) for l in range(nl)] if nl > 0 else []
        self.dt_years = None
        self.decay = np.zeros((kmax, nl))
        self.omd = np.ones((kmax, nl))
        # second-order placement (ramp) and the analytic rate read-out: 1 - phi(s dt) per mode, its
        # per-degree sums, and the rate weights a_k s_k (per year).  Used only when ramp is set
        # (commit) or when a rate is read (read_rate); the default path does not touch them.
        self.ramp = False
        self.phi = np.ones((kmax, nl))
        self.omphi = np.zeros((kmax, nl))
        self.w_sl = a_sl * s / 1000.0
        self.w_rr = a_rr * s / 1000.0
        self.set_dt(dt_years)

        self.R = np.zeros((kmax, nl, nl), dtype=np.complex128)
        self.p_load = np.zeros((nl, nl), dtype=np.complex128)
        self.p_ice = np.zeros((nl, nl), dtype=np.complex128)
        self.p_ocean = np.zeros((nl, nl), dtype=np.complex128)
        self.prev_ice = np.zeros((nl, nl), dtype=np.complex128)

        self.r_il = np.zeros((kmax, 3, 3))
        self.r_m = np.zeros((kmax, 3))
        self.r_lambda = np.zeros((kmax, 3), dtype=np.complex128)
        self.p_il = np.zeros((3, 3))
        self.p_m = np.zeros(3)
        self.p_lambda = np.zeros(3, dtype=np.complex128)

        self._viscous = np.zeros((nl, nl), dtype=np.complex128)
        self._radial = np.zeros((nl, nl), dtype=np.complex128)
        self._tmp = np.zeros((nl, nl), dtype=np.complex128)

        # ---- optional row-blocked, multi-threaded hereditary advance --------
        # The (k, l) update and the k-accumulation order are *unchanged* inside
        # every row block, so the result is identical to the serial loop;
        # blocking over l merely lets independent row ranges run concurrently
        # (NumPy's ufunc inner loops release the GIL).  nthreads = 1 keeps the
        # original serial code path exactly.
        self.nthreads = int(nthreads)
        self._ex = None
        if self.nthreads > 1:
            from concurrent.futures import ThreadPoolExecutor

            self._ex = ThreadPoolExecutor(max_workers=self.nthreads)
            nblk = self.nthreads
            b = np.linspace(0, nl, nblk + 1).astype(int)
            self._blk = [(int(b[i]), int(b[i + 1])) for i in range(nblk) if b[i + 1] > b[i]]
            self._blk_tmp = [np.empty((l1 - l0, nl), dtype=np.complex128) for l0, l1 in self._blk]

    def set_dt(self, dt_years):
        """Set the length of the next step(s) (years): d_k = exp(-s_k dt / 1000).

        Recomputed only when dt changes (cached), so a uniform run computes the factors once and a
        variable-step run once per distinct step length.  The arrays are updated in place."""
        dt_years = float(dt_years)
        if dt_years == self.dt_years:
            return
        decay = np.zeros_like(self.decay)
        for l in range(1, self.lmax + 1):
            nm = self._nmod[l]
            decay[:nm, l] = np.exp(-self._s[:nm, l] * dt_years / 1000.0)
        self.decay[...] = decay
        self.omd[...] = 1.0 - decay
        self.dt_years = dt_years
        self._set_ramp_coeffs(dt_years)

    def _set_ramp_coeffs(self, dt_years):
        """phi, 1 - phi of every mode for the step length dt_years and the derived per-degree sums:
        c_sl, c_rr = sum_k a_k (1 - phi_k)  (same-step viscous weight of the current increment),
        ws_sl, ws_rr = sum_k a_k s_k, wsphi_sl, wsphi_rr = sum_k a_k s_k phi_k (rate read-out),
        and at degree 2 the rotational c_il, c_m, c_tsl, c_trr."""
        phi, omp = ramp_factor(self._s * (dt_years / 1000.0))
        exists = self._s > 0
        self.phi[...] = np.where(exists, phi, 1.0)
        self.omphi[...] = np.where(exists, omp, 0.0)
        self.c_sl = (self.a_sl * self.omphi).sum(0)
        self.c_rr = (self.a_rr * self.omphi).sum(0)
        self.ws_sl = self.w_sl.sum(0)
        self.ws_rr = self.w_rr.sum(0)
        self.wsphi_sl = (self.w_sl * self.phi).sum(0)
        self.wsphi_rr = (self.w_rr * self.phi).sum(0)
        n2 = self.nmod2
        c2 = self.omphi[:n2, 2] if self.lmax >= 2 else np.zeros(0)
        self.c_il = float(c2 @ self.a_il[:n2])
        self.c_m = float(c2 @ self.a_m[:n2])
        self.c_tsl = float(c2 @ self.a_tsl[:n2])
        self.c_trr = float(c2 @ self.a_trr[:n2])

    def _advance_block(self, b):
        l0, l1 = self._blk[b]
        t = self._blk_tmp[b]
        P = self.p_load[l0:l1]
        visc = self._viscous[l0:l1]
        rad = self._radial[l0:l1]
        for k in range(self.kmax):
            Rk = self.R[k, l0:l1]
            Rk *= self.decay[k, l0:l1, None]
            np.multiply(P, self.omd[k, l0:l1, None], out=t)
            Rk += t
            np.multiply(Rk, self.a_sl[k, l0:l1, None], out=t)
            visc += t
            np.multiply(Rk, self.a_rr[k, l0:l1, None], out=t)
            rad += t

    def advance_and_read(self):
        """Advance every hereditary pole one timestep and return the response."""
        visc = self._viscous
        rad = self._radial
        tmp = self._tmp
        visc[:] = 0.0
        rad[:] = 0.0
        P = self.p_load
        if self._ex is not None:
            list(self._ex.map(self._advance_block, range(len(self._blk))))
        else:
            for k in range(self.kmax):
                Rk = self.R[k]
                Rk *= self.decay[k][:, None]
                np.multiply(P, self.omd[k][:, None], out=tmp)
                Rk += tmp
                np.multiply(Rk, self.a_sl[k][:, None], out=tmp)
                visc += tmp
                np.multiply(Rk, self.a_rr[k][:, None], out=tmp)
                rad += tmp

        sum_il = np.zeros((3, 3))
        sum_m = np.zeros(3)
        sum_tsl = np.zeros(3, dtype=np.complex128)
        sum_trr = np.zeros(3, dtype=np.complex128)
        if self.lmax >= 2:
            for k in range(self.nmod2):
                d = self.decay[k, 2]
                omd = 1.0 - d
                self.r_il[k] = d * self.r_il[k] + omd * self.p_il
                self.r_m[k] = d * self.r_m[k] + omd * self.p_m
                self.r_lambda[k] = d * self.r_lambda[k] + omd * self.p_lambda
                sum_il += self.a_il[k] * self.r_il[k]
                sum_m += self.a_m[k] * self.r_m[k]
                sum_tsl += self.a_tsl[k] * self.r_lambda[k]
                sum_trr += self.a_trr[k] * self.r_lambda[k]
        return visc, rad, sum_il, sum_m, sum_tsl, sum_trr

    def horizontal_weights(self, love):
        """a_hl[k, l] = l_k / s_k and the degree-2 tidal a_thl[k] (horizontal read-out)."""
        nl = self.lmax + 1
        a_hl = np.zeros((self.kmax, nl))
        a_thl = np.zeros(self.kmax)
        for l in range(1, self.lmax + 1):
            nm = int(love.nmod[l])
            a_hl[:nm, l] = love.rl[:nm, l] / love.s[:nm, l]
            if l == 2:
                a_thl[:nm] = love.rlT[:nm, l] / love.s[:nm, l]
        self.a_hl, self.a_thl = a_hl, a_thl

    def read_horizontal(self):
        """sum_k a_hl R_k over the CURRENT states (output only; state unchanged).

        Valid anywhere between this step's advance_and_read and the next one:
        commit only adds to p_load, it does not touch R_k."""
        hz = np.zeros((self.lmax + 1, self.lmax + 1), dtype=np.complex128)
        for k in range(self.kmax):
            hz += self.a_hl[k][:, None] * self.R[k]
        thl = np.zeros(3, dtype=np.complex128)
        if self.lmax >= 2:
            for k in range(self.nmod2):
                thl += self.a_thl[k] * self.r_lambda[k]
        return hz, thl

    def _ramp_commit_load(self, current_load):
        for k in range(self.kmax):
            self.R[k] += self.omphi[k][:, None] * current_load

    def _rate_sums(self):
        """sum_k a_k s_k R_k for sea level and radial displacement (per year; state unchanged)."""
        nl = self.lmax + 1
        out_sl = np.zeros((nl, nl), dtype=np.complex128)
        out_rr = np.zeros((nl, nl), dtype=np.complex128)
        for k in range(self.kmax):
            out_sl += self.w_sl[k][:, None] * self.R[k]
            out_rr += self.w_rr[k][:, None] * self.R[k]
        return out_sl, out_rr

    def read_rate(self, current_load):
        """Analytic left time derivative (per year) of the viscous sea-level and radial terms at the
        end of the current step, the current increment ``current_load`` taken to accrue linearly over
        the step: sum_k a_k s_k (P + dL - R_k - (1 - phi_k) dL) = ws P - sum_k a_k s_k R_k + wsphi dL,
        with P = p_load the load committed before this step.  Call between advance_and_read and
        commit of the step (state unchanged).  Each state obeys dR_k/dt = s_k (P(t) - R_k) inside
        the step, read with the ramp-consistent state R_k + (1 - phi_k) dL."""
        q_sl, q_rr = self._rate_sums()
        rs = self.ws_sl[:, None] * self.p_load - q_sl + self.wsphi_sl[:, None] * current_load
        rr = self.ws_rr[:, None] * self.p_load - q_rr + self.wsphi_rr[:, None] * current_load
        return rs, rr

    def read_rate_rotation(self, dil, dm, dlam):
        """Degree-2 rotational states, ramp-consistent at the end of the step, for the rate
        read-out: r_k + (1 - phi_k) d(.), mode by mode (k < nmod2).  State unchanged."""
        n2 = self.nmod2
        c = self.omphi[:n2, 2]
        r_il = self.r_il[:n2] + c[:, None, None] * dil[None]
        r_m = self.r_m[:n2] + c[:, None] * dm[None]
        r_la = self.r_lambda[:n2] + c[:, None] * dlam[None]
        return r_il, r_m, r_la, self._s[:n2, 2] / 1000.0

    def commit(self, current_load, dil=None, dm=None, dlambda=None):
        if self.ramp:
            # second-order placement: the increment of this step already carries (1 - phi_k) of
            # its response in every mode; one step later it contributes d c + (1 - d) = 1 - d phi
            self._ramp_commit_load(current_load)
            if self.lmax >= 2:
                n2 = self.nmod2
                c = self.omphi[:n2, 2]
                if dil is not None:
                    self.r_il[:n2] += c[:, None, None] * dil[None]
                if dm is not None:
                    self.r_m[:n2] += c[:, None] * dm[None]
                if dlambda is not None:
                    self.r_lambda[:n2] += c[:, None] * dlambda[None]
        self.p_load += current_load
        if dil is not None:
            self.p_il += dil
        if dm is not None:
            self.p_m += dm
        if dlambda is not None:
            self.p_lambda += dlambda


# ------------------------------------------------------------------------------ solver
def read_grid_text(path, nlat=512, nlon=1024):
    """Read one of the legacy text grid files (linear stream order -> [lat, lon])."""
    a = np.fromfile(path, sep=" ")
    if a.size != nlat * nlon:
        raise ValueError(f"{path}: got {a.size} values, expected {nlat*nlon}")
    return a.reshape(nlat, nlon)


class SeaLevelSolver:
    """Exact-memory sea-level solver, NumPy/BLAS implementation."""

    def __init__(
        self,
        sht: SHT,
        love: LoveNumbers,
        topo0: np.ndarray,
        ice0: np.ndarray,
        *,
        norder=512,
        dt=1,
        tpw=True,
        checkmarine=True,
        calcRG=True,
        epsilon1=1.0e-5,
        planet=None,
        mem_threads=1,
        ramp=False,
    ):
        self.sht = sht
        self.love = love
        self.norder = norder
        self.nl = norder + 1
        self.dt = dt
        self.tpw = tpw
        self.checkmarine = checkmarine
        self.calcRG = calcRG
        self.horizontal = False  # set True for the optional horizontal read-out
        self.epsilon1 = epsilon1
        self.p = planet or Planet()
        self.nlat, self.nlon = topo0.shape

        # ---------------- sl_solver_init ------------------------------------
        self.tinit_0 = np.ascontiguousarray(topo0, dtype=np.float64)
        self.tinit = self.tinit_0  # Travel == 0 for the standard chronology
        self.cxy0 = (self.tinit_0 < 0).astype(np.float64)

        # checkmarine helpers (tinit is constant for the standard chronology)
        self._marine_pos = self.tinit > 0
        self._marine_thr = np.abs(self.tinit) * self.p.rhow / self.p.rhoi

        # grounded-ice mask (Kendall et al. 2005, eqs 43-44): beta = 1 where I* = 0
        g0 = self._icestar(ice0) if self.checkmarine else ice0
        self.beta0 = (g0 < EPS0).astype(np.float64)
        self.cstar0 = self.cxy0 * self.beta0
        self.cxy = self.cxy0.copy()  # nmelt == 1 at the first marching step
        self._topo_prev = self.tinit.copy()  # flotation reference of the next step
        self.gr_counts = {"reguess": 0, "ice_analyses": 0}

        nl = self.nl
        self.clm0 = sht.analysis(self.cxy0, out=np.zeros((nl, nl), complex))
        self.ocean_area0 = self.clm0[0, 0].real * 4 * PI * self.p.radius**2
        cst = sht.analysis(self.cstar0, out=np.zeros((nl, nl), complex))
        self.oceanBeta_area0 = cst[0, 0].real * 4 * PI * self.p.radius**2
        if self.checkmarine:
            ic0 = self._icestar(ice0)
        else:
            ic0 = ice0
        self.ice_volume0 = (
            sht.analysis(ic0, out=np.zeros((nl, nl), complex))[0, 0].real * 4 * PI * self.p.radius**2
        )

        self.state = ExactMemory(love, norder, float(dt), tpw, nthreads=mem_threads)
        # optional second-order load placement (default off: step placement of Kendall et al. 2005)
        self.ramp = bool(ramp)
        self.state.ramp = self.ramp

        # rotational module variables (persist across timesteps, as in Fortran)
        self.il = np.zeros((3, 3))
        self.mm = np.zeros(3)
        self.lam = np.zeros((3, 3), dtype=np.complex128)

        # ---------------- cached Love combinations --------------------------
        l = np.arange(nl, dtype=np.float64)
        ttl = np.zeros(nl)
        ttl[1:] = (4.0 * PI * (self.p.radius**3)) / (self.p.mass * (2 * l[1:] + 1.0))
        ekhl = np.zeros(nl)
        ekhl[1:] = 1 + love.ke[1:] - love.he[1:]
        self._ttl = ttl[:, None]
        self._ttlekhl = (ttl * ekhl)[:, None]
        self._ttlhe = (ttl * np.concatenate(([0.0], love.he[1:])))[:, None]

        # ---------------- persistent work arrays ----------------------------
        sh = (self.nlat, self.nlon)
        self.icestarxy = np.zeros(sh)
        self.beta = np.zeros(sh)
        self.cstarxy = np.zeros(sh)
        self.tOxy = np.zeros(sh)
        self.rOxy = np.zeros(sh)
        self.dslxy = np.zeros(sh)
        self.deltaslxy = np.zeros(sh)
        self.topoxy = np.zeros(sh)
        self.rrxy = np.zeros(sh)
        self.tmpg = np.zeros(sh)

        z = lambda: np.zeros((nl, nl), dtype=np.complex128)
        self.icestarlm = z()
        self.cstarlm = z()
        self.tOlm = z()
        self.rOlm = z()
        self.load_inc_lm = z()
        self.old_load_inc_lm = z()
        self.sl_inc_lm = z()
        self.sl_tot_lm = z()
        self.rrlm = z()
        self.deltaS_cur = z()
        self.deltaicestarlm = z()
        self.dicestarlm = z()
        self.oldicestarlm = z()
        self.current_load = z()
        self.clm = z()
        self.slclm = z()
        self._sp1 = z()
        self._sp2 = z()

        self._gr_new = None
        self._gr_nacc = 0
        self.iter = 0
        self.icexy_prev = np.ascontiguousarray(ice0, dtype=np.float64)
        self.prof = {}

    # -------------------------------------------------------------- helpers
    def _icestar(self, ice, out=None):
        """checkmarine grounded-ice projection (eq. 44 region test)."""
        if out is None:
            out = np.zeros_like(ice)
        keep = self._marine_pos | (ice > self._marine_thr)
        np.multiply(ice, keep, out=out)
        return out

    def _icestar_on(self, ice, topo, out):
        """Grounded ice on the topography ``topo`` (flotation test of eq. 44)."""
        keep = (topo > 0) | (ice > np.abs(topo) * self.p.rhow / self.p.rhoi)
        np.multiply(ice, keep, out=out)
        return out

    def _grounded_reguess(self, icexy_cur, oldice):
        """checkmarine: grounded ice I*, beta on the topography iterate ``self.topoxy``.

        Called wherever the ocean function is re-guessed (after ``topoxy``/``cxy`` are set, before
        ``cstarxy``).  Only if the grounded flag of at least one cell changed are I*, beta, the
        ice spectrum (one analysis) and the ice-load increments updated.  Returns that count.
        """
        if self._gr_nacc >= GR_MAXACC:
            return 0
        if self._gr_new is None:
            self._gr_new = np.zeros_like(self.icestarxy)
        self._icestar_on(icexy_cur, self.topoxy, out=self._gr_new)
        bnew = self._gr_new < EPS0
        nd = int(np.count_nonzero(bnew != (self.beta != 0.0)))
        self.gr_counts["reguess"] += 1
        if nd:
            self.icestarxy[:] = self._gr_new
            self.beta[:] = bnew
            self.sht.analysis(self.icestarxy, out=self.icestarlm)
            np.subtract(self.icestarlm, oldice, out=self.dicestarlm)
            np.add(self.state.p_ice, self.dicestarlm, out=self.deltaicestarlm)
            self.gr_counts["ice_analyses"] += 1
            self._gr_nacc += 1
        return nd

    # ------------------------------------------------------------- one step
    def horizontal_lm(self, deltaicestar, deltaS_prev, load_inc_lm, rot_tot):
        """Optional horizontal read-out, returned in the spharmt (U) convention.

        Psi_lm = T_l [ l^E_l L_lm + sum_k (l_k/s_k) R_k,lm ] (+ degree-2 rotation
        with the tidal l numbers), u_h = grad_h Psi (theta = colatitude).  Uses
        the same modal states and the same instantaneous load as R.  Output
        only: nothing here feeds back into the solver.
        """
        p = self.p
        st = self.state
        if getattr(st, "a_hl", None) is None:
            st.horizontal_weights(self.love)
            le = np.concatenate(([0.0], self.love.le[1:]))
            self._ttlle = self._ttl * le[:, None]
        visc_h, thl = st.read_horizontal()
        A = p.rhoi * deltaicestar + p.rhow * deltaS_prev
        A += p.rhow * load_inc_lm
        hz = self._ttlle * A + self._ttl * visc_h
        if self.ramp:
            # same-step viscous weight of the current increment (second-order placement)
            c_hl = (st.a_hl * st.omphi).sum(0)
            hz += (self._ttl[:, 0] * c_hl)[:, None] * (p.rhoi * self.dicestarlm + p.rhow * load_inc_lm)
        hz[0, :] = 0.0
        sg = getattr(self, "_msign", 1.0)
        if self.tpw and rot_tot is not None and self.norder >= 2:
            if self.ramp:
                c_thl = float(st.omphi[:st.nmod2, 2] @ st.a_thl[:st.nmod2])
                thl = thl + c_thl * (rot_tot - st.p_lambda)
            rot = (self.love.lTE[2] * rot_tot + thl) / p.gacc
            rot[1] *= sg
            hz[2, 0:3] += rot
        if sg != 1.0:  # raw (sign-free) backend: V = (-1)^m U -> U
            hz[:, 1::2] *= -1.0
        return hz

    # ------------------------------------------------------- analytic rates
    def _rotation_rate(self, oldil, oldm, dlambda_prev2):
        """Analytic left derivative (per year) at the end of the step of the rotational sea-level and
        radial terms, with the converged il, m, Lambda of the step.  The m equation is
        differentiated explicitly: the current-increment term of the ramp becomes s_k (m - r_m,k);
        Lambda-dot by the chain rule.  Returns (rate_sl, rate_rr) (3, 3) in the U convention and
        the rotational rates (il-dot, m-dot, Lambda-dot)."""
        p = self.p
        st = self.state
        h = float(self.dt)
        dil = self.il - oldil
        dm = self.mm - oldm
        tot = self.lam[2, 0:3].copy()
        r_il, r_m, r_la, sk = st.read_rate_rotation(dil, dm, tot - dlambda_prev2)
        n2 = st.nmod2
        kTE2, hTE2, kE2 = self.love.kTE[2], self.love.hTE[2], self.love.ke[2]
        ild = dil / h
        sild = np.einsum("k,kij->ij", st.a_il[:n2] * sk, self.il[None] - r_il)
        smd = np.einsum("k,ki->i", st.a_m[:n2] * sk, self.mm[None] - r_m)
        alpha = 1 / (1 - kTE2 / p.kf)
        md = np.zeros(3)
        for i in (0, 1):
            md[i] = alpha * ((1 / (p.moiC - p.moiA)) * ((1 + kE2) * ild[i, 2] + sild[i, 2]) + smd[i] / p.kf)
        md[2] = (1.0 / p.moiC) * ((1 + kE2) * ild[2, 2] + sild[2, 2])
        m1, m2, m3 = self.mm
        d1, d2, d3 = md
        r2w2 = p.radius**2 * p.omega**2
        lad = np.zeros((3, 3), dtype=np.complex128)
        lad[0, 0] = (r2w2 / 3.0) * (2 * m1 * d1 + 2 * m2 * d2 + 2 * m3 * d3 + 2.0 * d3)
        lad[2, 0] = (r2w2 / (6.0 * math.sqrt(5.0))) * (2 * m1 * d1 + 2 * m2 * d2 - 4 * m3 * d3 - 4 * d3)
        lad[2, 1] = (r2w2 / math.sqrt(30.0)) * ((d1 - II * d2) * (1.0 + m3) + (m1 - II * m2) * d3)
        lad[2, 2] = (r2w2 / math.sqrt(5.0 * 24.0)) * (2 * m2 * d2 - 2 * m1 * d1 + 2.0 * II * (d1 * m2 + m1 * d2))
        tsld = np.einsum("k,ki->i", st.a_tsl[:n2] * sk, tot[None] - r_la)
        trrd = np.einsum("k,ki->i", st.a_trr[:n2] * sk, tot[None] - r_la)
        rate_sl = np.zeros((3, 3), dtype=np.complex128)
        rate_rr = np.zeros((3, 3), dtype=np.complex128)
        rate_sl[2, 0:3] = ((1 + kTE2 - hTE2) * lad[2, 0:3] + tsld) / p.gacc
        rate_rr[2, 0:3] = (hTE2 * lad[2, 0:3] + trrd) / p.gacc
        return rate_sl, rate_rr, {"ildot": ild, "mdot": md, "lambdadot": lad}

    def _rate_spectra(self, current_load, oldil, oldm, dlambda_prev2):
        """Rate spectra (per year) of sea level (degrees >= 1) and radial displacement at the end of
        the step; call before the commit.  The current increment is taken to accrue linearly over
        the step (elastic part dL / dt), with or without the ramp placement."""
        st = self.state
        h = float(self.dt)
        rv_sl, rv_rr = st.read_rate(current_load)
        dsl = self._ttlekhl * current_load / h + self._ttl * rv_sl
        drr = self._ttlhe * current_load / h + self._ttl * rv_rr
        dsl[0, :] = 0.0
        drr[0, :] = 0.0
        rot = None
        if self.tpw and self.norder >= 2:
            rs, rr, rot = self._rotation_rate(oldil, oldm, dlambda_prev2)
            sg = getattr(self, "_msign", 1.0)
            rs[2, 1] *= sg
            rr[2, 1] *= sg
            dsl[0:3, 0:3] += rs
            drr[0:3, 0:3] += rr
        return dsl, drr, rot

    def _rate_grids(self, dsl, drr, dice00, rot):
        """Synthesise the rate spectra (after the end-of-step ocean function is set).  Degree 0 of
        dSL/dt is the rate of the conservation constant on the ocean function C* = C beta of the
        end of the step: (-(rho_i / rho_w) dI*_00 / dt - mean(C* dSL/dt|l>=1)) / mean(C*)."""
        p = self.p
        sht = self.sht
        nl = self.nl
        g = np.zeros_like(self.topoxy)
        sht.synthesis(dsl, out=g)
        cst = self.cxy * self.beta
        a = sht.analysis(g * cst, out=np.zeros((nl, nl), complex))[0, 0].real
        c00 = sht.analysis(cst, out=np.zeros((nl, nl), complex))[0, 0].real
        dsl[0, 0] = (-(p.rhoi / p.rhow) * dice00.real / float(self.dt) - a) / c00
        sht.synthesis(dsl, out=g)
        gr = np.zeros_like(self.topoxy)
        sht.synthesis(drr, out=gr)
        out = {"dSLdt": g, "dRdt": gr, "dSLdt_lm": dsl, "dRdt_lm": drr}
        if rot is not None:
            out.update(rot)
        return out

    def set_step_dt(self, dt):
        """Length (years) of the step(s) that follow; dt=None keeps the current one."""
        if dt is not None:
            self.dt = float(dt)
            self.state.set_dt(self.dt)

    def step(self, icexy_cur, profile=False, dt=None, rates=False):
        """Advance one physical timestep.  Returns a dict of the step's outputs.

        dt: length of this step in years (variable steps); None = the constructor's dt.
        rates: also return the analytic rates dSL/dt, dR/dt (grids, m/yr) at the end of the step and
        the rotational rates (output only; result key "rates")."""
        self.set_step_dt(dt)
        p = self.p
        sht = self.sht
        nl = self.nl
        self.iter += 1
        it = self.iter
        tm = {}
        t_all = time.perf_counter()

        def mark(key, t0):
            if profile:
                tm[key] = tm.get(key, 0.0) + time.perf_counter() - t0

        # ---- exact hereditary advance (once per step, before the inner loop)
        t0 = time.perf_counter()
        viscous, radial, sum_il, sum_m, exact_tsl, exact_trr = self.state.advance_and_read()
        mark("memory", t0)
        deltaS_prev = self.state.p_ocean
        oldil = self.state.p_il.copy()
        oldm = self.state.p_m.copy()
        dlambda_prev2 = self.state.p_lambda.copy()  # deltalambda(2,0:2,prev)

        # ---- ice part -------------------------------------------------------
        t0 = time.perf_counter()
        gr_cm = self.checkmarine
        self._gr_nacc = 0
        if gr_cm:
            # flotation on the converged topography of the previous step (step 0: tinit)
            self._icestar_on(icexy_cur, self._topo_prev, out=self.icestarxy)
            src = self.icestarxy
        else:
            src = icexy_cur
        sht.analysis(src, out=self.icestarlm)
        if it == 1:
            if self.checkmarine:
                self._icestar(self.icexy_prev, out=self.tmpg)
                sht.analysis(self.tmpg, out=self.oldicestarlm)
            else:
                sht.analysis(self.icexy_prev, out=self.oldicestarlm)
        else:
            self.oldicestarlm[:] = self.state.prev_ice
        mark("sht", t0)

        np.subtract(self.icestarlm, self.oldicestarlm, out=self.dicestarlm)
        np.add(self.state.p_ice, self.dicestarlm, out=self.deltaicestarlm)
        dicestar = self.dicestarlm
        deltaicestar = self.deltaicestarlm

        oldice = self.oldicestarlm
        # ---- ocean part -----------------------------------------------------
        if gr_cm:
            self.beta[:] = self.icestarxy < EPS0  # grounded-ice mask (Kendall et al. 2005, eqs 43-44)
        else:
            self.beta[:] = icexy_cur < EPS0
        np.multiply(self.cxy, self.beta, out=self.cstarxy)
        t0 = time.perf_counter()
        sht.analysis(self.cstarxy, out=self.cstarlm)
        mark("sht", t0)
        np.subtract(self.cstarxy, self.cstar0, out=self.tOxy)
        np.multiply(self.tinit, self.tOxy, out=self.tOxy)
        t0 = time.perf_counter()
        sht.analysis(self.tOxy, out=self.tOlm)
        mark("sht", t0)

        # legacy grouping: (cstarlm / cstarlm(0,0)) * ((-rhoi/rhow) * dicestar(0,0))
        np.divide(self.cstarlm, self.cstarlm[0, 0], out=self.load_inc_lm)
        self.load_inc_lm *= (-p.rhoi / p.rhow) * dicestar[0, 0]
        self.deltaS_cur[:] = 0.0  # Fortran: deltaS(:,:,cur) is fresh-zero here

        dsl_rot = np.zeros((3, 3), dtype=np.complex128)
        rr_rot = np.zeros((3, 3), dtype=np.complex128)
        dil_cur = np.zeros((3, 3))
        dm_cur = np.zeros(3)
        dlam_cur2 = np.zeros(3, dtype=np.complex128)

        kTE2 = self.love.kTE[2]
        hTE2 = self.love.hTE[2]
        kE2 = self.love.ke[2]
        r2w2 = p.radius**2 * p.omega**2
        if self.ramp:
            ttl_csl = (self._ttl[:, 0] * self.state.c_sl)[:, None]

        ninner = 1
        while True:
            # ------------------------------ rotation ------------------------
            if self.tpw:
                t0 = time.perf_counter()
                z20 = p.rhow * self.deltaS_cur[2, 0] + p.rhoi * deltaicestar[2, 0]
                z21 = p.rhow * self.deltaS_cur[2, 1] + p.rhoi * deltaicestar[2, 1]
                self.il[2, 2] = (-8 * PI * p.radius**4 / (3 * math.sqrt(5.0))) * z20.real
                self.il[0, 2] = (8 * PI * p.radius**4 / (math.sqrt(30.0))) * z21.real
                self.il[1, 2] = (-8 * PI * p.radius**4 / (math.sqrt(30.0))) * (-1.0 * II * z21).real
                dil_cur = self.il - oldil

                if self.ramp:
                    # current increments carry (1 - phi) of their viscous response: implicit in m
                    st_ = self.state
                    sil_r = sum_il + st_.c_il * dil_cur
                    al_ = 1 / (1 - kTE2 / p.kf)
                    for i in (0, 1):
                        self.mm[i] = al_ * (
                            (1 / (p.moiC - p.moiA)) * ((1 + kE2) * self.il[i, 2] + sil_r[i, 2])
                            + (1 / p.kf) * (sum_m[i] - st_.c_m * oldm[i])
                        ) / (1 - al_ * st_.c_m / p.kf)
                    self.mm[2] = (1.0 / p.moiC) * ((1 + kE2) * self.il[2, 2] + sil_r[2, 2])
                else:
                    for i in (0, 1):
                        self.mm[i] = (1 / (1 - kTE2 / p.kf)) * (
                            (1 / (p.moiC - p.moiA)) * ((1 + kE2) * self.il[i, 2] + sum_il[i, 2])
                            + (1 / p.kf) * sum_m[i]
                        )
                    self.mm[2] = (1.0 / p.moiC) * ((1 + kE2) * self.il[2, 2] + sum_il[2, 2])
                dm_cur = self.mm - oldm

                m1, m2, m3 = self.mm[0], self.mm[1], self.mm[2]
                self.lam[0, 0] = (r2w2 / 3.0) * ((m1**2 + m2**2 + m3**2) + 2.0 * m3)
                self.lam[2, 0] = (r2w2 / (6.0 * math.sqrt(5.0))) * (m1**2 + m2**2 - 2.0 * m3**2 - 4.0 * m3)
                self.lam[2, 1] = (r2w2 / math.sqrt(30.0)) * ((m1 - II * m2) * (1.0 + m3))
                self.lam[2, 2] = (r2w2 / math.sqrt(5.0 * 24.0)) * (m2**2 - m1**2 + 2.0 * II * m1 * m2)
                dlam_cur2 = self.lam[2, 0:3] - dlambda_prev2

                ekhTE = 1 + kTE2 - hTE2
                tot = dlambda_prev2 + dlam_cur2
                tsl_v, trr_v = exact_tsl, exact_trr
                if self.ramp:
                    tsl_v = exact_tsl + self.state.c_tsl * dlam_cur2
                    trr_v = exact_trr + self.state.c_trr * dlam_cur2
                dsl_rot[:] = 0.0
                dsl_rot[2, 0:3] = (ekhTE * tot + tsl_v) / p.gacc
                if self.calcRG:
                    rr_rot[:] = 0.0
                    rr_rot[2, 0:3] = (hTE2 * tot + trr_v) / p.gacc
                mark("rot", t0)
            else:
                dsl_rot[:] = 0.0
                rr_rot[:] = 0.0

            # ------------------------------ sea level (eq. B18) -------------
            t0 = time.perf_counter()
            A = self._sp1
            np.multiply(deltaicestar, p.rhoi, out=A)
            np.multiply(deltaS_prev, p.rhow, out=self._sp2)
            A += self._sp2
            np.multiply(self.load_inc_lm, p.rhow, out=self._sp2)
            A += self._sp2
            np.multiply(A, self._ttlekhl, out=self.sl_inc_lm)
            np.multiply(viscous, self._ttl, out=self._sp2)
            self.sl_inc_lm += self._sp2
            if self.ramp:
                self.sl_inc_lm += ttl_csl * (p.rhoi * self.dicestarlm + p.rhow * self.load_inc_lm)
            self.sl_inc_lm[0:3, 0:3] += dsl_rot
            mark("spec", t0)

            t0 = time.perf_counter()
            sht.synthesis(self.sl_inc_lm, out=self.dslxy)
            mark("sht", t0)

            np.multiply(self.dslxy, self.cstarxy, out=self.rOxy)
            t0 = time.perf_counter()
            sht.analysis(self.rOxy, out=self.rOlm)
            mark("sht", t0)

            conserv = (
                (1 / self.cstarlm[0, 0])
                * (-1.0 * (p.rhoi / p.rhow) * deltaicestar[0, 0] - self.rOlm[0, 0] + self.tOlm[0, 0])
            ).real

            self.old_load_inc_lm, self.load_inc_lm = self.load_inc_lm, self.old_load_inc_lm
            old = self.old_load_inc_lm
            new = self.load_inc_lm
            np.multiply(deltaS_prev, -1.0, out=new)
            new += self.rOlm
            np.multiply(self.cstarlm, conserv, out=self._sp2)
            new += self._sp2
            new -= self.tOlm
            np.add(deltaS_prev, new, out=self.deltaS_cur)

            self.sl_tot_lm[:] = self.sl_inc_lm
            self.sl_tot_lm[0, 0] += conserv
            t0 = time.perf_counter()
            sht.synthesis(self.sl_tot_lm, out=self.deltaslxy)
            mark("sht", t0)

            sa = np.abs(new).sum()
            so = np.abs(old).sum()
            if abs(sa - so) < EPS0 and abs(so) < EPS0:
                xi = 0.0
            elif abs(so) < EPS0:
                xi = abs((sa - so) / (EPS0 * 10))
            else:
                xi = abs((sa - so) / so)

            if xi > self.epsilon1:
                np.subtract(self.tinit, self.deltaslxy, out=self.topoxy)
                self.cxy[:] = self.topoxy < 0.0
                if gr_cm:
                    self._grounded_reguess(icexy_cur, oldice)
                np.multiply(self.cxy, self.beta, out=self.cstarxy)
                t0 = time.perf_counter()
                sht.analysis(self.cstarxy, out=self.cstarlm)
                mark("sht", t0)
                np.subtract(self.cstarxy, self.cstar0, out=self.tOxy)
                np.multiply(self.tinit, self.tOxy, out=self.tOxy)
                t0 = time.perf_counter()
                sht.analysis(self.tOxy, out=self.tOlm)
                mark("sht", t0)

            if xi <= self.epsilon1:
                break
            if ninner == 9999:
                raise RuntimeError(f"inner loop failed to converge, xi={xi}")
            ninner += 1

        # ------------------------------------------------------------ R / G
        if self.calcRG:
            t0 = time.perf_counter()
            A = self._sp1
            np.multiply(deltaicestar, p.rhoi, out=A)
            np.multiply(deltaS_prev, p.rhow, out=self._sp2)
            A += self._sp2
            np.multiply(self.load_inc_lm, p.rhow, out=self._sp2)
            A += self._sp2
            np.multiply(A, self._ttlhe, out=self.rrlm)
            np.multiply(radial, self._ttl, out=self._sp2)
            self.rrlm += self._sp2
            if self.ramp:
                self.rrlm += (self._ttl[:, 0] * self.state.c_rr)[:, None] * (
                    p.rhoi * self.dicestarlm + p.rhow * self.load_inc_lm)
            self.rrlm[0:3, 0:3] += rr_rot
            mark("spec", t0)
            t0 = time.perf_counter()
            sht.synthesis(self.rrlm, out=self.rrxy)
            mark("sht", t0)

        hz = None
        if self.horizontal:
            hz = self.horizontal_lm(deltaicestar, deltaS_prev, self.load_inc_lm,
                                    dlambda_prev2 + dlam_cur2 if self.tpw else None)

        # -------------------------------------------------------- commit ----
        t0 = time.perf_counter()
        np.multiply(dicestar, p.rhoi, out=self.current_load)
        np.multiply(self.load_inc_lm, p.rhow, out=self._sp2)
        self.current_load += self._sp2
        self.current_load[0, :] = 0.0  # Fortran loop runs l = 1..norder
        rate_sp = None
        if rates:
            rate_sp = self._rate_spectra(self.current_load, oldil, oldm, dlambda_prev2)
            dice00 = complex(dicestar[0, 0])
        if self.tpw:
            self.state.commit(self.current_load, dil_cur, dm_cur, dlam_cur2)
        else:
            self.state.commit(self.current_load)
        self.state.p_ice[:] = deltaicestar
        self.state.p_ocean[:] = self.deltaS_cur
        self.state.prev_ice[:] = self.icestarlm
        mark("memory", t0)

        # ------------------------------------------------- topography/ocean -
        np.subtract(self.tinit, self.deltaslxy, out=self.topoxy)
        self.cxy[:] = self.topoxy < 0.0
        self._topo_prev[:] = self.topoxy

        t0 = time.perf_counter()
        np.multiply(self.deltaslxy, self.cxy, out=self.tmpg)
        sht.analysis(self.tmpg, out=self.slclm)
        sht.analysis(self.cxy, out=self.clm)
        gmslc_ocn = (self.slclm[0, 0] / self.clm[0, 0]).real
        np.multiply(self.tmpg, self.beta, out=self.tmpg)
        sht.analysis(self.tmpg, out=self.slclm)
        np.multiply(self.cxy, self.beta, out=self.cstarxy)
        sht.analysis(self.cstarxy, out=self.cstarlm)
        gmslc_ocnBeta = (self.slclm[0, 0] / self.cstarlm[0, 0]).real
        mark("sht", t0)

        ice_volume = self.icestarlm[0, 0].real * 4 * PI * p.radius**2
        self.icexy_prev = icexy_cur
        rate_out = self._rate_grids(rate_sp[0], rate_sp[1], dice00, rate_sp[2]) if rates else None

        wall = time.perf_counter() - t_all
        if profile:
            tm["total"] = wall
            self.prof = tm
        return {
            "iter": it,
            "numiter": ninner,
            "tgrid": self.topoxy,
            "ocean": self.cxy,
            "beta": self.beta,
            "R": self.rrxy if self.calcRG else None,
            "G": (self.deltaslxy + self.rrxy) if self.calcRG else None,
            "horiz": hz,
            "dS_converged": self.deltaS_cur,
            "ocean_area": self.clm[0, 0].real * 4 * PI * p.radius**2,
            "oceanBeta_area": self.cstarlm[0, 0].real * 4 * PI * p.radius**2,
            "gmslc_ocn": gmslc_ocn,
            "gmslc_ocnBeta": gmslc_ocnBeta,
            "ice_volume": ice_volume,
            "il": self.il.copy(),
            "mm": self.mm.copy(),
            "lambda": self.lam.copy(),
            "rates": rate_out,
            "wall": wall,
            "prof": tm if profile else None,
        }
