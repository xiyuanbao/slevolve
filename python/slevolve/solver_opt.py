"""
solver_opt.py -- optimised step for the exact-memory sea-level solver.

``OptSeaLevelSolver`` subclasses :class:`solver.SeaLevelSolver` and overrides
only ``step``; ``solver.py`` itself is untouched, so the default code path
(``run_ais.py`` with no new flags) still runs the reference implementation
exactly.  Every optimisation below is individually switchable through the
``opt`` set:

``raw``
    the backend returns/consumes ``V[l,m] = (-1)^m U[l,m]`` and folds the legacy
    ``1/(1.41421356 sqrt(2 pi))`` into ducc0's ``ringfactor``.  Removes one
    512x512-complex scaling pass per analysis and one full copy per synthesis.
    The ``(-1)^m`` is absorbed into the solver's own convention: every spectral
    statement in the solver is either elementwise with an *l*-dependent real
    weight, an ``abs()`` reduction, or an ``m = 0`` element, all of which are
    invariant; the only places that are not are the four degree-2, order-1
    rotational entries, which get an exact sign flip.

``diag``
    the four post-loop analyses (``slclm``, ``clm``, ``slclm`` again,
    ``cstarlm``) feed *only* ``[0,0]`` -- the ocean areas and the two GMSL
    numbers, all of them pure diagnostics that no later timestep reads.  The
    degree-0 coefficient of an analysis is a quadrature average, so they are
    replaced by one fused pass of row sums and a 512-long dot product.
    Saves 4 of the step's transforms.

``dsl``
    ``sl_tot_lm`` differs from ``sl_inc_lm`` only by ``conserv`` in the (0,0)
    coefficient, and ``Y_00`` is constant over the sphere, so
    ``deltaslxy == dslxy + conserv * synth(delta_00)`` *identically*.
    Saves one synthesis per inner iteration (``numiter`` per step).

``numba``
    the fused compiled kernels of ``kernels.py`` (hereditary advance, grid
    algebra, spectral algebra), each identical to the NumPy statements it
    replaces.

``sums``
    the convergence reduction.  (a) Carry ``sum(abs(load_inc_lm))`` from one inner
    iteration to the next instead of recomputing it -- ``old_load_inc_lm`` *is* the
    previous ``load_inc_lm`` buffer, so the carried value is the same IEEE number;
    this halves the number of reductions.  (b) Compute the remaining one with
    a four-way row split over a ``ThreadPoolExecutor`` (``np.abs`` releases the
    GIL), then sum the buffer serially: elementwise work split by rows is
    identical, and the pairwise summation of the float buffer is
    untouched.  0.41 ms -> 0.18 ms per reduction.  NumPy 2.5's complex ``abs``
    is a hand-vectorised kernel that agrees with neither ``hypot`` nor
    ``sqrt(x*x+y*y)`` exactly, so it is *not* reimplemented in numba.
"""

from __future__ import annotations

import math
import time

import numpy as np

from .solver import EPS0, PI, ExactMemory, SeaLevelSolver, axial_solve

II = 1j

ALL_OPTS = ("raw", "diag", "dsl", "numba", "sums")
INNER_MODES = ("off", "extrap", "anderson")


class _Anderson:
    """Type-II Anderson mixing (Walker & Ni 2011 form) on a complex array.

    Keeps the last ``depth + 1`` (input x_i, output g_i = G(x_i)) pairs, with
    residual f_i = g_i - x_i.  With dF_j = f_{j+1} - f_j, dG_j = g_{j+1} - g_j
    it solves  min_gamma || f_k - dF gamma ||_2  (real inner product on the
    real/imag parts of the spectral array) and returns
    x_{k+1} = g_k - dG gamma.  If ||f_k|| grew relative to ||f_{k-1}|| the
    history is dropped and the plain step x_{k+1} = g_k is returned (``mix``
    returns "fallback"; "plain" = fewer than two pairs; "singular").
    """

    def __init__(self, depth, shape):
        self.depth = depth
        self.x = [np.zeros(shape, complex) for _ in range(depth + 1)]
        self.g = [np.zeros(shape, complex) for _ in range(depth + 1)]
        self.f = [np.zeros(shape, complex) for _ in range(depth + 1)]
        self.reset()

    def reset(self):
        self.n = 0          # pairs stored
        self.head = -1      # slot of the newest pair
        self.fn = []        # residual norms of the stored pairs, oldest first

    @staticmethod
    def _dot(a, b):
        return float(np.dot(a.view(np.float64).ravel(), b.view(np.float64).ravel()))

    def push(self, x, g):
        self.head = (self.head + 1) % (self.depth + 1)
        h = self.head
        self.x[h][:] = x
        self.g[h][:] = g
        np.subtract(g, x, out=self.f[h])
        self.n = min(self.n + 1, self.depth + 1)
        self.fn.append(math.sqrt(self._dot(self.f[h], self.f[h])))
        self.fn = self.fn[-(self.depth + 1):]

    def _slot(self, i):
        """i = 0 newest, 1 previous, ..."""
        return (self.head - i) % (self.depth + 1)

    def mix(self, out):
        h = self.head
        if self.n >= 2 and self.fn[-1] > self.fn[-2]:
            # residual grew: fall back to plain iteration, restart history
            out[:] = self.g[h]
            gx, gg = self.x[h].copy(), self.g[h].copy()
            self.reset()
            self.push(gx, gg)
            return "fallback"
        if self.n < 2:
            out[:] = self.g[h]
            return "plain"
        m = self.n - 1
        dF = []
        dG = []
        for j in range(m):
            a, b = self._slot(j), self._slot(j + 1)
            dF.append(self.f[a] - self.f[b])
            dG.append(self.g[a] - self.g[b])
        A = np.empty((m, m))
        rhs = np.empty(m)
        for i in range(m):
            rhs[i] = self._dot(dF[i], self.f[h])
            for j in range(i, m):
                A[i, j] = A[j, i] = self._dot(dF[i], dF[j])
        try:
            gam = np.linalg.lstsq(A, rhs, rcond=1e-12)[0]
        except np.linalg.LinAlgError:
            out[:] = self.g[h]
            return "singular"
        if not np.all(np.isfinite(gam)):
            out[:] = self.g[h]
            return "singular"
        out[:] = self.g[h]
        for j in range(m):
            out -= gam[j] * dG[j]
        return "mixed"


def extrap_guess(h1, h2, dt0, dt1, dt2, out=None):
    """Extrapolated inner start of step j (paper Sect. 2.5):

        dS0_j = dt_j [ Sdot_{j-1} + (Sdot_{j-1} - Sdot_{j-2}) (dt_j + dt_{j-1}) / (dt_{j-1} + dt_{j-2}) ],
        Sdot_i = dS_i / dt_i,

    with h1 = dS_{j-1}, h2 = dS_{j-2}, dt0 = dt_j, dt1 = dt_{j-1}, dt2 = dt_{j-2}.  Equal steps give
    exactly 2 dS_{j-1} - dS_{j-2}."""
    if out is None:
        out = np.empty_like(h1)
    if dt0 == dt1 == dt2:
        np.multiply(h1, 2.0, out=out)
        out -= h2
    else:
        w = (dt0 + dt1) / (dt1 + dt2)
        np.multiply(h1, dt0 / dt1 * (1.0 + w), out=out)
        out -= (dt0 / dt2 * w) * h2
    return out


def parse_opts(spec: str) -> set:
    """``"all"`` / ``""`` / ``"raw,dsl"`` -> a set of option names."""
    if spec is None:
        return set()
    spec = spec.strip()
    if spec in ("", "none", "0"):
        return set()
    if spec == "all":
        return set(ALL_OPTS)
    out = set()
    for tok in spec.split(","):
        tok = tok.strip()
        if not tok:
            continue
        if tok not in ALL_OPTS:
            raise ValueError(f"unknown --opt {tok!r} (known: {','.join(ALL_OPTS)})")
        out.add(tok)
    return out


# ------------------------------------------------- hereditary memory, fused
class ExactMemoryFast(ExactMemory):
    """``ExactMemory`` with the advance done by one compiled kernel.

    The state is stored as ``R[l, k, m]`` instead of ``R[k, l, m]`` so that the
    ``kmax`` mode rows belonging to one degree are contiguous (208 kB), and the
    26 separate 4 MB streams of the NumPy version collapse into one pass.  The
    ``k``-accumulation order is unchanged, so the result is identical.
    """

    def __init__(self, love, lmax, dt_years, tpw, nthreads=1):
        super().__init__(love, lmax, dt_years, tpw, nthreads=1)
        nl = lmax + 1
        self.R = None  # replaced by the [l, k, m] layout below
        self.Rlk = np.zeros((nl, self.kmax, nl), dtype=np.complex128)
        self.nmod_l = np.minimum(np.asarray(love.nmod[: nl], dtype=np.int64), self.kmax)
        # contiguous per-mode weight rows (kernel reads them as [k, l])
        self._decay = np.ascontiguousarray(self.decay)
        self._omd = np.ascontiguousarray(self.omd)
        self._a_sl = np.ascontiguousarray(self.a_sl)
        # decay/omd are C-contiguous already: share them so that set_dt reaches the kernel
        assert self._decay is self.decay and self._omd is self.omd
        self._a_rr = np.ascontiguousarray(self.a_rr)

    def advance_and_read(self):
        from . import kernels

        kernels.advance_and_read(
            self.Rlk, self._decay, self._omd, self._a_sl, self._a_rr,
            self.nmod_l, self.p_load, self._viscous, self._radial,
        )
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
        return self._viscous, self._radial, sum_il, sum_m, sum_tsl, sum_trr

    def _ramp_commit_load(self, current_load):
        from . import kernels

        kernels.ramp_commit(self.Rlk, self.omphi, self.nmod_l, current_load)

    def _rate_sums(self):
        from . import kernels

        nl = self.lmax + 1
        out_sl = np.zeros((nl, nl), dtype=np.complex128)
        out_rr = np.zeros((nl, nl), dtype=np.complex128)
        kernels.rate_sums(self.Rlk, self.w_sl, self.w_rr, self.nmod_l, out_sl, out_rr)
        return out_sl, out_rr

    def read_horizontal(self):
        """sum_k a_hl R_k on the [l, k, m] layout (output only; state unchanged)."""
        nl = self.lmax + 1
        hz = np.zeros((nl, nl), dtype=np.complex128)
        for l in range(1, nl):
            nm = int(self.nmod_l[l])
            Rl = self.Rlk[l]
            for k in range(nm):
                hz[l] += self.a_hl[k, l] * Rl[k]
        thl = np.zeros(3, dtype=np.complex128)
        if self.lmax >= 2:
            for k in range(self.nmod2):
                thl += self.a_thl[k] * self.r_lambda[k]
        return hz, thl


# =========================================================== the fast solver
class OptSeaLevelSolver(SeaLevelSolver):
    def __init__(self, *args, opt=(), inner="off", anderson_depth=3,
                 max_mask_refresh=3, anderson_freeze=2, anderson_keep=False,
                 reuse="off", reuse_store=None, **kw):
        self.opt = set(opt)
        # ---- cross-outer-iteration reuse ----------------------------------------
        # reuse_store is a dict owned by the caller that outlives this solver;
        # keys are step numbers.  "warm": start step n's inner loop from the
        # stored converged load_inc_lm_n and ocean mask of the previous outer
        # iteration.  "warm+ice": also reuse the stored ice analysis when the
        # grounded-ice support (icestar != 0) is identical.
        # "warm+corr": warm start plus the pass-to-pass difference of the previous step,
        #   dS_n^(0) = dS_n^(k-1) + (dS_{n-1}^(k) - dS_{n-1}^(k-1));
        # "warm+corr2": the same difference extrapolated linearly, 2 D_{n-1} - D_{n-2}.
        # Initial guess only; the fixed point is unchanged.
        if reuse not in ("off", "warm", "warm+ice", "warm+corr", "warm+corr2"):
            raise ValueError(f"reuse must be off|warm|warm+ice|warm+corr|warm+corr2, got {reuse!r}")
        self.reuse = reuse
        self.reuse_store = reuse_store if reuse_store is not None else {}
        self.reuse_stats = {"warm_steps": 0, "ice_reused": 0, "ice_recomputed": 0}
        if inner not in INNER_MODES:
            raise ValueError(f"inner must be one of {INNER_MODES}, got {inner!r}")
        self.inner = inner
        self.anderson_depth = int(anderson_depth)
        self.max_mask_refresh = int(max_mask_refresh)
        self.anderson_freeze = int(anderson_freeze)
        self.anderson_keep = bool(anderson_keep)
        if "numba" in self.opt:
            kw["mem_threads"] = 1  # the compiled advance does its own threading
        super().__init__(*args, **kw)

        sht = self.sht
        raw = bool(getattr(sht, "raw", False))
        if ("raw" in self.opt) != raw:
            raise ValueError("--opt raw must be used with a raw-mode SHT backend")
        # V[l,m] = (-1)^m U[l,m]: an odd-order rotational entry changes sign
        self._msign = -1.0 if raw else 1.0

        if "numba" in self.opt:
            from . import kernels

            self.state = None  # free the [k, l, m] state before allocating [l, k, m]
            self.state = ExactMemoryFast(
                self.love, self.norder, float(self.dt), self.tpw
            )
            self.state.ramp = self.ramp
            self._k = kernels

        # 1-D copies of the l-weights for the compiled kernels
        self._ttl1 = np.ascontiguousarray(self._ttl[:, 0])
        self._ttlekhl1 = np.ascontiguousarray(self._ttlekhl[:, 0])
        self._ttlhe1 = np.ascontiguousarray(self._ttlhe[:, 0])

        # constant synthesised by a unit (0,0) coefficient -- measured through
        # the backend itself so that it is the number the backend would give
        nl = self.nl
        e00 = np.zeros((nl, nl), dtype=np.complex128)
        e00[0, 0] = 1.0
        self._K00 = float(sht.synthesis(e00, out=np.zeros_like(self.topoxy))[0, 0])

        # degree-0 quadrature row weights for the diagnostic block
        self._w00 = sht.mean00_weights() if hasattr(sht, "mean00_weights") else None
        if "diag" in self.opt and self._w00 is None:
            raise ValueError("--opt diag needs a backend exposing mean00_weights (ducc)")

        # persistent scratch
        self._absbuf = np.empty((nl, nl))
        self._absex = None
        if "sums" in self.opt:
            nthr = max(1, int(getattr(sht, "nthreads", 1)))
            if nthr > 1:
                from concurrent.futures import ThreadPoolExecutor

                self._absex = ThreadPoolExecutor(max_workers=nthr)
                b = np.linspace(0, nl, nthr + 1).astype(int)
                self._absblk = [(int(b[i]), int(b[i + 1])) for i in range(nthr)
                                if b[i + 1] > b[i]]
        self._rows = np.empty((4, self.nlat))
        self._dsl_rot = np.zeros((3, 3), dtype=np.complex128)
        self._rr_rot = np.zeros((3, 3), dtype=np.complex128)
        self._marine_pos_b = np.ascontiguousarray(self._marine_pos)
        self._marine_thr_c = np.ascontiguousarray(self._marine_thr)
        self._4piR2 = 4 * PI * self.p.radius**2

        # ---- inner-loop acceleration state (unused when inner=off) --
        if self.inner != "off":
            self._dS_hist = [z_ for z_ in (np.zeros((nl, nl), complex), np.zeros((nl, nl), complex))]
            self._dS_nhist = 0  # number of converged steps stored (max 2)
            self._dt_hist = [float(self.dt), float(self.dt)]  # dt of steps n-1, n-2
        if self.inner == "anderson":
            self._aa = _Anderson(self.anderson_depth, (nl, nl))
            self._mask_new = np.zeros_like(self.cxy)
        self.inner_stats = {"mixed": 0, "plain": 0, "fallback": 0, "singular": 0,
                            "mask_refresh": 0, "frozen_iters": 0}

    def _abs_rows(self, a):
        """``np.abs(a, out=self._absbuf)`` split by rows over a thread pool."""
        ab = self._absbuf
        list(self._absex.map(lambda b: np.abs(a[b[0]:b[1]], out=ab[b[0]:b[1]]),
                             self._absblk))

    # ---------------------------------------------------------------- step
    def step(self, icexy_cur, profile=False, dt=None, rates=False):
        self.set_step_dt(dt)
        p = self.p
        sht = self.sht
        opt = self.opt
        use_nb = "numba" in opt
        K = self._k if use_nb else None
        self.iter += 1
        it = self.iter
        tm = {}
        t_all = time.perf_counter()
        pc = time.perf_counter

        def mark(key, t0):
            if profile:
                tm[key] = tm.get(key, 0.0) + pc() - t0

        # ---- hereditary advance -----------------------------------------
        t0 = pc()
        viscous, radial, sum_il, sum_m, exact_tsl, exact_trr = self.state.advance_and_read()
        mark("memory", t0)
        deltaS_prev = self.state.p_ocean
        oldil = self.state.p_il.copy()
        oldm = self.state.p_m.copy()
        dlambda_prev2 = self.state.p_lambda.copy()
        ax_v = 0.0
        if self._ax_on and self.tpw and self.norder >= 2:
            # axial rotation option: viscous tidal response to Lambda_20, sum_k a_m,k r_lambda,k (fixed in the step)
            _n2 = self.state.nmod2
            ax_v = float(np.real(self.state.a_m[:_n2] @ self.state.r_lambda[:_n2, 0]))

        rs = self.reuse_store.get(it) if self.reuse != "off" else None
        if rs is not None:
            self.cxy[:] = rs["cxy"]  # converged mask of step n, previous outer

        # ---- ice + ocean grid block (one fused pass) ---------------------
        t0 = pc()
        gr_cm = self.checkmarine
        self._gr_nacc = 0
        if gr_cm:
            # grounded-ice mask, flotation on the converged topography of the previous step
            self._icestar_on(icexy_cur, self._topo_prev, out=self.icestarxy)
            src = self.icestarxy
            self.beta[:] = self.icestarxy < EPS0
            np.multiply(self.cxy, self.beta, out=self.cstarxy)
            np.subtract(self.cstarxy, self.cstar0, out=self.tOxy)
            np.multiply(self.tinit, self.tOxy, out=self.tOxy)
        elif use_nb:
            K.beta_cstar_tO(icexy_cur, self.cxy, self.cstar0, self.tinit, EPS0,
                            self.beta, self.cstarxy, self.tOxy)
            src = icexy_cur
        else:
            src = icexy_cur
            self.beta[:] = icexy_cur < EPS0
            np.multiply(self.cxy, self.beta, out=self.cstarxy)
            np.subtract(self.cstarxy, self.cstar0, out=self.tOxy)
            np.multiply(self.tinit, self.tOxy, out=self.tOxy)
        mark("grid", t0)

        t0 = pc()
        use_ice = self.reuse == "warm+ice"
        gm_cur = (src != 0.0) if use_ice else None
        if use_ice and rs is not None and np.array_equal(gm_cur, rs["gm"]):
            self.icestarlm[:] = rs["icelm"]
            self.reuse_stats["ice_reused"] += 1
        else:
            sht.analysis(src, out=self.icestarlm)
            if use_ice:
                self.reuse_stats["ice_recomputed"] += 1
        if it == 1:
            if self.checkmarine:
                if use_nb:
                    K.icestar(self.icexy_prev, self._marine_pos_b, self._marine_thr_c,
                              self.tmpg)
                else:
                    self._icestar(self.icexy_prev, out=self.tmpg)
                oldsrc = self.tmpg
            else:
                oldsrc = self.icexy_prev
            r0 = self.reuse_store.get(0) if use_ice else None
            gm0 = (oldsrc != 0.0) if use_ice else None
            if r0 is not None and np.array_equal(gm0, r0["gm"]):
                self.oldicestarlm[:] = r0["icelm"]
                self.reuse_stats["ice_reused"] += 1
            else:
                sht.analysis(oldsrc, out=self.oldicestarlm)
                if use_ice:
                    self.reuse_stats["ice_recomputed"] += 1
                    self.reuse_store[0] = {"gm": gm0, "icelm": self.oldicestarlm.copy()}
            oldice = self.oldicestarlm
        else:
            oldice = self.state.prev_ice
        mark("sht", t0)

        t0 = pc()
        if use_nb:
            K.dice_and_deltaice(self.icestarlm, oldice, self.state.p_ice,
                                self.dicestarlm, self.deltaicestarlm)
        else:
            np.subtract(self.icestarlm, oldice, out=self.dicestarlm)
            np.add(self.state.p_ice, self.dicestarlm, out=self.deltaicestarlm)
        mark("spec", t0)
        dicestar = self.dicestarlm
        deltaicestar = self.deltaicestarlm

        t0 = pc()
        sht.analysis(self.cstarxy, out=self.cstarlm)
        sht.analysis(self.tOxy, out=self.tOlm)
        mark("sht", t0)

        t0 = pc()
        # left to NumPy on purpose: complex division rounds differently in LLVM
        np.divide(self.cstarlm, self.cstarlm[0, 0], out=self.load_inc_lm)
        self.load_inc_lm *= (-p.rhoi / p.rhow) * dicestar[0, 0]
        self.deltaS_cur[:] = 0.0
        mark("spec", t0)

        dsl_rot = self._dsl_rot
        rr_rot = self._rr_rot
        dil_cur = np.zeros((3, 3))
        dm_cur = np.zeros(3)
        dlam_cur2 = np.zeros(3, dtype=np.complex128)

        kTE2 = self.love.kTE[2]
        hTE2 = self.love.hTE[2]
        kE2 = self.love.ke[2]
        r2w2 = p.radius**2 * p.omega**2
        sg = self._msign
        cm = self.checkmarine  # noqa: F841  (kept for symmetry with solver.py)
        if self.ramp:
            ttl_csl = np.ascontiguousarray(self._ttl1 * self.state.c_sl)

        # ---- extrapolated initial guess (inner = extrap|anderson) --
        acc = self.inner != "off"
        use_aa = self.inner == "anderson"
        if acc and self._dS_nhist >= 2:
            h1, h2 = self._dS_hist  # dS_{n-1}, dS_{n-2}
            # linear-in-time extrapolation of the increment: for uniform steps
            # 2 dS_{n-1} - dS_{n-2}; for variable steps the rate dS/dt is
            # extrapolated between step midpoints and scaled by dt_n
            extrap_guess(h1, h2, float(self.dt), self._dt_hist[0], self._dt_hist[1],
                         out=self.load_inc_lm)
            # the rotation block of iteration 1 reads deltaS_cur: give it the
            # same guess (plain iteration 1 reads 0 here)
            np.add(deltaS_prev, self.load_inc_lm, out=self.deltaS_cur)
        if rs is not None:
            # warm start: stored converged increment replaces eustatic/extrap guess
            self.load_inc_lm[:] = rs["dS"]
            if self.reuse in ("warm+corr", "warm+corr2"):
                p1 = self.reuse_store.get(it - 1)
                if p1 is not None and "dS_old" in p1:
                    d1 = p1["dS"] - p1["dS_old"]
                    p2 = self.reuse_store.get(it - 2) if self.reuse == "warm+corr2" else None
                    if p2 is not None and "dS_old" in p2:
                        self.load_inc_lm += 2.0 * d1 - (p2["dS"] - p2["dS_old"])
                    else:
                        self.load_inc_lm += d1
                    self.reuse_stats["corr_steps"] = self.reuse_stats.get("corr_steps", 0) + 1
            np.add(deltaS_prev, self.load_inc_lm, out=self.deltaS_cur)
            self.reuse_stats["warm_steps"] += 1
        frozen = False
        refreshing = False
        n_refresh = 0
        if use_aa:
            self._aa.reset()

        # the carried convergence sum: old_load_inc_lm at iteration 1 is this load_inc_lm
        so = np.abs(self.load_inc_lm).sum() if "sums" in opt else None

        ninner = 1
        while True:
            # ------------------------------ rotation ------------------------
            if self.tpw:
                t0 = pc()
                z20 = p.rhow * self.deltaS_cur[2, 0] + p.rhoi * deltaicestar[2, 0]
                z21 = sg * (p.rhow * self.deltaS_cur[2, 1] + p.rhoi * deltaicestar[2, 1])
                self.il[2, 2] = (-8 * PI * p.radius**4 / (3 * math.sqrt(5.0))) * z20.real
                # order 1: the U basis has no Condon-Shortley phase; Milne & Mitrovica (1998, GJI 133)
                # eqs. 33c, 43a,b are written with it (their eq. 33e), so il13, il23 and Lambda21 take the
                # opposite sign here: il_i3 = -int sigma x_i x_3 dA (Mitrovica & Wahr 2011, eq. 10)
                self.il[0, 2] = (-8 * PI * p.radius**4 / (math.sqrt(30.0))) * z21.real
                self.il[1, 2] = (8 * PI * p.radius**4 / (math.sqrt(30.0))) * (-1.0 * II * z21).real
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
                    # axial Liouville equation, sign of Mitrovica & Wahr (2011, AREPS 39, eq. 15): -1/C
                    self.mm[2] = -(1.0 / p.moiC) * ((1 + kE2) * self.il[2, 2] + sil_r[2, 2])
                    if self._ax_on:  # axial rotation option: + rotational I33 (and trace), implicit in m3
                        self.mm[2] = axial_solve(
                            p, self._ax, (1 + kE2) * self.il[2, 2] + sil_r[2, 2],
                            ax_v - st_.c_m * float(np.real(dlambda_prev2[0])), kTE2 + st_.c_m,
                            self.mm[0] ** 2 + self.mm[1] ** 2, self.mm[2])
                else:
                    for i in (0, 1):
                        self.mm[i] = (1 / (1 - kTE2 / p.kf)) * (
                            (1 / (p.moiC - p.moiA)) * ((1 + kE2) * self.il[i, 2] + sum_il[i, 2])
                            + (1 / p.kf) * sum_m[i]
                        )
                    # axial Liouville equation, sign of Mitrovica & Wahr (2011, AREPS 39, eq. 15): -1/C
                    self.mm[2] = -(1.0 / p.moiC) * ((1 + kE2) * self.il[2, 2] + sum_il[2, 2])
                    if self._ax_on:  # axial rotation option: + rotational I33 (and trace), implicit in m3
                        self.mm[2] = axial_solve(
                            p, self._ax, (1 + kE2) * self.il[2, 2] + sum_il[2, 2], ax_v, kTE2,
                            self.mm[0] ** 2 + self.mm[1] ** 2, self.mm[2])
                dm_cur = self.mm - oldm

                m1, m2, m3 = self.mm[0], self.mm[1], self.mm[2]
                self.lam[0, 0] = (r2w2 / 3.0) * ((m1**2 + m2**2 + m3**2) + 2.0 * m3)
                self.lam[2, 0] = (r2w2 / (6.0 * math.sqrt(5.0))) * (m1**2 + m2**2 - 2.0 * m3**2 - 4.0 * m3)
                # order 1 in the U basis (no Condon-Shortley phase; MM98 eq. 33c is written with it, their eq. 33e)
                self.lam[2, 1] = -(r2w2 / math.sqrt(30.0)) * ((m1 - II * m2) * (1.0 + m3))
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
                dsl_rot[2, 1] *= sg
                if self.calcRG:
                    rr_rot[:] = 0.0
                    rr_rot[2, 0:3] = (hTE2 * tot + trr_v) / p.gacc
                    rr_rot[2, 1] *= sg
                mark("rot", t0)
            else:
                dsl_rot[:] = 0.0
                rr_rot[:] = 0.0

            # ------------------------------ sea level (eq. B18) -------------
            t0 = pc()
            if use_nb:
                K.b18(deltaicestar, deltaS_prev, self.load_inc_lm, viscous, p.rhoi, p.rhow,
                      self._ttlekhl1, self._ttl1, self.sl_inc_lm)
            else:
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
                if use_nb:
                    K.ramp_add(self.dicestarlm, self.load_inc_lm, p.rhoi, p.rhow, ttl_csl, self.sl_inc_lm)
                else:
                    self.sl_inc_lm += ttl_csl[:, None] * (p.rhoi * self.dicestarlm + p.rhow * self.load_inc_lm)
            self.sl_inc_lm[0:3, 0:3] += dsl_rot
            mark("spec", t0)

            t0 = pc()
            sht.synthesis(self.sl_inc_lm, out=self.dslxy)
            mark("sht", t0)

            t0 = pc()
            if use_nb:
                K.rO_only(self.dslxy, self.cstarxy, self.rOxy)
            else:
                np.multiply(self.dslxy, self.cstarxy, out=self.rOxy)
            mark("grid", t0)
            t0 = pc()
            sht.analysis(self.rOxy, out=self.rOlm)
            mark("sht", t0)

            t0 = pc()
            conserv = (
                (1 / self.cstarlm[0, 0])
                * (-1.0 * (p.rhoi / p.rhow) * deltaicestar[0, 0] - self.rOlm[0, 0] + self.tOlm[0, 0])
            ).real

            self.old_load_inc_lm, self.load_inc_lm = self.load_inc_lm, self.old_load_inc_lm
            old = self.old_load_inc_lm
            new = self.load_inc_lm
            if use_nb:
                K.dS_update(deltaS_prev, self.rOlm, self.cstarlm, self.tOlm, conserv,
                            new, self.deltaS_cur)
            else:
                np.multiply(deltaS_prev, -1.0, out=new)
                new += self.rOlm
                np.multiply(self.cstarlm, conserv, out=self._sp2)
                new += self._sp2
                new -= self.tOlm
                np.add(deltaS_prev, new, out=self.deltaS_cur)
            mark("spec", t0)

            # ---------------------- deltaslxy -------------------------------
            if "dsl" in opt:
                t0 = pc()
                cadd = conserv * self._K00
                if use_nb:
                    K.add_const(self.dslxy, cadd, self.deltaslxy)
                else:
                    np.add(self.dslxy, cadd, out=self.deltaslxy)
                mark("grid", t0)
            else:
                t0 = pc()
                self.sl_tot_lm[:] = self.sl_inc_lm
                self.sl_tot_lm[0, 0] += conserv
                mark("spec", t0)
                t0 = pc()
                sht.synthesis(self.sl_tot_lm, out=self.deltaslxy)
                mark("sht", t0)

            t0 = pc()
            ab = self._absbuf
            if self._absex is not None:
                self._abs_rows(new)
            else:
                np.abs(new, out=ab)
            sa = ab.sum()
            if "sums" not in opt:
                np.abs(old, out=ab)
                so = ab.sum()
            if abs(sa - so) < EPS0 and abs(so) < EPS0:
                xi = 0.0
            elif abs(so) < EPS0:
                xi = abs((sa - so) / (EPS0 * 10))
            else:
                xi = abs((sa - so) / so)
            so = sa
            mark("spec", t0)

            if use_aa:
                # ---------- Anderson mixing on load_inc_lm, ocean mask frozen ------
                # Iterations 1 and 2 are plain evaluations.  At the end of
                # iteration 1 the plain mask update is done once (cxy,
                # cstarxy, tOxy and their analyses cstarlm, tOlm from
                # iteration 1's deltaslxy) and then FROZEN, i.e. iteration 2
                # onwards run on that mask.  The (input, output) pairs of
                # iterations 1 and 2 seed the history, so the first mixed
                # input is x_3 = AA(pairs 1, 2).  While frozen, cxy/cstarxy/tOxy/cstarlm/tOlm are not
                # recomputed, i.e. each evaluation costs 1 synthesis + 1
                # analysis instead of 1 + 3.  beta, tinit, cstar0 never change
                # inside a step anyway.  Everything else (load_inc_lm, deltaS_cur,
                # the rotation block, conserv, deltaslxy) is evaluated afresh.
                if xi <= self.epsilon1:
                    t0 = pc()
                    # converged under the frozen mask: accept only if the mask
                    # implied by this evaluation's deltaslxy equals it
                    if use_nb:
                        K.topo_ocean(self.tinit, self.deltaslxy, self.topoxy, self._mask_new)
                    else:
                        np.subtract(self.tinit, self.deltaslxy, out=self.topoxy)
                        self._mask_new[:] = self.topoxy < 0.0
                    same = np.array_equal(self._mask_new, self.cxy)
                    mark("grid", t0)
                    if same or n_refresh >= self.max_mask_refresh:
                        break
                    n_refresh += 1
                    self.inner_stats["mask_refresh"] += 1
                    refreshing = True
                    # refresh the mask below (plain update) and re-freeze it
                elif frozen:
                    t0 = pc()
                    self._aa.push(old, new)
                    res = self._aa.mix(out=self.load_inc_lm)  # load_inc_lm := next input x_{k+1}
                    self.inner_stats[res] = self.inner_stats.get(res, 0) + 1
                    self.inner_stats["frozen_iters"] += 1
                    np.add(deltaS_prev, self.load_inc_lm, out=self.deltaS_cur)
                    if "sums" in opt:
                        so = np.abs(self.load_inc_lm).sum()
                    mark("spec", t0)
                    if ninner == 9999:
                        raise RuntimeError(f"inner loop failed to converge, xi={xi}")
                    ninner += 1
                    continue
                # plain step (load_inc_lm already holds G(x_k)) + mask update, then freeze
                if refreshing:
                    if not self.anderson_keep:
                        self._aa.reset()
                    refreshing = False
                    frozen = True
                elif ninner >= self.anderson_freeze:
                    self._aa.reset()
                    self._aa.push(old, new)  # this pair seeds the history
                    frozen = True
                t0 = pc()
                if gr_cm:
                    np.subtract(self.tinit, self.deltaslxy, out=self.topoxy)
                    self.cxy[:] = self.topoxy < 0.0
                    self._grounded_reguess(icexy_cur, oldice)
                    np.multiply(self.cxy, self.beta, out=self.cstarxy)
                    np.subtract(self.cstarxy, self.cstar0, out=self.tOxy)
                    np.multiply(self.tinit, self.tOxy, out=self.tOxy)
                elif use_nb:
                    K.topo_ocean_cstar_tO(self.tinit, self.deltaslxy, self.beta,
                                          self.cstar0, self.topoxy, self.cxy,
                                          self.cstarxy, self.tOxy)
                else:
                    np.subtract(self.tinit, self.deltaslxy, out=self.topoxy)
                    self.cxy[:] = self.topoxy < 0.0
                    np.multiply(self.cxy, self.beta, out=self.cstarxy)
                    np.subtract(self.cstarxy, self.cstar0, out=self.tOxy)
                    np.multiply(self.tinit, self.tOxy, out=self.tOxy)
                mark("grid", t0)
                t0 = pc()
                sht.analysis(self.cstarxy, out=self.cstarlm)
                sht.analysis(self.tOxy, out=self.tOlm)
                mark("sht", t0)
                if ninner == 9999:
                    raise RuntimeError(f"inner loop failed to converge, xi={xi}")
                ninner += 1
                continue

            if xi > self.epsilon1:
                t0 = pc()
                if gr_cm:
                    np.subtract(self.tinit, self.deltaslxy, out=self.topoxy)
                    self.cxy[:] = self.topoxy < 0.0
                    self._grounded_reguess(icexy_cur, oldice)
                    np.multiply(self.cxy, self.beta, out=self.cstarxy)
                    np.subtract(self.cstarxy, self.cstar0, out=self.tOxy)
                    np.multiply(self.tinit, self.tOxy, out=self.tOxy)
                elif use_nb:
                    K.topo_ocean_cstar_tO(self.tinit, self.deltaslxy, self.beta,
                                          self.cstar0, self.topoxy, self.cxy,
                                          self.cstarxy, self.tOxy)
                else:
                    np.subtract(self.tinit, self.deltaslxy, out=self.topoxy)
                    self.cxy[:] = self.topoxy < 0.0
                    np.multiply(self.cxy, self.beta, out=self.cstarxy)
                    np.subtract(self.cstarxy, self.cstar0, out=self.tOxy)
                    np.multiply(self.tinit, self.tOxy, out=self.tOxy)
                mark("grid", t0)
                t0 = pc()
                sht.analysis(self.cstarxy, out=self.cstarlm)
                sht.analysis(self.tOxy, out=self.tOlm)
                mark("sht", t0)
            else:
                break
            if ninner == 9999:
                raise RuntimeError(f"inner loop failed to converge, xi={xi}")
            ninner += 1

        if acc:
            h1, h2 = self._dS_hist
            h2[:] = self.load_inc_lm  # becomes dS_{n-1}; old dS_{n-1} becomes dS_{n-2}
            self._dS_hist = [h2, h1]
            self._dt_hist = [float(self.dt), self._dt_hist[0]]
            self._dS_nhist = min(2, self._dS_nhist + 1)

        # ------------------------------------------------------------ R / G
        if self.calcRG:
            t0 = pc()
            if use_nb:
                K.b18(deltaicestar, deltaS_prev, self.load_inc_lm, radial, p.rhoi, p.rhow,
                      self._ttlhe1, self._ttl1, self.rrlm)
            else:
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
                ttl_crr = np.ascontiguousarray(self._ttl1 * self.state.c_rr)
                if use_nb:
                    K.ramp_add(self.dicestarlm, self.load_inc_lm, p.rhoi, p.rhow, ttl_crr, self.rrlm)
                else:
                    self.rrlm += ttl_crr[:, None] * (p.rhoi * self.dicestarlm + p.rhow * self.load_inc_lm)
            self.rrlm[0:3, 0:3] += rr_rot
            mark("spec", t0)
            t0 = pc()
            sht.synthesis(self.rrlm, out=self.rrxy)
            mark("sht", t0)

        hz = None
        if self.horizontal:
            hz = self.horizontal_lm(deltaicestar, deltaS_prev, self.load_inc_lm,
                                    dlambda_prev2 + dlam_cur2 if self.tpw else None)

        # -------------------------------------------------------- commit ----
        t0 = pc()
        if use_nb:
            K.current_load(dicestar, self.load_inc_lm, p.rhoi, p.rhow, self.current_load)
        else:
            np.multiply(dicestar, p.rhoi, out=self.current_load)
            np.multiply(self.load_inc_lm, p.rhow, out=self._sp2)
            self.current_load += self._sp2
            self.current_load[0, :] = 0.0
        rate_sp = None
        if rates:
            rate_sp = self._rate_spectra(self.current_load, oldil, oldm, dlambda_prev2)
            dice00 = complex(dicestar[0, 0])
        if self.tpw:
            self.state.commit(self.current_load, dil_cur, dm_cur, dlam_cur2)
        else:
            self.state.commit(self.current_load)
        if use_nb:
            K.copy3(deltaicestar, self.state.p_ice, self.deltaS_cur, self.state.p_ocean,
                    self.icestarlm, self.state.prev_ice)
        else:
            self.state.p_ice[:] = deltaicestar
            self.state.p_ocean[:] = self.deltaS_cur
            self.state.prev_ice[:] = self.icestarlm
        mark("memory", t0)

        # ------------------------------------------------- topography/ocean -
        t0 = pc()
        if use_nb:
            K.topo_ocean(self.tinit, self.deltaslxy, self.topoxy, self.cxy)
        else:
            np.subtract(self.tinit, self.deltaslxy, out=self.topoxy)
            self.cxy[:] = self.topoxy < 0.0
        self._topo_prev[:] = self.topoxy
        mark("grid", t0)

        # --------------------------------------------------- GMSL / areas ---
        if "diag" in opt:
            t0 = pc()
            r = self._rows
            if use_nb:
                K.gmsl_rowsums(self.deltaslxy, self.cxy, self.beta, r[0], r[1], r[2], r[3])
            else:
                np.multiply(self.deltaslxy, self.cxy, out=self.tmpg)
                r[0] = self.tmpg.sum(axis=1)
                r[1] = self.cxy.sum(axis=1)
                np.multiply(self.tmpg, self.beta, out=self.tmpg)
                r[2] = self.tmpg.sum(axis=1)
                np.multiply(self.cxy, self.beta, out=self.tmpg)
                r[3] = self.tmpg.sum(axis=1)
            w = self._w00
            slc_ocn00 = float(w @ r[0])
            c00 = float(w @ r[1])
            slc_beta00 = float(w @ r[2])
            cstar00 = float(w @ r[3])
            gmslc_ocn = slc_ocn00 / c00
            gmslc_ocnBeta = slc_beta00 / cstar00
            ocean_area = c00 * self._4piR2
            oceanBeta_area = cstar00 * self._4piR2
            mark("grid", t0)
        else:
            t0 = pc()
            np.multiply(self.deltaslxy, self.cxy, out=self.tmpg)
            sht.analysis(self.tmpg, out=self.slclm)
            sht.analysis(self.cxy, out=self.clm)
            gmslc_ocn = (self.slclm[0, 0] / self.clm[0, 0]).real
            np.multiply(self.tmpg, self.beta, out=self.tmpg)
            sht.analysis(self.tmpg, out=self.slclm)
            np.multiply(self.cxy, self.beta, out=self.cstarxy)
            sht.analysis(self.cstarxy, out=self.cstarlm)
            gmslc_ocnBeta = (self.slclm[0, 0] / self.cstarlm[0, 0]).real
            ocean_area = self.clm[0, 0].real * self._4piR2
            oceanBeta_area = self.cstarlm[0, 0].real * self._4piR2
            mark("sht", t0)

        if self.reuse != "off":
            ent = self.reuse_store.get(it)
            if ent is None:
                ent = self.reuse_store[it] = {"dS": np.empty_like(self.load_inc_lm),
                                              "cxy": np.empty(self.cxy.shape, bool)}
            elif self.reuse in ("warm+corr", "warm+corr2"):
                if "dS_old" not in ent:
                    ent["dS_old"] = np.empty_like(self.load_inc_lm)
                ent["dS_old"][:] = ent["dS"]  # converged dS of this step in the previous pass
            ent["dS"][:] = self.load_inc_lm
            np.not_equal(self.cxy, 0.0, out=ent["cxy"])  # 0/1 mask stored as bool
            if use_ice:
                ent["gm"] = (src != 0.0)  # I* may have changed in the inner loop
                if "icelm" not in ent:
                    ent["icelm"] = np.empty_like(self.icestarlm)
                ent["icelm"][:] = self.icestarlm

        ice_volume = self.icestarlm[0, 0].real * self._4piR2
        self.icexy_prev = icexy_cur
        rate_out = self._rate_grids(rate_sp[0], rate_sp[1], dice00, rate_sp[2]) if rates else None

        wall = pc() - t_all
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
            "ocean_area": ocean_area,
            "oceanBeta_area": oceanBeta_area,
            "gmslc_ocn": gmslc_ocn,
            "gmslc_ocnBeta": gmslc_ocnBeta,
            "ice_volume": ice_volume,
            "il": self.il.copy(),
            "mm": self.mm.copy(),
            "lambda": self.lam.copy(),
            "rates": rate_out,
            "wall": wall,
            "prof": tm if profile else None,
            "inner_stats": dict(self.inner_stats),
        }
