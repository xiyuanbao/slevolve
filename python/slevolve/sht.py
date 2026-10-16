"""
sht.py -- NumPy/BLAS spherical-harmonic transform that reproduces, to round-off,
the conventions of the Fortran ``spharmt`` module of the legacy code of Han et al. (2022)
(with the corrected ``fft99b``; ``spat2spec`` / ``spec2spat``).

Conventions (established empirically against the Fortran)
----------------------------------------------------------------------------
Grid
    ``grid[j, i]``, ``j = 0..nlat-1`` Gauss-Legendre colatitude rows running
    north -> south (``mu_0 = +0.99998...``), ``i = 0..nlon-1`` longitude,
    longitude varying fastest in memory.  This is exactly the linear element
    order of the Fortran arrays and of the one-value-per-line text files.

Legendre functions
    ``P~_lm(mu)`` orthonormal on [-1,1] (``int P~_lm^2 dmu = 1``), **no**
    Condon-Shortley phase:  ``P~_00 = 1/sqrt(2)``, ``P~_11 = sqrt(3)/2 sin(th)``.
    Built with the exact recursion of Fortran ``LEGEND``.

Analysis  (``spat2spec``)
    ``am[j, m] = rfft(grid, axis=1)[j, m] / nlon``          (Temperton fft991, isign=-1)
    ``U[l, m] = ( sum_j P~_lm(mu_j) w_j am[j, m] ) / 1.41421356``

Synthesis (``spec2spat``)
    ``am[j, m] = sum_l ( U[l, m] * 1.41421356 ) P~_lm(mu_j)``
    ``grid     = irfft(am * nlon, n=nlon, axis=1)``

Note ``1.41421356`` is the *truncated* literal that appears in the legacy source, not
``sqrt(2)``; it is reproduced verbatim (it cancels exactly in a round trip but
sets the absolute scale of the coefficients, which the physics uses).
Analysis of a constant field ``c`` therefore returns
``U[0,0] = c * sqrt(2)/1.41421356 = c * (1 + 1.7e-9)``.

Implementation
--------------
* Longitude: ``scipy.fft.rfft`` / ``irfft`` (pocketfft, ``workers=``).
* Legendre: BLAS.  The ``P~`` table is stored only for the northern half of the
  grid and split by the parity of ``l-m`` (``P~_lm(-mu) = (-1)^(l+m) P~_lm(mu)``),
  which halves both flops and memory traffic.  Blocks of ``chunk`` consecutive
  ``m`` are padded to a common length so the whole Legendre stage is a handful
  of *batched* ``numpy.matmul`` calls with a right-hand side of width 1
  (the shape for which this OpenBLAS build reaches full memory bandwidth).

Interface (fixed, so a GPU backend can be dropped in):
    ``SHT(nlat, nlon, lmax)``
    ``analysis(grid[nlat, nlon] float64) -> coeffs[lmax+1, lmax+1] complex128``  ([l, m], zero above the diagonal)
    ``synthesis(coeffs) -> grid``
"""

from __future__ import annotations

import numpy as np

try:
    import scipy.fft as _sfft

    _HAVE_SCIPY_FFT = True
except Exception:  # pragma: no cover
    _HAVE_SCIPY_FFT = False

# The literal that appears in the legacy spat2spec / spec2spat.  NOT sqrt(2).
SQRT2_LEGACY = 1.41421356


# ----------------------------------------------------------------------------- gaulw
def gaulw(nlat: int):
    """Gauss-Legendre nodes/weights, faithful port of ``spharmt::gaulw``.

    Returns ``(sinlats, wts)``; ``sinlats[0]`` is the northernmost node (+1 side).
    """
    nprec = 15  # precision(real64)
    converg = 0.1**nprec
    itermax = 10
    pi = 4.0 * np.arctan(1.0)

    z = np.cos(pi * (np.arange(1, nlat + 1, dtype=np.float64) - 0.25) / (nlat + 0.5))
    z_out = np.zeros(nlat)
    pp_out = np.zeros(nlat)
    done = np.zeros(nlat, dtype=bool)

    for _ in range(itermax):
        p1 = np.ones(nlat)
        p2 = np.zeros(nlat)
        for j in range(1, nlat + 1):
            p3 = p2
            p2 = p1
            p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j
        pp = nlat * (z * p1 - p2) / (z * z - 1.0)
        z1 = z
        z = z1 - p1 / pp
        conv = np.abs(z - z1) < converg
        newly = conv & ~done
        z_out[newly] = z[newly]
        pp_out[newly] = pp[newly]
        done |= conv
        if done.all():
            break
    if not done.all():
        raise RuntimeError("gaulw: abscissas failed to converge")

    wts = 2.0 / ((1.0 - z_out * z_out) * pp_out * pp_out)
    return z_out, wts


# ---------------------------------------------------------------------------- legend
def legend_columns(x: np.ndarray, ntrunc: int):
    """Yield ``(m, Pm)`` with ``Pm[l-m, j] = P~_lm(x_j)``, exact port of ``LEGEND``.

    ``x`` is a 1-D float64 array of ``mu`` values.
    """
    nj = x.size
    nmdim = (ntrunc + 1) * (ntrunc + 2) // 2
    sinsq = 1.0 - x * x
    A = 1.0
    B = 0.0
    prod = np.ones(nj)
    nmstrt = 0  # 0-based offset of column m in the packed vector

    for i in range(1, ntrunc + 2):
        m = i - 1
        nmax = ntrunc + 1 - m
        if m != 0:
            A = A + 2.0
            B = B + 2.0
            prod = prod * sinsq * A / B
        Pm = np.zeros((nmax, nj))
        pmnjm2 = np.sqrt(0.5 * prod)
        Pm[0] = pmnjm2
        pmnjm1 = np.sqrt(float(2 * m + 3)) * x * pmnjm2
        if nmstrt + 1 != nmdim:  # Fortran: IF (NM .NE. NMDIM) PMN(NM+1) = PMNJM1
            if nmax >= 2:
                Pm[1] = pmnjm1
        np1 = m + 1
        epsnmp1 = np.sqrt(float(np1 * np1 - m * m) / float(4 * np1 * np1 - 1))
        epspmn = x * pmnjm1 - epsnmp1 * pmnjm2
        for j in range(3, nmax + 1):
            n = m + j - 1
            eps = np.sqrt(float(n * n - m * m) / float(4 * n * n - 1))
            pmnj = epspmn / eps
            Pm[j - 1] = pmnj
            epspmn = x * pmnj - eps * pmnjm1
            pmnjm2 = pmnjm1
            pmnjm1 = pmnj
        nmstrt += nmax
        yield m, Pm


# ------------------------------------------------------------------------------- SHT
class SHT:
    """Pseudo-spectral transform on a Gauss-Legendre grid, legacy conventions."""

    def __init__(self, nlat: int, nlon: int, lmax: int, workers: int = 1, chunk: int = 16):
        if nlat % 2 != 0:
            raise ValueError("nlat must be even")
        if lmax > nlon // 2:
            raise ValueError("lmax must be <= nlon/2")
        self.nlat = int(nlat)
        self.nlon = int(nlon)
        self.lmax = int(lmax)
        self.nm = self.lmax + 1
        self.nh = self.nlat // 2
        self.workers = int(workers)
        self.chunk = int(chunk)
        self.nfreq = self.nlon // 2 + 1

        self.gaulats, self.weights = gaulw(self.nlat)

        nh, nm, nlat_, nlon_ = self.nh, self.nm, self.nlat, self.nlon
        # weights folded with the exact power-of-two 1/nlon of the FFT scaling
        self._wn = (self.weights[:nh] / nlon_).copy()
        self._ws = (self.weights[nlat_ - 1 : nh - 1 : -1] / nlon_).copy()

        xn = self.gaulats[:nh].copy()

        # ---- chunked, parity-split, padded Legendre tables ----------------------
        self._blocks = []  # (m0, cc, nemax, nomax, Pe, Po, PeT, PoT)
        off_e = 0
        off_o = 0
        offs = []
        for m0 in range(0, nm, self.chunk):
            cc = min(self.chunk, nm - m0)
            nl0 = nm - m0
            nemax = (nl0 + 1) // 2
            nomax = nl0 // 2
            offs.append((m0, cc, nemax, nomax, off_e, off_o))
            off_e += cc * nemax
            off_o += cc * nomax
        self._tot_e, self._tot_o = off_e, off_o

        Pe_all = [np.zeros((cc, nemax, nh)) for (_, cc, nemax, _, _, _) in offs]
        Po_all = [np.zeros((cc, nomax, nh)) for (_, cc, _, nomax, _, _) in offs]

        src_e, dst_e, src_o, dst_o = [], [], [], []
        for m, Pm in legend_columns(xn, self.lmax):
            ci = m // self.chunk
            i = m - offs[ci][0]
            nemax, nomax, oe, oo = offs[ci][2], offs[ci][3], offs[ci][4], offs[ci][5]
            ne = (Pm.shape[0] + 1) // 2
            no = Pm.shape[0] // 2
            Pe_all[ci][i, :ne] = Pm[0::2]
            Po_all[ci][i, :no] = Pm[1::2]
            for k in range(ne):
                src_e.append(oe + i * nemax + k)
                dst_e.append((m + 2 * k) * nm + m)
            for k in range(no):
                src_o.append(oo + i * nomax + k)
                dst_o.append((m + 2 * k + 1) * nm + m)

        self.src_e = np.asarray(src_e, dtype=np.intp)
        self.dst_e = np.asarray(dst_e, dtype=np.intp)
        self.src_o = np.asarray(src_o, dtype=np.intp)
        self.dst_o = np.asarray(dst_o, dtype=np.intp)

        for (m0, cc, nemax, nomax, oe, oo), Pe, Po in zip(offs, Pe_all, Po_all):
            PeT = np.ascontiguousarray(np.swapaxes(Pe, 1, 2))
            PoT = np.ascontiguousarray(np.swapaxes(Po, 1, 2))
            self._blocks.append((m0, cc, nemax, nomax, oe, oo, Pe, Po, PeT, PoT))

        # ---- persistent work buffers -------------------------------------------
        self._BEre = np.zeros(self._tot_e)
        self._BEim = np.zeros(self._tot_e)
        self._BOre = np.zeros(self._tot_o)
        self._BOim = np.zeros(self._tot_o)
        self._ve = [
            (
                self._BEre[oe : oe + cc * nemax].reshape(cc, nemax, 1),
                self._BEim[oe : oe + cc * nemax].reshape(cc, nemax, 1),
                self._BOre[oo : oo + cc * nomax].reshape(cc, nomax, 1),
                self._BOim[oo : oo + cc * nomax].reshape(cc, nomax, 1),
            )
            for (m0, cc, nemax, nomax, oe, oo, *_r) in self._blocks
        ]
        # analysis right-hand sides, (nm, nh)
        self._rs_re = np.zeros((nm, nh))
        self._rs_im = np.zeros((nm, nh))
        self._rd_re = np.zeros((nm, nh))
        self._rd_im = np.zeros((nm, nh))
        # synthesis per-m half-grid accumulators, (nm, nh)
        self._SEre = np.zeros((nm, nh))
        self._SEim = np.zeros((nm, nh))
        self._SOre = np.zeros((nm, nh))
        self._SOim = np.zeros((nm, nh))
        self._sv = [
            (
                self._SEre[m0 : m0 + cc].reshape(cc, nh, 1),
                self._SEim[m0 : m0 + cc].reshape(cc, nh, 1),
                self._SOre[m0 : m0 + cc].reshape(cc, nh, 1),
                self._SOim[m0 : m0 + cc].reshape(cc, nh, 1),
            )
            for (m0, cc, *_r) in self._blocks
        ]
        self._tmpA = np.zeros((nh, nm), dtype=np.complex128)
        self._tmpB = np.zeros((nh, nm), dtype=np.complex128)
        self._Fout = np.zeros((nlat_, self.nfreq), dtype=np.complex128)
        self._coef_out = np.zeros((nm, nm), dtype=np.complex128)
        self._grid_out = np.zeros((nlat_, nlon_))

    # -- byte size of the precomputed tables ---------------------------------
    @property
    def table_bytes(self) -> int:
        return sum(b[6].nbytes + b[7].nbytes + b[8].nbytes + b[9].nbytes for b in self._blocks)

    # ------------------------------------------------------------------ FFT
    def _rfft(self, grid):
        if _HAVE_SCIPY_FFT:
            return _sfft.rfft(grid, axis=1, workers=self.workers)
        return np.fft.rfft(grid, axis=1)

    def _irfft(self, spec, out=None):
        if _HAVE_SCIPY_FFT:
            return _sfft.irfft(spec, n=self.nlon, axis=1, workers=self.workers)
        return np.fft.irfft(spec, n=self.nlon, axis=1)

    # ------------------------------------------------------------- analysis
    def analysis(self, grid: np.ndarray, out: np.ndarray | None = None) -> np.ndarray:
        """grid[nlat, nlon] float64 -> coeffs[lmax+1, lmax+1] complex128 ([l, m])."""
        nh, nm, nlat = self.nh, self.nm, self.nlat
        F = self._rfft(grid)
        np.multiply(F[:nh, :nm], self._wn[:, None], out=self._tmpA)
        np.multiply(F[nlat - 1 : nh - 1 : -1, :nm], self._ws[:, None], out=self._tmpB)
        S = self._tmpA + self._tmpB
        D = self._tmpA
        np.subtract(self._tmpA, self._tmpB, out=D)

        self._rs_re[:] = S.real.T
        self._rs_im[:] = S.imag.T
        self._rd_re[:] = D.real.T
        self._rd_im[:] = D.imag.T

        for (m0, cc, nemax, nomax, oe, oo, Pe, Po, PeT, PoT), (be_r, be_i, bo_r, bo_i) in zip(
            self._blocks, self._ve
        ):
            np.matmul(Pe, self._rs_re[m0 : m0 + cc, :, None], out=be_r)
            np.matmul(Pe, self._rs_im[m0 : m0 + cc, :, None], out=be_i)
            np.matmul(Po, self._rd_re[m0 : m0 + cc, :, None], out=bo_r)
            np.matmul(Po, self._rd_im[m0 : m0 + cc, :, None], out=bo_i)

        if out is None:
            out = self._coef_out
        flat = out.reshape(-1)
        flat[self.dst_e] = (self._BEre[self.src_e] + 1j * self._BEim[self.src_e]) / SQRT2_LEGACY
        flat[self.dst_o] = (self._BOre[self.src_o] + 1j * self._BOim[self.src_o]) / SQRT2_LEGACY
        return out

    # ------------------------------------------------------------ synthesis
    def synthesis(self, coeffs: np.ndarray, out: np.ndarray | None = None) -> np.ndarray:
        """coeffs[lmax+1, lmax+1] complex128 ([l, m]) -> grid[nlat, nlon] float64."""
        nh, nm, nlat = self.nh, self.nm, self.nlat
        # SQRT2_LEGACY * nlon is exact (nlon a power of two), so folding the FFT
        # normalisation in here changes nothing in floating point.
        scale = SQRT2_LEGACY * self.nlon
        cf = coeffs.reshape(-1)
        ge = cf[self.dst_e] * scale
        go = cf[self.dst_o] * scale
        self._BEre[self.src_e] = ge.real
        self._BEim[self.src_e] = ge.imag
        self._BOre[self.src_o] = go.real
        self._BOim[self.src_o] = go.imag

        for (m0, cc, nemax, nomax, oe, oo, Pe, Po, PeT, PoT), (be_r, be_i, bo_r, bo_i), (
            se_r,
            se_i,
            so_r,
            so_i,
        ) in zip(self._blocks, self._ve, self._sv):
            np.matmul(PeT, be_r, out=se_r)
            np.matmul(PeT, be_i, out=se_i)
            np.matmul(PoT, bo_r, out=so_r)
            np.matmul(PoT, bo_i, out=so_i)

        E = self._SEre + 1j * self._SEim  # (nm, nh)
        O = self._SOre + 1j * self._SOim
        north = (E + O).T  # (nh, nm)
        south = (E - O).T
        Fo = self._Fout
        Fo[:nh, :nm] = north
        Fo[nlat - 1 : nh - 1 : -1, :nm] = south
        if nm < self.nfreq:
            Fo[:, nm:] = 0.0
        g = self._irfft(Fo)
        if out is not None:
            out[...] = g
            return out
        return g


__all__ = ["SHT", "gaulw", "legend_columns", "SQRT2_LEGACY"]
