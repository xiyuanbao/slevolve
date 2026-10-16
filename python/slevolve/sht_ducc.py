"""
sht_ducc.py -- ducc0-backed spherical-harmonic transform with *exactly* the
interface and the conventions of ``python/slevolve/sht.py`` (the legacy ``spharmt``).

Drop-in replacement for :class:`sht.SHT`:

    DuccSHT(nlat, nlon, lmax, workers=1)
        .analysis(grid[nlat, nlon] float64, out=None) -> coeffs[lmax+1, lmax+1] complex128
        .synthesis(coeffs, out=None)                  -> grid[nlat, nlon] float64

Convention mapping (derived empirically)
--------------------------------
``ducc0.sht.analysis_2d`` with ``geometry="GL"`` returns the coefficients of
the *standard, orthonormal, Condon-Shortley* complex spherical harmonics

    a_lm = \\int f(\\Omega) Y_lm^*(\\Omega) d\\Omega ,
    Y_lm = (-1)^m (1/sqrt(2 pi)) Ptilde_lm(mu) e^{i m phi} ,
    \\int |Y_lm|^2 d\\Omega = 1 ,

with ring 0 at the *northernmost* Gauss-Legendre colatitude -- i.e. the same
north->south row order as the legacy grid and as ``topo0_archive``.  the legacy
``spat2spec`` instead returns

    U[l,m] = (1/1.41421356) \\int dmu (1/2pi) \\int dphi  Ptilde_lm(mu) f e^{-i m phi}
           = (1/(1.41421356 sqrt(2 pi))) (-1)^m a_lm .

Hence the *only* per-coefficient factor is

    U[l,m] = LEGACY_OVER_DUCC * (-1)^m * a_lm,
    LEGACY_OVER_DUCC = 1 / (1.41421356 * sqrt(2 pi)) = 0.28209479224724215

(compare 1/sqrt(4 pi) = 0.28209479177387814; the 1.7e-9 offset is the legacy
truncated ``1.41421356`` literal, reproduced verbatim), and synthesis uses the
reciprocal.  Verified against the Fortran transform to <= 1e-13 relative.

Storage
-------
``ducc0`` addresses ``alm`` through ``mstart``/``lstride``, so the dense
``[l, m]`` array that the Python solver uses *is* a valid ducc0 alm layout:
``a_lm`` lives at flat index ``l*(lmax+1) + m``, i.e. ``mstart[m] = m`` and
``lstride = lmax+1``.  The backend therefore hands ducc0 the caller's own
buffer and needs **no** packing/unpacking gather-scatter; only a single
elementwise multiply by the per-``m`` factor row.

Restrictions (vs. the NumPy backend)
------------------------------------
* exact Gauss-Legendre analysis needs ``lmax <= nlat - 1`` -> ``lmax = 511``
  on the 512 x 1024 production grid (degree 512 is aliased anyway;
  the legacy grid does not resolve it);
* ``mmax <= (nlon-1)//2 = 511``, i.e. the Nyquist order ``m = 512`` that the
  NumPy backend retains is simply absent.  With ``lmax = mmax = 511`` both
  backends agree exactly on which orders exist.
"""

from __future__ import annotations

import numpy as np

try:  # keep the NumPy backend importable in environments without ducc0
    import ducc0.sht as _dsht
except ImportError:  # pragma: no cover
    _dsht = None

from .sht import SQRT2_LEGACY

#: U[l,m] = LEGACY_OVER_DUCC * (-1)**m * alm[l,m]
LEGACY_OVER_DUCC = 1.0 / (SQRT2_LEGACY * np.sqrt(2.0 * np.pi))


class DuccSHT:
    """ducc0 ``analysis_2d``/``synthesis_2d`` in the legacy ``spharmt`` conventions."""

    backend = "ducc"

    def __init__(self, nlat: int, nlon: int, lmax: int, workers: int = 1, chunk: int = 16,
                 raw: bool = False):
        if _dsht is None:
            raise ImportError("ducc0 is not installed in this interpreter (pip install ducc0)")
        if lmax > nlat - 1:
            raise ValueError(
                f"ducc0 GL analysis requires lmax <= nlat-1 ({nlat - 1}); got lmax={lmax}"
            )
        if lmax > (nlon - 1) // 2:
            raise ValueError(f"ducc0 requires mmax <= (nlon-1)//2 = {(nlon - 1) // 2}")
        self.nlat = int(nlat)
        self.nlon = int(nlon)
        self.lmax = int(lmax)
        self.nm = self.lmax + 1
        self.workers = int(workers)
        self.nthreads = int(workers)

        nl = self.nm
        # dense [l, m] addressing: a_lm at flat index l*nl + m
        self._mstart = np.arange(nl, dtype=np.uint64)
        self._lstride = nl

        sgn = np.where(np.arange(nl) % 2 == 0, 1.0, -1.0)
        self._fac = (LEGACY_OVER_DUCC * sgn)[None, :]  # (1, nl), applied along m
        self._ifac = (1.0 / (LEGACY_OVER_DUCC * sgn))[None, :]

        # ---- "raw" (sign-free) mode -------------------------------------
        # In raw mode the backend returns  V[l,m] = LEGACY_OVER_DUCC * a_lm
        #                                        = (-1)^m U[l,m]
        # and synthesises from that same convention.  The constant is applied
        # by ducc0 itself through ``ringfactor`` (fused into the ring loop, no
        # extra pass over the 512x512 spectral array) and the (-1)^m sign is
        # simply *not applied* -- it is absorbed into the caller's convention
        # This removes both the
        # 512x512-complex scale pass of ``analysis`` and the 512x512-complex
        # copy of ``synthesis``.
        self.raw = bool(raw)
        self._rf_a = np.full(self.nlat, LEGACY_OVER_DUCC, dtype=np.float64)
        self._rf_s = np.full(self.nlat, 1.0 / LEGACY_OVER_DUCC, dtype=np.float64)

        # persistent buffers
        self._coef_out = np.zeros((nl, nl), dtype=np.complex128)
        self._alm_s = np.zeros((nl, nl), dtype=np.complex128)  # synthesis scratch
        self._grid_out = np.zeros((self.nlat, self.nlon), dtype=np.float64)

        # Diagnostics only (north -> south).  ducc0 computes its own quadrature
        # internally; these are NumPy's GL nodes/weights, exposed for
        # verification against the NumPy backend.
        x, w = np.polynomial.legendre.leggauss(self.nlat)
        self.gaulats = x[::-1].copy()
        self.weights = w[::-1].copy()

    # ------------------------------------------------------------------ misc
    @property
    def table_bytes(self) -> int:
        """ducc0 builds its Legendre recurrence on the fly; no persistent table."""
        return 0

    # -------------------------------------------------------------- analysis
    def analysis(self, grid: np.ndarray, out: np.ndarray | None = None) -> np.ndarray:
        """grid[nlat, nlon] float64 -> coeffs[lmax+1, lmax+1] complex128 ([l, m])."""
        if out is None:
            out = self._coef_out
        g = np.ascontiguousarray(grid, dtype=np.float64).reshape(1, self.nlat, self.nlon)
        _dsht.analysis_2d(
            map=g,
            alm=out.reshape(1, -1),
            spin=0,
            lmax=self.lmax,
            mmax=self.lmax,
            geometry="GL",
            mstart=self._mstart,
            lstride=self._lstride,
            nthreads=self.nthreads,
            **({"ringfactor": self._rf_a} if self.raw else {}),
        )
        if not self.raw:
            np.multiply(out, self._fac, out=out)
        return out

    # ------------------------------------------------------------- synthesis
    def synthesis(self, coeffs: np.ndarray, out: np.ndarray | None = None) -> np.ndarray:
        """coeffs[lmax+1, lmax+1] complex128 ([l, m]) -> grid[nlat, nlon] float64."""
        if self.raw:
            # coefficients are already  V = LEGACY_OVER_DUCC * a_lm ; the 1/LEGACY
            # rescale rides on ducc0's ``ringfactor`` (applied to the output
            # map), so no copy of the spectral array is made at all.  Im a_l0
            # is discarded by the c2r ring transform itself (verified: the
            # result is exactly insensitive to it).
            a = coeffs
        else:
            a = self._alm_s
            np.multiply(coeffs, self._ifac, out=a)
            # the m = 0 column must be real: the legacy irfft (and ducc0's own real
            # synthesis) discard Im a_l0, so discard it here too, explicitly.
            a[:, 0].imag = 0.0
        g = self._grid_out if out is None else out
        gv = g.reshape(1, self.nlat, self.nlon)
        _dsht.synthesis_2d(
            alm=a.reshape(1, -1),
            map=gv,
            spin=0,
            lmax=self.lmax,
            mmax=self.lmax,
            geometry="GL",
            ntheta=self.nlat,
            nphi=self.nlon,
            mstart=self._mstart,
            lstride=self._lstride,
            nthreads=self.nthreads,
            **({"ringfactor": self._rf_s} if self.raw else {}),
        )
        return g

    # ------------------------------------------------------- degree-0 only
    def mean00_weights(self):
        """Row weights ``wrow`` with ``analysis(g)[0,0] == wrow @ g.sum(axis=1)``.

        The degree-0 coefficient of an analysis is nothing but a quadrature
        average of the field; the four post-loop GMSL / area analyses in the
        solver use *only* that one number, so they can be replaced by a single
        pass over the grid instead of a full Legendre sweep
        Uses ducc0's own Gauss-Legendre ring
        weights, so it agrees with a full ``analysis`` to ~1e-16 relative.
        """
        gw = _dsht.get_gridweights("GL", self.nlat)  # sums to 4*pi
        # a_00 = (1/sqrt(4 pi)) * sum_j gw_j * rowmean_j ; U[0,0] = LEGACY_OVER_DUCC * a_00
        # (in raw mode V[0,0] == U[0,0] because m = 0)
        return (LEGACY_OVER_DUCC / (np.sqrt(4.0 * np.pi) * self.nlon)) * gw


def make_sht(backend: str, nlat: int, nlon: int, lmax: int, workers: int = 1, chunk: int = 16,
             raw: bool = False):
    """Factory used by ``run_ais.py``, ``run_outer.py`` and ``init_topo_segmented.py``."""
    if backend == "numpy":
        from .sht import SHT

        if raw:
            raise ValueError("raw (sign-free) mode is only implemented for the ducc backend")
        return SHT(nlat, nlon, lmax, workers=workers, chunk=chunk)
    if backend == "ducc":
        return DuccSHT(nlat, nlon, lmax, workers=workers, chunk=chunk, raw=raw)
    raise ValueError(f"unknown SHT backend {backend!r} (expected 'numpy' or 'ducc')")


__all__ = ["DuccSHT", "LEGACY_OVER_DUCC", "make_sht"]
