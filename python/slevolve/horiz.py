"""
horiz.py -- optional horizontal (tangential) displacement read-out.

For a spherically symmetric Earth the surface horizontal displacement is the
surface gradient of one scalar,

    u_h = grad_h Psi ,   u_theta = dPsi/dtheta ,   u_phi = (1/sin theta) dPsi/dphi ,

    Psi_lm(t) = T_l [ l^E_l L_lm(t) + sum_k (l_k/s_k) R_k,lm(t) ]   (+ degree-2 rotation),

theta = colatitude (u_theta > 0 points south), phi = east longitude, Psi in metres
(T_l L_lm is the load potential divided by g, the same factor as in the radial
displacement R_lm = T_l [ h^E_l L_lm + sum_k (h_k/s_k) R_k,lm ]).  The solvers
write the coefficients Psi_lm in the spharmt (U) convention of the other spectral
outputs (``horiz<N>.npy``, '<c16' (L+1, L+1), a[l, m]); this module turns them
into grid components.

``gradient(U, nlat, nlon)`` uses ducc0's spin-1 gradient synthesis
(``mode="GRAD_ONLY"``) on the Gauss-Legendre grid (rows north -> south, first
column at phi = 0), with the coefficient mapping of ``sht_ducc.py``.
"""

from __future__ import annotations

import numpy as np

from .sht_ducc import LEGACY_OVER_DUCC


def legacy_to_ducc(U: np.ndarray) -> np.ndarray:
    """spharmt U[l, m] -> orthonormal Condon-Shortley a_lm (dense [l, m])."""
    nl = U.shape[0]
    sgn = np.where(np.arange(nl) % 2 == 0, 1.0, -1.0)
    a = np.asarray(U, dtype=np.complex128) * (sgn / LEGACY_OVER_DUCC)[None, :]
    a[:, 0].imag = 0.0
    return a


def gradient(U: np.ndarray, nlat: int, nlon: int, nthreads: int = 1):
    """(u_theta, u_phi) on the (nlat, nlon) GL grid from Psi coefficients U[l, m].

    U may hold degrees beyond what the grid resolves only up to
    lmax <= nlat - 1, mmax <= (nlon - 1)//2 (as for DuccSHT).
    """
    import ducc0.sht as dsht

    nl = U.shape[0]
    lmax = nl - 1
    a = legacy_to_ducc(U)
    # spin-1 E-mode coefficients of grad Psi in ducc0's GRAD_ONLY convention:
    # +sqrt(l(l+1)) a_lm (sign checked against analytic fields, see
    # tools/horiz_grad.py --selftest; component 0 = d/dtheta, 1 = d/dphi / sin)
    ell = np.arange(nl, dtype=np.float64)
    a = a * np.sqrt(ell * (ell + 1.0))[:, None]
    mstart = np.arange(nl, dtype=np.uint64)
    out = np.zeros((2, nlat, nlon), dtype=np.float64)
    dsht.synthesis_2d(alm=a.reshape(1, -1), map=out, spin=1, lmax=lmax, mmax=lmax,
                      geometry="GL", ntheta=nlat, nphi=nlon, mstart=mstart,
                      lstride=nl, nthreads=nthreads, mode="GRAD_ONLY")
    return out[0], out[1]


__all__ = ["gradient", "legacy_to_ducc"]
