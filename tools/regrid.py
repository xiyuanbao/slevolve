"""Regridding from the PSUICE3D 10-km polar-stereographic grid to the
512x1024 Gauss-Legendre grid.

The mapping of Han et al. (2022) was recovered from the data, not guessed: every non-zero
value in the archived `AISload_Run85_0_0_gd_0` is an exact float32 value of the
native PSUICE3D grounded thickness field at t = 0, so the mapping is a
nearest-source SAMPLING, and matching value -> source cell gives the actual
index map.  Fitting that map shows

    row  j  <->  Gauss-Legendre latitude  lat[j]        (north -> south)
    col  i  <->  longitude  (i + 1) * 360/1024          (note the one-cell
                                                         longitude offset)
    value   =   h(source cell nearest to that point on the sphere)

which reproduces 99.0 % of the 48 701 archived grounded cells cell-for-cell and the archived
ice volume to 4e-5.
"""
import os
import sys

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "python"))
from slevolve.sht import gaulw  # noqa: E402

NLAT, NLON = 512, 1024
RE = 6.371e6          # Earth radius of the default (Han et al. 2022) constants
RHOI = 920.0
RHOW = 1000.0
LON_SHIFT = 1         # recovered from the archived grounded files


def gl_grid(nlat=NLAT, nlon=NLON):
    mu, w = gaulw(nlat)                      # mu = sin(lat), decreasing (N -> S)
    lat = np.degrees(np.arcsin(mu))
    lon = np.arange(nlon) * (360.0 / nlon)
    edges = np.empty(nlat + 1)
    edges[0] = 1.0
    edges[1:] = 1.0 - np.cumsum(w)
    edges[-1] = -1.0
    area = (RE ** 2) * w[:, None] * np.full((1, nlon), 2 * np.pi / nlon)
    return lat, lon, mu, w, edges, area


def _unit_xyz(latd, lond):
    la = np.radians(latd)
    lo = np.radians(lond)
    c = np.cos(la)
    return np.stack([c * np.cos(lo), c * np.sin(lo), np.sin(la)], axis=-1)


class NearestSourceSampler:
    """Nearest-source sampling from any source grid (2-D arrays of latitude and longitude in degrees,
    alat, alon) onto the GL grid: every GL cell within maxdist_km of a source point takes the value of
    the nearest source point, all other cells are 0.  lon_shift = 1 (default) is the one-cell offset of
    the archived Case 1 files; use lon_shift=0 to sample onto the solver grid lon = i*360/nlon."""

    def __init__(self, alat, alon, nlat=NLAT, nlon=NLON, maxdist_km=12.0,
                 lon_shift=LON_SHIFT):
        lat, lon, mu, w, edges, area = gl_grid(nlat, nlon)
        self.nlat, self.nlon, self.area, self.lat = nlat, nlon, area, lat
        al = alat.ravel().astype(np.float64)
        ao = alon.ravel().astype(np.float64)
        latmin, latmax = al.min(), al.max()
        glat = np.repeat(lat, nlon)
        glon = np.tile((np.arange(nlon) + lon_shift) * (360.0 / nlon), nlat)
        band = (glat >= latmin - 1.0) & (glat <= latmax + 1.0)
        self.band = np.flatnonzero(band)
        tree = cKDTree(_unit_xyz(al, ao))
        dist, near = tree.query(_unit_xyz(glat[band], glon[band]), k=1)
        chord = 2.0 * np.sin(0.5 * (maxdist_km * 1e3 / RE))
        self.near = near
        self.ok = dist <= chord

    def sample(self, field):
        f = field.ravel().astype(np.float64)
        out = np.zeros(self.nlat * self.nlon)
        v = f[self.near]
        v[~self.ok] = 0.0
        out[self.band] = v
        return out.reshape(self.nlat, self.nlon)


def grounded_mask_from_flotation(h, hb, rhoi=RHOI, rhow=RHOW):
    """PSUICE3D flotation test: floating where h*rhoi/rhow < -hb."""
    return ~(h * (rhoi / rhow) < -hb)
