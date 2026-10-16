"""
npyio.py -- shared file conventions of the Python and Fortran
(fortran/src/npy_io_mod.f90) solvers.  See fortran/IO_FORMATS.md.

Grids (ice, topography, tgrid/ocean/beta/R/G snapshots)
    .npy, shape (nlat, nlon) = (512, 1024), fortran_order=False, dtype '<f8'
    (ice inputs may be '<f4' when lossless).  a[j, i]: j = Gauss-Legendre row
    (north -> south), i = longitude; the payload is the text-file stream order,
    identical to the legacy headerless .bin caches.
Spectral fields (dS_converged)
    '<c16', shape (N+1, N+1), fortran_order=True, a[l, m] = dS(l, m).
Masks: '|i1'.

Every reader returns float64 C-contiguous arrays.
"""
import os
import numpy as np

NLAT, NLON = 512, 1024

# name -> (printf format, delimiter, values per line; 0 = one per line)
TEXT_FORMATS = {
    "es16.9": ("%16.9E", "", 1),       # ES16.9E2 one-value-per-line output format
    "row-e16.7": ("%16.7e", "", 1024),   # OSF archive topo / tgrid files
    "row-e.18": ("%.18e", " ", 1024),  # np.savetxt default (OSF ice files)
    "row-E18.10": ("%18.10E", " ", 1024),  # nadi-case ice files
}


def read_text_grid(path, nlat=NLAT, nlon=NLON):
    """Legacy list-directed text grid (any values per line) -> (nlat, nlon) float64, stream order."""
    a = np.fromfile(path, sep=" ")
    if a.size != nlat * nlon:
        raise ValueError(f"{path}: got {a.size} values, expected {nlat * nlon}")
    return a.reshape(nlat, nlon)


def write_text_grid(path, a, fmt="es16.9", ncol=None):
    """Write a grid in one of the legacy text layouts (C-flat order)."""
    f, delim, per = TEXT_FORMATS[fmt]
    if ncol is not None and per != 1:
        per = ncol
    flat = np.ascontiguousarray(a, dtype=np.float64).ravel(order="C")
    with open(path, "w") as fh:
        if per == 1:
            fh.write("\n".join(f % v for v in flat.tolist()))
            fh.write("\n")
        else:
            vals = [f % v for v in flat.tolist()]
            for s in range(0, len(vals), per):
                fh.write(delim.join(vals[s : s + per]))
                fh.write("\n")


def f4_exact(a):
    a = np.asarray(a, dtype=np.float64)
    with np.errstate(over="ignore"):
        return bool(np.array_equal(a.astype(np.float32).astype(np.float64), a))


def save_grid(path, a, dtype="f8"):
    """Save a grid as .npy (C order).  dtype 'f4' only if lossless, else f8.
    Returns the dtype actually used."""
    a = np.ascontiguousarray(a, dtype=np.float64)
    if dtype == "f4" and f4_exact(a):
        np.save(path, a.astype("<f4"))
        return "f4"
    np.save(path, a.astype("<f8", copy=False))
    return "f8"


def save_spectral(path, z):
    """complex (N+1, N+1) spectral field, a[l, m]; written fortran_order=True
    so the payload equals the Fortran dS(0:N,0:N) storage order."""
    np.save(path, np.asfortranarray(np.asarray(z, dtype="<c16")))


def load_npy(path, mmap=False):
    return np.load(path, mmap_mode="r" if mmap else None)


def read_grid(path, nlat=NLAT, nlon=NLON):
    """Read a grid from .npy (any real dtype/order), legacy .bin (headerless LE
    f8 stream) or legacy text.  Without an extension, prefers <path>.npy, then
    <path>.bin, then the text file itself (the Fortran pp_read_grid order)."""
    if not path.endswith((".npy", ".bin")):
        for ext in (".npy", ".bin"):
            if os.path.exists(path + ext):
                return read_grid(path + ext, nlat, nlon)
        return read_text_grid(path, nlat, nlon)
    if path.endswith(".npy"):
        a = np.load(path)
        if a.shape != (nlat, nlon):
            if a.size != nlat * nlon:
                raise ValueError(f"{path}: shape {a.shape}, expected {(nlat, nlon)}")
            a = a.reshape(nlat, nlon)  # C-flat (stream) order, as the Fortran reader
        return np.ascontiguousarray(a, dtype=np.float64)
    if path.endswith(".bin"):
        n = os.path.getsize(path)
        if n != 8 * nlat * nlon:
            raise ValueError(f"{path}: {n} bytes, expected {8 * nlat * nlon} (LE f8)")
        return np.fromfile(path, dtype="<f8").reshape(nlat, nlon)
    return read_text_grid(path, nlat, nlon)


def find_grid(dirname, name, nlat=NLAT, nlon=NLON):
    """<dir>/<name>.npy | <dir>/<name>.bin | <dir>/<name> (text)."""
    return read_grid(os.path.join(dirname, name), nlat, nlon)


def grid_source(dirname, name):
    for ext in (".npy", ".bin", ""):
        p = os.path.join(dirname, name + ext)
        if os.path.exists(p):
            return p
    return None


class IceDir:
    """ice[k] read on demand from <dir>/<pfx><k>{.npy,.bin,} (npy preferred)."""

    def __init__(self, d, pfx, nlat=NLAT, nlon=NLON):
        self.d, self.pfx, self.nlat, self.nlon = d, pfx, nlat, nlon

    def __getitem__(self, k):
        return find_grid(self.d, f"{self.pfx}{k}", self.nlat, self.nlon)


# ------------------------------------------------------------ ice-cache manifest
def ice_manifest(dirname, pfx, nfiles):
    """Key of a binary ice cache: ice directory, prefix, file count and the
    (name, size, mtime_ns) of every source file (the one npyio would read)."""
    files = []
    for k in range(nfiles):
        src = grid_source(dirname, f"{pfx}{k}")
        if src is None:
            raise FileNotFoundError(os.path.join(dirname, f"{pfx}{k}"))
        st = os.stat(src)
        files.append([os.path.basename(src), int(st.st_size), int(st.st_mtime_ns)])
    return {"icedir": os.path.abspath(dirname), "prefix": pfx, "nfiles": int(nfiles), "files": files}


def manifest_path(cache_path):
    return cache_path + ".json"


def cache_matches(cache_path, man):
    """(True, '') if the cache at cache_path was built from the files in man
    (a cache built from more files is accepted when its first entries match)."""
    import json

    if not os.path.exists(cache_path):
        return False, "no cache"
    mp = manifest_path(cache_path)
    if not os.path.exists(mp):
        return False, "no manifest"
    try:
        with open(mp) as fh:
            old = json.load(fh)
    except (OSError, ValueError):
        return False, "unreadable manifest"
    if old.get("icedir") != man["icedir"]:
        return False, "ice directory changed"
    if old.get("prefix") != man["prefix"]:
        return False, "prefix changed"
    if int(old.get("nfiles", -1)) < man["nfiles"]:
        return False, "file count changed"
    if [list(x) for x in old.get("files", [])[: man["nfiles"]]] != man["files"]:
        return False, "ice file size or mtime changed"
    return True, ""


def write_manifest(cache_path, man):
    import json

    with open(manifest_path(cache_path), "w") as fh:
        json.dump(man, fh, indent=0)


def read_times(path, nsteps=None):
    """Read an optional times file: one time t_k (years) per line for ice files k = 0..N, strictly
    increasing (t_0 any value).  Returns a float64 array; with nsteps given, at least nsteps + 1
    values are required and the first nsteps + 1 are returned.  Step j has length t_j - t_{j-1}."""
    t = np.loadtxt(path, dtype=np.float64, ndmin=1)
    if t.ndim != 1 or t.size < 2:
        raise ValueError(f"{path}: need one time per line, at least 2 values")
    if nsteps is not None:
        if t.size < nsteps + 1:
            raise ValueError(f"{path}: {t.size} times, need {nsteps + 1} (files 0..{nsteps})")
        t = t[: nsteps + 1]
    if not np.all(np.isfinite(t)) or not np.all(np.diff(t) > 0):
        raise ValueError(f"{path}: times must be finite and strictly increasing")
    return t
