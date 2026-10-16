"""Default locations used by the drivers.

In a repository checkout the package lives in <repo>/python/slevolve; the ice caches then go to
<repo>/python/cache and the development root is <repo> (both as in the script layout).  For an
installed package the caches go to $SLEVOLVE_CACHE_DIR, else ~/.cache/slevolve, and the root is the
current directory.  $SLEVOLVE_CACHE_DIR overrides the cache location in both cases.
"""

import os

PKG = os.path.dirname(os.path.abspath(__file__))
_PYDIR = os.path.dirname(PKG)
IN_CHECKOUT = os.path.isfile(os.path.join(_PYDIR, "run_ais.py")) and os.path.basename(_PYDIR) == "python"
ROOT = os.path.dirname(_PYDIR) if IN_CHECKOUT else os.getcwd()
CACHE = os.environ.get("SLEVOLVE_CACHE_DIR") or (
    os.path.join(_PYDIR, "cache") if IN_CHECKOUT else os.path.join(os.path.expanduser("~"), ".cache", "slevolve"))
