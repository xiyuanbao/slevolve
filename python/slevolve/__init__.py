"""SLEvolve: exact-memory sea-level solver (Python implementation).

The main classes are exposed lazily, so ``import slevolve`` does not import NumPy; the drivers
(slevolve.run_ais, slevolve.run_outer) pin the thread counts before NumPy is first imported.

    from slevolve import SeaLevelSolver, OptSeaLevelSolver, LoveNumbers, Planet, SHT, make_sht
"""

__version__ = "1.0.0"

_EXPORTS = {
    "SeaLevelSolver": "solver", "LoveNumbers": "solver", "Planet": "solver", "ExactMemory": "solver",
    "read_grid_text": "solver", "OptSeaLevelSolver": "solver_opt", "parse_opts": "solver_opt",
    "SHT": "sht", "make_sht": "sht_ducc",
}

__all__ = ["__version__", *_EXPORTS]


def __getattr__(name):
    if name in _EXPORTS:
        import importlib
        return getattr(importlib.import_module(f".{_EXPORTS[name]}", __name__), name)
    raise AttributeError(f"module 'slevolve' has no attribute {name!r}")
