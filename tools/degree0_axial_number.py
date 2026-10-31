#!/usr/bin/env python3
"""degree0_axial_number.py -- elastic axial number tau_c of an Earth model (axial rotation option `full`).

    python tools/degree0_axial_number.py                                  # lovenum/prem_continental.deck
    python tools/degree0_axial_number.py --convergence                    # step sizes 40, 20, 10, 5 km
    python tools/degree0_axial_number.py my_model.deck --surface 6371e3
    python tools/degree0_axial_number.py --prem-polynomial ti --surface 6368e3   # check model

What it computes.  With axial rotation `full` the axial Liouville equation carries Delta T^R / 3, a third of the
change in the trace A + B + C of the inertia tensor caused by the degree-0 part of the rotational potential.
A change m3 of the spin rate changes that part by Lambda_00(r) = (2/3) Omega^2 r^2 m3 (Y_00 = 1; the solver's
Lambda_00 = K00 (q + m3^2 + 2 m3), K00 = a^2 Omega^2 / 3).  Let U_0(r) be the elastic radial displacement of the
Earth model under the degree-0 potential (r/a)^2 (body force rho grad of it).  The displacement for m3 is then
U = (2/3) Omega^2 a^2 m3 U_0, and since a mass element at radius r adds 2 r^2 to the trace,
    delta(A + B + C) = 2 int rho delta(r^2) dV = 16 pi int_0^a rho U r^3 dr,
    Delta T^R / 3 = tau_c C m3,   tau_c = (32 pi / 9) (Omega^2 a^2 / C) int_0^a rho U_0 r^3 dr,
with a the surface radius of the Earth model (6371 km for the shipped deck) and, by default, the solvers' Omega =
7.292e-5 rad/s and C = 0.3307007 M a^2, M = 5.976e24 kg, a = 6.371e6 m (Planet, planets_mod).
In the solver's notation T = Delta T^R / 3 per unit Lambda_00(a) = (16 pi / 3) int rho U_0 r^3 dr = (3/2) tau_c C /
(a^2 Omega^2).  The tool also prints m_c = -sqrt(4 pi) int_0^a rho U_0 r^3 dr [kg s^2], the same response in
the response per unit degree-0 potential (r/a)^2 Y00; tau_c = (16 sqrt(pi) / 9)
(-m_c) Omega^2 a^2 / C.

Equations (static, elastic, compressible, self-gravitating, hydrostatic prestress; Takeuchi & Saito 1972 and
Dahlen & Tromp 1998, radial modes at zero frequency).  The displacement is radial, u = U(r) r_hat; sigma is the
radial incremental (Lagrangian) stress; for a transversely isotropic medium with radial symmetry axis (moduli
A, C, F, N; isotropic: A = C = lambda + 2 mu, F = lambda, N = mu)
    sigma        = C U' + 2 F U / r
    sigma_tt     = F U' + 2 (A - N) U / r                       (tangential stress)
    sigma'       = -(2 / r) (sigma - sigma_tt) - 4 rho g U / r - rho dPhi/dr,   Phi = (r/a)^2.
The degree-0 gravity perturbation is fixed by mass conservation (the mass inside r changes by -4 pi r^2 rho U),
and its 4 pi G rho^2 U term cancels the one from g' = 4 pi G rho - 2 g / r: only g enters.  A fluid layer
(N = 0, A = C = F = kappa) needs no special treatment at degree 0 (no tangential motion).  Conditions: U regular
at r = 0 (U ~ r), U and sigma continuous across every interface (solid-fluid included), sigma(a) = 0.
Solution: the regular homogeneous solution and a particular solution are integrated outward together with
int rho U r^3 dr (classical RK4, steps <= h, geometrically refined towards r = 0, step / r <= h / 1000 km),
combined so that sigma(a) = 0.  Fourth order in h.

Earth models.
  deck (default lovenum/prem_continental.deck): columns r [m], rho, vp, vs [, g]; centre outward; rho, lambda,
       mu and g linear between nodes (as in the code that computed the shipped Love numbers); two nodes at one
       radius mark a discontinuity.  Isotropic.  g from the deck, or from the density with --g density.
  --prem-polynomial ti|iso: PREM of Dziewonski & Anderson (1981, Table I) as exact polynomials (x = r / 6371 km),
       the ocean replaced by upper crust up to --surface (6371e3: the variant of the shipped Love numbers;
       6368e3: PREM without its ocean).  `ti`: transversely isotropic 6151-6346.6 km; `iso`: the isotropic
       version.  g from the density.  Used for the check of the solver.
Needs numpy only.
"""
import argparse
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_DECK = os.path.join(HERE, "..", "lovenum", "prem_continental.deck")
G_NEWTON = 6.6743e-11
# SLEvolve's planet (python/slevolve/solver.py Planet; Fortran planets_mod): a, Omega, M, C / (M a^2)
A_SLE, OMEGA_SLE, MASS_SLE, CFAC_SLE = 6.371e6, 7.292e-5, 5.976e24, 0.3307007

_GL_X, _GL_W = np.polynomial.legendre.leggauss(6)


# ----------------------------------------------------------------------------- Earth models
class Layer:
    """One radial interval [r0, r1] with smooth properties; props(r) -> (rho, A, C, F, N, g)."""

    def __init__(self, r0, r1, props):
        self.r0, self.r1, self.props = r0, r1, props


def read_deck(path):
    """Deck columns r [m], rho [kg/m^3], vp, vs [m/s] (, g [m/s^2]); '#' lines are comments."""
    d = np.loadtxt(path, comments="#")
    if d.ndim != 2 or d.shape[1] < 4:
        raise ValueError(f"{path}: need at least 4 columns r rho vp vs")
    if np.any(np.diff(d[:, 0]) < 0) or d[0, 0] != 0.0:
        raise ValueError(f"{path}: radii must start at 0 and be non-decreasing")
    return d


def _mass_profile(layers, G):
    """Mass below the bottom of each layer (exact for the polynomial densities used here)."""
    m, out = 0.0, []
    for L in layers:
        out.append(m)
        rr = 0.5 * (L.r1 - L.r0) * _GL_X + 0.5 * (L.r1 + L.r0)
        m += 0.5 * (L.r1 - L.r0) * sum(w * 4 * math.pi * L.props(r)[0] * r * r for r, w in zip(rr, _GL_W))
    return out, m


def _with_density_gravity(layers, G):
    """Replace g in every layer by G M(r) / r^2 from the density."""
    mbelow, _ = _mass_profile(layers, G)
    out = []
    for L, mb in zip(layers, mbelow):
        def props(r, L=L, mb=mb):
            rho, A, C, F, N, _g = L.props(r)
            if r <= 0.0:
                return rho, A, C, F, N, 0.0
            x = 0.5 * (r - L.r0) * _GL_X + 0.5 * (r + L.r0)
            m = mb + 0.5 * (r - L.r0) * sum(w * 4 * math.pi * L.props(s)[0] * s * s for s, w in zip(x, _GL_W))
            return rho, A, C, F, N, G * m / (r * r)
        out.append(Layer(L.r0, L.r1, props))
    return out


def deck_model(path, gravity="deck", G=G_NEWTON):
    """Layers of a deck: rho, lambda, mu, g linear between nodes."""
    d = read_deck(path)
    r, rho, vp, vs = d[:, 0], d[:, 1], d[:, 2], d[:, 3]
    mu = rho * vs**2
    lam = rho * vp**2 - 2 * mu
    gd = d[:, 4] if d.shape[1] >= 5 else None
    if gravity == "deck" and gd is None:
        raise ValueError(f"{path}: no g column; use --g density")
    layers = []
    for i in range(len(r) - 1):
        if r[i + 1] <= r[i]:
            continue  # discontinuity
        def props(x, i=i):
            f = (x - r[i]) / (r[i + 1] - r[i])
            ro = rho[i] + f * (rho[i + 1] - rho[i])
            m_ = mu[i] + f * (mu[i + 1] - mu[i])
            la = lam[i] + f * (lam[i + 1] - lam[i])
            g = 0.0 if gd is None else gd[i] + f * (gd[i + 1] - gd[i])
            return ro, la + 2 * m_, la + 2 * m_, la, m_, g
        layers.append(Layer(r[i], r[i + 1], props))
    if gravity == "density":
        layers = _with_density_gravity(layers, G)
    return layers, r[-1]


# PREM, Dziewonski & Anderson (1981) Table I: (r0, r1) km, coefficients in x = r / 6371 km (g/cm^3, km/s)
_PREM = [
    (0.0, 1221.5, dict(rho=(13.0885, 0, -8.8381), vp=(11.2622, 0, -6.3640), vs=(3.6678, 0, -4.4475))),
    (1221.5, 3480.0, dict(rho=(12.5815, -1.2638, -3.6426, -5.5281), vp=(11.0487, -4.0362, 4.8023, -13.5732),
                          vs=(0.0,))),
    (3480.0, 3630.0, dict(rho=(7.9565, -6.4761, 5.5283, -3.0807), vp=(15.3891, -5.3181, 5.5242, -2.5514),
                          vs=(6.9254, 1.4672, -2.0834, 0.9783))),
    (3630.0, 5600.0, dict(rho=(7.9565, -6.4761, 5.5283, -3.0807), vp=(24.9520, -40.4673, 51.4832, -26.6419),
                          vs=(11.1671, -13.7818, 17.4575, -9.2777))),
    (5600.0, 5701.0, dict(rho=(7.9565, -6.4761, 5.5283, -3.0807), vp=(29.2766, -23.6027, 5.5242, -2.5514),
                          vs=(22.3459, -17.2473, -2.0834, 0.9783))),
    (5701.0, 5771.0, dict(rho=(5.3197, -1.4836), vp=(19.0957, -9.8672), vs=(9.9839, -4.9324))),
    (5771.0, 5971.0, dict(rho=(11.2494, -8.0298), vp=(39.7027, -32.6166), vs=(22.3512, -18.5856))),
    (5971.0, 6151.0, dict(rho=(7.1089, -3.8045), vp=(20.3926, -12.2569), vs=(8.9496, -4.4597))),
    (6151.0, 6346.6, dict(rho=(2.6910, 0.6924), vp=(4.1875, 3.9382), vs=(2.1519, 2.3481),
                          vpv=(0.8317, 7.2180), vph=(3.5908, 4.6172), vsv=(5.8582, -1.4678),
                          vsh=(-1.0839, 5.7176), eta=(3.3687, -2.4778))),
    (6346.6, 6356.0, dict(rho=(2.900,), vp=(6.800,), vs=(3.900,))),
    (6356.0, 6371.0, dict(rho=(2.600,), vp=(5.800,), vs=(3.200,))),  # upper crust, up to the chosen surface
]


def prem_polynomial(surface=6371e3, variant="ti", G=G_NEWTON):
    """PREM layers (ocean replaced by upper crust up to `surface`), exact polynomials; g from the density."""
    if not 6356e3 < surface <= 6371e3:
        raise ValueError("surface must lie in the upper crust (6356-6371 km)")
    layers = []
    for r0, r1, c in _PREM:
        r0, r1 = r0 * 1e3, min(r1 * 1e3, surface)
        ti = variant == "ti" and "vpv" in c

        def props(r, c=c, ti=ti):
            x = r / 6371e3
            p = lambda k: sum(a * x**n for n, a in enumerate(c[k]))  # noqa: E731
            rho = p("rho") * 1e3
            if ti:
                A, C = rho * (p("vph") * 1e3)**2, rho * (p("vpv") * 1e3)**2
                L_, N = rho * (p("vsv") * 1e3)**2, rho * (p("vsh") * 1e3)**2
                return rho, A, C, p("eta") * (A - 2 * L_), N, 0.0
            vp, vs = p("vp") * 1e3, p("vs") * 1e3
            M = rho * vp * vp
            N = rho * vs * vs
            return rho, M, M, M - 2 * N, N, 0.0
        layers.append(Layer(r0, r1, props))
    return _with_density_gravity(layers, G), surface


# ----------------------------------------------------------------------------- solver
def _rhs(L, r, y, a):
    """d/dr of (U_h, s_h, U_p, s_p, I_h, I_p); h = homogeneous, p = forced by Phi = (r/a)^2."""
    rho, A, C, F, N, g = L.props(r)
    out = []
    for U, s, forced in ((y[0], y[1], False), (y[2], y[3], True)):
        dU = (s - 2 * F * U / r) / C
        stt = F * dU + 2 * (A - N) * U / r
        ds = -2 * (s - stt) / r - 4 * rho * g * U / r - (rho * 2 * r / (a * a) if forced else 0.0)
        out += [dU, ds]
    r3 = rho * r**3
    return out + [r3 * y[0], r3 * y[2]]


def _mesh(L, h, first):
    """Step end points in [L.r0, L.r1]: steps <= h; in the first layer also step <= r h / 1000 km."""
    if not first:
        n = max(1, int(math.ceil((L.r1 - L.r0) / h - 1e-9)))
        return list(np.linspace(L.r0, L.r1, n + 1)[1:])
    q = h / 1.0e6
    pts, r = [], 1e-6 * h  # start of the outward integration
    while r < L.r1:
        r = min(r + min(h, q * r), L.r1)
        if L.r1 - r < 1e-9 * L.r1:
            r = L.r1
        pts.append(r)
    return pts


def solve(layers, a, h=10e3):
    """int_0^a rho U_0 r^3 dr for the potential (r/a)^2, and U_0(a).  h: maximum radial step [m]."""
    L0 = layers[0]
    r = 1e-6 * h
    rho, A, C, F, N, g = L0.props(0.0)
    y = [r, (C + 2 * F) * r, 0.0, 0.0, 0.0, 0.0]  # regular solution U ~ r at the centre
    for k, L in enumerate(layers):
        for r1 in _mesh(L, h, k == 0):
            dr = r1 - r
            k1 = _rhs(L, r, y, a)
            k2 = _rhs(L, r + dr / 2, [yi + dr / 2 * ki for yi, ki in zip(y, k1)], a)
            k3 = _rhs(L, r + dr / 2, [yi + dr / 2 * ki for yi, ki in zip(y, k2)], a)
            k4 = _rhs(L, r1, [yi + dr * ki for yi, ki in zip(y, k3)], a)
            y = [yi + dr / 6 * (p + 2 * q + 2 * s + t) for yi, p, q, s, t in zip(y, k1, k2, k3, k4)]
            r = r1
    c = -y[3] / y[1]  # sigma(a) = 0
    return {"I0": y[5] + c * y[4], "U0_surface": y[2] + c * y[0]}


def axial_numbers(I0, a=A_SLE, omega=OMEGA_SLE, C=CFAC_SLE * MASS_SLE * A_SLE**2):
    """m_c [kg s^2], T (Delta T^R / 3 per unit Lambda_00(a)) [kg s^2] and tau_c; a = surface radius of the model."""
    return {"m_c": -2 * math.sqrt(math.pi) * I0, "T": 16 * math.pi / 3 * I0,
            "tau_c": 32 * math.pi / 9 * omega**2 * a**2 / C * I0}


def tau_c_of_deck(path=DEFAULT_DECK, h=10e3, gravity="deck"):
    """tau_c of a deck with SLEvolve's a, Omega, C (the value used by the solvers)."""
    layers, a = deck_model(path, gravity)
    return axial_numbers(solve(layers, a, h)["I0"], a)["tau_c"]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("deck", nargs="?", default=DEFAULT_DECK, help="Earth-model deck (r rho vp vs [g])")
    ap.add_argument("--prem-polynomial", choices=("ti", "iso"), help="use the PREM polynomials instead of a deck")
    ap.add_argument("--surface", type=float, default=6371e3, help="surface radius [m] for --prem-polynomial")
    ap.add_argument("--g", choices=("deck", "density"), default="deck", help="gravity of a deck (default: deck)")
    ap.add_argument("--h", type=float, default=10.0, help="maximum radial step [km] (default 10)")
    ap.add_argument("--convergence", action="store_true", help="also h = 4h, 2h, h/2 and the observed order")
    ap.add_argument("--omega", type=float, default=OMEGA_SLE)
    ap.add_argument("--mass", type=float, default=MASS_SLE)
    ap.add_argument("--moi-factor", type=float, default=CFAC_SLE, help="C / (M a^2)")
    ap.add_argument("--radius", type=float, default=A_SLE, help="a in C = f M a^2 (default: SLEvolve's 6.371e6)")
    a_ = ap.parse_args(argv)
    if a_.prem_polynomial:
        layers, rs = prem_polynomial(a_.surface, a_.prem_polynomial)
        name = f"PREM polynomials ({a_.prem_polynomial}), surface {rs / 1e3:g} km"
    else:
        layers, rs = deck_model(a_.deck, a_.g)
        name = f"{os.path.basename(a_.deck)} (g: {a_.g}), surface {rs / 1e3:g} km"
    C = a_.moi_factor * a_.mass * a_.radius**2
    hs = [4 * a_.h, 2 * a_.h, a_.h, a_.h / 2] if a_.convergence else [a_.h]
    print(f"model: {name}; Omega = {a_.omega:g} rad/s, C = {a_.moi_factor} M a^2 with M = {a_.mass:g} kg, "
          f"a = {a_.radius:g} m")
    res = []
    for h in hs:
        s = solve(layers, rs, h * 1e3)
        n = axial_numbers(s["I0"], rs, a_.omega, C)
        res.append(n)
        print(f"h = {h:7.3f} km   m_c = {n['m_c']:.15e} kg s^2   tau_c = {n['tau_c']:.15e}   "
              f"U_0(a) = {s['U0_surface']:.6e} s^2/m")
    if len(res) == 4:
        e = [abs(res[i]["m_c"] - res[i + 1]["m_c"]) for i in range(3)]
        orders = [math.log2(e[i] / e[i + 1]) if e[i + 1] > 0 else float("inf") for i in range(2)]
        print("successive differences of m_c (relative): " + ", ".join(f"{x / abs(res[-1]['m_c']):.2e}" for x in e)
              + "; observed order " + ", ".join(f"{o:.2f}" for o in orders))
    return res


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
