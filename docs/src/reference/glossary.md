# Glossary

This page defines the terms used throughout the manual.

* *Topography* T (`tgrid`): bedrock elevation relative to the sea surface, m, positive up; the
  ocean is where T < 0. Relative sea-level change is minus the change of T. All output fields are
  defined in [Output fields](../outputs.md).
* *Ocean function* C (`ocean`): 1 where T < 0, else 0. *beta* (`beta`): 1 where the cell is free of
  ice, else 0. The water load sits on C * beta.
* *Marine check* (`checkmarine`): ice on a cell with initial topography T0 <= 0 counts as load only if
  it is grounded, `rhoi * H > rhow * |T0|` (920 and 1000 kg m^-3).
* *Inner iteration*: at each step, the ocean-load increment is updated until the ocean load, the
  deformation, the gravity field and the ocean function are consistent. *epsilon1* is its stopping
  tolerance, a dimensionless ratio (Kendall et al. 2005, eq. 83): with `S_i` the sum of the absolute
  values of the spherical-harmonic coefficients of the ocean-load increment at iteration i, the loop
  stops when `|S_i - S_(i-1)| / S_(i-1) <= epsilon1` (default 1e-5). It is not a length in metres.
* *Outer iteration* (ice age runs): repeats the whole run with an updated initial topography until the
  predicted present topography matches the observed one. The number of passes is set by the user
  (`run_outer.py --iters`, `./exe LAST K0 K1`); the Fortran parameter `epsilon2` is defined but not
  used by the code.
* *calcRG*: also compute the radial displacement R and the sea-surface change G (`R<N>`, `G<N>`).
* *TPW*: rotational feedback (true polar wander).

**Next:** [Method overview](../method.md).
