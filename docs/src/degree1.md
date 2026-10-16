# Degree-1 Love numbers and reference frames

This page explains the degree-1 frame of the shipped Love-number files and the degree-1 collocation initialisation defect (a copy of `lovenum/DEGREE1_NOTES.md`).

## 1. Degree-1 frames
At degree 1 a rigid translation of the Earth is a free mode, so h_1, l_1 and k_1 are defined only once an origin is chosen. Moving the origin adds the same constant to h_1, l_1 and k_1. Common conventions:
- **CM** (centre of mass of Earth + load): k_1 = -1;
- **CE** (centre of mass of the solid Earth): k_1 = 0;
- **CF** (centre of surface figure): h_1 + 2 l_1 = 0.

All Love-number files in `lovenum/` use the CE frame (k_1 = 0) for their elastic and fluid limits (see Sect. 2 for the degree-1 samples of the Case 1 file), and SLEvolve applies no frame conversion: degree-1 displacement and geoid outputs are in the frame of the Love-number file. Differences such as h_1 - k_1 and l_1 - k_1 do not depend on the frame, and sea level, which depends on 1 + k - h, is frame-independent.

## 2. The degree-1 collocation initialisation defect
The collocation method (Peltier 1985; Mitrovica and Peltier 1992) fits sum_j a_j/(s_i + s_j) = X(s_i) - X_el at fixed poles s_i, where X_el is the elastic limit (s -> infinity), computed separately. By construction the left side is the relaxation part only: it vanishes as s -> infinity. The fit therefore assumes that the samples X(s_i) and the limit X_el are in the same degree-1 frame.

In earlier library files (the Love-number library of Austermann and Mitrovica 2015, from which the Case 2 file was previously taken), the degree-1 samples were initialised in CM while the limits were in CE. In plain terms: every sample then carries a constant offset of -1 in h, l and k, and a sum of decaying terms a_j/(s + s_j) cannot represent a constant, so the fit compensates with large alternating residues; the degree-1 relaxation curves h_1(t), l_1(t) and k_1(t) are then wrong by the same time-dependent amount. The Case 2 and Case 3 files in `lovenum/` were computed with the samples and limits in the same frame, which removes the defect (their degree-1 k residues are exactly 0). The Case 1 file (`prem_512.l60K2C.sum18p6.dum19p2.tz19p4.lm22`, distributed with Han et al. 2022 and included unchanged) carries the same degree-1 defect, with no effect on sea level and topography.

Because the same error enters h, l and k, it cancels in h - k: sea level and topography are unaffected (degree-1 h - k, a dimensionless Love number, unchanged to 3e-5). The separately output degree-1 fields (radial displacement, geoid and horizontal displacement) differ.

## References
- Austermann, J. and Mitrovica, J. X.: Calculating gravitationally self-consistent sea level changes driven by dynamic topography, Geophys. J. Int., 203, 1909-1922, 2015.
- Han, H. K., Gomez, N., and Wan, J. X. W.: Geosci. Model Dev., 15, 1355-1373, 2022; data: doi:10.17605/OSF.IO/8PTFM.
- Mitrovica, J. X. and Peltier, W. R.: A comparison of methods for the inversion of viscoelastic relaxation spectra, Geophys. J. Int., 108, 410-414, 1992.
- Peltier, W. R.: The LAGEOS constraint on deep mantle viscosity: Results from a new normal mode method for the inversion of viscoelastic relaxation spectra, J. Geophys. Res., 90, 9411-9421, 1985.

**Next:** [Release notes and known limitations](release_notes.md).
