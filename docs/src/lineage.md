# Lineage and references

This page describes the code lineage of SLEvolve and lists the references of that lineage.

The pseudo-spectral solution of the sea-level equation was developed by Mitrovica and Peltier (1991)
and extended to a rotating Earth by Milne and Mitrovica (1998); its numerical formulation with
migrating shorelines and the outer iteration on the initial topography is that of Kendall et al.
(2005), which SLEvolve retains in full. Only the evaluation of the viscous history term is changed.

The ice age sea-level code implementing this formulation was written at Harvard by S. Goldberg
(2015-2016, senior thesis supervised by J. X. Mitrovica; see also Goldberg et al., 2016), as recorded
in the headers of his source and of the released legacy code. It was modified by E. N.-H. Chan
(McGill University, 2017-2019) for text input and output and a Mars option, and then extended by
Han et al. (2022) into 1DSeaLevelModel_FWTW with the forward algorithm of Gomez et al. (2010), an
ice-sheet coupling interface following Gomez et al. (2012) and the time window algorithm (MIT
license). The Fortran version of SLEvolve is derived from that code; the Python version was written
from scratch following the same formulation. SLcode (Austermann and Mitrovica, 2015) is an
independent implementation of the same formulation.

## References

- Mitrovica, J. X. and Peltier, W. R.: On postglacial geoid subsidence over the equatorial oceans,
  J. Geophys. Res. Solid Earth, 96, 20053-20071, https://doi.org/10.1029/91jb01284, 1991.
- Milne, G. A. and Mitrovica, J. X.: Postglacial sea-level change on a rotating Earth, Geophys. J.
  Int., 133, 1-19, https://doi.org/10.1046/j.1365-246x.1998.1331455.x, 1998.
- Kendall, R. A., Mitrovica, J. X., and Milne, G. A.: On post-glacial sea level - II. Numerical
  formulation and comparative results on spherically symmetric models, Geophys. J. Int., 161,
  679-706, https://doi.org/10.1111/j.1365-246x.2005.02553.x, 2005.
- Gomez, N., Mitrovica, J. X., Tamisiea, M. E., and Clark, P. U.: A new projection of sea level change
  in response to collapse of marine sectors of the Antarctic Ice Sheet, Geophys. J. Int., 180, 623-634,
  https://doi.org/10.1111/j.1365-246X.2009.04419.x, 2010.
- Gomez, N., Pollard, D., Mitrovica, J. X., Huybers, P., and Clark, P. U.: Evolution of a coupled marine
  ice sheet-sea level model, J. Geophys. Res. Earth Surf., 117, F01013,
  https://doi.org/10.1029/2011JF002128, 2012.
- Austermann, J. and Mitrovica, J. X.: Calculating gravitationally self-consistent sea level changes
  driven by dynamic topography, Geophys. J. Int., 203, 1909-1922,
  https://doi.org/10.1093/gji/ggv371, 2015.
- Goldberg, S. L., Lau, H. C. P., Mitrovica, J. X., and Latychev, K.: The timing of the Black Sea flood
  event: Insights from modeling of glacial isostatic adjustment, Earth Planet. Sci. Lett., 452,
  178-184, https://doi.org/10.1016/j.epsl.2016.06.016, 2016.
- Han, H. K., Gomez, N., and Wan, J. X. W.: Capturing the interactions between ice sheets, sea level
  and the solid Earth on a range of timescales: a new "time window" algorithm, Geosci. Model Dev.,
  15, 1355-1373, https://doi.org/10.5194/gmd-15-1355-2022, 2022.

References for the optional self-consistent axial rotation ([Method](method.md#self-consistent-axial-rotation-optional)):

- Mitrovica, J. X., Wahr, J., Matsuyama, I., and Paulson, A.: The rotational stability of an ice-age
  earth, Geophys. J. Int., 161, 491-506, https://doi.org/10.1111/j.1365-246X.2005.02609.x, 2005.
- Mitrovica, J. X. and Wahr, J.: Ice age Earth rotation, Annu. Rev. Earth Planet. Sci., 39, 577-616,
  https://doi.org/10.1146/annurev-earth-040610-133404, 2011.

Data and model references cited on [Data provenance](provenance.md) and [Running the benchmark cases](benchmarks.md):

- Argus, D. F., Peltier, W. R., Drummond, R., and Moore, A. W.: The Antarctica component of postglacial
  rebound model ICE-6G_C (VM5a) based on GPS positioning, exposure age dating of ice thicknesses, and
  relative sea level histories, Geophys. J. Int., 198, 537-563, https://doi.org/10.1093/gji/ggu140, 2014.
- DeConto, R. M., Pollard, D., Alley, R. B., Velicogna, I., Gasson, E., Gomez, N., Sadai, S., Condron, A.,
  Gilford, D. M., Ashe, E. L., Kopp, R. E., Li, D., and Dutton, A.: The Paris Climate Agreement and future
  sea-level rise from Antarctica, Nature, 593, 83-89, https://doi.org/10.1038/s41586-021-03427-0, 2021.
- Han, H. K., Gomez, N., Pollard, D., and DeConto, R.: Modeling Northern Hemispheric ice sheet dynamics,
  sea level change, and solid Earth deformation through the last glacial cycle, J. Geophys. Res. Earth
  Surf., 126, e2020JF006040, https://doi.org/10.1029/2020JF006040, 2021.
- Lucas, E., Kodama, S., Pico, T., Stoker, B. J., Lesnek, A., and Dalton, A. S.: Data from: ICE-NADI-1:
  A new deglacial ice sheet thickness history (25-1 ka) for North America using empirically-constrained
  ice margins with implications for global sea-level contributions, relative sea level, and vertical
  land motion, Zenodo [data set], https://doi.org/10.5281/zenodo.21464139, 2026.
- Peltier, W. R., Argus, D. F., and Drummond, R.: Space geodesy constrains ice age terminal deglaciation:
  The global ICE-6G_C (VM5a) model, J. Geophys. Res. Solid Earth, 120, 450-487,
  https://doi.org/10.1002/2014JB011176, 2015.
- Yousefi, M., Wan, J., Pan, L., Gomez, N., Latychev, K., Mitrovica, J. X., Pollard, D., and DeConto, R. M.:
  The Influence of the Solid Earth on the Contribution of Marine Sections of the Antarctic Ice Sheet to
  Future Sea-Level Change, Geophys. Res. Lett., 49, e2021GL097525, https://doi.org/10.1029/2021gl097525, 2022.
- Yousefi, M., Wan, J., Pan, L., Gomez, N., Latychev, K., Mitrovica, J. X., Pollard, D., and DeConto, R. M.:
  Data for: The influence of the solid Earth on the contribution of marine sections of the Antarctic Ice
  Sheet to future sea-level change, Dryad [data set],
  https://datadryad.org/share/Zgxv_caT-dpJFY4mS9rFutjiQI-MQDefD5AvgZwlHc4, 2022.

**Next:** [Home](index.md).
