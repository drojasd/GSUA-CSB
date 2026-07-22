# GSUA-CSB

Global Sensitivity and Uncertainty Analysis - Confidence Subcontour Box (GSUA-CSB) is a MATLAB toolbox for validating mathematical models implemented with Symbolic Math Toolbox or Simulink.

The toolbox supports:

- variance-based sensitivity analysis;
- uncertainty analysis;
- parameter estimation;
- confidence sub-contour box estimation for fitted parameters;
- visualization workflows for model-validation studies.

GSUA-CSB was developed at Universidad EAFIT and builds on previous work by Carlos Mario Velez on GSUA for dynamical systems using variance-based methods.

## Why It Matters

Mathematical models used in epidemiology, public health, engineering, and biological systems often depend on uncertain parameters. GSUA-CSB helps researchers understand which parameters matter, how uncertainty propagates through the model, and how identifiable fitted parameters are under available data.

## Links

- User guide and examples: [GSUA-CSB documentation](https://drojasd.github.io/GSUA-CSB/gsua_userguide)
- MATLAB File Exchange: [View GSUA-CSB on File Exchange](https://www.mathworks.com/matlabcentral/fileexchange/72637-gsua-csb)
- Latest GitHub release: [v1.6](https://github.com/drojasd/GSUA-CSB/releases/tag/v1.6)
- DOI: [Zenodo DOI](https://zenodo.org/badge/latestdoi/205731654)

[![DOI](https://zenodo.org/badge/205731654.svg)](https://zenodo.org/badge/latestdoi/205731654)

## Methods Implemented

Sensitivity index estimators implemented in this toolbox are based on:

1. Saltelli, A., Annoni, P., Azzini, I., Campolongo, F., Ratto, M., and Tarantola, S. (2010). Variance based sensitivity analysis of model output: design and estimator for the total sensitivity index. *Computer Physics Communications*, 181(2), 259-270.
2. Xiao, S., Lu, Z., and Wang, P. (2018). Multivariate global sensitivity analysis based on distance components decomposition. *Risk Analysis*, 38(12), 2703-2721.

## Citation

If you use this toolbox, cite:

Rojas-Diaz, Daniel and Velez-Sanchez, Carlos Mario (2019). GSUA-CSB (https://www.github.com/drojasd/GSUA-CSB), GitHub. doi:10.5755281/zenodo.3383316.

## License

This project is distributed under the MIT License.

