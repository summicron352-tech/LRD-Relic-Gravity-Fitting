# LRD Relic Gravity Fitting Pipeline

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.xxxxxxx.svg)](https://doi.org/10.5281/zenodo.xxxxxxx)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

This repository contains the data and code used in the paper: **"A Redshift-Dependent Spectral Distortion in 260 Little Red Dots: Constraints on Compactness Evolution at Cosmic Dawn"** by Xin Tan.

## Overview
This pipeline systematically evaluates the "Gravitational Relic Hypothesis" against standard dust attenuation models across 260 confirmed James Webb Space Telescope (JWST) "Little Red Dots" (LRDs). By incorporating a time-dependent effective gravitational metric ($G_{\rm eff}$), it effectively bypasses the classical mid-infrared and X-ray radiation paradoxes associated with these extreme high-redshift objects.

## Repository Contents

- `lrd_fitting_pipeline.py`: The core Python script performing the 5-parameter phenomenological SED fitting. This script compares both standard chromatic dust pathways ($A_V$) and the achromatic gravitational relic displacement ($z_{dist}$).
- `Bulletproof_Results.csv`: The complete statistical catalog containing best-fit parameters ($A_V$, $z_{dist}$, $\alpha_{\rm AGN}$, $\Delta$BIC, $\Delta \chi^2_\nu$) for all 260 LRDs analyzed in this study.
- `Table1_Representative_Sources.csv`: A subset of 20 highly representative sources showcasing the statistical dichotomy between relic supporters and conventional starbursts.

## Requirements

To run this pipeline locally, you will need Python 3.8+ and standard scientific computing libraries.

```bash
pip install numpy pandas scipy matplotlib astropy
```

## Reproducibility

The baseline empirical photometric catalog is derived from Kokorev et al. (2024).

1. Clone this repository:
   ```bash
   git clone https://github.com/xxxx/lrd-relic-gravity.git
   cd lrd-relic-gravity
   ```
2. Execute the fitting pipeline:
   ```bash
   python lrd_fitting_pipeline.py
   ```

Running the provided script will reproduce the quantitative multi-panel SED fits (as seen in Figure 3 of the paper) and generate the statistical metrics used to construct the cumulative distribution functions.

## Citation

If you use this code or data in your research, please cite the corresponding paper and the Zenodo archive:

```bibtex
@article{Tan2026LRD,
  title={A Redshift-Dependent Spectral Distortion in 260 Little Red Dots: Constraints on Compactness Evolution at Cosmic Dawn},
  author={Tan, Xin},
  journal={arXiv preprint (or proper journal after acceptance)},
  year={2026}
}
```

*Note to users: The Zenodo DOI linking to the persistent release of this specific code version can be found in the badge at the top of this document.*
