# Surface Density as the Organizing Principle of Little Red Dot Spectral Properties

**Author:** Tan Xin (Independent Researcher, IAIP)
**Status:** Submitted to *The Astrophysical Journal* (ApJ)
**Date:** May 2026

## Overview

This repository contains the complete data, code, and figures supporting the paper **"Surface Density as the Organizing Principle of Little Red Dot Spectral Properties"**, which presents a density-dependent effective gravity framework $G_{\rm eff}(\Sigma)$ as the organizing principle behind the spectral properties of JWST-observed Little Red Dots (LRDs) at $z \sim 4$--$9$.

### Key Findings

- **5.6$\sigma$ partial-correlation** between stellar surface density $\Sigma$ and $F444W - F150W$ color in 260 JWST LRDs, controlling for redshift and luminosity
- **34.6% color-reversal rate** in 1,626 control pairs: denser sources are systematically redder
- **KS test** yields $D = 0.574$ (6.2$\sigma$) between density-separated color distributions
- **JADES DR5 cross-validation** on 1,467 candidates confirms directional consistency
- **SPARC local-universe test** yields $\chi^2_{\rm red} = 122.3$ (zero detection), validating the high-redshift specificity prediction

## Repository Structure

```
.
├── manuscript/                    # Full paper and cover letter
│   ├── main.tex                   # ApJ manuscript source (aastex701)
│   ├── main.pdf                   # Compiled manuscript
│   ├── cover_letter.tex           # Submission cover letter
│   ├── cover_letter.pdf           # Compiled cover letter
│   └── aastex701.cls              # AAS LaTeX class file
│
├── code/                          # Analysis scripts (Python 3)
│   ├── Figure1_pubquality_pathA.py          # Fig 1: Sample overview
│   ├── control_pair_analysis.py             # Fig 3: Control-pair reversal analysis
│   ├── jades_cross_validation_v2.py         # Fig 6: JADES DR5 cross-validation
│   ├── sparc_GSigma_fitting.py              # Fig 7: SPARC local-universe test
│   ├── mcmc_GSigma_corner.py                # Fig D1: MCMC posterior
│   ├── magnitude_deepdive.py                # Fig E1: Gravitational redshift calculation
│   ├── vms_spurs_figure_v3.py               # Fig F1: SPURS VMS comparison
│   ├── plot_figure3_cosmos_pair.py          # Control pair visualization
│   ├── plot_figure5_phase_transition.py     # Phase transition diagram
│   ├── plot_master_paper.py                 # Master plotting utilities
│   ├── plot_rubies_validation.py            # RUBIES spectroscopic validation
│   ├── R1_R10_Robustness_Tests.py           # Robustness tests (R1-R10)
│   └── uid_validate.py                      # UID formula 5-layer validation
│
├── figures/
│   ├── main/                        # Paper figures (10 figures)
│   │   ├── Fig1_SampleOverview.png
│   │   ├── Fig2_Color_vs_Sigma.png
│   │   ├── Fig3_ControlPair_Reversal.png
│   │   ├── Fig4_FWHM_vs_Sigma.png
│   │   ├── Fig6_JADES_vs_COSMOS.png
│   │   ├── Fig7_SPARC_Chi2_ZeroDetection.png
│   │   ├── FigD1_MCMC_CornerPlot.png
│   │   ├── FigE1_Magnitude_DeepDive.png
│   │   └── FigF1_SPURS_IMF.png
│   └── robustness/                  # Supplementary robustness analysis
│       ├── Fig_mcmc_corner_GSigma.png
│       ├── Fig_mcmc_model_fit.png
│       ├── Fig_mcmc_spearman_posterior.png
│       ├── Fig_reversal_magnitude_histogram.png
│       ├── Fig_deltamag_mc_distribution.png
│       ├── Fig_R4_PerField_BarChart.png
│       ├── Fig_R9_LOWESS_Slope.png
│       ├── mcmc_results.json
│       ├── dmag_mc_results.json
│       ├── presubmission_analysis_results.json
│       └── reversal_histogram_results.json
│
└── data/
    ├── master_catalog/              # Primary LRD sample
    │   ├── LRD_Master_Combined_AllParams.csv   # 260-source master catalog
    │   ├── Kokorev_LRDs_Full.csv                # Original Kokorev et al. catalog
    │   ├── PlotData_SampleA_260.csv             # Sample A plotting data
    │   ├── PlotData_SampleB_36.csv              # Sample B plotting data
    │   ├── lrd_table_v1.1.fits                  # FITS binary table
    │   ├── LRD_density_analysis_results.csv     # Density analysis summary
    │   ├── LRD_StellarMass_Estimates.csv        # Stellar mass estimates
    │   ├── LRD_XRay_UpperLimits.csv             # Chandra X-ray upper limits
    │   └── lrd_magnitude_analysis_v2.csv        # Magnitude deep-dive analysis
    ├── control_pairs/               # Control pair analysis
    │   ├── FULL_ControlPairs_AllCandidates.csv  # 1,626 control pairs
    │   ├── control_pair_results.csv             # Pair matching results
    │   ├── MasterPaper_ResultsTable.csv         # Master paper results
    │   ├── Table1_Representative_Sources.csv    # Table 1 representative sources
    │   ├── KS_Test_Results_AllParams.csv        # KS test (all parameters)
    │   ├── KS_Test_Results_ThreeGroups.csv      # KS test (density groups)
    │   └── GSigma_Parameters.csv                # $G_{\rm eff}(\Sigma)$ fit parameters
    ├── rubies_validation/           # RUBIES spectroscopic validation
    │   ├── RUBIES_BroadBalmer_Catalog.csv
    │   ├── RUBIES_LRD_Catalog.csv
    │   ├── RUBIES_LRD_Summary.csv
    │   ├── RUBIES_LRD_xDJA_full.csv
    │   ├── RUBIES_LRD_xDJA_matched.csv
    │   └── RUBIES_Validation_Statistics.csv
    ├── jades_cross_validation/      # JADES DR5 independent test
    │   ├── JADES_LRD_candidates_v2.csv          # 1,467 JADES candidates
    │   ├── cross_validation_summary_v2.csv      # Cross-validation summary
    │   └── selection_cuts_v2.csv                # Applied selection criteria
    └── robustness_results/          # Additional robustness tests
        └── Triple_Scan_Results.csv              # Triple scan analysis results
```

## Data Sources

- **LRD catalog:** Kokorev et al. (2024, 2025), via the [DJA (Data Release of JADES-Webb)](https://dja.webb.community/)
- **Spectroscopy:** RUBIES survey (Goulding et al. 2024)
- **Cross-validation:** JADES DR5 (Eisenstein et al. 2025)
- **Local-universe test:** SPARC database (Lelli et al. 2016)
- **X-ray limits:** Chandra COSMOS-Legacy (Civano et al. 2016)

## How to Compile the Manuscript

```bash
cd manuscript/
pdflatex main.tex
pdflatex main.tex
pdflatex main.tex   # Triple pass for cross-references
```

## Requirements

- Python >= 3.8 with `numpy`, `scipy`, `matplotlib`, `astropy`, `pandas`, `emcee`
- LaTeX with `aastex701` class (included in `manuscript/`)

## Citation

If you use this code or data, please cite:

```bibtex
@article{Tan2026,
  author = {Tan, Xin},
  title = {Surface Density as the Organizing Principle of Little Red Dot Spectral Properties},
  journal = {The Astrophysical Journal},
  year = {2026},
  note = {Submitted}
}
```

## License

This repository is released under the [MIT License](LICENSE).

## Contact

Tan Xin | summicron352@gmail.com | [IAIP](https://iaip.org)
