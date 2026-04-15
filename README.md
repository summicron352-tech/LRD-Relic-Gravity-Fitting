# LRD Relic Gravity Fitting — MNRAS Submission

## A Redshift-Dependent Spectral Distortion in 260 Little Red Dots:
## Constraints on Compactness Evolution at Cosmic Dawn

**Status:** Submitted to MNRAS (2026-04-15)  
**Authors:** Xin Tan (Independent Researcher, IAIP)  
**arXiv:** [to be added]  
**GitHub Repository:** https://github.com/summicron352-tech/LRD-Relic-Gravity-Fitting/

---

### Quick Links
- 📄 [Paper PDF](paper/LRDs_Paper_Draft_v2.pdf)
- 📊 [Statistical Over-Constraint Figure](figures/Figure_StatisticalOverConstraint.png) *(new)*
- 📝 [Cover Letter (MNRAS)](submission/CoverLetter_MNRAS.pdf)

### Directory Structure
```
├── paper/          # LaTeX source, BibTeX, compiled PDF, FITS table
├── figures/        # All publication-quality figures (PNG/PDF)
├── data/
│   ├── csv/        # All analysis results & input catalogs
│   └── fits/       # Source catalog in FITS format
├── code/           # Python scripts for reproducibility
└── submission/     # Cover letter
```

### Key Findings
1. **83.5% of LRDs favor or accept** the gravitational redshift model ($z_{\rm dist}=0.17$) over dust attenuation
2. **Single mechanism explains six independent puzzles**: V-shaped SED, MIR silence, X-ray darkness, low dust mass, FWHM paradox, [N II] weakness
3. **Statistical over-constraint**: $P_{\rm joint}=1.8\times10^{-6}$ → predicts ~0 sources out of 260; observed **217** ($>100\sigma$)

### Reproducibility
All fitting codes are provided in `code/`. The main analysis pipeline:
```bash
python3 agn_three_comp_final.py          # AGN 3-component SED fitting
python3 gen_figure_overconstraint.py      # Statistical over-constraint figure
python3 triple_scan_classification.py     # Triple-scan classification
```

Input data: `data/csv/Kokorev_LRDs_Full.csv` (260 LRDs from Kokorev et al. 2024)

### Citation
If you use this work, please cite:
> Tan, X.\ 2026, MNRAS, submitted  
> *A Redshift-Dependent Spectral Distortion in 260 Little Red Dots: Constraints on Compactness Evolution at Cosmic Dawn*

### License
This repository is released under the MIT License.
