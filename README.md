# LRD Relic Gravity Fitting

**Two complementary papers on JWST Little Red Dots (LRDs) by Xin Tan (Independent Researcher, IAIP)**

---

## 📄 Paper I — MNRAS (Dual-Probe Paper)

### *A Redshift-Dependent Spectral Distortion in 260 Little Red Dots: Constraints on Compactness Evolution at Cosmic Dawn*

**Status:** Submitted to MNRAS (MN-26-1077-P, 2026-04-15) · Awaiting Reviewer Selection  
**GitHub Repo:** [Dense-Environments-Evidence-from-260-High-Redshift-Little-Red-Dots](https://github.com/summicron352-tech/Dense-Environments-Evidence-from-260-High-Redshift-Little-Red-Dots)

**Quick Links:** [Paper PDF](paper/LRDs_Paper_Draft_v2.pdf) · [Cover Letter](submission/CoverLetter_MNRAS.pdf)

### Key Findings
1. **83.5% of LRDs favor or accept** the gravitational redshift model ($z_{\rm dist}=0.17$) over dust attenuation
2. **Single mechanism explains six independent puzzles**: V-shaped SED, MIR silence, X-ray darkness, low dust mass, FWHM paradox, [N II] weakness
3. **Statistical over-constraint**: $P_{\rm joint}=1.8\times10^{-6}$ → predicts ~0 sources out of 260; observed **217** ($>100\sigma$)

### Reproducibility (MNRAS Paper)
```bash
python3 agn_three_comp_final.py          # AGN 3-component SED fitting
python3 gen_figure_overconstraint.py     # Statistical over-constraint figure
python3 triple_scan_classification.py     # Triple-scan classification
```

---

## 📄 Paper II — APJL (Companion Letter)

### *Density-Dependent Spectral Distortion in Extreme Dense Environments: Evidence from 260 High-Redshift Little Red Dots*

**Status:** Submitted to *The Astrophysical Journal Letters* (APJL)  
**GitHub Repo:** [Dense-Environments-Evidence-from-260-High-Redshift-Little-Red-Dots](https://github.com/summicron352-tech/Dense-Environments-Evidence-from-260-High-Redshift-Little-Red-Dots)

**Quick Links:** [Paper PDF](paper_apjl/PAPER_DRAFT_EN.pdf) · [Outline](paper_apjl/APJL_OUTLINE_v3.md) · [中文说明](paper_apjl/PAPER_DRAFT_ZH.md)

### Key Results

| Metric | Value |
|--------|-------|
| Partial correlation (ρₚ) | +0.341 |
| Significance | **5.6σ** |
| KS statistic (D) | 0.574 |
| KS significance | **6.2σ** |
| Control (F444W/F356W) | < 1σ ✓ |

### Reproducibility (APJL Paper)
```bash
python3 code/apjl_pathA_flux_ratio_analysis.py  # Partial correlation analysis
python3 code/apjl_pathB_free_zdist_fitting.py    # Free z_dist fitting
python3 code/apjl_Figure1_pubquality.py         # Generate Figure 1
```

---

## 🗂️ Directory Structure

```
LRD-Relic-Gravity-Fitting/
├── paper/              # MNRAS paper (LaTeX + PDF)
├── paper_apjl/         # APJL letter (LaTeX + PDF)
│   ├── PAPER_DRAFT_EN.tex
│   ├── PAPER_DRAFT_EN.pdf
│   └── ...
├── figures/            # All figures
│   ├── *               # MNRAS paper figures
│   └── apjl_*         # APJL letter figures
├── data/
│   ├── csv/           # MNRAS analysis results
│   ├── fits/          # FITS catalog
│   └── Kokorev_LRDs_Full.csv  # LRD catalog
├── code/
│   ├── *              # MNRAS analysis scripts
│   └── apjl_*         # APJL analysis scripts
└── submission/        # MNRAS cover letter
```

---

## 📦 Data

All analysis uses the public LRD catalog from **Kokorev et al. (2024)**, via [GitHub](https://github.com/Vlas-Sokolov/Kokorev-LRD-catalog).

---

## 📚 Citation

**MNRAS Paper:**
> Tan, X. 2026, MNRAS, submitted  
> *A Redshift-Dependent Spectral Distortion in 260 Little Red Dots: Constraints on Compactness Evolution at Cosmic Dawn*

**APJL Letter:**
> Tan, X. 2026, ApJL, submitted  
> *Density-Dependent Spectral Distortion in Extreme Dense Environments: Evidence from 260 High-Redshift Little Red Dots*

---

## ⚠️ Disclaimer

This repository is archived on **Zenodo**:
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.19587085-blue)](https://doi.org/10.5281/zenodo.19587085)

---

## 📜 License

MIT License.
