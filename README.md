# Three-Way Convergence in Little Red Dots

**Compactness-Correlated Spectroscopic Anomalies in Little Red Dots: Convergent Evidence from Six JWST Samples and an Engine-Degenerate Physical Framework**

Code, data, and figures accompanying the paper.

## Repository Structure

```
.
├── code/                   # All Python analysis scripts
│   ├── lrd_dyson_gradient_model.py      # Dyson dense cluster gradient model
│   ├── test_engine_degeneracy.py        # Engine degeneracy radiative transfer
│   ├── glimmir_variability_model.py     # Variability collision model
│   ├── dyson_no_imbh_explore.py         # Pure cluster core (no IMBH)
│   ├── dyson_parameter_scan.py          # Dyson parameter space scan
│   ├── cosmos_lrd_crossmatch.py         # COSMOSWeb LRD cross-matching
│   ├── generate_paper_figures.py        # Paper figure generation
│   ├── diagnose_sigma_definitions.py    # Sigma definition validation
│   ├── lrd_batch_geff_validation.py     # G_eff batch validation
│   ├── bird_stress_test.py              # BiRD Dyson stress test
│   ├── bird_stress_test_v2.py           # BiRD epsilon_g parameter scan
│   ├── loiacono_5sources_dyson.py       # 5-source Loiacono Dyson analysis
│   ├── loiacono_eps_threshold.py        # Epsilon_g threshold analysis
│   ├── weibel2026_analysis.py           # Weibel+ BH* analysis
│   ├── withers_analysis.py              # Withers+ REG analysis
│   ├── xiao_massive_sigma_estimate.py   # Xiao+ massive galaxy Sigma
│   └── nikopoulos_crossmatch.py         # Nikopoulos+ cross-match
│
├── data/                   # Processed data tables
│   ├── path1_merged_38sources.csv       # Path 1 cross-matched sample
│   ├── kokorev_260_sb.csv              # Kokorev 260 LRD catalog
│   ├── yanagisawa_overlap_15sources.csv # Yanagisawa overlap (15 sources)
│   ├── withers_9_regs.csv              # Withers+ 9 REGs
│   ├── lord_compactness.json           # LRD compactness catalog
│   ├── deGraaff2025_lrds.fits         # de Graaff+ LRD catalog
│   └── cosmos_lrd_crossmatch.json      # COSMOSWeb cross-match
│
├── figures/                # Key paper figures
├── results/                # Analysis results (JSON)
├── paper/                  # Manuscript and references
│   ├── main_paper.tex
│   ├── main_paper.pdf
│   └── references.bib
├── references/             # Reference papers
└── 三路汇聚论文大纲.md       # Paper outline (Chinese)
```

## Key Results

- **N=291 SB-Σ correlation**: ρ = −0.478, p = 5.3 × 10⁻¹⁸ — largest LRD compactness-spectroscopic sample
- **COSMOSWeb 664k blind verification**: Σ-color partial correlation ρ = +0.516 at z = 7–9 (>5.6σ)
- **Engine degeneracy**: Dense stellar cluster cores produce identical NIRCam color sequences to AGN engines
- **Dyson gradient model**: Viscous velocity cascade naturally produces observed FWHM without SMBH

## Citation

If you use this code or data, please cite the accompanying paper.

## License

MIT License
