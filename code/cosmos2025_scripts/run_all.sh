#!/bin/bash
# Run all COSMOS2025 analysis scripts in order
cd "$(dirname "$0")"
echo "=== 01: Extract Data ===" && python3 01_extract_data.py
echo "=== 02: ρ(z) Evolution ===" && python3 02_rho_evolution.py
echo "=== 03: Paired Analysis ===" && python3 03_paired_analysis.py
echo "=== 04: Environment Density ===" && python3 04_env_density.py
echo "=== 05: Tully-Fisher + Orientation ===" && python3 05_tully_fisher_orientation.py
echo "=== 06: Bulge+Disk Decomposition ===" && python3 06_bulge_disk_decomp.py
echo "=== 07: CIGALE SFH ===" && python3 07_cigale_sfh.py
echo "=== All done ==="
