import numpy as np

G_N = 4.30091e-3
r_c, r_b = 0.003, 0.3
alpha_in, alpha_out = 2.0, 1.2
Sigma0, beta_v = 1e3, 0.7
r_grid = np.logspace(-4, 3, 2000)

sources = {
    'BiRD':    {'z': 2.33, 'FWHM': 2000, 'MBH_vir': 24e7,  'cls': 'LRD', 'b_opt': 2.71},
    'ID 1017': {'z': 2.50, 'FWHM': 1073, 'MBH_vir': 1.7e7, 'cls': 'mid', 'b_opt': 0.00},
    'ID 3646': {'z': 2.41, 'FWHM': 1237, 'MBH_vir': 3.5e7, 'cls': 'LBD', 'b_opt': -1.79},
    'ID 8511': {'z': 4.50, 'FWHM': 1668, 'MBH_vir': 1.6e7, 'cls': 'LRD', 'b_opt': 2.46},
    'ID 9008': {'z': 4.50, 'FWHM': 1067, 'MBH_vir': 0.83e7,'cls': 'LRD', 'b_opt': 1.88},
}

def find_eps_threshold(Mstar, r_eff, z, FWHM_target):
    v_target = FWHM_target / 2.355
    rho_c_g = 1e6
    rho = np.zeros_like(r_grid)
    inner = r_grid < r_b
    rho[inner] = rho_c_g * (r_grid[inner] / r_c)**(-alpha_in)
    outer = r_grid >= r_b
    rho_at_b = rho_c_g * (r_b / r_c)**(-alpha_in)
    rho[outer] = rho_at_b * (r_grid[outer] / r_b)**(-alpha_out)
    rho[r_grid < r_c*0.1] = rho_c_g * (0.1)**(-alpha_in)
    dr = np.diff(np.concatenate([[0], r_grid]))
    M_enc_init = np.cumsum(rho * dr * 4 * np.pi * r_grid**2)
    idx_eff = np.argmin(np.abs(r_grid - r_eff))
    rho *= Mstar / M_enc_init[idx_eff]
    M_enc = np.cumsum(rho * dr * 4 * np.pi * r_grid**2)
    Sigma_r = M_enc / (np.pi * r_grid**2)
    idx_blr = np.argmin(np.abs(r_grid - 0.03))
    Menc_blr = M_enc[idx_blr]
    Sig_blr = Sigma_r[idx_blr]
    v0_sq = G_N * Menc_blr / 0.03
    Sig_factor = (Sig_blr / Sigma0)**beta_v
    z_factor = (1+z)**3
    if v0_sq >= v_target**2:
        eps_threshold = 0.0
    else:
        eps_threshold = (v_target**2 / v0_sq - 1) / (z_factor * Sig_factor)
    return eps_threshold, Sig_blr, np.sqrt(v0_sq)

scenarios = [(5e9, 100), (1e10, 150), (2e10, 50)]

hdr = f"{'Source':>10s} {'cls':>6s} {'z':>5s} {'FWHM':>6s} {'b_opt':>8s} "
hdr += f"{'M*[e9]':>8s} {'re[pc]':>8s} {'Sig@0.03':>10s} {'eps_th':>10s} {'v_newt':>8s} {'OK?':>6s}"
print(hdr)
print('-' * len(hdr))

for name, s in sources.items():
    for Mstar, reff in scenarios:
        eps_th, Sig_blr, v_newt = find_eps_threshold(Mstar, reff, s['z'], s['FWHM'])
        ok = 'YES' if eps_th < 0.15 else ('maybe' if eps_th < 0.5 else 'no')
        print(f'{name:>10s} {s["cls"]:>6s} {s["z"]:>5.2f} {s["FWHM"]:>6d} {s["b_opt"]:+8.2f} '
              f'{Mstar/1e9:>8.1f} {reff:>8.0f} {Sig_blr:>10.1e} {eps_th:>10.5f} {v_newt:>8.0f} {ok:>6s}')
