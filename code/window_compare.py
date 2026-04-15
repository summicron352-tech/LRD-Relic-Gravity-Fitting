#!/usr/bin/env python3
"""
v10: Window Function Comparison Experiment
============================================
Test 5 different window forms against 260 LRDs.
Uses scipy.optimize.curve_fit for speed (vs slow grid search).
"""

import sys, csv, os, warnings
from math import log, exp, sqrt, isnan, isinf, isfinite
warnings.filterwarnings("ignore")

# ============================================================
# Paths & Constants
# ============================================================
BASE = "/Users/tanxin/Desktop/数据处理"
CSV_IN = os.path.join(BASE, "Kokorev_LRDs_Full.csv")
OUT_CSV = os.path.join(BASE, "WindowCompare_Results.csv")
OUT_TXT = os.path.join(BASE, "WindowCompare_Summary.txt")
OUT_PNG = os.path.join(BASE, "WindowCompare_Functions.png")

AG_OPT = 0.189
GAMMA_OPT = 0.50
BETA_OPT = 1.58
ZON_OPT = 2.0

FILTERS = [
    (0.90, "f090w"), (1.15, "f115w"), (1.50, "f150w"),
    (2.00, "f200w"), (2.77, "f277w"), (3.56, "f356w"), (4.44, "f444w"),
]
FILT_WAVE = [fw for fw, _ in FILTERS]
FILT_KEYS = [fk for _, fk in FILTERS]


def my_tanh(x):
    if x > 20: return 1.0
    if x < -20: return -1.0
    ex = exp(x)
    return (ex - 1/ex) / (ex + 1/ex)


def w1_tanh(z, AG=AG_OPT, beta=BETA_OPT, z_on=ZON_OPT):
    return AG * 0.5 * (my_tanh(beta * (z - z_on)) + 1)

def w2_gauss_single(z, AG=AG_OPT, z0=5.0, sigma=2.0):
    return AG * exp(-0.5 * ((z - z0) / sigma) ** 2)

def w3_gauss_double(z, AG=AG_OPT):
    g1 = 0.6 * exp(-0.5 * ((z - 4.4) / 1.0) ** 2)
    g2 = 0.4 * exp(-0.5 * ((z - 6.4) / 0.8) ** 2)
    return AG * (g1 + g2)

def w4_tanh_wide(z, AG=AG_OPT, beta=0.8, z_on=1.5):
    return AG * 0.5 * (my_tanh(beta * (z - z_on)) + 1)

def w5_tanh_decay(z, AG=AG_OPT, beta=BETA_OPT, z_on=ZON_OPT, z_cut=7.0):
    base = 0.5 * (my_tanh(beta * (z - z_on)) + 1)
    decay = exp(-0.5 * max(0, (z - z_cut) / 1.0) ** 2)
    return AG * base * decay


WINDOWS = [
    ("Tanh (baseline)",       w1_tanh,   "#2ecc71", "-"),
    ("Gaussian single",      w2_gauss_single, "#3498db", "--"),
    ("Gaussian double",      w3_gauss_double, "#e74c3c", "--"),
    ("Wide tanh",            w4_tanh_wide, "#9b59b6", ":"),
    ("Tanh + high-z decay",  w5_tanh_decay, "#f39c12", "-."),
]

# ============================================================
# SED model functions
# ============================================================

T_DUST = 60.0

def sed_uv(wave, Cb):
    """UV power law F ~ lambda^-0.5"""
    return Cb * (wave / 0.15) ** (-0.5)

def sed_ir(wave, Cr, zd=0.0):
    """Dust blackbody redshifted by zd"""
    lam_r = wave / (1.0 + zd)
    if lam_r <= 0:
        return 0.0
    x_hc_kt = 0.014387769 / lam_r / T_DUST
    if x_hc_kt > 100:
        return 0.0
    try:
        val = Cr * wave ** (-4.0) / (exp(x_hc_kt) - 1.0)
        if not isfinite(val):
            return 0.0
        return val
    except:
        return 0.0

def sed_total(params, waves, zd):
    Cb, Cr = params
    return [sed_uv(w, Cb) + sed_ir(w, Cr, zd) for w in waves]


def fit_one_source(fluxes, errors, waves, zd_pred):
    """Fit Null (zd=0) vs G_eff (zd=zd_pred), return stats."""
    from scipy.optimize import curve_fit
    
    nf = len(waves)
    p0 = [1.0, 10.0]  # [Cb, Cr]
    
    result_dict = {"chi2nu_null": 0, "chi2nu_geff": 0,
                   "delta_chi2_nu": 0, "delta_BIC": 0, "verdict": "error"}
    
    try:
        # --- Null model ---
        popt_n, _ = curve_fit(
            lambda w, Cb, Cr: sed_total([Cb, Cr], w, zd=0.0),
            waves, fluxes, p0=p0, sigma=errors,
            absolute_sigma=True, maxfev=3000
        )
        model_n = sed_total(popt_n, waves, zd=0.0)
        chi2_n = sum((fluxes[i] - model_n[i])**2 / errors[i]**2 for i in range(nf))
        
        # --- G_eff model ---
        popt_g, _ = curve_fit(
            lambda w, Cb, Cr: sed_total([Cb, Cr], w, zd=zd_pred),
            waves, fluxes, p0=p0, sigma=errors,
            absolute_sigma=True, maxfev=3000
        )
        model_g = sed_total(popt_g, waves, zd=zd_pred)
        chi2_g = sum((fluxes[i] - model_g[i])**2 / errors[i]**2 for i in range(nf))
        
        dof = max(nf - 2, 1)
        chi2nu_n = chi2_n / dof
        chi2nu_g = chi2_g / dof
        
        dc2nu = chi2nu_n - chi2nu_g
        dbic = (chi2_g + 2*log(nf)) - (chi2_n + 2*log(nf))
        
        if dc2nu >= 10:
            verdict = "Support"
        elif dc2nu > 2:
            verdict = "Weak"
        elif dc2nu >= -2:
            verdict = "Neutral"
        else:
            verdict = "Null wins"
        
        result_dict = {
            "chi2nu_null": round(chi2nu_n, 4),
            "chi2nu_geff": round(chi2nu_g, 4),
            "delta_chi2_nu": round(dc2nu, 4),
            "delta_BIC": round(dbic, 2),
            "verdict": verdict,
        }
    except Exception as e:
        pass
    
    return result_dict


# ============================================================
# Main
# ============================================================

print("=" * 70)
print("  v10: Window Comparison — 260 sources x 5 windows (curve_fit)")
print("=" * 70)

all_data = []
with open(CSV_IN) as f:
    reader = csv.DictReader(f)
    for row in reader:
        all_data.append(row)
print("  Read %d sources" % len(all_data))

all_results = {}
for wname, wfunc, color, ls in WINDOWS:
    print("\n  >> Computing: %s ..." % wname)
    results = []
    
    for ri, row in enumerate(all_data):
        rid = row.get("id", "?")
        rfld = row.get("field", "?")
        
        try:
            z_phot = float(row["z_phot"])
            
            # Collect photometry
            waves, fluxes, errors = [], [], []
            for fw, fkey in FILTERS:
                fk = fkey + "_flux"
                ek = fkey + "_fluxerr"
                fv = row.get(fk, "")
                ev = row.get(ek, "")
                if fv and fv.strip() != "" and fv.strip() != "nan":
                    fval = float(fv)
                    e_val = float(ev) if ev and ev.strip() != "" else max(fval * 0.1, 1e-30)
                    if fval > 0 and isfinite(fval) and e_val > 0:
                        obs_wave = fw * (1.0 + z_phot)
                        waves.append(obs_wave)
                        fluxes.append(fval)
                        errors.append(e_val)
            
            nw = len(waves)
            if nw < 4:
                results.append({"id": rid, "field": rfld, "z_phot": z_phot,
                               "n_bands": 0, "verdict": "skip",
                               "delta_chi2_nu": 0, "delta_BIC": 0})
                continue
            
            dG = wfunc(z_phot)
            zd = GAMMA_OPT * dG
            res = fit_one_source(fluxes, errors, waves, zd)
            res["id"] = rid
            res["field"] = rfld
            res["z_phot"] = z_phot
            res["n_bands"] = nw
            res["z_dist"] = round(zd, 6)
            res["delta_G"] = round(dG, 6)
            results.append(res)
        
        except Exception as e:
            print("    [ERR] %s: %s" % (rid, e))
            results.append({"id": rid, "field": rfld, "z_phot": -1,
                           "verdict": "error", "delta_chi2_nu": 0})
        
        if (ri + 1) % 50 == 0:
            print("    ... %d/%d" % (ri+1, len(all_data)))
    
    all_results[wname] = results
    valid = [r for r in results if r.get("verdict") in ("Support","Weak","Neutral","Null wins")]
    n_pos = sum(1 for r in valid if r["verdict"] in ("Support","Weak"))
    dcs = [r["delta_chi2_nu"] for r in valid]
    mn = sum(dcs)/len(dcs) if dcs else 0
    print("    -> N=%d, support=%d/%d(%d%%), Mean_dC2=%+.1f" % (
          len(valid), n_pos, len(valid),
          n_pos*100//max(len(valid),1), mn))

# ============================================================
# Output CSV
# ============================================================
print("\n  Writing %s" % OUT_CSV)
with open(OUT_CSV, "w") as f:
    header = ["id","field","z_phot"]
    safe_names = []
    for wn, _, _, _ in WINDOWS:
        sn = wn.replace(" ","_").replace("+","p").replace("-","m").replace("(","").replace(")","")
        safe_names.append(sn)
        header.extend([sn + "_dC2", sn + "_verdict"])
    f.write(",".join(header) + "\n")
    
    for i, row in enumerate(all_data):
        vals = [row.get("id","?"), row.get("field","?"), str(row.get("z_phot","0"))]
        for wn, _, _, _ in WINDOWS:
            r = all_results[wn][i]
            vals.append(str(r.get("delta_chi2_nu", 0)))
            vals.append(r.get("verdict", "?"))
        f.write(",".join(vals) + "\n")

# ============================================================
# Output text summary  
# ============================================================
print("  Writing %s" % OUT_TXT)
with open(OUT_TXT, "w") as f:
    f.write("=" * 72 + "\n")
    f.write("  Window Function Comparison Results\n")
    f.write("  Sample: 260 CEERS/UDS/GDS LRDs | gamma=0.50 fixed\n")
    f.write("=" * 72 + "\n\n")
    
    hdr = "%-28s %4s %5s %6s %6s %6s %6s %8s %7s\n" % (
           "Window", "N", "Sup", "Weak", "Neut", "Null", "%+", "MeanC2", "MedC2")
    f.write(hdr)
    f.write("-" * 88 + "\n")
    
    for wname, _, _, _ in WINDOWS:
        valid = [r for r in all_results[wname] if r.get("verdict") in 
                 ("Support","Weak","Neutral","Null wins")]
        ns = len(valid)
        counts = {}
        for v in ["Support", "Weak", "Neutral", "Null wins"]:
            counts[v] = sum(1 for r in valid if r["verdict"] == v)
        pct_pos = (counts["Support"]+counts["Weak"])*100/ns if ns > 0 else 0
        dcs = sorted([r["delta_chi2_nu"] for r in valid])
        mn = sum(dcs)/len(dcs) if dcs else 0
        med = dcs[len(dcs)//2] if dcs else 0
        f.write("%-28s %4d %5d %6d %6d %6d %4.1f%% %+8.1f %+7.1f\n" % (
               wname[:28], ns, counts["Support"], counts["Weak"],
               counts["Neutral"], counts["Null wins"], pct_pos, mn, med))
    
    # Redshift bins for baseline
    f.write("\n--- By redshift bin (Tanh baseline) ---\n")
    valid_base = [r for r in all_results["Tanh (baseline)"] 
                  if r.get("verdict") in ("Support","Weak","Neutral","Null wins")]
    bins = [(3.5,4.5,"z=3.5-4.5"),(4.5,5.5,"z=4.5-5.5"),(5.5,6.5,"z=5.5-6.5"),
            (6.5,7.5,"z=6.5-7.5"),(7.5,9.5,"z>7.5")]
    for lo, hi, label in bins:
        sub = [r for r in valid_base if lo <= float(r["z_phot"]) < hi]
        if sub:
            dcs = [r["delta_chi2_nu"] for r in sub]
            pos = sum(1 for d in dcs if d > 0)
            f.write("  %-12s N=%3d +:%d/%d(%d%%) mean=%+.1f\n" % (
                   label, len(sub), pos, len(sub),
                   pos*100//max(len(sub),1), sum(dcs)/len(dcs)))
    
    # Per-source table
    f.write("\n=== Per-source delta-chi2 across windows ===\n")
    f.write("%5s %5s" % ("ID", "z"))
    for sn in safe_names:
        f.write(" %14s" % sn[:14])
    f.write("\n")
    for i, row in enumerate(all_data):
        rid = row.get("id","?")
        rz = row.get("z_phot","0")
        f.write("%5s %5s" % (rid, rz))
        for wn, _, _, _ in WINDOWS:
            r = all_results[wn][i]
            d = r.get("delta_chi2_nu", 0)
            v = r.get("verdict", "?")[0]
            f.write(" %+7.1f[%s]" % (d, v))
        f.write("\n")


# ============================================================
# Plotting
# ============================================================
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec

    fig = plt.figure(figsize=(20, 18))
    gs = GridSpec(3, 3, figure=fig, hspace=0.35, wspace=0.3)

    # Panel A: Window shapes
    ax_a = fig.add_subplot(gs[0, :])
    z_grid = [0.1*i for i in range(10, 95)]
    for wname, wfunc, color, ls in WINDOWS:
        vals = [wfunc(z) for z in z_grid]
        ax_a.plot(z_grid, vals, color=color, linestyle=ls, linewidth=2.5,
                  label=wname, alpha=0.85)
    ax_a.set_xlabel("Redshift $z$", fontsize=13)
    ax_a.set_ylabel("$\\delta G / G_N$", fontsize=13)
    ax_a.set_title("Panel A: Five Window Function Forms Compared", fontsize=14, fontweight="bold")
    ax_a.legend(loc="upper left", fontsize=10, ncol=2)
    ax_a.axhline(0, color="gray", linewidth=0.5)
    ax_a.set_xlim(3, 9.5)
    ax_a.grid(True, alpha=0.3)

    # Panels B-F: Delta-chi2 per source
    positions = [(1,0), (1,1), (1,2), (2,0), (2,1)]
    for idx, (wname, wfunc, color, _) in enumerate(WINDOWS):
        pr, pc = positions[idx]
        ax = fig.add_subplot(gs[pr, pc])
        valid = [r for r in all_results[wname] if r.get("verdict") in 
                 ("Support","Weak","Neutral","Null wins")]
        dcs = [r["delta_chi2_nu"] for r in valid]
        bar_colors = []
        for d in dcs:
            if d >= 10: bar_colors.append("#27ae60")
            elif d > 2: bar_colors.append("#2ecc71")
            elif d >= -2: bar_colors.append("#f39c12")
            else: bar_colors.append("#e74c3c")
        ax.bar(range(len(valid)), dcs, color=bar_colors, alpha=0.7, edgecolor="none")
        ax.axhline(0, color="black", linewidth=1)
        ax.axhline(2, color="gray", linewidth=0.5, ls="--", alpha=0.5)
        ax.axhline(-2, color="gray", linewidth=0.5, ls="--", alpha=0.5)
        ax.set_xlabel("Source Index", fontsize=11)
        ax.set_ylabel("$\Delta\\chi^2_\\nu$", fontsize=12)
        n_pos = sum(1 for d in dcs if d > 0)
        pct = n_pos*100/len(valid) if valid else 0
        mn = sum(dcs)/len(dcs) if dcs else 0
        ax.title.set_text("%s) %s\nP+=%d%% Mean=%+.1f" % (chr(66+idx), wname, pct, mn))
        max_val = max(dcs) if dcs else 0
        min_val = min(dcs) if dcs else 0
        if max_val - min_val < 5:
            ax.set_ylim(-5, 5)
        else:
            ylo = min(-5, min_val - 5)
            yhi = max(5, max_val * 1.2)
            ax.set_ylim(ylo, yhi)
        ax.tick_params(labelsize=8)

    # Panel G: Summary bar chart
    ax_g = fig.add_subplot(gs[2, 2])
    names, pcts_list, means = [], [], []
    for wname, _, color, _ in WINDOWS:
        valid = [r for r in all_results[wname] if r.get("verdict") in 
                 ("Support","Weak","Neutral","Null wins")]
        pos = sum(1 for r in valid if r["verdict"] in ("Support", "Weak"))
        pcts_list.append(pos*100/max(len(valid),1))
        means.append(sum(r["delta_chi2_nu"] for r in valid)/max(len(valid),1))
        names.append(wname.replace(" ", "\n"))

    xp = range(len(names))
    bars = ax_g.bar(xp, pcts_list, color=[w[2] for w in WINDOWS], alpha=0.75, edgecolor="black")
    ax_g.set_xticks(xp)
    ax_g.set_xticklabels(names, fontsize=9)
    ax_g.set_ylabel("G_eff Support Rate (%)", fontsize=12)
    ax_g.title.set_text("G) Support Rate Comparison")
    ymax = max(pcts_list)*1.25 if max(pcts_list) > 0 else 20
    ax_g.set_ylim(0, ymax)
    ax_g.axhline(pcts_list[0], color="#2ecc71", ls=":", lw=1.5, alpha=0.7,
                 label="Tanh baseline: %.0f%%" % pcts_list[0])
    ax_g.legend(fontsize=9)
    for i, (p, m) in enumerate(zip(pcts_list, means)):
        ax_g.text(i, p+1, "%.0f%%\nm=%+.1f" % (p, m), ha="center", va="bottom", 
                 fontsize=8, fontweight="bold")

    plt.suptitle("Window Function Comparison\n260 LRDs x 5 Window Forms",
                 fontsize=16, fontweight="bold", y=0.95)
    plt.savefig(OUT_PNG, dpi=180, bbox_inches="tight")
    plt.close()
    print("  Saved: %s" % OUT_PNG)
except ImportError:
    print("  [!] No matplotlib")

print("\n  Done!")
