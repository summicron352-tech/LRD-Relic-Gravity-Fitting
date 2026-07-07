"""
COSMOSWeb 660k × LRD 样本交叉匹配
====================================
目的: 证明LRD样本的Σ-color相关性不是挑选偏倚导致的

策略:
  1. 从COSMOSWeb 664k全样本中提取与LRD同Σ/z/M★参数空间的子样本
  2. 对比LRD子样本 vs COSMOS匹配子样本的Σ-color偏相关
  3. 如果两者一致 → 数据本身在说话，不是挑数据
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
# 中文字体配置
plt.rcParams['font.sans-serif'] = ['Hiragino Sans GB', 'Arial Unicode MS', 'Heiti TC', 'PingFang HK']
plt.rcParams['axes.unicode_minus'] = False
from scipy.stats import spearmanr
from scipy.stats import rankdata
import os, json

# ===========================
# 路径
# ===========================
data_dir = '/Users/tanxin/Desktop/数据处理/14_ThreeWay_Convergence_Paper/data'
fig_dir = '/Users/tanxin/Desktop/数据处理/14_ThreeWay_Convergence_Paper/figures'
cosmos_dir = '/Users/tanxin/Desktop/数据处理/COSMOS2025'
os.makedirs(fig_dir, exist_ok=True)

print("=" * 70)
print("COSMOSWeb 660k × LRD 交叉匹配")
print("=" * 70)

# ===========================
# 1. 加载 COSMOSWeb 664k
# ===========================
print("\n1. 加载 COSMOSWeb 664k 数据...")
cosmos = np.load(os.path.join(cosmos_dir, 'cosmos2025_extracted.npz'))
z_cosmos = cosmos['z']
logM_cosmos = cosmos['logM']
logSigma_cosmos = cosmos['logSigma']
color_cosmos = cosmos['color']
reff_cosmos = cosmos['reff_pc']
m444_cosmos = cosmos['m444']
N_cosmos = len(z_cosmos)

print(f"   COSMOSWeb: {N_cosmos} 源")
print(f"   z: [{z_cosmos.min():.2f}, {z_cosmos.max():.2f}]")
print(f"   logM★: [{logM_cosmos.min():.2f}, {logM_cosmos.max():.2f}]")
print(f"   logΣ: [{logSigma_cosmos.min():.2f}, {logSigma_cosmos.max():.2f}]")
print(f"   color: [{color_cosmos.min():.2f}, {color_cosmos.max():.2f}]")

# ===========================
# 2. 偏Spearman相关函数
# ===========================
def psc(x, y, ctrls):
    """偏Spearman秩相关，控制ctrls变量"""
    rx = rankdata(x)
    ry = rankdata(y)
    for c in ctrls:
        rc = rankdata(c)
        b = np.cov(rx, rc)[0, 1] / np.var(rc)
        rx = rx - b * rc
        b = np.cov(ry, rc)[0, 1] / np.var(rc)
        ry = ry - b * rc
    return spearmanr(rx, ry)[0]

def bootstrap_psc(x, y, ctrls, n_boot=500):
    """Bootstrap偏相关"""
    n = len(x)
    rho_boot = []
    for _ in range(n_boot):
        idx = np.random.choice(n, n, replace=True)
        try:
            r = psc(x[idx], y[idx], [c[idx] for c in ctrls])
            rho_boot.append(r)
        except:
            pass
    rho_boot = np.array(rho_boot)
    return np.mean(rho_boot), np.std(rho_boot)

# ===========================
# 3. COSMOS全样本 Σ-color 按红移箱
# ===========================
print("\n2. COSMOS全样本 Σ-color 偏相关 (控制z, logM★)")
print("-" * 50)

z_bins = [(0,1), (1,3), (3,5), (5,7), (7,9), (9,15)]
z_bin_labels = ['0-1', '1-3', '3-5', '5-7', '7-9', '9-15']
full_results = []

for (z_lo, z_hi), label in zip(z_bins, z_bin_labels):
    mask = (z_cosmos >= z_lo) & (z_cosmos < z_hi)
    n_src = mask.sum()
    
    if n_src < 10:
        continue
    
    rho_psc = psc(logSigma_cosmos[mask], color_cosmos[mask], 
                  [z_cosmos[mask], logM_cosmos[mask]])
    
    full_results.append({
        'z_range': label, 'N': n_src,
        'rho': rho_psc,
        'logSigma_med': np.median(logSigma_cosmos[mask]),
        'logM_med': np.median(logM_cosmos[mask]),
        'z_med': np.median(z_cosmos[mask])
    })
    
    print(f"   z={label:>5s}: N={n_src:>6d}  ρ={rho_psc:>+7.4f}  "
          f"logΣ_med={np.median(logSigma_cosmos[mask]):.2f}  "
          f"z_med={np.median(z_cosmos[mask]):.2f}")

# ===========================
# 4. 加载 LRD 样本并匹配
# ===========================
print("\n3. LRD样本加载与COSMOSWeb参数空间匹配")
print("-" * 50)

# 加载各LRD样本
lrd_samples = {}

# Path 1: 38 sources
try:
    path1 = np.genfromtxt(
        os.path.join(data_dir, 'path1_merged_38sources.csv'),
        delimiter=',', names=True, dtype=None, encoding='utf-8')
    path1_keys = path1.dtype.names
    lrd_samples['Path1_38'] = path1
    print(f"   ✓ Path1 (38源): 加载成功, {len(path1_keys)}列")
except Exception as e:
    print(f"   ✗ Path1: 加载失败 - {e}")

# 尝试用deGraaff FITS
try:
    from astropy.io import fits
    fits_file = os.path.join(data_dir, 'deGraaff2025_lrds.fits')
    if os.path.exists(fits_file):
        hdul = fits.open(fits_file)
        deg_data = hdul[1].data
        lrd_samples['deGraaff_LRDs'] = deg_data
        print(f"   ✓ deGraaff LRDs: {len(deg_data)}源, {len(deg_data.columns)}列")
        hdul.close()
except Exception as e:
    print(f"   ✗ deGraaff: 加载失败 - {e}")

# 尝试Kokorev 260源
try:
    kokorev = np.genfromtxt(
        os.path.join(data_dir, 'kokorev_260_sb.csv'),
        delimiter=',', names=True, dtype=None, encoding='utf-8')
    lrd_samples['Kokorev_260'] = kokorev
    print(f"   ✓ Kokorev (260源): 加载成功, {len(kokorev.dtype.names)}列")
except Exception as e:
    print(f"   ✗ Kokorev: 加载失败 - {e}")

# 如果都没有，用已有的CSV
if len(lrd_samples) == 0:
    print("\n   尝试备用加载方案...")
    for fname in ['apjl_38sources_sb.csv', 'yanagisawa_overlap_15sources.csv']:
        try:
            data = np.genfromtxt(
                os.path.join(data_dir, fname),
                delimiter=',', names=True, dtype=None, encoding='utf-8')
            lrd_samples[fname] = data
            print(f"   ✓ {fname}: 加载成功, {len(data)}源")
        except Exception as e:
            print(f"   ✗ {fname}: {e}")

# ===========================
# 5. 匹配分析
# ===========================
if len(lrd_samples) > 0:
    print(f"\n4. LRD vs COSMOS匹配子样本 Σ-color 对比")
    print("-" * 60)
    
    # 创建一个表格放结果
    all_matches = []
    
    for lrd_name, lrd_data in lrd_samples.items():
        print(f"\n   === {lrd_name} ===")
        
        # 提取LRD的z, logM, logSigma, color
        z_lrd = None
        logM_lrd = None
        logSigma_lrd = None
        color_lrd = None
        
        # 探测列名
        if hasattr(lrd_data, 'dtype') and lrd_data.dtype.names is not None:
            col_names = list(lrd_data.dtype.names)
        elif isinstance(lrd_data, np.ndarray) and hasattr(lrd_data, 'dtype') and lrd_data.dtype.names:
            col_names = list(lrd_data.dtype.names)
        else:
            col_names = []
        
        # col_names as list for .fits
        cn = list(col_names) if not isinstance(col_names, list) else col_names
        cn_lower = [c.lower() for c in cn]
        
        # 显式列名映射（按优先级）
        # z 映射
        for z_key in ['z', 'zfinal', 'zspec', 'z_phot', 'zphot', 'redshift']:
            try:
                idx = cn_lower.index(z_key)
                z_lrd = np.array(lrd_data[cn[idx]], dtype=float)
                break
            except (ValueError, IndexError):
                pass
        
        # logM 映射
        for m_key in ['logm', 'logmstar', 'logMstar_best', 'logmstar_best', 'mass_med', 'mass']:
            try:
                idx = cn_lower.index(m_key)
                logM_lrd = np.array(lrd_data[cn[idx]], dtype=float)
                break
            except (ValueError, IndexError):
                pass
        
        # logSigma 映射
        for s_key in ['logsigma', 'logsigma_mstar', 'logSigma_Mstar', 'logsigma_Mstar']:
            try:
                idx = cn_lower.index(s_key)
                logSigma_lrd = np.array(lrd_data[cn[idx]], dtype=float)
                print(f"     找到 logSigma: {s_key} → 值范围 [{logSigma_lrd.min():.2f}, {logSigma_lrd.max():.2f}]")
                break
            except (ValueError, IndexError):
                pass
        
        # 如果没有直接的logSigma，用logM和r_eff计算
        if logSigma_lrd is None:
            r_eff_lrd = None
            for r_key in ['r_eff_50_phys', 'reff_pc', 'r_eff', 'radius']:
                try:
                    idx = cn_lower.index(r_key)
                    r_eff_lrd = np.array(lrd_data[cn[idx]], dtype=float)
                    print(f"     找到 r_eff: {r_key} → 用于计算 logSigma")
                    break
                except (ValueError, IndexError):
                    pass
            if logM_lrd is not None and r_eff_lrd is not None:
                logSigma_lrd = logM_lrd - np.log10(np.pi * r_eff_lrd**2)
                print(f"     计算 logSigma = logM - log10(π r^2)")
        
        # color 映射
        for c_key in ['color', 'f150w_f444w']:
            try:
                idx = cn_lower.index(c_key)
                color_lrd = np.array(lrd_data[cn[idx]], dtype=float)
                break
            except (ValueError, IndexError):
                pass
        
        # 如果color不存在，用flux计算
        if color_lrd is None:
            f150w = None
            f444w = None
            for k in ['mag_auto_f150w', 'f150w_mag', 'f150w']:
                if k in cn_lower:
                    idx = cn_lower.index(k)
                    f150w = np.array(lrd_data[cn[idx]], dtype=float)
                    break
            for k in ['mag_auto_f444w', 'f444w_mag', 'f444w']:
                if k in cn_lower:
                    idx = cn_lower.index(k)
                    f444w = np.array(lrd_data[cn[idx]], dtype=float)
                    break
            if f150w is not None and f444w is not None:
                color_lrd = f150w - f444w
                print(f"     计算 color = F150W - F444W")
        
        if z_lrd is None or logM_lrd is None:
            print(f"   ⚠ 无法从列名识别 z/logM: {cn[:15]}...")
            continue
        
        # 转换为float数组
        z_lrd = np.array(z_lrd, dtype=float)
        logM_lrd = np.array(logM_lrd, dtype=float)
        
        # 过滤: NaN + 物理上合理的范围
        valid = np.isfinite(z_lrd) & np.isfinite(logM_lrd)
        valid &= (z_lrd > 0) & (z_lrd < 20)
        valid &= (logM_lrd > 6) & (logM_lrd < 14)  # 排除 SED 拟合失败
        
        if logSigma_lrd is not None:
            valid &= np.isfinite(logSigma_lrd)
            valid &= (logSigma_lrd > -5) & (logSigma_lrd < 15)
        
        if valid.sum() < len(z_lrd) * 0.8:
            print(f"     已过滤 {len(z_lrd) - valid.sum()}/{len(z_lrd)} 个异常源")
        
        z_lrd = z_lrd[valid]
        logM_lrd = logM_lrd[valid]
        if logSigma_lrd is not None:
            logSigma_lrd = logSigma_lrd[valid]
        
        if len(z_lrd) < 5:
            print(f"   ⚠ 仅有{len(z_lrd)}个有效源，跳过")
            continue
        
        # LRD的红移箱分布
        for (z_lo, z_hi), label in zip(z_bins, z_bin_labels):
            mask_lrd = (z_lrd >= z_lo) & (z_lrd < z_hi)
            n_lrd = mask_lrd.sum()
            if n_lrd < 3:
                continue
            
            z_lrd_bin = z_lrd[mask_lrd]
            logM_lrd_bin = logM_lrd[mask_lrd]
            
            z_med = np.median(z_lrd_bin)
            logM_med = np.median(logM_lrd_bin)
            
            # 从COSMOS 660k中, 匹配相同z/M范围
            # 放宽匹配边界(50%更宽以保证统计量)
            z_range = max(0.5, (z_lrd_bin.max() - z_lrd_bin.min()) * 1.5)
            z_center = z_med
            M_range = max(0.5, (logM_lrd_bin.max() - logM_lrd_bin.min()) * 1.5)
            M_center = logM_med
            
            match_mask = (
                (z_cosmos >= z_center - z_range/2) & 
                (z_cosmos < z_center + z_range/2) &
                (logM_cosmos >= M_center - M_range/2) & 
                (logM_cosmos < M_center + M_range/2)
            )
            n_match = match_mask.sum()
            
            if n_match < 10:
                continue
            
            # COSMOS匹配子样本的偏相关
            rho_match = psc(logSigma_cosmos[match_mask], color_cosmos[match_mask],
                           [z_cosmos[match_mask], logM_cosmos[match_mask]])
            
            # COSMOS全红移箱的偏相关（参考）
            rho_full = None
            for fr in full_results:
                if fr['z_range'] == label:
                    rho_full = fr['rho']
                    break
            
            # LRD样本自身的Σ-color相关
            # 如果没有颜色数据，只报告匹配结果
            if color_lrd is not None and logSigma_lrd is not None:
                color_lrd_arr = np.array(color_lrd, dtype=float)
                logSigma_lrd_arr = np.array(logSigma_lrd, dtype=float)
                valid_lrd = np.isfinite(color_lrd_arr) & np.isfinite(logSigma_lrd_arr)
                if valid_lrd.sum() >= 3:
                    rho_lrd_direct = spearmanr(logSigma_lrd_arr[valid_lrd], 
                                              color_lrd_arr[valid_lrd])[0]
                else:
                    rho_lrd_direct = None
            else:
                rho_lrd_direct = None
            
            print(f"   z={label:>5s} | LRD({n_lrd}) vs COSMOS匹配({n_match})")
            print(f"     COSMOS全ρ={rho_full:>+.4f} | 匹配ρ={rho_match:>+.4f} | "
                  f"LRD直接ρ={rho_lrd_direct if rho_lrd_direct else 'N/A':>7}")
            
            all_matches.append({
                'sample': lrd_name,
                'z_bin': label,
                'n_lrd': n_lrd,
                'n_match': n_match,
                'rho_cosmos_full': rho_full,
                'rho_cosmos_match': rho_match,
                'rho_lrd_direct': rho_lrd_direct,
                'z_med': z_med,
                'logM_med': logM_med
            })
    
    # ===========================
    # 6. 可视化
    # ===========================
    print(f"\n5. 生成可视化...")
    
    # 图1: Σ-color散点图 — COSMOS全样本 + LRD叠加
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    axes = axes.flatten()
    
    for i, ((z_lo, z_hi), label) in enumerate(zip(z_bins, z_bin_labels)):
        if i >= 6:
            break
        ax = axes[i]
        
        mask = (z_cosmos >= z_lo) & (z_cosmos < z_hi)
        
        # 抽取子样本（太多点画不下，随机取2000）
        sub_idx = np.random.choice(np.where(mask)[0], min(2000, mask.sum()), replace=False)
        ax.scatter(logSigma_cosmos[sub_idx], color_cosmos[sub_idx], 
                  c='lightgray', s=1, alpha=0.3, label=f'COSMOS ({mask.sum()})')
        
        # 添加拟合线
        ax.set_xlabel('log Σ [Msun/pc^2]', fontsize=11)
        ax.set_ylabel('F150W-F444W', fontsize=11)
        ax.set_title(f'z={label}  (N={mask.sum()})', fontsize=12)
        ax.grid(True, alpha=0.2)
    
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, 'cosmos_vs_lrd_scatter.png'), dpi=150)
    plt.close()
    print("  ✓ cosmos_vs_lrd_scatter.png")
    
    # 图2: 偏相关随红移演化 — COSMOS全 vs 匹配
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # COSMOS全样本
    z_centers = [(fr['z_med']) for fr in full_results]
    rho_vals = [fr['rho'] for fr in full_results]
    n_vals = [fr['N'] for fr in full_results]
    z_labels = [fr['z_range'] for fr in full_results]
    
    ax.plot(z_centers, rho_vals, 'o-', color='blue', linewidth=2, 
            markersize=8, label='COSMOS全样本', zorder=5)
    
    # 标注N
    for zc, rv, nv in zip(z_centers, rho_vals, n_vals):
        ax.annotate(f'N={nv}', (zc, rv), textcoords='offset points', 
                   xytext=(0, 10), fontsize=8, ha='center')
    
    # LRD匹配子样本
    if len(all_matches) > 0:
        match_z = []
        match_rho = []
        match_sample = []
        for m in all_matches:
            if m['rho_cosmos_match'] is not None:
                match_z.append(m['z_med'])
                match_rho.append(m['rho_cosmos_match'])
                match_sample.append(m['sample'])
        
        if len(match_z) > 0:
            ax.scatter(match_z, match_rho, marker='s', s=80, 
                      c='red', edgecolors='darkred', 
                      label='COSMOS匹配LRD子样本', zorder=6)
    
    ax.set_xlabel('红移 z', fontsize=13)
    ax.set_ylabel('偏Spearman ρ(Σ, color | z, M★)', fontsize=13)
    ax.set_title('Σ-color偏相关随红移演化: COSMOS全样本 vs LRD匹配子样本', fontsize=14)
    ax.axhline(0, color='gray', linestyle='--', alpha=0.5)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, 'cosmos_vs_lrd_rho_evolution.png'), dpi=150)
    plt.close()
    print("  ✓ cosmos_vs_lrd_rho_evolution.png")

    # ===========================
    # 保存结果
    # ===========================
    results_out = {
        'full_cosmos': full_results,
        'lrd_matches': all_matches
    }
    
    # 改为标准json可序列化
    def make_serializable(obj):
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj) if np.isfinite(obj) else None
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return obj
    
    # 手动构建可序列化结构
    serializable = {
        'full_cosmos': [{k: make_serializable(v) for k, v in fr.items()} for fr in full_results],
        'lrd_matches': [{k: make_serializable(v) for k, v in m.items()} for m in all_matches]
    }
    
    with open(os.path.join(data_dir, 'cosmos_lrd_crossmatch.json'), 'w') as f:
        json.dump(serializable, f, indent=2)
    print(f"\n   ✓ 结果保存至 cosmos_lrd_crossmatch.json")

else:
    print("\n⚠ 无有效LRD样本加载，跳过匹配分析")
    print("  请检查数据路径和格式")

print("\n" + "=" * 70)
print("交叉匹配完成")
print("=" * 70)
