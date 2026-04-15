#!/usr/bin/env python3
"""
FINAL approach: 3-phase MIRI F770W coverage check.
Phase 1: query_region → obs table with positions (PROVEN)
Phase 2: batch get_product_list for MIRI obs only → filter info  
Phase 3: source-to-obs position match using Phase 1 data
"""
import os
for k in ['HTTP_PROXY','HTTPS_PROXY','http_proxy','https_proxy','ALL_PROXY']:
    os.environ.pop(k, None)

import pandas as pd
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table
from astroquery.mast import Observations
import time, warnings, sys
warnings.filterwarnings('ignore')

DATA_DIR = "/Users/tanxin/Desktop/数据处理"
catalog = pd.read_csv(f"{DATA_DIR}/Kokorev_LRDs_Full.csv")
fields = catalog['field'].unique()

# ============================================================
# PHASE 1: Get ALL obs per field + store positions
# ============================================================
print("="*70)
print("PHASE 1: Querying MAST regions...")
print("="*70)

all_obs_data = []  # {field, obsid, s_ra, s_dec, instrument}

for fi, field in enumerate(fields):
    fs = catalog[catalog['field'] == field]
    ra_c, dec_c = fs['ra'].mean(), fs['dec'].mean()
    
    print(f"\n[{fi+1}/{len(fields)}] {field} ({len(fs)} sources)", end='', flush=True)
    
    try:
        obs = Observations.query_region(
            coordinates=f"{ra_c} {dec_c}", radius="8 arcmin"
        )
        
        # Filter MIRI locally
        instr = [str(s).upper() for s in obs['instrument_name']]
        for i, row in enumerate(obs):
            if 'MIRI' in instr[i]:
                all_obs_data.append({
                    'field': field,
                    'obsid': int(row['obsid']),
                    's_ra': float(row['s_ra']),
                    's_dec': float(row['s_dec']),
                })
        
        n_miri = sum(1 for d in all_obs_data if d['field'] == field)
        print(f" → total={len(obs)}, MIRI={n_miri}")
        
    except Exception as e:
        print(f" → ERROR: {str(e)[:50]}")
    
    time.sleep(1)

miri_df = pd.DataFrame(all_obs_data)
miri_df.drop_duplicates(subset=['field','obsid'], inplace=True)
miri_df.to_csv(f"{DATA_DIR}/MAST_MIRI_Obs_Positions.csv", index=False)
n_total_miri = len(miri_df)
print(f"\n✅ Total unique MIRI observations: {n_total_miri}")

if n_total_miri == 0:
    sys.exit("No MIRI data found at all!")

# ============================================================
# PHASE 2: Get products for MIRI obs → find F770W ones
# ============================================================
print(f"\n{'='*70}")
print(f"PHASE 2: Querying products for {n_total_miri} MIRI obs...")
print("="*70)

f770w_obsids_per_field = {}  # field -> set of obsids with F770W
BATCH = 20
n_batches = (n_total_miri // BATCH) + 1

all_f770w_products = []

for bi in range(0, n_total_miri, BATCH):
    batch_df = miri_df.iloc[bi:bi+BATCH]
    obsid_list = [str(int(x)) for x in batch_df['obsid'].tolist()]
    
    bnum = bi // BATCH + 1
    
    try:
        prods = Observations.get_product_list(obsid_list)
        
        for row in prods:
            flt = ''
            if 'filters' in prods.colnames:
                fv = str(row.get('filters', '') or '')
                if fv.strip() not in ('', '--', '--', 'nan'):
                    flt = fv.strip()
            
            poid = int(row.get('parent_obsid', 0) or 0)
            
            if flt == 'F770W':
                # Find which field this belongs to
                field_match = miri_df[miri_df['obsid'] == poid]['field']
                if len(field_match) > 0:
                    fld = field_match.values[0]
                    
                    if fld not in f770w_obsids_per_field:
                        f770w_obsids_per_field[fld] = set()
                    f770w_obsids_per_field[fld].add(poid)
                    
                    all_f770w_products.append({
                        'field': fld,
                        'parent_obsid': poid,
                        'filename': str(row.get('productFilename', ''))[:60],
                        'proposal_id': str(row.get('proposal_id', '') or ''),
                    })
        
        nf770w_new = sum(len(v) for v in f770w_obsids_per_field.values())
        if bnum % 5 == 0 or bnum == n_batches:
            print(f"  Batch {bnum}/{n_batches}: {len(prods)} prods, "
                  f"F770W files so far={len(all_f770w_products)}, "
                  f"unique F770W obs={nf770w_new}")
        
    except Exception as e:
        print(f"  Batch {bnum} error: {str(e)[:60]}")
    
    time.sleep(2)

# ============================================================
# PHASE 3: Source matching
# ============================================================
print(f"\n{'='*70}")
print("PHASE 3: Source-to-F770W matching")
print("="*70)

coverage_summary = []
covered_all = []

for field in fields:
    fs = catalog[catalog['field'] == field].copy()
    n_fs = len(fs)
    
    f770w_ids = f770w_obsids_per_field.get(field, set())
    f770w_obs_in_field = miri_df[
        (miri_df['field'] == field) & 
        (miri_df['obsid'].isin(f770w_ids))
    ].copy()
    
    print(f"\n{'─'*55}\n📍 {field.upper()} ({n_fs} sources)")
    print(f"   F770W observations available: {len(f770w_obs_in_field)}")
    
    if len(f770w_obs_in_field) == 0:
        print(f"   ❌ No F770W coverage")
        coverage_summary.append({
            'Field': field, 'N_sources': n_fs, 'N_covered': 0,
            'Pct': 0.0, 'z_range': '-', 'N_f770w_obs': 0
        })
        continue
    
    src_coord = SkyCoord(ra=fs['ra']*u.deg, dec=fs['dec']*u.deg)
    covered_this_field = []
    
    for _, orow in f770w_obs_in_field.iterrows():
        oc = SkyCoord(orow['s_ra'], orow['s_dec'], unit=u.deg)
        sep = src_coord.separation(oc).arcsec
        
        close_idx = np.where(sep < 30)[0]  # 30 arcsec
        
        for ci in close_idx:
            sid = int(fs.iloc[ci]['id'])
            z_val = fs.iloc[ci]['z_phot']
            covered_this_field.append({
                'ID': sid, 'z_phot': z_val, 'sep_arcsec': sep[ci],
                'obs_id': int(orow['obsid']), 'field': field
            })
    
    if covered_this_field:
        cov_df = (pd.DataFrame(covered_this_field)
                   .drop_duplicates(subset='ID')
                   .sort_values('sep_arcsec'))
        
        ncov = len(cov_df)
        z_min = cov_df['z_phot'].min()
        z_max = cov_df['z_phot'].max()
        
        print(f"   ✅ COVERED: {ncov}/{n_fs} ({ncov/n_fs*100:.0f}%), "
              f"z=[{z_min:.2f}-{z_max:.2f}]")
        
        for _, r in cov_df.iterrows():
            print(f"     ID={r['ID']:>6d}, z={r['z_phot']:.2f}, "
                  f"sep={r['sep_arcsec']:>5.1f}\", obs_id={int(r['obs_id'])}")
        
        coverage_summary.append({
            'Field': field, 'N_sources': n_fs, 'N_covered': ncov,
            'Pct': ncov/n_fs*100, 'z_range': f'{z_min:.2f}-{z_max:.2f}',
            'N_f770w_obs': len(f770w_obs_in_field)
        })
        covered_all.extend(cov_df.to_dict('records'))
    else:
        print(f"   ⚠️ F770W exists but no source <30\" match")
        coverage_summary.append({
            'Field': field, 'N_sources': n_fs, 'N_covered': 0,
            'Pct': 0.0, 'z_range': '-', 'N_f770w_obs': len(f770w_obs_in_field)
        })

# ============================================================
# FINAL REPORT
# ============================================================
print(f"\n\n{'='*70}")
print("📊 FINAL REPORT: JWST MIRI F770W COVERAGE OF LRD SAMPLE")
print("="*70)

summary_df = pd.DataFrame(coverage_summary)
print("\n" + summary_df.to_string(index=False))

total_cov = sum(r['N_covered'] for r in coverage_summary)
total_src = sum(r['N_sources'] for r in coverage_summary)

print(f"\n🏆 GRAND TOTAL: {total_cov}/{total_src} LRD sources "
      f"have MIRI F770W data available ({total_cov/total_src*100:.1f}%)")

# Save everything
summary_df.to_csv(f"{DATA_DIR}/MIRI_F770W_Coverage_Summary.csv", index=False)

if covered_all:
    full_cov = pd.DataFrame(covered_all)
    full_cov.to_csv(f"{DATA_DIR}/MIRI_F770W_Source_Coverage.csv", index=False)
    print(f"\n💾 Files saved:")
    print(f"   {DATA_DIR}/MIRI_F770W_Coverage_Summary.csv")
    print(f"   {DATA_DIR}/MIRI_F770W_Source_Coverage.csv")

if all_f770w_products:
    prod_df = pd.DataFrame(all_f770w_products)
    prod_df.to_csv(f"{DATA_DIR}/MIRI_F770W_AllProducts.csv", index=False)
    print(f"   {DATA_DIR}/MIRI_F770W_AllProducts.csv "
          f"({len(prod_df)} F770W files)")
