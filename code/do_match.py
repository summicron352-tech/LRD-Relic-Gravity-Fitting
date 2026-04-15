import pandas as pd, numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u

cat = pd.read_csv('Kokorev_LRDs_Full.csv')
obs = pd.read_csv('MAST_MIRI_Obs_Positions.csv')
print('LRD: %d, MIRI obs: %d' % (len(cat), len(obs)))

all_results = []

for field in cat['field'].unique():
    fs = cat[cat['field']==field]
    fobs = obs[obs['field']==field]
    n_fs = len(fs)
    
    if len(fobs) == 0:
        print('\n%s (%d): NO MIRI OBSERVATIONS' % (field, n_fs))
        continue
    
    src = SkyCoord(ra=fs['ra'].values*u.deg, dec=fs['dec'].values*u.deg)
    oc = SkyCoord(ra=fobs['s_ra'].values*u.deg, dec=fobs['s_dec'].values*u.deg)
    
    matches = []
    for si in range(len(fs)):
        seps = src[si].separation(oc).arcsec
        idx_min = int(np.argmin(seps))
        matches.append({
            'id': int(fs.iloc[si]['id']),
            'z_phot': float(fs.iloc[si]['z_phot']),
            'sep_arcsec': float(seps[idx_min]),
            'closest_obsid': int(fobs.iloc[idx_min]['obsid']),
            'field': str(field),
        })
    
    bm = pd.DataFrame(matches)
    
    print('\n' + '='*60)
    print('%s (%d sources, %d MIRI obs)' % (field.upper(), n_fs, len(fobs)))
    print('='*60)
    
    for rad in [10, 20, 30]:
        n_close = int((bm['sep_arcsec'] < rad).sum())
        
        if rad == 30 and n_close > 0:
            close = bm[bm['sep_arcsec'] < rad].sort_values('sep_arcsec')
            print('  F770W COVERAGE (<%d\"): %d/%d (%.0f%%)' % (rad, n_close, n_fs, 100.*n_close/n_fs))
            
            z_vals = close['z_phot'].values
            print('  z_range: [%.2f, %.2f]' % (z_vals.min(), z_vals.max()))
            
            for _, r in close.iterrows():
                print('    ID=%d, z=%.2f, sep=%.1f", obs_id=%d' % (
                    r['id'], r['z_phot'], r['sep_arcsec'], r['closest_obsid']))
            
            all_results.extend(close.to_dict('records'))
        else:
            print('  Within %d": %d' % (rad, n_close))
    
    if (bm['sep_arcsec'] < 30).sum() == 0:
        all_results.append({'field': field, 'n_sources': n_fs, 'n_covered': 0})

# Grand total
print('\n' + '='*60)
print('GRAND TOTAL')

if all_results:
    has_id = [r for r in all_results if 'id' in r]
    unique_ids = set(r['id'] for r in has_id) if has_id else set()
    print('Total LRD with MIRI data within 30": %d / %d' % (len(unique_ids), len(cat)))
    
    # Save
    if has_id:
        cov_df = pd.DataFrame(has_id).drop_duplicates(subset='id').sort_values(['field','id'])
        cov_df.to_csv('MIRI_F770W_Source_Coverage.csv', index=False)
        print('\nSaved: MIRI_F770W_Source_Coverage.csv (%d sources covered)' % len(cov_df))
else:
    print('No coverage found!')
