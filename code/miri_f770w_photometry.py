#!/usr/bin/env python3
"""
MIRI F770W Aperture Photometry for LRD Sample
==============================================
Performs aperture photometry on pre-downloaded MIRI F770W Level3 images.
Uses existing files in MIRI_F770W_Images/ directory.
"""

import os
import time
import warnings
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.wcs import WCS
from astropy.stats import sigma_clipped_stats
from photutils.aperture import (
    CircularAperture, CircularAnnulus, aperture_photometry
)

warnings.filterwarnings('ignore', category=UserWarning)
warnings.filterwarnings('ignore', category=RuntimeWarning)  # FITS warnings

# ============================================================
# CONFIGURATION
# ============================================================
WORK_DIR = '/Users/tanxin/Desktop/数据处理'
DOWNLOAD_DIR = os.path.join(WORK_DIR, 'MIRI_F770W_Images')
OUTPUT_FILE = os.path.join(WORK_DIR, 'LRD_MIRI_F770W_Photometry.csv')

# MIRI F770W parameters
MIRI_FWHM_ARCSEC = 0.68   # FWHM at 7.7 microns (diffraction limited)
APERTURE_RADIUS = 0.35     # arcsec (~0.5 * FWHM for optimal S/N)
ANNULUS_INNER = 0.45       # arcsec  
ANNULUS_OUTER = 1.2        # arcsec

# Aperture correction: 0.35" -> infinity for MIRI F770W
# From MIRI instrument handbook Table 10.5 or Labiano et al. 2021
APERTURE_CORRECTION = 1.28


def find_miri_images():
    """Scan download directory and return {obs_id_str: filepath}."""
    print('='*70)
    print('STEP 1: SCANNING FOR MIRI F770W IMAGES')
    print('='*70)
    
    image_map = {}
    
    if not os.path.exists(DOWNLOAD_DIR):
        print(f'ERROR: Download dir not found: {DOWNLOAD_DIR}')
        return image_map
    
    # Find all .fits files recursively
    for root, dirs, files in os.walk(DOWNLOAD_DIR):
        for f in files:
            if not f.endswith('.fits'):
                continue
            
            fpath = os.path.join(root, f)
            
            # Extract obs_id from filename like jw03794-o025_t025_miri_f770w_i2d.fits
            # Format: jw<proposal>-o<obsid>_t<visit>_<detector>_<filter>_type.fits
            base = f.replace('.fits','').split('_')[0]  # e.g., "jw03794-o025"
            if '-' in base:
                prop_obs = base.split('-')[-1]  # e.g., "o025"
                if prop_obs.startswith('o'):
                    obs_num = prop_obs[1:].lstrip('0') or '0'
                    image_map[obs_num] = fpath
    
    print(f'Found {len(image_map)} MIRI F770W FITS images:\n')
    for oid, path in sorted(image_map.items(), key=lambda x: str(x[0])):
        size_mb = os.path.getsize(path) / 1024 / 1024
        print(f'  obs_id={oid:>6s}: {os.path.basename(path)} ({size_mb:.0f} MB)')
    
    return image_map


def do_aperture_photometry(image_data, wcs, ra, dec,
                          aper_radius=APERTURE_RADIUS,
                          ann_in=ANNULUS_INNER, ann_out=ANNULUS_OUTER):
    """
    Perform aperture photometry at given sky coordinates.
    Returns dict with flux, error, background, SNR.
    """
    coord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
    
    # Check if coordinate is within WCS bounds first
    try:
        x_pix, y_pix = wcs.world_to_pixel(coord)
        x_pix = float(np.atleast_1d(x_pix)[0])
        y_pix = float(np.atleast_1d(y_pix)[0])
    except Exception:
        return {'status': 'WCS_ERROR', 'flux': np.nan, 'flux_err': np.nan,
                'background': np.nan, 'snr': np.nan,
                'x_pix': np.nan, 'y_pix': np.nan}
    
    if not (np.isfinite(x_pix) and np.isfinite(y_pix)):
        return {'status': 'OUT_OF_WCS', 'flux': np.nan, 'flux_err': np.nan,
                'background': np.nan, 'snr': np.nan,
                'x_pix': x_pix, 'y_pix': y_pix}
    
    ny, nx = image_data.shape
    margin = ann_out + 5
    
    if x_pix < margin or x_pix > nx - margin or y_pix < margin or y_pix > ny - margin:
        return {'status': 'EDGE', 'flux': np.nan, 'flux_err': np.nan,
                'background': np.nan, 'snr': np.nan, 'x_pix': x_pix, 'y_pix': y_pix}
    
    # Handle NaN/Inf in data
    clean_data = np.where(np.isfinite(image_data), image_data, 0.0)
    
    # Apertures
    aperture = CircularAperture((x_pix, y_pix), r=aper_radius)
    annulus = CircularAnnulus((x_pix, y_pix), r_in=ann_in, r_out=ann_out)
    
    # Background from annulus
    ann_mask = annulus.to_mask(method='center')
    bkg_values = []
    ann_data = ann_mask.multiply(clean_data)
    ann_1d = ann_data[ann_mask.data > 0]
    if len(ann_1d) > 10:
        _, median_bkg, _ = sigma_clipped_stats(ann_1d, sigma=3.0)
        bkg_values.append(median_bkg)
    
    if len(bkg_values) == 0:
        return {'status': 'NO_BKG', 'flux': np.nan, 'flux_err': np.nan,
                'background': np.nan, 'snr': np.nan, 'x_pix': x_pix, 'y_pix': y_pix}
    
    background = float(np.median(bkg_values))
    
    # Photometry with Poisson errors
    sub_data = np.maximum(clean_data - background, 0.0)
    error_data = np.sqrt(np.maximum(clean_data, 0.0))
    
    try:
        phot_table = aperture_photometry(sub_data, aperture, error=error_data)
        raw_flux = phot_table['aperture_sum'][0]
        raw_err = phot_table['aperture_sum_err'][0]
    except Exception as e:
        return {'status': 'PHOT_ERR', 'flux': np.nan, 'flux_err': np.nan,
                'background': background, 'snr': np.nan, 'x_pix': x_pix, 'y_pix': y_pix}
    
    # Apply aperture correction
    flux_corr = raw_flux * APERTURE_CORRECTION
    flux_err_corr = raw_err * APERTURE_CORRECTION
    
    snr = flux_corr / flux_err_corr if flux_err_corr > 0 else np.nan
    status = 'DET' if snr >= 3.0 else ('LOW_SNR' if snr >= 2.0 else 'NONDET')
    
    return {
        'status': status, 'flux': float(flux_corr), 'flux_err': float(flux_err_corr),
        'background': float(background), 'snr': float(snr),
        'x_pix': float(x_pix), 'y_pix': float(y_pix)
    }


def get_image_info(fpath):
    """Get basic info about a MIRI image."""
    with fits.open(fpath) as hdul:
        # Find science extension
        sci_idx = None
        for i, hdu in enumerate(hdul):
            if hdu.data is not None and hdu.data.ndim == 2:
                hdr = hdu.header
                extname = str(hdr.get('EXTNAME', ''))
                if extname == 'SCI' or 'BUNIT' in hdr or 'MJY' in str(hdr.get('BUNIT', '')):
                    sci_idx = i
                    break
        
        if sci_idx is None:
            # First 2D extension
            for i, hdu in enumerate(hdul):
                if hdu.data is not None and hdu.data.ndim == 2:
                    sci_idx = i
                    break
        
        if sci_idx is None:
            return None, None, None
        
        header = hdul[sci_idx].header
        data = hdul[sci_idx].data
        wcs = WCS(header)
        
        # Get BUNIT
        bunit = header.get('BUNIT', 'UNKNOWN')
        
        # Get pixel scale
        cdelt1 = abs(header.get('CD1_1', header.get('CDELT1', 0.11)))
        pixel_scale_arcsec = cdelt1 * 3600  # deg -> arcsec
        
        # Get error extension
        err_data = None
        for i, hdu in enumerate(hdul):
            ename = str(hdu.header.get('EXTNAME', ''))
            if ename == 'ERR':
                err_data = hdu.data
                break
        
        return data, wcs, {'bunit': bunit, 'pixel_scale': pixel_scale_arcsec, 
                           'err': err_data}


def find_image_for_obs(obs_id_str, image_map):
    """Find matching image for an observation ID (try multiple formats)."""
    # Try direct match
    if obs_id_str in image_map:
        return image_map[obs_id_str]
    
    # Try partial match: obs_id may contain the short obs number
    for oid_key, fpath in image_map.items():
        if oid_key in obs_id_str or obs_id_str in oid_key:
            return fpath
    
    # Try matching by filename containing obs_id substring
    # e.g., coverage has 213002751, filename is jw03794-o025_... -> look for "025" in both
    for oid_key, fpath in image_map.items():
        fname = os.path.basename(fpath)
        if oid_key in fname:
            return fpath
        # Extract numeric parts of filename and check
        import re
        nums_in_fname = re.findall(r'\d+', fname)
        for num in nums_in_fname:
            if num in obs_id_str or obs_id_str in num:
                return fpath
    
    return None


def run_photometry(coverage_df, master_df, image_map):
    """Run photometry pipeline on all covered sources."""
    print('\n' + '='*70)
    print('STEP 2: APERTURE PHOTOMETRY ON MIRI F770W IMAGES')
    print('='*70)
    
    results = []
    obs_col = 'closest_obsid' if 'closest_obsid' in coverage_df.columns else \
              ('obs_id' if 'obs_id' in coverage_df.columns else None)
    
    # Merge: coverage + master coordinates
    merged = coverage_df.merge(master_df[['id','ra','dec','z_phot']], on='id', how='left', 
                                suffixes=('', '_master'))
    
    # Use correct column name for z_phot (may have suffix after merge)
    z_col = 'z_phot' if 'z_phot' in merged.columns else ('z_phot_master' if 'z_phot_master' in merged.columns else None)
    
    n_total = len(merged)
    
    # Cache loaded images to avoid re-reading same file
    img_cache = {}
    
    for idx, row in merged.iterrows():
        src_id = row['id']
        z = row[z_col] if z_col else np.nan
        ra = row.get('ra', row.get('ra_master', np.nan))
        dec = row.get('dec', row.get('dec_master', np.nan))
        
        sep_col_name = 'sep_arcsec' if 'sep_arcsec' in merged.columns else \
                       ('separation_arcsec' if 'separation_arcsec' in merged.columns else None)
        sep_val = row[sep_col_name] if sep_col_name else 99
        
        obs_id_str = str(int(row[obs_col])) if obs_col else '?'
        
        # Find matching image
        img_path = find_image_for_obs(obs_id_str, image_map)
        
        if img_path is None:
            result = {'id': src_id, 'obs_id': obs_id_str, 'status': 'NO_IMAGE',
                     'flux_njy': np.nan, 'flux_err_njy': np.nan,
                     'snr': np.nan, 'f_abmag': np.nan,
                     'x_pix': np.nan, 'y_pix': np.nan, 'z_phot': z, 'sep_arcsec': sep_val}
            results.append(result)
            print(f'[{idx+1:2d}/{n_total}] ID={str(src_id):>6s} z={z:.2f}: NO_IMAGE (obs={obs_id_str})')
            continue
        
        
        # Load image (with caching)
        if img_path not in img_cache:
            data, wcs, info = get_image_info(img_path)
            if data is None:
                result = {'id': src_id, 'obs_id': obs_id_str, 'status': 'READ_ERR',
                         'flux_njy': np.nan, 'flux_err_njy': np.nan,
                         'snr': np.nan, 'f_abmag': np.nan,
                         'x_pix': np.nan, 'y_pix': np.nan, 'z_phot': z, 'sep_arcsec': sep_val}
                results.append(result)
                print(f'[{idx+1:2d}/{n_total}] ID={str(src_id):>6s} z={z:.2f}: READ_ERROR')
                continue
            img_cache[img_path] = (data, wcs, info)
        else:
            data, wcs, info = img_cache[img_path]
        
        # Do photometry
        photo = do_aperture_photometry(data, wcs, ra, dec)
        
        # Convert to physical units based on BUNIT
        flux_raw = photo['flux'] if not np.isnan(photo['flux']) else np.nan
        flux_err_raw = photo['flux_err'] if not np.isnan(photo['flux_err']) else np.nan
        
        bunit = info['bunit'].upper() if info else ''
        
        # Convert to nJy
        flux_njy = np.nan
        flux_err_njy = np.nan
        
        if not np.isnan(flux_raw) and bunit:
            if 'MJY/SR' in bunit or 'MJY STERADIAN' in bunit:
                # Surface brightness: need to multiply by pixel area in steradians
                ps = info['pixel_scale'] / 3600  # arcsec -> deg
                pix_area_sr = (ps * np.pi/180)**2  # deg^2 -> sr
                flux_njy = flux_raw * 1e9 * pix_area_sr  # MJy/sr/pix * sr -> Jy -> nJy
                flux_err_njy = flux_err_raw * 1e9 * pix_area_sr
            elif 'MJY' in bunit:
                flux_njy = flux_raw * 1e9  # MJy -> nJy
                flux_err_njy = flux_err_raw * 1e9
            elif 'JY' in bunit or 'JANSKY' in bunit:
                flux_njy = flux_raw * 1e9  # Jy -> nJy
                flux_err_njy = flux_err_raw * 1e9
            elif 'DN' in bunit or 'ELECTRON' in bunit:
                # DN/s or count rate — need PHOTMJSR conversion
                flux_njy = flux_raw  # placeholder, needs calibration
                flux_err_njy = flux_err_raw
            else:
                flux_njy = flux_raw
                flux_err_njy = flux_err_raw
        
        # AB magnitude
        ab_mag = (-2.5 * np.log10(flux_njy) + 31.4) if (not np.isnan(flux_njy) and flux_njy > 0) else np.nan
        ab_mag_err = (2.5/np.log(10) * flux_err_njy/flux_njy) if (not np.isnan(flux_err_njy) and flux_err_njy > 0 and flux_njy > 0) else np.nan
        
        status_icon = {'DET':'✓', 'LOW_SNR':'~', 'NONDET':'○', 'EDGE':'△', 'NO_BKG':'✗', 'PHOT_ERR':'?', 'READ_ERR':'?'}
        icon = status_icon.get(photo['status'], '?')
        
        result = {
            'id': src_id, 'obs_id': obs_id_str, 'status': photo['status'],
            'flux_njy': round(float(flux_njy), 3) if not np.isnan(flux_njy) else np.nan,
            'flux_err_njy': round(float(flux_err_njy), 3) if not np.isnan(flux_err_njy) else np.nan,
            'snr': round(float(photo['snr']), 2) if not np.isnan(photo['snr']) else np.nan,
            'f_abmag': round(float(ab_mag), 2) if not np.isnan(ab_mag) else np.nan,
            'f_abmag_err': round(float(ab_mag_err), 2) if not np.isnan(ab_mag_err) else np.nan,
            'x_pix': round(float(np.atleast_1d(photo['x_pix'])[0]), 2),
            'y_pix': round(float(np.atleast_1d(photo['y_pix'])[0]), 2),
            'z_phot': z, 'sep_arcsec': sep_val,
            'bunit': bunit
        }
        results.append(result)
        
        print(f'  [{idx+1:2d}/{n_total}] ID={str(src_id):>6s} z={z:.2f} sep={sep_val:5.1f}" | '
              f'{icon} {photo["status"]:>8s} | SNR={photo["snr"]:+.1f} | '
              f'f={flux_njy:+.1f} nJy | m_AB={ab_mag:.1f}')
    
    return pd.DataFrame(results)


def summarize(df_photo):
    """Print summary statistics."""
    print('\n' + '='*70)
    print('MIRI F770W PHOTOMETRY SUMMARY')
    print('='*70)
    
    n = len(df_photo)
    
    for status in ['DET', 'LOW_SNR', 'NONDET', 'EDGE', 'NO_BKG', 'PHOT_ERR', 'READ_ERR', 'NO_IMAGE']:
        cnt = len(df_photo[df_photo['status']==status])
        if cnt > 0:
            pct = cnt/n*100
            print(f'  {status:>8s}: {cnt:3d} ({pct:.0f}%)')
    
    det = df_photo[df_photo['status']=='DET']
    if len(det) > 0:
        print(f'\n--- Detected sources ({len(det)}) ---')
        print(f'  Flux range:     {det["flux_njy"].min():.1f} – {det["flux_njy"].max():.1f} nJy')
        print(f'  AB mag range:   {det["f_abmag"].min():.1f} – {det["f_abmag"].max():.1f}')
        print(f'  Median SNR:     {det["snr"].median():.1f}')
        print(f'  Median sep:     {det["sep_arcsec"].median():.1f}"')
    
    nondet = df_photo[df_photo['status'].isin(['NONDET','LOW_SNR'])]
    if len(nondet) > 0:
        ulimits = (3 * nondet['flux_err_njy']).dropna()
        if len(ulimits) > 0:
            print(f'\n--- Non-detections ({len(nondet)}) ---')
            print(f'  Median 3σ upper limit: {ulimits.median():.1f} nJy')


if __name__ == '__main__':
    t0 = time.time()
    
    # Load coverage list
    cov_path = os.path.join(WORK_DIR, 'MIRI_F770W_Source_Coverage.csv')
    coverage = pd.read_csv(cov_path)
    print(f'Loaded {len(coverage)} covered sources\n')
    
    # Load master catalog
    master = pd.read_csv(os.path.join(WORK_DIR, 'Kokorev_LRDs_Full.csv'))
    print(f'Loaded master: {len(master)} sources\n')
    
    # Step 1: Find images
    image_map = find_miri_images()
    
    # Step 2: Photometry
    if len(image_map) > 0:
        photo_df = run_photometry(coverage, master, image_map)
        
        summarize(photo_df)
        
        # Save results
        photo_df.to_csv(OUTPUT_FILE, index=False)
        print(f'\n✅ Results saved: {OUTPUT_FILE}')
        
        # Quick reference table
        cols = ['id','z_phot','status','snr','f_abmag','flux_njy','sep_arcsec','obs_id']
        print(f'\n--- Full Results ---\n{photo_df[cols].to_string(index=False)}')
    else:
        print('\n⚠ No MIRI images found! Cannot run photometry.')
        print(f'Please ensure F770W .fits files exist in: {DOWNLOAD_DIR}')
    
    elapsed = time.time() - t0
    print(f'\nTotal time: {elapsed:.0f}s ({elapsed/60:.1f} min)')
