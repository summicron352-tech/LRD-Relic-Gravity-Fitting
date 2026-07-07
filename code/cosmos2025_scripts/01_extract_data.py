#!/usr/bin/env python3
"""Extract key columns from COSMOSWeb_mastercatalog_v1.1.fits and compute Σ, color."""
from astropy.io import fits
from astropy.cosmology import Planck18
import astropy.units as u
import numpy as np
import time

t0 = time.time()
hdul = fits.open('../COSMOSWeb_mastercatalog_v1.1.fits', memmap=True)

phot = hdul[1].data; lephare = hdul[2].data

zfinal = lephare['zfinal']; mass_med = lephare['mass_med']
obj_type = lephare['Type'] if 'Type' in lephare.columns.names else np.zeros(len(zfinal))
radius_deg = phot['radius_sersic']; mag_f444w = phot['mag_auto_f444w']
mag_f150w = phot['mag_auto_f150w']; warn_flag = phot['warn_flag']

mask = ((obj_type==0) & (zfinal>0) & (zfinal<99) & np.isfinite(zfinal) &
        (mass_med>-90) & np.isfinite(mass_med) & (radius_deg>0) & np.isfinite(radius_deg) &
        (warn_flag==0) & np.isfinite(mag_f444w) & np.isfinite(mag_f150w))

z = zfinal[mask]; logM = mass_med[mask]; r_deg = radius_deg[mask]
m444 = mag_f444w[mask]; m150 = mag_f150w[mask]

r_rad = np.deg2rad(r_deg); DA = Planck18.angular_diameter_distance(z).to(u.pc).value
r_pc = r_rad * DA; Sigma = 10**logM / (np.pi * r_pc**2)
logSigma = np.log10(Sigma); color = 0.4*(m150-m444)

valid = np.isfinite(logSigma) & np.isfinite(color)
z=z[valid]; logM=logM[valid]; logSigma=logSigma[valid]; color=color[valid]; r_pc=r_pc[valid]

np.savez('../cosmos2025_extracted.npz', z=z, logM=logM, logSigma=logSigma,
         color=color, reff_pc=r_pc, n=len(z))
print(f"Extracted {len(z)} galaxies in {time.time()-t0:.1f}s")
hdul.close()
