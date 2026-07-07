***

Catalog of the photometric sample of BH* dominated LRDs presented in Weibel et al. 2026. 
Please cite this paper when using the catalog.

Columns are described below, <fil> stands for the HST or NIRCam filter name in lowercase letters (e.g., f444w):

id                          n/a                 identification number from the original SourceExtractor run
ra                          deg                 right ascension (J2000)
dec                         deg                 declination (J2000)
field			    n/a			field from which the source was selected
x_image                     pixel               centroid x-position of the object (SourceExtractor)
y_image                     pixel               centroid y-position of the object (SourceExtractor)
a_image                     pixel               isophotal image major axis (SourceExtractor)
b_image                     pixel               isophotal image minor axis (SourceExtractor)
theta_image                 deg                 isophotal image position angle (SourceExtractor)
kron_radius                 n/a                 Kron radius in units of A or B (SourceExtractor)
sn_lw_det_img	    	    n/a			S/N in the LW detection image (inverse variance weighted 
						F277W+F356W+F444W stack)
aper_corr                   n/a                 aperture correction factor derived in F444W based on the Kron
						ellipse derived in the detection image
enclosed_flux_fraction      n/a			fraction of energy encircled by the Kron-ellipse in the F444W PSF. 
						The total correction applied to the aperture fluxes is 
						= use aper_corr / enclosed_flux_fraction

faper_<fil>		    nJy			aperture flux in a 0.16" radius circular aperture (PSF-matched)
eaper_<fil>		    nJy		        uncertainty in the aperture flux
f_<fil>		            nJy		        total flux (aperture-corrected and PSF-matched)
e_<fil>			    nJy			error in the total flux, with a 5% error floor
sn_<fil>		    n/a			signal-to-noise ratio before applying the 5% error floor
		
c_f444w			    n/a			compactness parameter f(D=0.4") / f(D=0.2") in F444W, 
						used to pre-select point sources
bhs_template		    n/a                 BH*-template producing the lowest chi^2 fit with the corresponding 
						sfhz+BH* template set
bhs_template_contribution   n/a	                Integrated contribution of the BH*-template to the best-fitting SED 
						at rest-frame wavelengths between 0.4 and 1 microns, or the longest 
						rest-wavelength covered by NIRCam

is_vshaped	            binary	  	whether the photometry satisfies the V-shaped color-selection 
						following Kokorev+24
						(https://ui.adsabs.harvard.edu/abs/2024ApJ...968...38K/abstract)

z_phot_eazy                 n/a                 best-fitting photometric redshift from eazy (peak of the PDF)
z_phot_160_eazy             n/a                 16th percentile of the redshift PDF
z_phot_840_eazy             n/a                 84th percentile of the redshift PDF
z_phot_chi2_eazy            n/a                 chi^2 of the fit at z_phot_eazy

delta_chi2_eazy_stars       n/a                 chi^2 difference between the best fit with the sfhz+BH* template set, 
						and the best fit with a large grid of stellar templates
bayes_factor_eazy_stars	    n/a			Bayes factor based on delta_chi2_eazy_stars, capped at 9e279
delta_chi2_sfhz		    n/a			chi^2 difference between the best fit with the sfhz+BH* template set, 
						and the best fit with the sfhz template set only, == -99 if eazy did 
						not find a fit with the sfhz template set
bayes_factor_sfhz	    n/a			Bayes factor based on delta_chi2_sfhz, set to 9e279 if 
						delta_chi2_sfhz == -99

z_phot_160_zls10_eazy	    n/a			16th percentile of the redshift PDF, ignoring PDF at z > 10
z_phot_500_zls10_eazy	    n/a			50th percentile of the redshift PDF, ignoring PDF at z > 10
z_phot_840_zls10_eazy	    n/a			84th percentile of the redshift PDF, ignoring PDF at z > 10

nirspec_gratings	    n/a			list of available NIRSpec gratings on the DJA
z_spec_prism		    n/a			spectroscopic redshift from PRISM spectrum (DJA)
z_spec_g395m		    n/a			spectroscopic redshift from G395M spectrum (DJA)
z_spec_g235m		    n/a			spectroscopic redshift from G235M spectrum (DJA)
spec_grade_dja		    n/a			grade for spectroscopic redshift from DJA 
						(best grade in case of multiple spectra)
z_spec_grism		    n/a			spectroscopic redshift from NIRCam/grism 
						(only if no NIRSpec redshift is available)
z_spec_other		    n/a			spectroscopic redshift from other instruments 
						(only if no JWST redshift is available)
spec_program 		    n/a			list of programs that obtained spectroscopy for this source 
						(usually JWST)
z_spec			    n/a			spectroscopic redshift used in paper, prioritizing 
						1) prism, 2) grating, 3) grism, 4) other

mu_lensing		    n/a			lensing magnification from UNCOVER DR4 
						(https://jwst-uncover.github.io/DR4.html#LensingMaps)
L5100			    erg/s		monochromatic luminosity at 5100 Angstrom
L5100_err		    erg/s		uncertainty in monochromatic luminosity at 5100 Angstrom, 
						derived from uncertainty in the rest-frame V-band (eazy)
balmer_break	            n/a			Balmer break strength measured from the best-fitting SED in f_nu 
						as the ratio between median fluxes in the bands 
						[4060,4140] Angstrom and [3560,3640] Angstrom
balmer_break_err	    n/a			uncertainty in the Balmer break strength, derived from 
						uncertainties in the rest-frame B- and U-bands (eazy)
Lbol       		    erg/s		bolometric luminosity derived from integrating over the BH*-template 
						contribution to the best-fitting SED
Lbol_err		    erg/s		uncertainty in bolometric luminosity, derived from 
						uncertainty in the rest-frame V-band (eazy)

***