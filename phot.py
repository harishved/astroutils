# SIMPLE APERTURE PHOTOMETRY WITH BACKGROUND SUBTRACTION
# USAGE:  phot.py imname ra_deg dec_deg aper_rad_asec
# AUTHOR: HARISH VEDANTHAM
#
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord
from photutils import SkyCircularAperture, SkyCircularAnnulus
from photutils import aperture_photometry
from astropy.io import fits as pf
from astropy.wcs import WCS
import sys


if len(sys.argv)!=5:
   print ("USAGE: phot.py imname ra_deg dec_deg aper_rad_asec")
 
im_name = sys.argv[1]

#
ra = float(sys.argv[2])
dec = float(sys.argv[3])
aper_rad = float(sys.argv[4])

pos = SkyCoord(ra*u.deg,dec*u.deg,frame='icrs')
aperture = SkyCircularAperture(pos, r=aper_rad * u.arcsec)
annulus = SkyCircularAnnulus(pos, r_in=5* u.arcsec, r_out=7*u.arcsec)

data = pf.getdata(im_name)
hdr = pf.getheader(im_name)
wcs = WCS(hdr)

pix_aperture = aperture.to_pixel(wcs)
pix_annulus = annulus.to_pixel(wcs)
annulus_masks = pix_annulus.to_mask(method='center')
annulus_data = annulus_masks.multiply(data)
mask = annulus_masks.data
annulus_data_1d = annulus_data[mask > 0]
bkg = np.nanmedian(annulus_data_1d)
sig = np.nanmedian(np.absolute(annulus_data_1d-bkg)) * 1.4826

print ("Annulus background = %.1e, sigma = %.1e"%(bkg,sig))

rapix,decpix = pix_aperture.positions
print ("RA/DEc pixel = %.1f %.1f"%(rapix,decpix))
phot_table = aperture_photometry(data-bkg, pix_aperture, error=sig*np.ones(data.shape))

# Make an image of the aperture
plt.imshow(pix_aperture.to_mask())
plt.colorbar()
plt.xlabel("X pixel"); plt.ylabel("Y pixel"); plt.title("Aperture mask %.1f arcsec radius"%aper_rad)
plt.tight_layout(); plt.savefig("eraseme.png")
#
# Print the photmetry table
print (phot_table)
flux_adu = phot_table['aperture_sum']
flux_ad_err = phot_table['aperture_sum_err'] 

#mag = -2.5*np.log10(np.absolute(phot_table['aperture_sum']))+25+2.5*np.log10(actual_exposure)
#mag_1sig = -2.5*np.log10(phot_table['aperture_sum_err'])+25+2.5*np.log10(actual_exposure)

print ("FROM HERE ON THINGS ARE TELESCOPE SCPECIFIC. CHECK CODE TO SEE WHAT IS GOIMG ON") 
flux_uJy = 640*1e6*10**(-24.473/2.5) * flux_adu 
flux_uJy_err=  640*1e6*10**(-24.473/2.5) * flux_ad_err

mag = -2.5*np.log10(flux_uJy*1e-6/3631)
mag_max = -2.5*np.log10((flux_uJy+flux_uJy_err)*1e-6/3631)
mag_min = -2.5*np.log10((flux_uJy-flux_uJy_err)*1e-6/3631)
mag_2sig_ul = -2.5*np.log10((flux_uJy_err)*1e-6/3631)

print("flux [uJy] = %.1e, err = %.1e"%(flux_uJy,flux_uJy_err))
print("mag = %.1f; mag 1sig range = %.1f - %.1f"%(mag,mag_min,mag_max))
print("mag 1sign UL = %.1f"%mag_2sig_ul)
