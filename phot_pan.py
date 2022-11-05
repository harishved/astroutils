# SIMPLE PHOTOMETRY OF PANSTARRS IMAGES
# Make sure you download not just image.fits but also weight and exp images
# USAGE: pan_phot.py imname ra_deg dec_deg aper_rad_asec
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

cov_length = 0.75 # Width of the pixel error covariance matrix in arcsec (corresponding to 1arcsec)

if len(sys.argv)!=5:
   print ("USAGE: pan_phot.py imname ra_deg dec_deg aper_rad_asec")
 
prefix = sys.argv[1]
im_name = prefix+".fits"
exp_name = prefix+".exp.fits"
wt_name = prefix+".wt.fits"
#
ra = float(sys.argv[2])
dec = float(sys.argv[3])
aper_rad = float(sys.argv[4])

Nind = (aper_rad/cov_length)**2
if Nind<1:
   Nind = 1
print ("Num of independent noise pixels = %.1f"%Nind)

pos = SkyCoord(ra*u.deg,dec*u.deg,frame='icrs')
aperture = SkyCircularAperture(pos, r=aper_rad * u.arcsec)
annulus = SkyCircularAnnulus(pos, r_in=4* u.arcsec, r_out=7*u.arcsec)

data = pf.getdata(im_name)
wt = pf.getdata(wt_name)
act_exp = pf.getdata(exp_name)

hdr = pf.getheader(im_name)
wcs = WCS(hdr)
exptime = hdr['EXPTIME'] # Avg EXP time
print ("Total exp time = %.1f s"%exptime)

#pix = wcs.all_world2pix([ra],[dec],[0])
#print (pix)
#exit(1)


pix_aperture = aperture.to_pixel(wcs)
pix_annulus = annulus.to_pixel(wcs)
annulus_masks = pix_annulus.to_mask(method='center')
annulus_data = annulus_masks.multiply(data)
mask = annulus_masks.data
annulus_data_1d = annulus_data[mask > 0]
bkg = np.nanmedian(annulus_data_1d)
print ("Annulus background = %.1e"%bkg)

rapix,decpix = pix_aperture.positions
print ("RA/DEc pixel = %.1f %.1f"%(rapix,decpix))
actual_exposure = act_exp[int(decpix),int(rapix)]
print ("Actual exposure at target location = %.1f"%actual_exposure)
phot_table = aperture_photometry(data-bkg, pix_aperture, error=wt**0.5)

samp_table = aperture_photometry(np.ones(data.shape), pix_aperture, error=wt**0.5)
Npix = samp_table['aperture_sum']

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

print("TELESCOPE SPECIFIC HARDCODED VALUES FROM HERE ON. CHECK CODE BEFORE BELIEVING ANYTHING")
flux_uJy = 0.3631 * flux_adu/actual_exposure
flux_uJy_err= 0.3631 * flux_ad_err/actual_exposure * (Npix/Nind)**0.5

mag = -2.5*np.log10(flux_uJy*1e-6/3631)
mag_max = -2.5*np.log10((flux_uJy+flux_uJy_err)*1e-6/3631)
mag_min = -2.5*np.log10((flux_uJy-flux_uJy_err)*1e-6/3631)
mag_2sig_ul = -2.5*np.log10((1*flux_uJy_err)*1e-6/3631)

print("flux [uJy] = %.1e, err = %.1e"%(flux_uJy,flux_uJy_err))
print("mag = %.1f; mag 1sig range = %.1f - %.1f"%(mag,mag_min,mag_max))
print("mag 1sign UL = %.1f"%mag_2sig_ul)
