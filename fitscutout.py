# MAKE A SMALL FITS CUTOUT IMAGE AROUND A GIVEN POSIUTION
# AUTHOR HARISH VEDANTHAM
# USAGE fitscutout.py ip_fname ra_deg dec_deg size_deg op_fname
#
import numpy as np
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.utils.data import download_file
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u
import sys

ip_f = sys.argv[1]
op_f = sys.argv[5]
ra = float(sys.argv[2])*u.deg
dec = float(sys.argv[3])*u.deg
size = float(sys.argv[4])*u.deg
pos = SkyCoord(ra, dec, frame='icrs')
# MIGHT HAVE TO CHANGE THE LAST ARGUMENT FROM 1 TO 0 DEPENING ON FITS FILE INPUTTED
hdu = fits.open(ip_f)[0]
wcs = WCS(hdu.header)
cutout = Cutout2D(hdu.data, position=pos, size=size, wcs=wcs)
new_data = cutout.data
new_header = cutout.wcs.to_header()

fits.writeto(op_f,new_data,new_header,overwrite=True)

#EOF