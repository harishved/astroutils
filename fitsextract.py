# Program to extract a subimage of fits file
# The given position is at the cntre pixel of resultant file
# Assumes that fits has 4 axes (typical for radio astronomy images)
# I am not sure that this is the correct way to do this in terms of WCS away from the extract centre
# It is mainly to be used as a quick tool to make thumbnails of very small sky areas
# AUTHOR: HARISH VEDANTHAM
#
from astropy.io.fits import getdata
from astropy.io.fits import getheader
from astropy.io.fits import writeto
import sys
from astropy.wcs import WCS 
from  astropy.wcs import WCS 
import numpy as np

if len(sys.argv)!=6:
    print("USAGE: fitsextract.py in_fname ra_deg dec_deg imsize_deg op_fname")
    exit(1)
header=getheader(sys.argv[1])
print (header['OBJECT'])
ra = float(sys.argv[2])
dec = float(sys.argv[3])
imsize = float(sys.argv[4])

h_wcs=WCS(header=header)
coords=np.array([[ra,dec,0,0],])
tp = h_wcs.all_world2pix(coords,0)

ra_pix=tp[0][0]
dec_pix=tp[0][1]

ra_pix_round=int(np.floor(tp[0][0]))
dec_pix_round = int(np.floor(tp[0][1]))

print(ra_pix_round,dec_pix_round)
ra_cnt = header['CRPIX1']
dec_cnt = header['CRPIX2']

dra=header['CDELT2']
npix_out = int(2*np.floor(imsize/dra/2)+1)

print (dec_pix_round-(npix_out-1)/2,dec_pix_round+(npix_out-1)/2+1,ra_pix_round-(npix_out-1)/2,ra_pix_round + (npix_out-1)/2+1)
imsub = getdata(sys.argv[1])[:,:,int(dec_pix_round-(npix_out-1)/2):int(dec_pix_round+(npix_out-1)/2+1),int(ra_pix_round-(npix_out-1)/2):int(ra_pix_round + (npix_out-1)/2+1)]

header['CRPIX1']=(npix_out+1)/2 + ra_pix-ra_pix_round
header['CRPIX2']=(npix_out+1)/2 + dec_pix-dec_pix_round

header['CRVAL1']=ra
header['CRVAL2']=dec

header['NAXIS1']=npix_out
header['NACXIS2']=npix_out

writeto(sys.argv[5],data=imsub,header=header,clobber=True)

