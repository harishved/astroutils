# GET CLOSEST WISE SOURCE AND ITS PHOTOMETRY
# AUTHOR: HARISH VEDANTHAM
#
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import sys
import astropy.io.fits as pf
import wget
import astropy.wcs as wcs
from PIL import Image
import matplotlib.pyplot as plt
from os import path
from astropy.table import Table
import requests
from io import BytesIO
import aplpy as ap
#
#
#
def get_name(ra,dec):
   rah = int(np.floor(ra/15))
   ram = int(np.floor(60*(ra/15-rah)))
   decd = int(np.floor(dec))
   decm = int(np.floor(60*(dec-decd)))
   return "BDC-J%02d%02d+%02d%02d"%(rah,ram,decd,decm)
#
#
#
ra = float(sys.argv[1])
dec = float(sys.argv[2])
offset = 2.0
length = 10.0
name = get_name(ra,dec)
#
# WISE PART
d = pf.getdata("tiles.fits")
tile_ra = d['ra']
tile_dec = d['dec']
tile_dis = ( (np.cos(np.pi/2*(tile_dec+dec/2) ))**2*(tile_ra-ra)**2 +\
             (tile_dec-dec)**2 )**0.5
tp = np.argsort(tile_dis)
I = tp[1]

filepath = "https://faun.rc.fas.harvard.edu/unwise/release/band-merged/%s.cat.fits"%d[I][0]
filename = "%s.cat.fits"%d[I][0]
if not path.exists("/Users/vedantham/Desktop/op_find_bds/%s"%filename):
   filename = wget.download(filepath)
#
d = pf.getdata(filename)
ra1 = d['ra12']
dec1 = d['dec12']
cd  = np.cos(0.5*np.pi/180*(dec+dec1))
dis =  (cd**2*(ra1-ra)**2 + (dec1-dec)**2 ) **0.5 * 3600
I = np.where(dis<6)[0]
print ("I am here")
print (np.argmin(dis))
for i in I:
   w1 = 22.5-2.5*np.log10(d['fluxlbs'][i][0])
   w2 = 22.5-2.5*np.log10(d['fluxlbs'][i][1])
   w1e = d['fluxlbs'][i][0]/d['dfluxlbs'][i][0]
   w2e = d['fluxlbs'][i][1]/d['dfluxlbs'][i][1]
   print (get_name(ra,dec),dis[i][0],w1,w1e,dis[i][1],w2,w2e)
