# GET PANSTARRS CUT OUT IMAGES IN DIFF FILTERS 
# USAGE python3 get_panstarrs.py <ra> <dec>
#
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import aplpy as ap
import matplotlib.pyplot as plt
import sys
import wget
import astropy.io.fits as pf
from os import path
import astropy.wcs as wcs
#
def get_name(ra,dec):
   rah = int(np.floor(ra/15))
   ram = int(np.floor(60*(ra/15-rah)))
   decd = int(np.floor(dec))
   decm = int(np.floor(60*(dec-decd)))
   return "BDC-J%02d%02d+%02d%02d"%(rah,ram,decd,decm)
#
#
ra   = float(sys.argv[1])
dec  = float(sys.argv[2])
offset = 2.0
length = 10.0
cd = np.cos(dec*np.pi/10.0)
name = get_name(ra,dec)
#

fname = wget.download("http://ps1images.stsci.edu/cgi-bin/ps1filenames.py?ra=%f&dec=%f"%(ra,dec))
f = open(fname,"r")
lines = f.readlines()
nlines = len(lines)
fil = []
img = []
for i in range(1,nlines):
   w = lines[i].strip("\n").split(" ")
   fil.append(w[4])
   img.append(w[7])

down_img_name = []

for imgpath in img:
   f = wget.download("http://ps1images.stsci.edu/cgi-bin/fitscut.cgi?ra=%f&dec=%f&red=%s&wcs=true&format=fits&size=256"%(ra,dec,imgpath))
   down_img_name.append(f)

'''
fil = ['g','i','r','y','z']
down_img_name = ['cutout_rings.v3.skycell.1980.090.stk.g.unconv.fits_sci.fits',\
'cutout_rings.v3.skycell.1980.090.stk.i.unconv.fits_sci.fits',\
'cutout_rings.v3.skycell.1980.090.stk.r.unconv.fits_sci.fits',\
'cutout_rings.v3.skycell.1980.090.stk.y.unconv.fits_sci.fits',\
'cutout_rings.v3.skycell.1980.090.stk.z.unconv.fits_sci.fits']
'''

for i in range(len(fil)):
   filename = down_img_name[i]
   #   
   d = pf.getdata(filename)
   med = np.nanmedian(d)
   mad = np.nanmedian(np.absolute(d-med))
   gc = ap.FITSFigure(filename,figsize=(4,4))
   gc.recenter(ra,dec,0.5/60)
   cd = np.cos(dec*np.pi/180)
   gc.show_lines([np.array([[ra+offset/3600/cd,ra+(offset + length)/3600/cd],[dec,dec]])],color="blue")
   gc.show_lines([np.array([[ra,ra],[dec+offset/3600,dec+(offset + length)/3600]])],color="blue")
   gc.show_colorscale(vmin=med-15*mad,vmax=med+15*mad,interpolation='none',cmap=mpl.cm.RdBu_r)
   gc.add_label(0.2,0.95,"PanSTARRS %s"%fil[i],relative=True,color="black",fontsize=12)
   gc.tick_labels.hide()
   gc.axis_labels.hide()
   gc.save("%s_%s.png"%(name,fil[i]),dpi=150)
   gc.close()
#
