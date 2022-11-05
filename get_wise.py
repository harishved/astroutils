# GET UNWISE IMAGE FOR A GIVEN RA DEC LOCATION ON SKY
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
I = tp[0]

path_w1 = "http://unwise.me/data/neo6/unwise-coadds/fulldepth/"+\
        str(d[I][0][:3])+"/"+str(d[I][0])+"/"+\
        "unwise-%s-w1-img-m.fits"%d[I][0]
path_w2 = "http://unwise.me/data/neo6/unwise-coadds/fulldepth/"+\
        str(d[I][0][:3])+"/"+str(d[I][0])+"/"+\
        "unwise-%s-w2-img-m.fits"%d[I][0]

filename_w1 = "unwise-%s-w1-img-m.fits"%d[I][0]
filename_w2 = "unwise-%s-w2-img-m.fits"%d[I][0]


print (filename_w1, filename_w2)


if not path.exists("/Users/vedantham/Desktop/op_find_bds/%s"%filename_w1):
    filename_w1 = wget.download(path_w1)
if not path.exists("/Users/vedantham/Desktop/op_find_bds/%s"%filename_w2):
   filename_w2 = wget.download(path_w2)

#
d = pf.getdata(filename_w1)
med = np.median(d)
mad = np.median(np.absolute(d-med))
gc = ap.FITSFigure(filename_w1,figsize=(4,4))
gc.recenter(ra,dec,0.5/60)
cd = np.cos(dec*np.pi/180)
gc.show_lines([np.array([[ra+offset/3600/cd,ra+(offset + length)/3600/cd],[dec,dec]])],color="blue")
gc.show_lines([np.array([[ra,ra],[dec+offset/3600,dec+(offset + length)/3600]])],color="blue")
gc.show_colorscale(vmin=med-15*mad,vmax=med+15*mad,interpolation='none',cmap=mpl.cm.RdBu_r)
gc.add_label(0.15,0.92,"WISE W1",relative=True,color="black",fontsize=12)
gc.tick_labels.hide()
gc.axis_labels.hide()
gc.save("%s_W1.png"%name,dpi=150)
gc.close()
#
#   
d = pf.getdata(filename_w2)
med = np.median(d)
mad = np.median(np.absolute(d-med))
gc = ap.FITSFigure(filename_w2,figsize=(4,4))
gc.recenter(ra,dec,0.5/60)
cd = np.cos(dec*np.pi/180)
gc.show_lines([np.array([[ra+offset/3600/cd,ra+(offset + length)/3600/cd],[dec,dec]])],color="blue")
gc.show_lines([np.array([[ra,ra],[dec+offset/3600,dec+(offset + length)/3600]])],color="blue")
gc.show_colorscale(vmin=med-15*mad,vmax=med+15*mad,interpolation='none',cmap=mpl.cm.RdBu_r)
gc.add_label(0.15,0.92,"WISE W2",relative=True,color="black",fontsize=12)
gc.tick_labels.hide()
gc.axis_labels.hide()
gc.save("%s_W2.png"%name,dpi=150)
gc.close()
#
