# TAKE IN A LIST OF SOURCE COMMON NAMES AND MAKE A LIST THAT CAN BE UPLOADED
# TO THE LOFAR OBS TOOL. 
#
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
import numpy as np

list = ["SIG CRB", "FK COM", "BF LYN", "DM UMA", "EV DRA", "OU AND", "WX UMA", "LP212-62"]
run = 0
for name in list:
   src = SkyCoord.from_name(name)
   ra = src.ra.value/15
   dec = src.dec.value
   rah = np.floor(ra)
   ram = np.floor(60*(ra-rah))
   ras = 3600*(ra-rah-ram/60.0)
   sgn = np.sign(dec)
   dec = np.absolute(dec)
   decd = np.floor(dec)
   decm = np.floor(60*(dec-decd))
   decs = 3600*(dec-decd-decm/60.0)
   print ("%s,%d:%d:%.2f,%d:%d:%.1f,j2000,%d,N,15600s,0.001,144e6,0.0,102-347,A,A, Obs. start LST %02d:%02d "%(name.replace(" ","-"),rah,ram,ras,decd*sgn,decm,decs, run, rah-2,ram))
   run+=1
