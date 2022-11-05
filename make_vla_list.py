# TAKE IN A LIST OF SOURCE COMMON NAMES AND OUTPUT A FILE THAT CAN BE READ IN
# BY THE VLA OBS PLANNING TOOL
#
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
import numpy as np

list = ["SIG CRB", "FK COM", "II PEG", "DM UMA", "EV DRA", "OU AND", "WX UMA", "LP212-62"]

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
   decm = 60*(dec-decd)
   decs = 3600*(dec-decd-decm/60.0)
   print ("%s; %s; Equatorial; J2000; %d:%d:%.2f; %d:%d:%.2f; Barycentric; Optical; 0; N;"%(name.replace(" ","-"),name.replace(" ","-"),rah,ram,ras,decd*sgn,decm,decs))
