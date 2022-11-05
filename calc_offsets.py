# CALCULATE ROUGH OFFSETS FOR OPTICAL SPECTROSCOPY
# THIS IS QUICK AND DIRTY BUT SHOULD BE ACCURATE ENOUGH FOR MOST AOPPLICATIONS
# AUTHOR: HARISH VEDANTHAM
#
import numpy as np

def fmt(ra,dec):
   h = np.floor(ra/15)
   m = np.floor((ra/15-h)*60)
   s = (ra/15-h-m/60.0)*3600
   print ("%02d %02d %.4f"%(h,m,s))

   dec_sgn = np.sign(dec)
   deca = np.absolute(dec)
   h = np.floor(deca)
   m = np.floor((deca-h)*60)
   s = (deca-h-m/60.0)*3600
   print ("%02d %02d %.4f"%(h*dec_sgn,m,s))
#


ra2 = (18.+ 48./60+ 29.22/3600)*15.0
dec2 = -(02.+ 10./60 + 01.32/3600)

rao = 282.11740788627
deco = -2.16667854905


cd = np.cos(0.5*np.pi/180.0*(deco+dec2))
print ("Offset star")
fmt(rao,deco)

dra = (ra2-rao)*cd*3600
ddec = (dec2-deco)*3600
print ("Offsets = %.2f, %.2f"%(dra,ddec))

