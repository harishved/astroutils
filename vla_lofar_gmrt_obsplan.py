# ELEVATION PLOT AT JVLA, LOFAR AND GMRT THAT CAN BE USED TO PLAN
# SIMULTANEOUS OBSERVATIONS
# USAGE: python3 vla_lofar_gmrt_obsplan.py <SRC NAME> <DATE>
# E.G. python3 vla_lofar_gmrt_obsplan.py "WX UMA" 2022-12-12
# 
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
import sys
import numpy as np
import matplotlib.pyplot as plt
from colors import lc

src_name = sys.argv[1]
date = sys.argv[2]
print ("Src name"+src_name)
src = SkyCoord.from_name(src_name)
vla = EarthLocation(lat=34.0784*u.deg, lon=-107.6184*u.deg, height=2124.0*u.m)
lofar= EarthLocation(lat=52.9153*u.deg, lon=6.8698*u.deg, height=0.0*u.m)
gmrt = EarthLocation(lat=19.0919*u.deg, lon=74.0506*u.deg, height=680*u.m)


utvec = np.linspace(0,1.0,100) # in fraction of days
t0 = Time("%sT00:00:00.0"%date, format="isot",scale="utc")
tvec = [Time(x,format="mjd") for x in t0.mjd+utvec]

plt.figure(figsize=(10/1.5,6/1.5))
vla_altaz = src.transform_to(AltAz(obstime=tvec,location=vla)).alt.value
plt.plot(utvec*24,vla_altaz,label="VLA",color=lc[0],linewidth=2)
lofar_altaz = src.transform_to(AltAz(obstime=tvec,location=lofar)).alt.value
plt.plot(utvec*24,lofar_altaz,label="LOFAR",color=lc[2],linewidth=2)
gmrt_altaz = src.transform_to(AltAz(obstime=tvec,location=gmrt)).alt.value
plt.plot(utvec*24,gmrt_altaz,label="GMRT",color=lc[4],linewidth=2)
plt.title("Src: %s, Date: %s"%(src_name,date))
plt.legend()
plt.plot([0,24],[30,30],'k--')
#for i in range(40,90,10):
#   plt.plot([0,24],[i,i],'0.5')

plt.plot([7.5,9.5],[60,60],linewidth=4,alpha=0.6,color=lc[4])
plt.plot([11.5,13.5],[60,60],linewidth=4,alpha=0.6,color=lc[0])
plt.plot([8.5,12.5],[54,54],linewidth=4,alpha=0.6,color=lc[2])

plt.text(s="LOFAR",x=10.5,y=54,color=lc[2],fontsize=14)
plt.text(s="VLA",x=12.5,y=60,color=lc[0],fontsize=14)
plt.text(s="GMRT",x=8.5,y=60,color=lc[4],fontsize=14)

plt.xlabel("UT / hours",fontsize=14)
plt.ylabel("Altitude / degree",fontsize=14)
#plt.fill_between([7.5,9.5],[0,0],[90,90],color=lc[0],alpha=0.5)
#plt.fill_between([11.5,13.5],[0,0],[90,90],color=None, edgecolor=lc[4],alpha=0.5,hatch="/")
#plt.fill_between([8.5,12.5],[0,0],[90,90],color=None, edgecolor=lc[2],alpha=0.5,hatch='x')

plt.minorticks_on()
plt.xlim([6,15])
plt.ylim([0,90])
plt.tight_layout()   
plt.savefig("obsplan.pdf")
print("Plot saved to obsplan.pdf")
plt.show()

plt.close()
