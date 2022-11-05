import numpy as np
def rastr2deg(inp):
   w = inp.split(":")
   out = 0.0
   place_value = 1.0
   for w1 in w:
      out+=float(w1)*place_value
      place_value/=60.0
   return out*15.0

def decstr2deg(inp):
   w = inp.split(":")
   sign = np.sign(float(w[0]))
   place_value = 1.0
   out = 0.0
   for w1 in w:
      out+=np.absolute(float(w1))*place_value
      place_value/=60.0
   return sign*out
#
#
name = []
rastr = []
decstr = []
radeg = []
decdeg = []
per = []
#
f = open("pulsars.txt")
lines = f.readlines()
for line in lines:
   w = line.strip("\n").split(",")
   name.append(w[0])
   rastr.append(w[1])
   decstr.append(w[2])
   radeg.append(rastr2deg(w[1]))
   decdeg.append(decstr2deg(w[2]))
   print (w)
   if "*" in w[3]:
      per.append(np.nan)
   else:
      per.append(float(w[3]))

np.savez("atnf_pulsars.npz",name=name,radeg=radeg,decdeg=decdeg,rastr=rastr,decstr=decstr,period = per)
