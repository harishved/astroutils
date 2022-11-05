# SOME IMPORTANT PLASMA PARAMETERS COMMONLY USED IN CALCULATIONS
# AUTHOR: HARISH VEDANTHAM
#
from numpy import pi as PI
from const import *
def cycl_freq(B):
	return C_e*B/(C_me*C_c) / (2.0*PI)
def plasma_freq(ne):
	return (4*PI*ne*C_e**2/C_me)**0.5 / (2*PI)
def culog(T,nu):
	return 20.0+np.log(T/1e6/(nu/100e6))
def alfven_speed(B,ne):
	return C_c*B/(B**2+4*PI*C_mp*ne*C_c**2)**0.5
def plasma_beta(B,T,ne):
	return B**2/(8*PI)/(ne*C_k*T)
def electron_thermal_speed(T):
	return (C_k*T/C_me)**0.5
def ion_thermal_speed(T):
	return (C_k*T/C_mp)**0.5

