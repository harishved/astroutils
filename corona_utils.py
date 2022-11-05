#  Routines important for stellar coronal structure calculations
#  Author : Harish Vedantham, ASTRON


import numpy as np
import matplotlib.pyplot as plt
from const import *


def parker_wind_r2v(r):	
    # Radius to valocity mapping for stellar wind solution
    # Radius is in units of sonic radius
    # Oupyt velocity is in units of sound speed
    # Takes parker.npz look-up table as input (same directory)
    # parker.npz has (alist,sol) pairs where sol = ln(v^2/vs^2_
    # a is the R.H.S of the wind equation given below.
    # (v^2/vs^2-ln(v^2/vs^2) = fln(r/rs) + 4rs/r -3)

    a = 4*np.log(r) + 4/r -3	

    tp = np.load("parker.npz")	# Load look-up table
    alist = tp['alist']
    sol_l = tp['sol_l']
    sol_u = tp['sol_u']
    a0 = alist[0]
    # Assumes lookup table is linearly spaced
    da = alist[1]-alist[0]
    nmax = len(sol_l)	
    vlist_l = np.zeros(r.shape)	# Lower branch velocity
    vlist_u = np.zeros(r.shape)	# Upper branch velocity
    for i in range(len(r)):		# For each requested radius
        I = int(np.floor((a[i]-a0)/da))	
        # Find nearest left neighbor (for interpolation)
        if I <= 0: # If outside the support then extrapolate lienarly
            i1 = 0; 
            i2 = 1;
        elif I >= nmax-1: # If outside the support then 
                          # extrapolate linearly
            i1 = nmax-2; 
            i2 = nmax-1;
        else:	# Else lienarly interpolate
            i1 = I; 
            i2 = I+1
		
        y1 = sol_l[i1]; 
        y2 = sol_l[i2]; 
        x1 = alist[i1]; 
        x2 = alist[i2]

        y = (y1*(a[i]-x2) - y2*(a[i]-x1)) / (x1-x2)
        vlist_l[i] = np.sqrt(np.exp(y))

        y1 = sol_u[i1]; 
        y2 = sol_u[i2]; 
	x1 = alist[i1]; 
	x2 = alist[i2]
	y = (y1*(a[i]-x2) - y2*(a[i]-x1)) / (x1-x2)
	vlist_u[i] = np.sqrt(np.exp(y))

    return vlist_l,vlist_u
    # For Parker wind for r/rs<1 take the lower branch 
    # and for r/rs>1 take the upper branch
    # Example usage follows:
    '''
	r=np.logspace(-1,1,100)
	vl,vu = wind_r2v(r)
	v=np.zeros(r.shape)
	v[r<1]=vl[r<1]
	v[r>1]=vu[r>1]
	#plt.plot(r,vl)
	#plt.plot(r,vu)
	plt.semilogy(r,1.0/v/r**2,'k')
	plt.show()
	plt.close()
    '''



def v_escape(Mstar,Rstar): # Escape velocity at stellar surface
    return (2*C_G*C_Msun*Mstar/C_Rsun/Rstar)**0.5


def density_scale_height(Mstar,Rstar,T):  # Hydrostatic isothermal
    # T is coronal temp; Mstar and Rstar are stellar parameters 
    # in units of solar values
    return 2*C_k*T*(C_Rsun*Rstar)**2 / (C_G*C_Msun*Mstar*C_mp)


def sonic_point(Mstar,Rstar,T):  # rs/r* = vesc^2/(4vs^2)
    return (v_escape(Mstar,Rstar)/v_sound(T)/2)**2 * Rstar*C_Rsun


def v_sound(T):	# Sound speed vs^2 = kT/(mu*mp)
    # ASsume mean molucular weight of 0,6 (for pure H it is 1/2)
    return (C_k*T/(0.6*C_mp))**0.5


def density_profile(Mstar,Rstar,T,rvec):  
    # Rvec is the radius in units of R* (rvec=1 is at the surface)
    rs = sonic_point(Mstar,Rstar,T)/(Rstar*C_Rsun)
    # Sonic point in units of R*	
    hp = density_scale_height(Mstar,Rstar,T)/(Rstar*C_Rsun)	
    # Scale height in units of R*
    density_norm = np.zeros(rvec.shape)	 # Initialize o/p vec
    I = np.where(rvec<rs) # DOmain of hydrostatic solution
    density_norm[I] = np.exp(-(1-1/rvec[I])*1.0/hp)	
    # Hydrostartic atmopshere
    I = np.argmin(np.absolute(rvec-rs))	
    # Closest point to the sonic radius (take this as Rs hereon 
    # to ensure smooth transition)
    vl,vu = parker_wind_r2v(rvec[I:])  # Parker wind velocity
    density_norm[I:] = (rvec[I:]/rvec[I])**-2 *vu**-1 \
                       * np.exp(-(1-1/rvec[I])*1/hp)  # Wind solution
    return density_norm 		

def emission_measure(Mstar,Rstar,T):
    # Calculate the emission measure for hydrostatic + wind solution
    # Assume hydrostatic atmopshere until rs and wind sol. for >rs
    # This returns the normalized EM.
    # It has to be multiplied by n0^2 where n0 = base density at r*

    hp = density_scale_height(Mstar,Rstar,T) / (Rstar*C_Rsun)
    rs = sonic_point(Mstar,Rstar,T) / (Rstar*C_Rsun)

    # First compute the EM from the hydrostatic part
    if rs <= 1:
        # Sonic point is inside the stellar radius
        # There is no hydrostatic part
        em_hydrostatic = 0.0
    else:
        rvec = np.arange(1.0,rs,min(hp/100.0,(rs-1)/100.0))
        d_rvec = rvec[1]-rvec[0]
        integrand = 2.0*np.pi*(Rstar*C_Rsun)**3 \
                    * rvec**2*np.exp(-1/hp*(1-1/rvec))
        em_hydrostatic = d_rvec * np.sum(integrand)

    # Next compute the EM from wind part
    rvec = np.arange(rs,5*rs,rs/100)
    d_rvec = rvec[1]-rvec[0]
    dummy,vr = parker_wind_r2v(rvec)
    nr = (rvec/rs)**-2 /vr *np.exp(-(1-1/rs)*1/hp)
    integrand = 2.0*np.pi*(Rstar*C_Rsun)**3 * nr**2
    em_wind = d_rvec * np.sum(integrand)

    return em_wind+em_hydrostatic

def Xray_lum(Emin,Emax,T,em):
    # Calculate the X-ray free-free luminosity
    # T = Coronal temperature
    # Emin = Min energy for calculationb (in eV)
    # Emax = Max ebergy for in eV
    # em = Volume averaged emission measure in cgs units 

    he_abundance = 6.08e-2 # From doi:10.1088/0004-637X/755/1/33

    evec = np.arange(Emin,Emax,Emin/10)
    d_evec = evec[1]-evec[0]
    uvec = evec*C_ev2ergs/(C_k*T)

    f_k = 6.8e-38	# See eqn 5.14a of Bybicki & Lightman
    Ry_by_kT = 13.6*C_ev2ergs/(C_k*T)
    
    emissivity_vec = f_k*em*np.exp(-u) *T**-0.5 \
                   * (he_abundance*4*get_gff(Ry_by_kT*4,uvec) \
                   + (1-he_abundance)*get_gff(Ry_by_kT,uvec))
    # This is volume avg emissivity in ergs/s/Hz
    # This has to be averaged over frequency (or energy)
    # dnu = dE/h
    luminosity = d_evec/C_h*np.sum(emissivity_vec)
    return luminosity
    
def get_gff(gamma2,u):
    # Free-free Gaunt factor from Sutherland 1998
    # This is the gff averaged over a maxwell sidtribution
    # gamma2 = Z^2Ry/(kT)
    # u = hnu/(kT)

    tp = np.loadtxt("gff_maxwell.txt",comments="#")
    # File format us u gamma2 gff
    # The stride is 81 samples
    # Data are log-spaced

   log_u = np.log(tp[:,0])
   log_g2 = np.log(tp[:,1])
   gff = tp[:,2]

   log_u_grid = log_u[::81]
   log_g2_grid = log_g2[:81]

   n_u_grid = len(log_u_grid)
   n_g2_grid = len(log_g2_grid)

   Iu = np.argmin(np,absolute(log_u_grid-np.log(u)))
   Ig2 = np.argmin(np,absolute(log_g2_grid-np.log(g2)))

   return gff[Iu*81+Ig2]
