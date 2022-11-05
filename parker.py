# Make a lookup table for solution of Parler wind
# Author: Harish Vedantham, ASTRON
# 
import numpy as np

niter=100	# Number of inerarion for Newton methid to find roots
alist = np.linspace(1.0,100.0,10000)	# a = R.H.S of wind solution = 4ln(r/rs) + 4rs/r -3 
sol_l=np.zeros(alist.shape)	# Lower solution branch (here sol = ln(v^2/vs^2)
sol_u=np.zeros(alist.shape)	# Upper solution branch
icnt=0				# Dummy increment variable for looping
for a in alist:			# For each R.H.S value
	x0=-1.0			# Initialize lower branch solution
	for i in range(niter):	#
		x1 = x0 - (np.exp(x0)-x0-a)/(np.exp(x0)-1)	# Newtown method solution update
		x0=x1;	
	sol_l[icnt]=x1
	icnt+=1
icnt=0		# Now reset counter and solve upper branch
for a in alist:
	x0=1.0	# Upper branch initialization
	for i in range(niter):
		x1 = x0 - (np.exp(x0)-x0-a)/(np.exp(x0)-1)
		x0=x1;
	sol_u[icnt]=x1
	icnt+=1
# Save results in npz file
np.savez("parker.npz",alist=alist,sol_l=sol_l,sol_u=sol_u)
