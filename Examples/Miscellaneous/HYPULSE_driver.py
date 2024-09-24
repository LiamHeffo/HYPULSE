"""
Shock and Detonation Toolbox Demo Program

Description of Program. What does it do? 
################################################################################
Theory, numerical methods and applications are described in the following report:

SDToolbox - Numerical Tools for Shock and Detonation Wave Modeling,
Explosion Dynamics Laboratory, Contributors: S. Browne, J. Ziegler,
N. Bitter, B. Schmidt, J. Lawson and J. E. Shepherd, GALCIT
Technical Report FM2018.001 Revised January 2021.

Please refer to LICENCE.txt or the above report for copyright and disclaimers.
http://shepherd.caltech.edu/EDL/PublicResources/sdt/

Written by Liam Heffernan using routines from the SDToolbox.
################################################################################ 
"""
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

from sdtoolbox.thermo import soundspeed_fr,soundspeed_eq
from sdtoolbox.postshock import CJspeed,PostShock_eq

# Initial state specification:
# P1 = Initial Pressure  
# T1 = Initial Temperature 
# U = Shock Speed 
# q = Initial Composition 
# mech = Cantera mechanism File name
# plots = True


# set cold driver state:
#P4 = 34e6
#T4 = 300
#gamma4 = 1.667
#R4 = 2077.1
#a4 = np.sqrt(gamma4*R4*T4)


#set initial detoantion driver fill state
#P4_P100 = 200
#P100 = 34e6/P4_P100
#T100 = 300
#q100 = 'O2:15 H2:30 HE:55'

def HYPULSE_driver(P4, T4, q4, P100, T100, q100):
		"""
		Purpose, What it does, whatever.
		"""

		mech = 'sandiego20161214_H2only.cti'

		# State 400 gas object
		gas4 = ct.Solution(mech)
		gas4.TPX = T4, P4, q4
		a4 = soundspeed_fr(gas4)
		D4 = gas4.density
		gamma4 = a4*a4*D4/P4 #check this

		# State 100 gas object
		gas100 = ct.Solution(mech)
		gas = ct.Solution(mech)
		gas100.TPX = T100, P100, q100
		a100_fr = soundspeed_fr(gas100)
		D100 = gas100.density
		gamma100_fr =  a100_fr*a100_fr*D100/P100;

		print('\nDetonation Driver Fill State:');
		print(' Composition ' + q100)
		print(' Pressure %.2f (MPa) ' % (P100*10**(-6)) )
		print(' Temperature %.2f (K) ' % (T100) )

		print('Cold Driver Fill State:')
		print(' Composition: He: 100 ' )
		print(' Pressure %.2f (MPa) ' % (P4*10**(-6)) )
		print(' Temperature %.2f (K) ' % (T4) )
		print(' a4 (frozen) %.2f (m/s)' % (a4) )
		print(' gamma4 (frozen) %.4f '% (gamma4) )

		print('\nComputing CJ state and isentrope for '+q100+' using '+mech)
		
		# compute CJ speed
		cj_speed = CJspeed(P100, T100, q100, mech);
		
		# compute equilibrium CJ state
		gas = PostShock_eq(cj_speed, P100, T100, q100, mech)

		T2 = gas.T
		P2 = gas.P
		D2 = gas.density
		V2 = 1./D2
		S2  = gas.entropy_mass
		w2 = D100*cj_speed/D2
		u2 = cj_speed-w2
		a2_eq = soundspeed_eq(gas)
		a2_fr = soundspeed_fr(gas)
		gamma2_fr =  a2_fr*a2_fr*D2/P2;
		gamma2_eq =  a2_eq*a2_eq*D2/P2;

		print('CJ speed = %.2f (m/s)' % (cj_speed))
		print ('CJ State')
		print(' Pressure %.2f (MPa) ' % (P2*10**(-6)) )
		print(' Temperature %.2f (K) ' % (T2) )
		#print(' Density %.3f (kg/m3) ' % (D2) )
		#print(' Entropy %.3f (J/kg-K) ' % (S2) )
		#print(' w2 (frozen) %.2f (m/s)' % (w2) )
		print(' u2 (frozen) %.2f (m/s)' % (u2) )
		#print(' a2 (frozen) %.2f (m/s)' % (a2_fr) )
		print(' a2 (equilibrium) %.2f (m/s)' % (a2_eq) )
		#print(' gamma2 (frozen) %.4f '% (gamma2_fr) )   
		print(' gamma2 (equilibrium) %.4f '% (gamma2_eq) )   

		#compute the effective value of q based on two-gamma model
		M1 = cj_speed/a100_fr;
		eparam = a100_fr**2*(M1**(-2)*(gamma2_eq/gamma100_fr)**2 \
*(1+gamma100_fr*M1**2)**2/(2*(gamma2_eq**2-1))-1/(gamma100_fr-1)-M1**2/2)
		print('Detonation CJ Mach number %.2f' % (M1))
		print('2-gamma energy parameter q %.3f (J/kg)' % (eparam))

		# Find points on the isentrope connected to the CJ state, evaluate velocity
		# in Taylor wave using trapezoidal rule to evaluate the Riemann function
		npoints=50000
		vv = V2
		V = np.zeros(npoints,float)
		P = np.zeros(npoints,float)
		D = np.zeros(npoints,float)
		a = np.zeros(npoints,float)
		u = np.zeros(npoints,float)
		T = np.zeros(npoints,float)
		V[1] = V2
		P[1] = P2
		D[1] = D2
		a[1] = a2_eq
		u[1] = u2
		T[1] = T2
		print('\nGenerating points on isentrope and computing Taylor wave velocity')

		i = 1
		has_converged = False

		while u[i] > 0 and i < npoints:
			i = i+1
			vv = vv*1.0001
			x = gas.X
			gas.SVX = S2, vv, x
			gas.equilibrate('SV')
			P[i] = gas.P
			D[i] = gas.density
			V[i] = 1/D[i]
			T[i] = gas.T
			a[i] = soundspeed_eq(gas)
			u[i] = u[i-1] + 0.5*(P[i]-P[i-1])*(1./(D[i]*a[i]) + 1./(D[i-1]*a[i-1]))

			# expanding state 4 to state 400 conditions:
			Re_tol = 1e-3
			P3 = P4*(1 - 0.5*(gamma4 - 1)*(u[i]/a4))**(2*gamma4/(gamma4 -1))
			if abs((P3 - P[i])/P3) < Re_tol:
				has_converged = True
				break

		nfinal = i
		P400 = P[nfinal]
		u400 = u[nfinal]
		a400 = a[nfinal]
		V400 = V[nfinal]
		T400 = T[nfinal]

		# State 400: 
		print('State 400 ')
		print(' Pressure %.2f (MPa)' % (P400*10**(-6)))
		print(' Temperature %.2f (K)' % (T400))
		#print(' Volume %.2f (m3/kg)' % (V400))
		print(' Velocity %.2f (m/s)' % (u400))
		
		# evaluate final state 400 to get sound speeds and effective gammas
		x = gas.X
		gas.SVX = S2, V400, x
		gas.equilibrate('SV')
		a400_fr = soundspeed_fr(gas);
		gamma400_fr =  a400_fr**2/(P400*V400);
		gamma400_eq =  a400**2/(V400*P400);
		#print(' Sound speed (frozen) %.2f (m/s)' % (a400_fr))
		print(' Sound speed (equilibrium) %.2f (m/s)' % (a400))
		#print(' Gamma (frozen) %.4f' % (gamma400_fr))
		print(' Gamma (equilibrium) %.4f' % (gamma400_eq))

		# State 3:
		u3 = u400
		T3 = T4*(P3/P4)**((gamma4-1)/gamma4)
		print('State 3:')
		print(' Pressure %.2f (MPa)' % (P3*10**(-6)))
		print(' Velocity %.2f (m/s)' % u3)
		print(' Temperature %.2f (K)' % (T3))

		# Calculating cold-driver-equivalent numbers from Lu:
		pe = P400*(1 + (gamma400_eq - 1)/2*(u400/a400_fr))**(2*gamma400_eq/(gamma400_eq - 1))
		ae = a400*(1 + (gamma400_eq - 1)/2*(u400/a400_fr))
		print('State 400 Equivalent')
		print(' pe = %.2f (MPa)' % (pe*10**(-6)))
		print(' ae = %.2f (m/s)' % (ae) )
		
		return P400, u400, a400, V400, T400


def HYPULSE_driver_equivalentconditions(P4, T4, q4, P100, T100, q100):
		"""
		Purpose, What it does, whatever.
		"""

		mech = 'sandiego20161214_H2only.cti'

		# State 400 gas object
		gas4 = ct.Solution(mech)
		gas4.TPX = T4, P4, q4
		a4 = soundspeed_fr(gas4)
		D4 = gas4.density
		gamma4 = a4*a4*D4/P4 #check this

		# State 100 gas object
		gas100 = ct.Solution(mech)
		gas = ct.Solution(mech)
		gas100.TPX = T100, P100, q100
		a100_fr = soundspeed_fr(gas100)
		D100 = gas100.density
		gamma100_fr =  a100_fr*a100_fr*D100/P100;

		print('\nDetonation Driver Fill State:');
		print(' Composition ' + q100)
		print(' Pressure %.2f (MPa) ' % (P100*10**(-6)) )
		print(' Temperature %.2f (K) ' % (T100) )

		print('Cold Driver Fill State:')
		print(' Composition: He: 100 ' )
		print(' Pressure %.2f (MPa) ' % (P4*10**(-6)) )
		print(' Temperature %.2f (K) ' % (T4) )
		print(' a4 (frozen) %.2f (m/s)' % (a4) )
		print(' gamma4 (frozen) %.4f '% (gamma4) )

		print('\nComputing CJ state and isentrope for '+q100+' using '+mech)
		
		# compute CJ speed
		cj_speed = CJspeed(P100, T100, q100, mech);
		
		# compute equilibrium CJ state
		gas = PostShock_eq(cj_speed, P100, T100, q100, mech)

		T2 = gas.T
		P2 = gas.P
		D2 = gas.density
		V2 = 1./D2
		S2  = gas.entropy_mass
		w2 = D100*cj_speed/D2
		u2 = cj_speed-w2
		a2_eq = soundspeed_eq(gas)
		a2_fr = soundspeed_fr(gas)
		gamma2_fr =  a2_fr*a2_fr*D2/P2;
		gamma2_eq =  a2_eq*a2_eq*D2/P2;

		print('CJ speed = %.2f (m/s)' % (cj_speed))
		print ('CJ State')
		print(' Pressure %.2f (MPa) ' % (P2*10**(-6)) )
		print(' Temperature %.2f (K) ' % (T2) )
		#print(' Density %.3f (kg/m3) ' % (D2) )
		#print(' Entropy %.3f (J/kg-K) ' % (S2) )
		#print(' w2 (frozen) %.2f (m/s)' % (w2) )
		print(' u2 (frozen) %.2f (m/s)' % (u2) )
		#print(' a2 (frozen) %.2f (m/s)' % (a2_fr) )
		print(' a2 (equilibrium) %.2f (m/s)' % (a2_eq) )
		#print(' gamma2 (frozen) %.4f '% (gamma2_fr) )   
		print(' gamma2 (equilibrium) %.4f '% (gamma2_eq) )   

		#compute the effective value of q based on two-gamma model
		M1 = cj_speed/a100_fr;
		eparam = a100_fr**2*(M1**(-2)*(gamma2_eq/gamma100_fr)**2 \
*(1+gamma100_fr*M1**2)**2/(2*(gamma2_eq**2-1))-1/(gamma100_fr-1)-M1**2/2)
		print('Detonation CJ Mach number %.2f' % (M1))
		print('2-gamma energy parameter q %.3f (J/kg)' % (eparam))

		# Find points on the isentrope connected to the CJ state, evaluate velocity
		# in Taylor wave using trapezoidal rule to evaluate the Riemann function
		npoints=50000
		vv = V2
		V = np.zeros(npoints,float)
		P = np.zeros(npoints,float)
		D = np.zeros(npoints,float)
		a = np.zeros(npoints,float)
		u = np.zeros(npoints,float)
		T = np.zeros(npoints,float)
		V[1] = V2
		P[1] = P2
		D[1] = D2
		a[1] = a2_eq
		u[1] = u2
		T[1] = T2
		print('\nGenerating points on isentrope and computing Taylor wave velocity')

		i = 1
		has_converged = False

		while u[i] > 0 and i < npoints:
			i = i+1
			vv = vv*1.0001
			x = gas.X
			gas.SVX = S2, vv, x
			gas.equilibrate('SV')
			P[i] = gas.P
			D[i] = gas.density
			V[i] = 1/D[i]
			T[i] = gas.T
			a[i] = soundspeed_eq(gas)
			u[i] = u[i-1] + 0.5*(P[i]-P[i-1])*(1./(D[i]*a[i]) + 1./(D[i-1]*a[i-1]))

			# expanding state 4 to state 400 conditions:
			Re_tol = 1e-3
			P3 = P4*(1 - 0.5*(gamma4 - 1)*(u[i]/a4))**(2*gamma4/(gamma4 -1))
			if abs((P3 - P[i])/P3) < Re_tol:
				has_converged = True
				break

		nfinal = i
		P400 = P[nfinal]
		u400 = u[nfinal]
		a400 = a[nfinal]
		V400 = V[nfinal]
		T400 = T[nfinal]

		# State 400: 
		print('State 400 ')
		print(' Pressure %.2f (MPa)' % (P400*10**(-6)))
		print(' Temperature %.2f (K)' % (T400))
		#print(' Volume %.2f (m3/kg)' % (V400))
		print(' Velocity %.2f (m/s)' % (u400))
		
		# evaluate final state 400 to get sound speeds and effective gammas
		x = gas.X
		gas.SVX = S2, V400, x
		gas.equilibrate('SV')
		a400_fr = soundspeed_fr(gas);
		gamma400_fr =  a400_fr**2/(P400*V400);
		gamma400_eq =  a400**2/(V400*P400);
		#print(' Sound speed (frozen) %.2f (m/s)' % (a400_fr))
		print(' Sound speed (equilibrium) %.2f (m/s)' % (a400))
		#print(' Gamma (frozen) %.4f' % (gamma400_fr))
		print(' Gamma (equilibrium) %.4f' % (gamma400_eq))

		# State 3:
		u3 = u400
		T3 = T4*(P3/P4)**((gamma4-1)/gamma4)
		print('State 3:')
		print(' Pressure %.2f (MPa)' % (P3*10**(-6)))
		print(' Velocity %.2f (m/s)' % u3)
		print(' Temperature %.2f (K)' % (T3))

		# Calculating cold-driver-equivalent numbers from Lu:
		pe = P400*(1 + (gamma400_eq - 1)/2*(u400/a400_fr))**(2*gamma400_eq/(gamma400_eq - 1))
		ae = a400*(1 + (gamma400_eq - 1)/2*(u400/a400_fr))
		print('State 400 Equivalent')
		print(' pe = %.2f (MPa)' % (pe*10**(-6)))
		print(' ae = %.2f (m/s)' % (ae) )
		
		return pe, ae

def JF16_driver(P100, T100, q100):
		"""
		Purpose, What it does, whatever.
		"""

		mech = 'sandiego20161214_H2only.cti'

		# State 100 gas object
		gas100 = ct.Solution(mech)
		gas = ct.Solution(mech)
		gas100.TPX = T100, P100, q100
		a100_fr = soundspeed_fr(gas100)
		D100 = gas100.density
		gamma100_fr =  a100_fr*a100_fr*D100/P100;

		print('\nDetonation Driver Fill State:');
		print(' Composition ' + q100)
		print(' Pressure %.2f (MPa) ' % (P100*10**(-6)) )
		print(' Temperature %.2f (K) ' % (T100) )

		print('\nComputing CJ state and isentrope for '+q100+' using '+mech)
		
		# compute CJ speed
		cj_speed = CJspeed(P100, T100, q100, mech);
		
		# compute equilibrium CJ state
		gas = PostShock_eq(cj_speed, P100, T100, q100, mech)

		T2 = gas.T
		P2 = gas.P
		D2 = gas.density
		V2 = 1./D2
		S2  = gas.entropy_mass
		w2 = D100*cj_speed/D2
		u2 = cj_speed-w2
		a2_eq = soundspeed_eq(gas)
		a2_fr = soundspeed_fr(gas)
		gamma2_fr =  a2_fr*a2_fr*D2/P2;
		gamma2_eq =  a2_eq*a2_eq*D2/P2;

		print('CJ speed = %.2f (m/s)' % (cj_speed))
		print ('CJ State')
		print(' Pressure %.2f (MPa) ' % (P2*10**(-6)) )
		print(' Temperature %.2f (K) ' % (T2) )
		#print(' Density %.3f (kg/m3) ' % (D2) )
		#print(' Entropy %.3f (J/kg-K) ' % (S2) )
		#print(' w2 (frozen) %.2f (m/s)' % (w2) )
		print(' u2 (frozen) %.2f (m/s)' % (u2) )
		#print(' a2 (frozen) %.2f (m/s)' % (a2_fr) )
		print(' a2 (equilibrium) %.2f (m/s)' % (a2_eq) )
		#print(' gamma2 (frozen) %.4f '% (gamma2_fr) )   
		print(' gamma2 (equilibrium) %.4f '% (gamma2_eq) )   

		#compute the effective value of q based on two-gamma model
		M1 = cj_speed/a100_fr;
		eparam = a100_fr**2*(M1**(-2)*(gamma2_eq/gamma100_fr)**2 \
*(1+gamma100_fr*M1**2)**2/(2*(gamma2_eq**2-1))-1/(gamma100_fr-1)-M1**2/2)
		print('Detonation CJ Mach number %.2f' % (M1))
		print('2-gamma energy parameter q %.3f (J/kg)' % (eparam))

		# Find points on the isentrope connected to the CJ state, evaluate velocity
		# in Taylor wave using trapezoidal rule to evaluate the Riemann function
		npoints=50000
		vv = V2
		V = np.zeros(npoints,float)
		P = np.zeros(npoints,float)
		D = np.zeros(npoints,float)
		a = np.zeros(npoints,float)
		u = np.zeros(npoints,float)
		T = np.zeros(npoints,float)
		V[1] = V2
		P[1] = P2
		D[1] = D2
		a[1] = a2_eq
		u[1] = u2
		T[1] = T2
		print('\nGenerating points on isentrope and computing Taylor wave velocity')

		i = 1
		has_converged = False

		while u[i] > 0 and i < npoints:
			i = i+1
			vv = vv*1.0001
			x = gas.X
			gas.SVX = S2, vv, x
			gas.equilibrate('SV')
			P[i] = gas.P
			D[i] = gas.density
			V[i] = 1/D[i]
			T[i] = gas.T
			a[i] = soundspeed_eq(gas)
			u[i] = u[i-1] + 0.5*(P[i]-P[i-1])*(1./(D[i]*a[i]) + 1./(D[i-1]*a[i-1]))


		nfinal = i
		P400 = P[nfinal]
		u400 = u[nfinal]
		a400 = a[nfinal]
		V400 = V[nfinal]
		T400 = T[nfinal]

		# State 400: 
		print('State 400 ')
		print(' Pressure %.2f (MPa)' % (P400*10**(-6)))
		print(' Temperature %.2f (K)' % (T400))
		#print(' Volume %.2f (m3/kg)' % (V400))
		print(' Velocity %.2f (m/s)' % (u400))
		
		# evaluate final state 400 to get sound speeds and effective gammas
		x = gas.X
		gas.SVX = S2, V400, x
		gas.equilibrate('SV')
		a400_fr = soundspeed_fr(gas);
		gamma400_fr =  a400_fr**2/(P400*V400);
		gamma400_eq =  a400**2/(V400*P400);
		#print(' Sound speed (frozen) %.2f (m/s)' % (a400_fr))
		print(' Sound speed (equilibrium) %.2f (m/s)' % (a400))
		#print(' Gamma (frozen) %.4f' % (gamma400_fr))
		print(' Gamma (equilibrium) %.4f' % (gamma400_eq))

		# Calculating cold-driver-equivalent numbers from Lu:
		pe = P400*(1 + (gamma400_eq - 1)/2*(u400/a400_fr))**(2*gamma400_eq/(gamma400_eq - 1))
		ae = a400*(1 + (gamma400_eq - 1)/2*(u400/a400_fr))
		print('State 400 Equivalent')
		print(' pe = %.2f (MPa)' % (pe*10**(-6)))
		print(' ae = %.2f (m/s)' % (ae) )
		
		return P400, u400, a400, V400, T400, ae, pe


def HYPULSE_SET_driver(P4, T4, q4, P100, T100, q100):
		"""
		Purpose, What it does, whatever.
		"""

		mech = 'sandiego20161214_H2only.cti'

		# State 100 gas object
		gas100 = ct.Solution(mech)
		gas = ct.Solution(mech)
		gas100.TPX = T100, P100, q100
		a100_fr = soundspeed_fr(gas100)
		D100 = gas100.density
		gamma100_fr =  a100_fr*a100_fr*D100/P100;

		# State 400 gas object
		gas4 = ct.Solution(mech)
		gas4.TPX = T4, P4, q4
		a4 = soundspeed_fr(gas4)
		D4 = gas4.density
		gamma4 = a4*a4*D4/P4 #check this


		print('\nDetonation Driver Fill State:');
		print(' Composition ' + q100)
		print(' Pressure %.2f (MPa) ' % (P100*10**(-6)) )
		print(' Temperature %.2f (K) ' % (T100) )

		print('Cold Driver Fill State:')
		print(' Composition: He: 100 ' )
		print(' Pressure %.2f (MPa) ' % (P4*10**(-6)) )
		print(' Temperature %.2f (K) ' % (T4) )
		print(' a4 (frozen) %.2f (m/s)' % (a4) )
		print(' gamma4 (frozen) %.4f '% (gamma4) )

		print('\nComputing CJ state and isentrope for '+q100+' using '+mech)
		
		# compute CJ speed
		cj_speed = CJspeed(P100, T100, q100, mech);
		
		# compute equilibrium CJ state
		gas = PostShock_eq(cj_speed, P100, T100, q100, mech)

		T400 = gas.T
		P400 = gas.P
		D400 = gas.density
		V400 = 1./D400
		S400  = gas.entropy_mass
		w400 = D100*cj_speed/D400
		u400 = cj_speed-w400
		a400 = soundspeed_eq(gas)
		a400_eq = a400
		a400_fr = soundspeed_fr(gas)
		gamma400_fr =  a400_fr*a400_fr*D400/P400;
		gamma400_eq =  a400_eq*a400_eq*D400/P400;

		print('CJ speed = %.2f (m/s)' % (cj_speed))
		print ('CJ State')
		print(' Pressure %.2f (MPa) ' % (P400*10**(-6)) )
		print(' Temperature %.2f (K) ' % (T400) )
		print(' Density %.3f (kg/m3) ' % (D400) )
		#print(' Entropy %.3f (J/kg-K) ' % (S400) )
		print(' w400 (frozen) %.2f (m/s)' % (w400) )
		print(' u400 (frozen) %.2f (m/s)' % (u400) )
		print(' a400 (frozen) %.2f (m/s)' % (a400_fr) )
		print(' a400 (equilibrium) %.2f (m/s)' % (a400_eq) )
		#print(' gamma400 (frozen) %.4f '% (gamma2_fr) )   
		print(' gamma400 (equilibrium) %.4f '% (gamma400_eq) )   

		pe = P400*(1 + (gamma400_eq - 1)/2*(u400/a400_fr))**(2*gamma400_eq/(gamma400_eq - 1))
		ae = a400*(1 + (gamma400_eq - 1)/2*(u400/a400_fr))

		#compute the effective value of q based on two-gamma model
		M1 = cj_speed/a100_fr;
		eparam = a100_fr**2*(M1**(-2)*(gamma400_eq/gamma100_fr)**2 \
*(1+gamma100_fr*M1**2)**2/(2*(gamma400_eq**2-1))-1/(gamma100_fr-1)-M1**2/2)
		print('Detonation CJ Mach number %.2f' % (M1))
		print('2-gamma energy parameter q %.3f (J/kg)' % (eparam))
	
		return P400, u400, a400, V400, T400, ae, pe
		

