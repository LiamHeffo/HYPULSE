"""
Theory, numerical methods and applications are described in the following report:

SDToolbox - Numerical Tools for Shock and Detonation Wave Modeling,
Explosion Dynamics Laboratory, Contributors: S. Browne, J. Ziegler,
N. Bitter, B. Schmidt, J. Lawson and J. E. Shepherd, GALCIT
Technical Report FM2018.001 Revised January 2021.

Please refer to LICENCE.txt or the above report for copyright and disclaimers.
http://shepherd.caltech.edu/EDL/PublicResources/sdt/

Written by Liam Heffernan using routines from the SDToolbox.
"""
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import newton

from sdtoolbox.thermo import soundspeed_fr,soundspeed_eq
from sdtoolbox.postshock import CJspeed,PostShock_eq
from matplotlib.ticker import MaxNLocator

# Initial state specification:
# P1 = Initial Pressure  
# T1 = Initial Temperature 
# U = Shock Speed 
# q = Initial Composition 
# mech = Cantera mechanism File name
# plots = True

def HYPULSE_driver(P4, T4, q4, P100, T100, q100):
	"""
    This function models a HyPulse detonation driver experiment using a specified equation of state and
    generates a stage-wise sequence of thermodynamic and gas dynamical states for the detonation and
    its subsequent expansion, while comparing the model's state predictions with experimental observations.
    
    Parameters
    -----------
    P4, P100: float 
        The initial pressures (Pa) at states '4' and '100'. 
    
    T4, T100: float 
        The initial temperatures (Kelvin) in states '4' and '100'.
    
    q4, q100: str
        The initial mole fractions at states '4' and '100'.

    Returns
    --------
    P400: float 
        The pressure (Pa) at final state '400'.
    
    u400: float 
        The velocity (m/s) at final state '400'.
    
    a400: float 
        The sound speed (m/s) at final state '400'.
    
    V400: float 
        The specific volume (m^3/kg) at final state '400'.
    
    T400: float
        The temperature (K) at final state '400'.

    Raises
    ------
    CanteraError
        If any input parameter is not valid, or if Cantera encounters an error during the simulation.

    Notes
    -----
    This function prints out many intermediate values and stages of the experiment. It might 
    be preferable to comment out some of the print statements when using this function in a larger script.

    """

	mech = 'sandiego20161214_H2only.cti'

	# State 400 gas object
	gas4 = ct.Solution(mech)
	gas4.TPX = T4, P4, q4
	a4 = soundspeed_fr(gas4)
	D4 = gas4.density
	gamma4 = a4*a4*D4/P4

	# State 100 gas object
	gas100 = ct.Solution(mech)
	gas100.TPX = T100, P100, q100
	a100_fr = soundspeed_fr(gas100)
	D100 = gas100.density
	gamma100_fr =  a100_fr*a100_fr*D100/P100

	print('\nDetonation Driver Fill State:')
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
	cj_speed = CJspeed(P100, T100, q100, mech)
	
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
	gamma2_eq =  a2_eq*a2_eq*D2/P2

	print('CJ speed = %.2f (m/s)' % (cj_speed))
	print ('CJ State')
	print(' Pressure %.2f (MPa) ' % (P2*10**(-6)) )
	print(' Temperature %.2f (K) ' % (T2) )
	print(' u2 (frozen) %.2f (m/s)' % (u2) )
	print(' a2 (equilibrium) %.2f (m/s)' % (a2_eq) )
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
	print(' Velocity %.2f (m/s)' % (u400))

	# evaluate final state 400 to get sound speeds and effective gammas
	x = gas.X
	gas.SVX = S2, V400, x
	gas.equilibrate('SV')
	a400_fr = soundspeed_fr(gas);
	gamma400_eq =  a400**2/(V400*P400);

	# plt.figure(1)
	# plt.plot(a[1:nfinal], P[1:nfinal] * 10 **(-6))
	# plt.plot(a[1:nfinal], 10510*a[1:nfinal]*10**(-6), color='red')
	# plt.title('Isentropic expansion after detonation, CJ speed = %.1f m/s' % (cj_speed))
	# plt.xlabel(r'Equilibrium Sound Speed (m/s)')
	# plt.ylabel('Pressure (MPa)')
	# plt.gca().invert_xaxis()
	# plt.show()

	print(' Sound speed (equilibrium) %.2f (m/s)' % (a400))
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

	# plt.figure(dpi=500)
	# plt.scatter(a[1], P[1] * 10 ** (-6), label='CJ State', color='red')
	# plt.scatter(a[nfinal], P[nfinal] * 10 ** (-6), label='Plateau Condition', color='green')
	# plt.plot(a[1:nfinal], P[1:nfinal] * 10 ** (-6), label='Isentrope')
	# # plt.plot(ae_array[1:nfinal], 19451 * ae_array[1:nfinal] * 10 ** (-6), color='red')
	# # plt.title('Isentropic expansion after detonation, CJ speed = %.1f m/s' % (cj_speed))
	# plt.gca().xaxis.set_major_locator(MaxNLocator(6))
	# plt.xlabel(r'Equilibrium Sound Speed (m/s)')
	# plt.ylabel('Pressure (MPa)')
	# plt.ylim(bottom=0, top=25)
	# plt.gca().invert_xaxis()
	# plt.legend()
	# # plt.show()
	# plt.savefig('Mach8Isentrope.png', bbox_inches='tight')

	return P400, u400, a400, V400, T400


def HYPULSE_driver_equivalentconditions(P4, T4, q4, P100, T100, q100):
	"""
	Models a HyPulse detonation driver experiment to find equivalent conditions using a tuned equation of state.

	This function creates two Cantera Solution objects for the initial and final states, computes 
	thermodynamic properties like sound speed, density, and heat capacity ratio, and uses them to find 
	equivalent conditions of pressure and temperature as per Lu's suggestion.

	Parameters
	-----------
	P4, P100 : float 
		Initial pressure (Pa) at states 4 and 100.

	T4, T100 : float 
		Initial temperature (K) at states 4 and 100.

	q4, q100 : str 
		Initial mole fractions at states 4 and 100.

	Returns
	--------
	Pe, Te : float 
		Equivalent pressure (Pa) and temperature (K) in the initial state.

	pe2, Te2 : float 
		Equivalent pressure (Pa) and temperature (K) in the final CJ state.

	Raises
	------
	CanteraError
		If any input parameter is not valid, or if Cantera encounters an error during the simulation.

	Note
	----
	This function prints out many intermediate values and stages of the experiment.

	"""

	mech = 'sandiego20161214_H2only.cti'

	# State 400 gas object
	gas4 = ct.Solution(mech)
	gas4.TPX = T4, P4, q4
	a4 = soundspeed_fr(gas4)
	D4 = gas4.density
	gamma4 = a4*a4*D4/P4

	# State 100 gas object
	gas100 = ct.Solution(mech)
	gas = ct.Solution(mech)
	gas100.TPX = T100, P100, q100
	a100_fr = soundspeed_fr(gas100)
	D100 = gas100.density
	gamma100_fr =  a100_fr*a100_fr*D100/P100

	print('\nDetonation Driver Fill State:')
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
	gamma2_eq =  a2_eq*a2_eq*D2/P2;

	print('CJ speed = %.2f (m/s)' % (cj_speed))
	print ('CJ State')
	print(' Pressure %.2f (MPa) ' % (P2*10**(-6)) )
	print(' Temperature %.2f (K) ' % (T2) )

	print(' u2 (frozen) %.2f (m/s)' % (u2) )
	print(' a2 (equilibrium) %.2f (m/s)' % (a2_eq) )
	print(' gamma2 (equilibrium) %.4f '% (gamma2_eq) )

	pe2 = P2 * (1 + (gamma2_eq - 1) / 2 * (u2 / a2_eq)) ** (2 * gamma2_eq / (gamma2_eq - 1))
	ae2 = a2_eq * (1 + (gamma2_eq - 1) / 2 * (u2 / a2_eq))
	print('State 2 Equivalent')
	print(' pe = %.2f (MPa)' % (pe2 * 10 ** (-6)))
	print(' ae = %.2f (m/s)' % (ae2))

	def calculate_a_eqn(T2):
		x = gas.X
		gas.TPX = T2, pe2, x
		gas.equilibrate('TP')
		a_guess = soundspeed_eq(gas)
		return a_guess

	def T2e_eqn(T2):
		# print('-' * 60)
		# print(f"guessed T4 = {T4} K")
		# print(f"related guessed a4 = {calculate_a4_eqn(T4)} m/s")
		# print(f"actual a4 - guessed a4 = {ae - calculate_a4_eqn(T4)}")
		return ae2 - calculate_a_eqn(T2)

	Te2 = newton(T2e_eqn, 3000)
	print(f"Final T2e = {Te2:.2f} K")

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
	pe_array = np.zeros(npoints,float)
	ae_array = np.zeros(npoints,float)
	V[1] = V2
	P[1] = P2
	D[1] = D2
	a[1] = a2_eq
	u[1] = u2
	T[1] = T2
	pe_array[1] = pe2
	ae_array[1] = ae2
	print('\nGenerating points on isentrope and computing Taylor wave velocity')

	i = 1

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
		gamma_i = a[i]**2 * D[i] / P[i];

		# Calculating cold-driver-equivalent numbers from Lu:
		pe_array[i] = P[i] * (1 + (gamma_i - 1) / 2 * (u[i] / a[i])) ** (2 * gamma_i / (gamma_i - 1))
		ae_array[i] = a[i] * (1 + (gamma_i - 1) / 2 * (u[i] / a[i]))

		# expanding state 4 to state 400 conditions:
		Re_tol = 1e-3
		P3 = P4*(1 - 0.5*(gamma4 - 1)*(u[i]/a4))**(2*gamma4/(gamma4 -1))
		if abs((P3 - P[i])/P3) < Re_tol:
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
	print(' Velocity %.2f (m/s)' % (u400))
	
	# evaluate final state 400 to get sound speeds and effective gammas
	x = gas.X
	gas.SVX = S2, V400, x
	gas.equilibrate('SV')
	gamma400_eq =  a400**2/(V400*P400);
	print(' Sound speed (equilibrium) %.2f (m/s)' % (a400))
	print(' Gamma (equilibrium) %.4f' % (gamma400_eq))

	# State 3:
	u3 = u400
	T3 = T4*(P3/P4)**((gamma4-1)/gamma4)
	print('State 3:')
	print(' Pressure %.2f (MPa)' % (P3*10**(-6)))
	print(' Velocity %.2f (m/s)' % u3)
	print(' Temperature %.2f (K)' % (T3))

	# Calculating cold-driver-equivalent numbers from Lu:
	pe = P400*(1 + (gamma400_eq - 1)/2*(u400/a400))**(2*gamma400_eq/(gamma400_eq - 1))
	ae = a400*(1 + (gamma400_eq - 1)/2*(u400/a400))
	print('State 400 Equivalent')
	print(' pe = %.2f (MPa)' % (pe*10**(-6)))
	print(' ae = %.2f (m/s)' % (ae) )

	def calculate_a4_eqn(T4):
		x = gas.X
		gas.TPX = T4, pe, x
		gas.equilibrate('TP')
		a4_guess = soundspeed_eq(gas)
		return a4_guess

	def T4_eqn(T4):
		# print('-' * 60)
		# print(f"guessed T4 = {T4} K")
		# print(f"related guessed a4 = {calculate_a4_eqn(T4)} m/s")
		# print(f"actual a4 - guessed a4 = {ae - calculate_a4_eqn(T4)}")
		return ae - calculate_a4_eqn(T4)

	Te = newton(T4_eqn, 3000)
	print(f"Final T4 = {Te:.2f} K")

	# plt.figure(dpi=500)
	# plt.scatter(ae_array[1], pe_array[1] * 10 ** (-6), label='CJ State', color='red')
	# plt.scatter(ae_array[nfinal], pe_array[nfinal] * 10 ** (-6), label='Plateau Condition', color='green')
	# plt.plot(ae_array[1:nfinal], pe_array[1:nfinal] * 10 ** (-6), label='Isentrope')
	# plt.plot(ae_array[1:nfinal], 19451 * ae_array[1:nfinal] * 10 ** (-6), color='red', label='Experimental Locus')
	# # plt.title('Isentropic expansion after detonation, CJ speed = %.1f m/s' % (cj_speed))
	# plt.xlabel(r'Equilibrium Sound Speed (m/s)')
	# plt.ylabel('Pressure (MPa)')
	# plt.gca().invert_xaxis()
	# plt.legend()
	# # plt.show()
	# plt.savefig('Shot68Isentrope.png', bbox_inches='tight')

	# k = 1
	# fn = 'cj_isentrope.txt'
	# outfile = open(fn, 'w')
	# outfile.write('# CJ state isentrope\n')
	# outfile.write('# INITIAL CONDITIONS\n')
	# outfile.write('# TEMPERATURE (K) %.1f\n' % T2)
	# outfile.write('# PRESSURE (ATM) %.1f\n' % (P2 * 10**(-6)))
	# outfile.write('# DENSITY (KG/M^3) %.4f\n' % D2)
	# outfile.write('# SPECIES MOLE FRACTIONS: ' + q2 + '\n')
	# outfile.write('# MECHANISM: ' + mech + '\n')
	# outfile.write('# CJ speed (M/S) %.2f\n\n' % cj_speed)
	# outfile.write('# THE OUTPUT DATA COLUMNS ARE:\n')
	# outfile.write(
	# 	'Variables = "Specific Volume (m3/kg)", "Temperature (K)", "Pressure (Pa)", "sound speed (eq)(m/s)", "fluid speed (m/s)"\n')
	# while k <= nfinal:
	# 	outfile.write('%.4E \t %.1f \t %.4E \t %.1f \t %.1f\n' % (V[k], T[k], pe_array[k], ae_array[k], u[k]))
	# 	k = k + 1
	# outfile.close()
	
	return Pe, Te, pe2, Te2


def HYPULSE_driver_exp_tuned(P4, T4, q4, P100, T100, q100):
	"""
    Models a HyPulse detonation driver experiment using the tuned parameters from Heffernan, L., James, C.M. and Jewell, J.S. (2024). HYPULSE Modelling Using Quasi-0D Analytical Framework. AIAA AVIATION FORUM AND ASCEND 2024. doi:https://doi.org/10.2514/6.2024-4566.

    This function creates two Cantera Solution objects for the initial and final states and uses them to compute
    various thermodynamic properties while comparing the model's predictions with experimental observations.

    Parameters
    -----------
    P4, P100 : float
        Initial pressure (Pa) at states 4 and 100.

    T4, T100 : float
        Initial temperature (K) at states 4 and 100.

    q4, q100 : str
        Initial mole fractions at states 4 and 100 as a string.

    Returns
    --------
    P400 : float
        Pressure at the final state (Pa).

    u400 : float
        Velocity at the final state (m/s).

    a400 : float
        Sound speed at the final state (m/s).

    V400 : float
        Specific volume at the final state (m^3/kg).

    T400 : float
        Temperature at the final state (K).

    Raises
    ------
    CanteraError
        If any input parameter is not valid, or if Cantera encounters an error during the simulation.

    Note
    ----
    This function prints out many intermediate values and stages of the experiment.

    """
	mech = 'sandiego20161214_H2only.cti'

	# State 400 gas object
	gas4 = ct.Solution(mech)
	gas4.TPX = T4, P4, q4
	a4 = soundspeed_fr(gas4)
	D4 = gas4.density
	gamma4 = a4 * a4 * D4 / P4

	# State 100 gas object
	gas100 = ct.Solution(mech)
	gas100.TPX = T100, P100, q100
	a100_fr = soundspeed_fr(gas100)
	D100 = gas100.density
	gamma100_fr = a100_fr * a100_fr * D100 / P100

	print('\nDetonation Driver Fill State:')
	print(' Pressure %.2f (MPa) ' % (P100 * 10 ** (-6)))
	print(' Temperature %.2f (K) ' % (T100))

	print('Cold Driver Fill State:')
	print(' Composition: He: 100 ')
	print(' Pressure %.2f (MPa) ' % (P4 * 10 ** (-6)))
	print(' Temperature %.2f (K) ' % (T4))
	print(' a4 (frozen) %.2f (m/s)' % (a4))
	print(' gamma4 (frozen) %.4f ' % (gamma4))

	print('\nComputing CJ state and isentrope for ' + q100 + ' using ' + mech)

	# compute CJ speed
	cj_speed = CJspeed(P100, T100, q100, mech)

	# compute equilibrium CJ state
	gas = PostShock_eq(cj_speed, P100, T100, q100, mech)

	T2 = gas.T
	P2 = gas.P
	D2 = gas.density
	V2 = 1. / D2
	S2 = gas.entropy_mass
	w2 = D100 * cj_speed / D2
	u2 = cj_speed - w2
	a2_eq = soundspeed_eq(gas)
	a2_fr = soundspeed_fr(gas)
	gamma2_eq = a2_eq * a2_eq * D2 / P2

	print('CJ speed = %.2f (m/s)' % (cj_speed))
	print('CJ State')
	print(' Pressure %.2f (MPa) ' % (P2 * 10 ** (-6)))
	print(' Temperature %.2f (K) ' % (T2))
	print(' u2 (frozen) %.2f (m/s)' % (u2))
	print(' a2 (equilibrium) %.2f (m/s)' % (a2_eq))
	print(' gamma2 (equilibrium) %.4f ' % (gamma2_eq))

	# compute the effective value of q based on two-gamma model
	M1 = cj_speed / a100_fr;
	eparam = a100_fr ** 2 * (M1 ** (-2) * (gamma2_eq / gamma100_fr) ** 2 \
							 * (1 + gamma100_fr * M1 ** 2) ** 2 / (2 * (gamma2_eq ** 2 - 1)) - 1 / (
										 gamma100_fr - 1) - M1 ** 2 / 2)
	print('Detonation CJ Mach number %.2f' % (M1))
	print('2-gamma energy parameter q %.3f (J/kg)' % (eparam))

	# Find points on the isentrope connected to the CJ state, evaluate velocity
	# in Taylor wave using trapezoidal rule to evaluate the Riemann function
	npoints = 50000
	vv = V2
	V = np.zeros(npoints, float)
	P = np.zeros(npoints, float)
	D = np.zeros(npoints, float)
	a = np.zeros(npoints, float)
	u = np.zeros(npoints, float)
	T = np.zeros(npoints, float)
	V[1] = V2
	P[1] = P2
	D[1] = D2
	a[1] = a2_eq
	u[1] = u2
	T[1] = T2
	print('\nGenerating points on isentrope and computing Taylor wave velocity')

	i = 1
	while u[i] > 0 and i < npoints:
		i = i + 1
		vv = vv * 1.0001
		x = gas.X
		gas.SVX = S2, vv, x
		gas.equilibrate('SV')
		P[i] = gas.P
		D[i] = gas.density
		V[i] = 1 / D[i]
		T[i] = gas.T
		a[i] = soundspeed_eq(gas)
		u[i] = u[i - 1] + 0.5 * (P[i] - P[i - 1]) * (1. / (D[i] * a[i]) + 1. / (D[i - 1] * a[i - 1]))

		# expanding state 4 to state 400 conditions:
		Re_tol = 1e-3
		P3 = P4 * (1 - 0.5 * (gamma4 - 1) * (u[i] / a4)) ** (2 * gamma4 / (gamma4 - 1))
		if abs((P3 - P[i]) / P3) < Re_tol:
			break

	nfinal = i
	u400 = u[nfinal]
	a400 = a[nfinal]
	V400 = V[nfinal]
	T400 = T[nfinal]
	P400 = 19451 * a400

	# State 400:
	print('State 400 ')
	print(' Pressure %.2f (MPa)' % (P400 * 10 ** (-6)))
	print(' Temperature %.2f (K)' % (T400))
	print(' Velocity %.2f (m/s)' % (u400))

	# evaluate final state 400 to get sound speeds and effective gammas
	x = gas.X
	gas.SVX = S2, V400, x
	gas.equilibrate('SV')
	gamma400_eq = a400 ** 2 / (V400 * P400);

	print(' Sound speed (equilibrium) %.2f (m/s)' % (a400))
	print(' Gamma (equilibrium) %.4f' % (gamma400_eq))

	# State 3:
	u3 = u400
	T3 = T4 * (P3 / P4) ** ((gamma4 - 1) / gamma4)
	print('State 3:')
	print(' Pressure %.2f (MPa)' % (P3 * 10 ** (-6)))
	print(' Velocity %.2f (m/s)' % u3)
	print(' Temperature %.2f (K)' % (T3))

	return P400, u400, a400, V400, T400


		

