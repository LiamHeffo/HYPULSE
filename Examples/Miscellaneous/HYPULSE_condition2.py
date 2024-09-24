"""
Shock and Detonation Toolbox Demo Program

Calculates  points on the isentrope and Taylor-Zeldovich expansion 
behind a CJ detonation.  Prints output file and makes plots.
 
################################################################################
Theory, numerical methods and applications are described in the following report:

SDToolbox - Numerical Tools for Shock and Detonation Wave Modeling,
Explosion Dynamics Laboratory, Contributors: S. Browne, J. Ziegler,
N. Bitter, B. Schmidt, J. Lawson and J. E. Shepherd, GALCIT
Technical Report FM2018.001 Revised January 2021.

Please cite this report and the website if you use these routines. 

Please refer to LICENCE.txt or the above report for copyright and disclaimers.

http://shepherd.caltech.edu/EDL/PublicResources/sdt/

################################################################################ 
Updated September 2018
Tested with: 
    Python 3.5 and 3.6, Cantera 2.3 and 2.4
Under these operating systems:
    Windows 8.1, Windows 10, Linux (Debian 9)
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
plots = True

#set initial detoantion driver fill state
P1 = 34/30.1*10**6
T1= 295

# set cold driver state:
P4 = 34*10**6
T4 = 295
gamma4 = 5/3
R4 = 2077.1
a4 = np.sqrt(gamma4*R4*T4)

q= 'O2:15 H2:30 He:55'
mech = 'Mevel2017.cti'
gas1 = ct.Solution(mech)
gas = ct.Solution(mech)
gas1.TPX = T1, P1, q
a1_fr = soundspeed_fr(gas1)
D1 = gas1.density
gamma1_fr =  a1_fr*a1_fr*D1/P1;

print('\nDetonation Driver Fill State:');
print(' Composition ' + q)
print(' Pressure %.2f (MPa) ' % (P1*10**(-6)) )
print(' Temperature %.2f (K) ' % (T1) )
print(' Density %.3f (kg/m3) ' % (D1) )
print(' a1 (frozen) %.2f (m/s)' % (a1_fr) )
print(' gamma1 (frozen) %.4f '% (gamma1_fr) )

print('Cold Driver Fill State:')
print(' Composition: He: 100 ' )
print(' Pressure %.2f (MPa) ' % (P4*10**(-6)) )
print(' Temperature %.2f (K) ' % (T4) )
print(' a4 (frozen) %.2f (m/s)' % (a4) )
print(' gamma4 (frozen) %.4f '% (gamma4) )



print('\nComputing CJ state and isentrope for '+q+' using '+mech)
# compute CJ speed
cj_speed = CJspeed(P1, T1, q, mech);   

# compute equilibrium CJ state
gas = PostShock_eq(cj_speed, P1, T1, q, mech)

T2 = gas.T
P2 = gas.P
D2 = gas.density
V2 = 1./D2
S2  = gas.entropy_mass
w2 = D1*cj_speed/D2
u2 = cj_speed-w2
a2_eq = soundspeed_eq(gas)
a2_fr = soundspeed_fr(gas)
gamma2_fr =  a2_fr*a2_fr*D2/P2;
gamma2_eq =  a2_eq*a2_eq*D2/P2;

print('CJ speed = %.2f (m/s)' % (cj_speed))
print ('CJ State')
print(' Pressure %.2f (Pa) ' % (P2) )
print(' Temperature %.2f (K) ' % (T2) )
print(' Density %.3f (kg/m3) ' % (D2) )
print(' Entropy %.3f (J/kg-K) ' % (S2) )
print(' w2 (frozen) %.2f (m/s)' % (w2) )
print(' u2 (frozen) %.2f (m/s)' % (u2) )
print(' a2 (frozen) %.2f (m/s)' % (a2_fr) )
print(' a2 (equilibrium) %.2f (m/s)' % (a2_eq) )
print(' gamma2 (frozen) %.4f '% (gamma2_fr) )   
print(' gamma2 (equilibrium) %.4f '% (gamma2_eq) )   

#compute the effective value of q based on two-gamma model
M1 = cj_speed/a1_fr;
eparam = a1_fr**2*(M1**(-2)*(gamma2_eq/gamma1_fr)**2 \
                   *(1+gamma1_fr*M1**2)**2/(2*(gamma2_eq**2-1))-1/(gamma1_fr-1)-M1**2/2)
print('Detonation CJ Mach number %.2f' % (M1))
print('2-gamma energy parameter q %.3f (J/kg)' % (eparam))

# Find points on the isentrope connected to the CJ state, evaluate velocity
# in Taylor wave using trapezoidal rule to evaluate the Riemann function
npoints=10000
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
    vv = vv*1.001
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
    Re_tol = 10000 #in Pa remember
    P3 = P4*(1 - (gamma4 - 1)/2*( u[i] / a4 ))**(2*gamma4/(gamma4 - 1))
    if abs(P3 - P[i]) < Re_tol:
        has_converged = True	
        break
print('number of points = ', i)
#if has_converged == True:
    #print('this worked')
#else:
 #   print("this hasn't worked")

# pulling last value as state 400 conditions 
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
print(' Volume %.2f (m3/kg)' % (V400))
print(' Velocity %.2f (m/s)' % (u400))
 
# evaluate final state 400 to get sound speeds and effective gammas
x = gas.X
gas.SVX = S2, V400, x
gas.equilibrate('SV')
a400_fr = soundspeed_fr(gas);
gamma400_fr =  a400_fr**2/(P400*V400);
gamma400_eq =  a400**2/(V400*P400);
print(' Sound speed (frozen) %.2f (m/s)' % (a400_fr))
print(' Sound speed (equilibrium) %.2f (m/s)' % (a400))
print(' Gamma (frozen) %.4f' % (gamma400_fr))
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


# output file
import datetime
d = datetime.date.today()
k = 1
fn = 'cj_isentrope.txt' 
outfile = open(fn, 'w')
outfile.write('# CJ state isentrope\n')
outfile.write('# CALCULATION RUN ON %s\n\n' % d)
outfile.write('# INITIAL CONDITIONS\n')
outfile.write('# TEMPERATURE (K) %.1f\n' % T1)
outfile.write('# PRESSURE (ATM) %.1f\n' % (P1/ct.one_atm))
outfile.write('# DENSITY (KG/M^3) %.4f\n' % D1)
outfile.write('# SPECIES MOLE FRACTIONS: ' + q + '\n')
outfile.write('# MECHANISM: ' + mech + '\n')
outfile.write('# CJ speed (M/S) %.2f\n\n' % cj_speed)
outfile.write('# THE OUTPUT DATA COLUMNS ARE:\n')
outfile.write('Variables = "Specific Volume (m3/kg)", "Temperature (K)", "Pressure (Pa)", "sound speed (eq)(m/s)", "fluid speed (m/s)"\n')
while k <= nfinal:
    outfile.write('%.4E \t %.1f \t %.4E \t %.1f \t %.1f\n' % (V[k], T[k], P[k], a[k], u[k]))
    k = k + 1
outfile.close()


# plotting results
if plots:
    plt.figure(1)
    plt.plot(V[1:nfinal],P[1:nfinal]/ct.one_atm)
    plt.title('Isentropic expansion after detonation, CJ speed = %.1f m/s' % (cj_speed))
    plt.xlabel(r'Volume (m$^3$/kg)')
    plt.ylabel('Pressure (atm)')
    
    plt.figure(2)
    plt.plot(V[1:nfinal],a[1:nfinal])
    plt.title('Isentropic expansion after detonation,  CJ speed = %.1f m/s' % (cj_speed))
    plt.xlabel(r'Volume (m$^3$/kg)')
    plt.ylabel('Equilibrium sound speed (m/s)')

    plt.figure(3)
    plt.plot(V[1:nfinal],u[1:nfinal])
    plt.title('Isentropic expansion after detonation,  CJ speed = %.1f m/s' % (cj_speed))
    plt.xlabel(r'Volume (m$^3$/kg)')
    plt.ylabel('Fluid speed (m/s)')

    plt.figure(4)
    plt.plot(V[1:nfinal],T[1:nfinal])
    plt.title('Isentropic expansion after detonation,  CJ speed = %.1f m/s' % (cj_speed))
    plt.xlabel(r'Volume (m$^3$/kg)')
    plt.ylabel('Temperature (K)')
    plt.show()
