"""
liam_HYPULSE_test_5.py
$ python3 liam_HYPULSE_test_5.py
Liam Heffernan (lheffern@purdue.edu) - 26/09/23
Chris James (c.james4@uq.edu.au) - 26/09/23

"""

#----------------------------------------------------------------------------------------------

# the code that runs PITOT3

from HYPULSE_driver import HYPULSE_driver
from pitot3 import run_pitot3


# conditions from Table 3 in Lu and Marren HYPULSE chapter (except temp whcih we are messing with) M15 condition
# Temperature is to match the measured primary shock speed of 2240 m/s

# FACILITY CONDITIONS:

	# DRIVERS:
driver_gas_model = 'CEAGas-preset'
driver_fill_gas_name = '30-h2-15-o2-55-he'

		# COLD DRIVER:
q4 = 'N2:78.07 O2:20.99 Ar:0.94'
p4 = 18.483e6 # Pa, 
T4 = 300 # K

		# DETONATION DRIVER:
q100 = 'O2:33 H2:66 HE:0'
p100 = 0.20301e6 # Pa, this is the detonation driver pressure
T100 = 300 # K

p400, u400, a400, v400, T400 = HYPULSE_driver(p4, T4, q4, p100, T100, q100)

