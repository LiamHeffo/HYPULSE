"""
TsaiBakos_RST_M10.py

Modelling the M10 Nozzle AR 175 RST condition from 'Mach 7-21 Flight Simulation in the HYPULSE Shock Tunnel'. C.-Y. Tsai and R.J.
Bakos, 2002.

$ python3 TsaiBakos_RST_M10.py

Liam Heffernan (liam.heffernan@uqconnect.edu.au)

"""

#----------------------------------------------------------------------------------------------
# the code that runs PITOT3

from pitot3 import run_pitot3
from HYPULSE_driver import HYPULSE_driver, HYPULSE_driver_equivalentconditions


driver_gas_model = 'CEAGas-preset'
driver_fill_gas_name = '467-h2-233-o2-30-ar'

facility_name = 'HYPULSE_rst_Mach_8_to_10_nozzle_area_ratio_175'

# conditions from Table 3 in Lu and Marren HYPULSE chapter (except temp whcih we are messing with) M15 condition
# Temperature is to match the measured primary shock speed of 2240 m/s

		# COLD DRIVER:
q4 = 'HE:100'
p4 = 83.3e6 # Pa, 
T4 = 300 # K
p4_p100 = 51.2

		# DETONATION DRIVER:
q100 = 'O2:23.3 H2:46.7 Ar:30'
p100 = p4/p4_p100 # Pa, this is the detonation driver pressure
T100 = 300 # K

# p400, u400, a400, v400, T400 = HYPULSE_driver(p4, T4, q4, p100, T100, q100)
p400, T400, p4CJ, T4CJ = HYPULSE_driver_equivalentconditions(p4, T4, q4, p100, T100, q100)  # Equivalent Conditions

p100_p1 = 23.1
p1 = p100/p100_p1 # Pa

### I'm using the state 400 values, but in PITOT3 the driver is state 4...
driver_dict = {'driver_condition_name' :'HYPULSE test', 'driver_condition_type' : 'empirical',
               'driver_gas_model' : driver_gas_model, 'driver_fill_gas_name' : driver_fill_gas_name,
               'p4' : p400, 'T4' : T400, 'v4': 0, 'M_throat' : 0.0}

config_dict = {'mode':'fully_theoretical','output_filename':'TsaiBakos_Shot_82',
               'facility':facility_name,
               'driver_condition':'custom_from_dict', 'driver_dict':driver_dict,
               'test_gas_gas_model':'CEAGas', 'test_gas_name':'n2-o2-with-ions',
               'p1':p1, 'p5':None,
                'vs1_tolerance' : 1.0e-4}

config_data, gas_path, object_dict, states_dict = run_pitot3(config_dict = config_dict)

#----------------------------------------------------------------------------------------------
# the code that lets us pull some useful things out of the output of the code

# the object dictionary has all of the states we need
# 'driver' is the driver
# 'shock_tube' is shock tube
# 'acceleration_tube' is the acceleration tube
# 'nozzle' is the nozzle
# 'test_section' is the test section
# there are also states for the diaphragms etc. as well

driver_object = object_dict['driver']
shock_tube_object = object_dict['shock_tube']
test_section_object = object_dict['test_section']

# lets us start by getting the shock speeds:

vs1 = shock_tube_object.get_shock_speed()

print('-'*60)
print("The calculated shock speeds are:")
print(f"vs1 = {vs1:.2f} m/s")

# I'm using the state 400 values, but in PITOT3 the driver is state 4...
driver_dict = {'driver_condition_name' :'HYPULSE test', 'driver_condition_type' : 'empirical',
               'driver_gas_model' : driver_gas_model, 'driver_fill_gas_name' : driver_fill_gas_name,
               'p4' : p4CJ, 'T4' : T4CJ, 'v4': 0, 'M_throat' : 0.0}

config_dict = {'mode':'fully_theoretical','output_filename':'TsaiBakos_Shot_70',
               'facility':facility_name,
               'driver_condition':'custom_from_dict', 'driver_dict':driver_dict,
               'test_gas_gas_model':'CEAGas', 'test_gas_name':'n2-o2-with-ions',
               'p1':p1, 'p5':0,
                'vs1_tolerance' : 1.0e-4}

config_data, gas_path, object_dict, states_dict = run_pitot3(config_dict = config_dict)

#----------------------------------------------------------------------------------------------
# the code that lets us pull some useful things out of the output of the code

# the object dictionary has all of the states we need
# 'driver' is the driver
# 'shock_tube' is shock tube
# 'acceleration_tube' is the acceleration tube
# 'nozzle' is the nozzle
# 'test_section' is the test section
# there are also states for the diaphragms etc. as well

driver_object = object_dict['driver']
shock_tube_object = object_dict['shock_tube']
test_section_object = object_dict['test_section']

# lets us start by getting the shock speeds:

vs1 = shock_tube_object.get_shock_speed()

print('-'*60)
print("The calculated shock speeds are:")
print(f"vs1 = {vs1:.2f} m/s")

















