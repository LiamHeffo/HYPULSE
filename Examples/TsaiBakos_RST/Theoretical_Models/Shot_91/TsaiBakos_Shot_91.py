"""
TsaiBakos_Shot_91.py

Testing more HYPULSE GASL Shots

$ python3 TsaiBakos_Shot_91.py

Liam Heffernan 14/01/24

"""

#----------------------------------------------------------------------------------------------
# the code that runs PITOT3

from pitot3 import run_pitot3
from HYPULSE_driver import HYPULSE_driver


driver_gas_model = 'CEAGas-preset'
driver_fill_gas_name = '30-h2-15-o2-55-ar-room-temperature-only-gas-model'

facility_name = 'HYPULSE_rst_Mach_8_to_10_nozzle_area_ratio_175'

		# COLD DRIVER:
q4 = 'HE:100'
p4 = 13.79e6 # Pa,
T4 = 298.15 # K
p4_p100 = 32.89

		# DETONATION DRIVER:
q100 = 'O2:15 H2:30 Ar:55'
p100 = p4/p4_p100 # Pa, this is the detonation driver pressure
T100 = 298.15 # K

p400, u400, a400, v400, T400 = HYPULSE_driver(p4, T4, q4, p100, T100, q100)

p100_p1 = 20.93
p1 = p100/p100_p1 # Pa

# I'm using the state 400 values, but in PITOT3 the driver is state 4...
driver_dict = {'driver_condition_name' :'HYPULSE test', 'driver_condition_type' : 'empirical',
               'driver_gas_model' : driver_gas_model, 'driver_fill_gas_name' : driver_fill_gas_name,
               'p4' : p400, 'T4' : T400, 'v4': u400, 'M_throat' : 0.0}

config_dict = {'mode':'fully_theoretical','output_filename':'TsaiBakos_Shot_91',
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

# it is useful to show how to get a tube fill state, shocked state and unsteadily expanded state
# so I'll get the shock tube fill pressure (state1), shock tube post-shock state (state2),
# and the unsteadily expanded state in the acceleration tube (state 7)

state1_facility_state = shock_tube_object.get_fill_state()
state1_gas_state = state1_facility_state.get_gas_state()
state1_fill_pressure = state1_gas_state.p
state1_fill_temperature = state1_gas_state.T

print('-'*60)
print(f"Shock tube fill pressure (p1) is {state1_fill_pressure:.2f} Pa")
print(f"Shock tube fill temperature (T1) is {state1_fill_temperature:.2f} K")

state2_facility_state = shock_tube_object.get_shocked_state()
state2_gas_state = state2_facility_state.get_gas_state()
state2_pressure = state2_gas_state.p
state2_temperature = state2_gas_state.T
state2_composition, state2_composition_units = state2_facility_state.get_reduced_composition_with_units()

print('-'*60)
print(f"Shock tube post-shock pressure (p2) is {state2_pressure:.2f} Pa")
print(f"Shock tube post-shock temperature (T2) is {state2_temperature:.2f} K")
print(f"Shock tube post-shock composition (by {state2_composition_units}):")
print(state2_composition)

state10f_facility_state = test_section_object.get_post_normal_shock_ideal_gas_state()
state10f_gas_state = state10f_facility_state.get_gas_state()
state10f_pressure = state10f_gas_state.p
state10f_temperature = state10f_gas_state.T

print('-'*60)
print(f"Test section frozen post-normal shock pressure (p10f) is {state10f_pressure:.2f} Pa")
print(f"Test section frozen post-normal shock temperature (T10f) is {state10f_temperature:.2f} K")

state10e_facility_state = test_section_object.get_post_normal_shock_state()
state10e_gas_state = state10e_facility_state.get_gas_state()
state10e_pressure = state10e_gas_state.p
state10e_temperature = state10e_gas_state.T
state10e_composition, state10e_composition_units = state10e_facility_state.get_reduced_composition_with_units()

print('-'*60)
print(f"Test section equilibrium post-normal shock pressure (p10e) is {state10e_pressure:.2f} Pa")
print(f"Test section equilibrium post-normal shock temperature (T10e) is {state10e_temperature:.2f} K")
print(f"Test section equilibrium post-normal shock composition (by {state10e_composition_units}):")
print(state10e_composition)
















