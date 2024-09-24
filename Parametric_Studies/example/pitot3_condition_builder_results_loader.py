"""
Just a program to load the condition builder result and make some plots...

Chris James (c.james4@uq.edu.au) - 01/11/23

"""

import json
import matplotlib.pyplot as plt

json_result_filename = 'HYPULSE_RST_M8_condition_builder_final_result_dict_output.json'

split_driver_conditions = True

with open(json_result_filename, "r") as json_result_file:
    json_result_dict = json.load(json_result_file)

print(json_result_dict.keys())

# I got this from the old pitot 2 condition builder plot...
label_dictionary = {}
label_dictionary['Ht'] = 'Stagnation enthalpy (Ht MJ/kg)'
label_dictionary['Ue'] = 'Flight equivalent velocity (m/s)'
label_dictionary['rho7'] = 'Freestream density (rho7, kg/m^3)'
label_dictionary['rho8'] = 'Freestream density (rho8, kg/m^3)'
label_dictionary['p7'] = 'Freestream pressure (p7, Pa)'
label_dictionary['p8'] = 'Freestream pressure (p8, Pa)'
label_dictionary['T7'] = 'Freestream temperature (T7, K)'
label_dictionary['T8'] = 'Freestream temperature (T8, K)'
label_dictionary['M7'] = 'Freestream Mach number (M7)'
label_dictionary['M8'] = 'Freestream Mach number (M8)'
label_dictionary['V8'] = 'Freestream velocity (m/s)'
label_dictionary['pitot_p'] = 'Pitot pressure (Pa)'

x_variable = 'Ue'
y_variable = 'rho8'

target_x_value = 4500 # Pa
target_y_value = 0.0037 # MJ/kg

x_values = json_result_dict[x_variable]
y_values = json_result_dict[y_variable]

fig, ax = plt.subplots(figsize = (10,10))

if not split_driver_conditions:

    ax.scatter(x_values, y_values)

else:
    # we need to go through and split the driver conditions...

    driver_condition_values = json_result_dict['driver_condition']

    driver_condition_based_results_dict = {}

    for driver_condition, x_value, y_value in zip(driver_condition_values, x_values, y_values):

        if driver_condition not in driver_condition_based_results_dict:
            # add it!

            driver_condition_based_results_dict[driver_condition] = {}

            driver_condition_based_results_dict[driver_condition]['x_values'] = []
            driver_condition_based_results_dict[driver_condition]['y_values'] = []

        driver_condition_based_results_dict[driver_condition]['x_values'].append(x_value)
        driver_condition_based_results_dict[driver_condition]['y_values'].append(y_value)

    for driver_condition in driver_condition_based_results_dict:

        x_values = driver_condition_based_results_dict[driver_condition]['x_values']
        y_values = driver_condition_based_results_dict[driver_condition]['y_values']

        ax.scatter(x_values, y_values, label = driver_condition)

ax.set_xlabel(label_dictionary[x_variable])
ax.set_ylabel(label_dictionary[y_variable])

if target_x_value:
    ax.axvline(target_x_value, label = 'x target')
if target_y_value:
    ax.axhline(target_y_value, label = 'y target')


if split_driver_conditions:
    ax.legend(loc = 'best')

output_filename = f'envelope_{x_variable}_versus_{y_variable}'

plt.savefig(output_filename + '.pdf', bbox_inches='tight')
plt.savefig(output_filename + '.png', bbox_inches='tight')

plt.show()

# now let us try to refine things a bit...
# we will instead make a dictionary where things are sorted by sim and then we can pull out sims that match certain requirements

test_name_results_dict = {}

for i, test_name in enumerate(json_result_dict['test_name']):

    # we make a dictionary for this test

    current_test_results_dictionary = {}

    for variable in json_result_dict:
        current_test_results_dictionary[variable] = json_result_dict[variable][i]

    test_name_results_dict[test_name] = current_test_results_dictionary

refined_test_name_results_dict = {}

variables_to_check = ['Ue', 'rho8']

variables_to_check_setpoint_dict = {'Ue':1000, 'Ht': 2e-2}
variables_to_check_percentage_difference_dict = {'pitot_p':0.01, 'Ht': 0.01} # 1 pc for both...

for test_name in test_name_results_dict.keys():

    current_test_results_dictionary = test_name_results_dict[test_name]

    # assume we will and tehn remo
    add_test_to_refined_results = False

    for variable in variables_to_check:

        variable_setpoint = variables_to_check_setpoint_dict[variable]
        variable_minimum_value = variable_setpoint - variable_setpoint*variables_to_check_percentage_difference_dict[variable]
        variable_maximum_value = variable_setpoint + variable_setpoint*variables_to_check_percentage_difference_dict[variable]


        variable_value = current_test_results_dictionary[variable]

        if variable_value and variable_minimum_value <= variable_value <= variable_maximum_value:
            add_test_to_refined_results = True
        else:
            add_test_to_refined_results = False

    if add_test_to_refined_results:
        refined_test_name_results_dict[test_name] = current_test_results_dictionary

#print(len(refined_test_name_results_dict))

#print(refined_test_name_results_dict)

variables_to_print = ['vs1', 'p1', 'p5', 'M8', 'T8', 'driver_condition']

for test_name in refined_test_name_results_dict.keys():
    current_test_results_dictionary =  refined_test_name_results_dict[test_name]

    output_line = f'test_name = {test_name}'

    for variable in variables_to_print:

        if '_over_' in variable:
            # we need to divide two values...

            variable_list = variable.split('_over_')

            variable_1 = variable_list[0]
            variable_2 = variable_list[1]

            variable_value = current_test_results_dictionary[variable_1]/current_test_results_dictionary[variable_2]

        else:
            variable_value = current_test_results_dictionary[variable]

        if isinstance(variable_value, float):
            output_line += f', {variable} = {variable_value:.2f}'
        else:
            output_line += f', {variable} = {variable_value}'

    print(output_line)









