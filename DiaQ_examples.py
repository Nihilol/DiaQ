import DiaQ.diagnostics as dg
import numpy as np

# Read the .yml file that contains the gate names.
# This file should be in the same order as the eventual leakage matrix.
# If not, the output will be non-sensical.

name_dict, diagnostics_values, gds_coordinates, gds_file = dg.read_yaml("DiaQ")

# cut the inner structure out of a larger .gds file, and save it as a new .gds file.

dg.cut_gds(gds_coordinates, gds_file)

# Populate the gates, meaning that the gates are either named manually, or 
# read from the Inter-language files directory
# naming should be set to True if this is first time running the script. 
# It will prompt you to name the gates in the terminal.
# The gate is an object class defined in the DiaQ library. The source code can be seen there.
# Most notable attributes are the leakage, neighbours, name, and overlap. 

gates = dg.populate_gates(max_distance = diagnostics_values['max-distance'], naming = False)

# Initialise a leakage matrix. This would be the output of the QCodes

leakage_matrix_init = np.random.random((29, 29))

# Make it a Leakage_Matrix object, with a threshold of 0.9, and the names of the gates.
# The threshold is the value at which the value of the matrix is actually considered a leakage,
# and not just, say, noise.

leakage_matrix = dg.Leakage_Matrix(leakage_matrix_init, 0.9, name_dict)

# Convert the leakage matrix to a binary matrix, and perform the failure modes analysis

leakage_matrix.failure_modes(gates, show_plot = True)

# Run through each objects, and print the leakage and cause of each gate.

for key in gates:
    print("The leakage of", key, "is: ", gates[key].leakage)
