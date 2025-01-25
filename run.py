import sys
import os

# Add the directory containing Binding_calculator.py to the Python path
sys.path.append(os.path.abspath('/home/openkim/test-drivers/BindingEnergyvsHydrostaticPressure/test_driver'))

from test_driver.test_driver import TestDriver
kim_model_name = "EAM_Dynamo_ErcolessiAdams_1994_Al__MO_123629422045_005"
test_driver = TestDriver(kim_model_name)

############################################################################
 
from ase.lattice.cubic import FaceCenteredCubic
a0 = 4.0
atoms = FaceCenteredCubic("Al", latticeconstant=a0)

print ('\nRUNNING TEST DRIVER ON WURTZITE ATOMS OBJECT\n')
test_driver(atoms,optimize=True,max_pressure_scale=1e-2,temperature_K=0)
############################################################################

print('\n--------------------------------------')

print('Lattice constant a: %f %s'%(
    test_driver.property_instances[0]['a']['source-value'],
    test_driver.property_instances[0]['a']['source-unit']
))

print('--------------------------------------\n')
###############################################################################
# You can also dump the instances to a file, by default ``output/results.edn``.
test_driver.write_property_instances_to_file()

"""
from kim_tools import query_crystal_genome_structures
list_of_queried_structures = query_crystal_genome_structures(
    kim_model_name=kim_model_name,
    stoichiometric_species=['Zn','S'],
    prototype_label='AB_hP4_186_b_b')

for queried_structure in list_of_queried_structures:
    print ('\nRUNNING TEST DRIVER ON QUERIED STRUCTURE\n')
    test_driver(**queried_structure,max_pressure_scale=0.01)


###############################################################################
# Remember that this will append the results to ``test_driver.property_instances``.
"""