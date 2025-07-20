import sys
import os

# Add the directory containing Binding_calculator.py to the Python path
sys.path.append(os.path.abspath('/home/openkim/test-drivers/BindingEnergyvsHydrostaticPressure/test_driver'))

# For Mo, the potential used is: "SNAP_LiHuChen_2018_NiMo__MO_468686727341_000"
# For Al, the potential used is: MEAM_LAMMPS_JeongParkDo_2018_PdAl__MO_616482358807_002
# for Fe: the potnetial used is: MEAM_LAMMPS_KoJimLee_2012_FeP__MO_179420363944_002
from test_driver.test_driver import TestDriver
kim_model_name = "MEAM_LAMMPS_KoJimLee_2012_FeP__MO_179420363944_002"
test_driver = TestDriver(kim_model_name)

############################################################################
 # Fe and Mo are BCC; Al and Pd are FCC
#from ase.lattice.cubic import BaseCenteredCubic
from ase.build import bulk
a0 = 4.0
#atoms = BaseCenteredCubic("Fe", latticeconstant=a0)
atoms = bulk("Fe", crystalstructure="bcc", a=a0)

#pressure is in units of ev/Angstrom^3
print ('\nRUNNING TEST DRIVER ON WURTZITE ATOMS OBJECT\n')
test_driver(atoms,optimize=True,max_pressure_scale=6e-2,num_steps=30,temperature_K=0)
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