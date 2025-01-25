"""
==========================================
kim-tools package. See https://kim-tools.readthedocs.io for more information.
"""
from kim_tools import CrystalGenomeTestDriver
from kim_tools import get_stoich_reduced_list_from_prototype
from ase.calculators.kim import KIM
from .helper_functions import binding_potential_energy

class TestDriver(CrystalGenomeTestDriver):
    def _calculate(self, max_pressure_scale: float = 1e-2, num_steps: int = 10, **kwargs):
        """
        Computes the energy vs. pressure curve for isotropic expansion and compression. 

        Args:
            max_pressure_scale:
                Maximum fractional change in pressure to investigate
            num_steps:
                Number of steps to take in each direction
        """
        
        # The base class provides self.atoms, an ase.Atoms object representing the initial configuration of the crystal.
        # Use this configuration to evaluate the material property of interest
        original_cell = self.atoms.get_cell()
        num_atoms = len(self.atoms)
        
        #Create temporary atoms object for negative pressure
        atoms_tmp = self.atoms.copy()
        atoms_tmp.calc = self._calc
        # Besides temperature, stress, and atoms, you may wish to access other attributes of the base class for information about 
        # the material, such as its symmetry-reduced AFLOW prototype label. Here we use it to get information about the stoichiometry of the crystal.
        # See the API documentation for CrystalGenomeTestDriver for more information:
        # https://kim-tools.readthedocs.io/en/latest/kim_tools.html#kim_tools.CrystalGenomeTestDriver
        num_atoms_in_formula = sum(get_stoich_reduced_list_from_prototype(self.prototype_label))
        original_step_size = max_pressure_scale/num_steps
        disclaimer = None

        print('\nPerforming energy scan...\n')
        #Calculation of binding energy for negative pressure
        step_size = -original_step_size
        neg_results = binding_potential_energy(num_steps, step_size, num_atoms, atoms_tmp, num_atoms_in_formula)
        volume_per_atom_1 = neg_results[0]
        volume_per_formula_1 = neg_results[1]
        binding_potential_energy_per_atom_1 = neg_results[2]
        binding_potential_energy_per_formula_1 = neg_results[3]
        pressure_array_neg = neg_results[4]
        
        #-------------------------------------------------------------------------------------------------------
        #redefine the atom object for positive pressure calculation
        atoms_tmp_1 = self.atoms.copy()
        atoms_tmp_1.calc = self._calc
        #Calculation of binding energy for positive pressure
        step_size = original_step_size
        pos_results = binding_potential_energy(num_steps, step_size, num_atoms, atoms_tmp_1, num_atoms_in_formula)
        volume_per_atom_2 = pos_results[0]
        volume_per_formula_2 = pos_results[1]
        binding_potential_energy_per_atom_2 = pos_results[2]
        binding_potential_energy_per_formula_2 = pos_results[3]
        pressure_array_pos = pos_results[4]
        #-------------------------------------------------------------------------------------------------------
        #sort the -ve array in ascending order first and then remove the last element
        volume_per_atom_1.sort()
        volume_per_atom_1 = volume_per_atom_1[:-1]
        
        volume_per_formula_1.sort()
        volume_per_formula_1 = volume_per_formula_1[:-1]

        binding_potential_energy_per_atom_1.sort()
        binding_potential_energy_per_atom_1 = binding_potential_energy_per_atom_1[:-1]

        binding_potential_energy_per_formula_1.sort()
        binding_potential_energy_per_formula_1 = binding_potential_energy_per_formula_1[:-1]

        pressure_array_neg.sort()
        pressure_array_neg = pressure_array_neg[:-1]
        
        #concatenate arrays for -ve and +ve pressure
        volume_per_atom = volume_per_atom_1 + volume_per_atom_2
        volume_per_formula = volume_per_formula_1 + volume_per_formula_2
        binding_potential_energy_per_atom = binding_potential_energy_per_atom_1 + binding_potential_energy_per_atom_2
        binding_potential_energy_per_formula = binding_potential_energy_per_formula_1 + binding_potential_energy_per_formula_2
        pressure_array = pressure_array_neg + pressure_array_pos
        
        print('pressure_array:',pressure_array)
        # Now it is time to write the output in the format you created in your Property Definition. The base class provides utility methods
        # to facilitate this process.
        
        # This method initializes the Property Instance and adds the keys common to all Crystal Genome properties.
        # property_name can be the full "property-id" field in your Property Definition, or the "Property Name",
        # which is just the short name after the slash, as used here. You can also specify whether your property
        # includes stress and temperature (no by default), and have the option to specify a disclaimer.
        self._add_property_instance_and_common_crystal_genome_keys("energy-vs-pressure-isotropic-crystal",
                                                                   write_stress=False, write_temp=False, disclaimer=disclaimer)

        # This method adds additional fields to your property instance by specifying the key names you defined
        # in your property definition and providing units if necessary.
        self._add_key_to_current_property_instance("volume-per-atom",volume_per_atom,
                                                   units="angstrom^3")
        self._add_key_to_current_property_instance("volume-per-formula",volume_per_formula,
                                                   units="angstrom^3")

        # You may also provide a dictionary supplying uncertainty information. It is optional, and normally
        # would not be reported for a deterministic calculation like this, only one involving some statistics,
        # such as molecular dynamics or fitting. There are many possible ways to report uncertainty information,
        # detailed at https://openkim.org/doc/schema/properties-framework/
        uncertainty_info = {
            "source-std-uncert-value": [0]*len(binding_potential_energy_per_atom)
        }

        self._add_key_to_current_property_instance("binding-potential-energy-per-atom",binding_potential_energy_per_atom,
                                                   units="eV",uncertainty_info=uncertainty_info)
        self._add_key_to_current_property_instance("binding-potential-energy-per-formula",binding_potential_energy_per_formula,
                                                   units="eV",uncertainty_info=uncertainty_info)
        
        #Adding the pressure in the property definition
        self._add_key_to_current_property_instance("pressure", pressure_array, 
                                                   units="GPa")
        # If your Test Driver reports multiple Property Instances, repeat the process above for each one.

