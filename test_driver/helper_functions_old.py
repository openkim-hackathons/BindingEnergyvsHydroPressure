'''
# function to calculate the binding energy and volume per atom to dump in the output 
def binding_potential_energy(no_of_steps, step_size, num_atoms, atoms_tmp, num_atoms_in_formula):
    #import minimize_wrapper from the kim tools
    from kim_tools import minimize_wrapper

    #initialize the arrays to store the volume and binding energy
    binding_potential_energy_per_atom = []
    binding_potential_energy_per_formula = []
    volume_per_atom = []
    volume_per_formula = []
    pressure_array = []

    #loop to calculate the binding energy and volume per atom
    for i in range(0, no_of_steps+2):
        if i == 0:
            pressure_scale = 0 # pressure scale
        elif i == 1:
            k = 0.001
            pressure_scale = k * step_size # pressure scale
        else:
            pressure_scale = step_size*(i-1)  # pressure scale
            #volume_scale = 1 + step_size*i
            #linear_scale = volume_scale ** (1/3)
            #self.atoms.set_cell(original_cell*linear_scale,scale_atoms=True)
            volume = atoms_tmp.get_volume()
            current_volume_per_atom = volume/num_atoms
            print('no. of atoms in the crystal structure: %d'%num_atoms)
            minimize_wrapper(atoms_tmp, variable_cell=True, fix_symmetry=True, flt_kwargs={'scalar_pressure':pressure_scale})
            problem_occurred = False
            try:
                # The self.atoms object comes pre-initialized with the calculator set to the interatomic model
                # the Test Driver was called with. If you need to access the name of the KIM model (for example,
                # if you are exporting the atomic configuration to run an external simulator like LAMMPS), it can
                # be accessed at self.kim_model_name
                potential_energy = atoms_tmp.get_potential_energy()                
                print('Volume: %5.5f Energy: %5.5f'%(volume,potential_energy))
            except:
                print('Failed to get energy at volume %f'%volume)
                problem_occurred = True
                disclaimer = "At least one of the requested deformations of the unit cell failed to compute a potential energy."                
            
            if not problem_occurred:
                current_binding_potential_energy_per_atom = potential_energy/num_atoms
                volume_per_atom.append(current_volume_per_atom)
                volume_per_formula.append(current_volume_per_atom*num_atoms_in_formula)
                binding_potential_energy_per_atom.append(current_binding_potential_energy_per_atom)
                binding_potential_energy_per_formula.append(current_binding_potential_energy_per_atom*num_atoms_in_formula)
                pressure_array.append(pressure_scale)

    return volume_per_atom, volume_per_formula, binding_potential_energy_per_atom, binding_potential_energy_per_formula, pressure_array

'''
# function to calculate the binding energy and volume per atom to dump in the output 
def binding_potential_energy(no_of_steps, step_size, num_atoms, atoms_tmp, num_atoms_in_formula):
    #import minimize_wrapper from the kim tools
    from kim_tools import minimize_wrapper

    #initialize the arrays to store the volume and binding energy
    binding_potential_energy_per_atom = []
    binding_potential_energy_per_formula = []
    volume_per_atom = []
    volume_per_formula = []
    pressure_array = []

    #loop to calculate the binding energy and volume per atom
    for i in range(0, no_of_steps + 1):
        pressure_scale = 0 + step_size * i  # pressure scale
        # Call minimize_wrapper for all iterations, including i = 0 and i = 1
        minimize_wrapper(atoms_tmp, variable_cell=True, fix_symmetry=True, flt_kwargs={'scalar_pressure': pressure_scale})
        problem_occurred = False

        try:
            # Get the potential energy and volume
            potential_energy = atoms_tmp.get_potential_energy()
            volume = atoms_tmp.get_volume()
            current_volume_per_atom = volume / num_atoms
            print('Step: %d, Pressure: %5.5f, Volume: %5.5f, Energy: %5.5f' % (i, pressure_scale, volume, potential_energy))
        except:
            print('Failed to get energy at step %d with pressure %f' % (i, pressure_scale))
            problem_occurred = True
            disclaimer = "At least one of the requested deformations of the unit cell failed to compute a potential energy."

        if not problem_occurred:
            # Calculate binding energy per atom and append values
            current_binding_potential_energy_per_atom = potential_energy / num_atoms
            volume_per_atom.append(current_volume_per_atom)
            volume_per_formula.append(current_volume_per_atom * num_atoms_in_formula)
            binding_potential_energy_per_atom.append(current_binding_potential_energy_per_atom)
            binding_potential_energy_per_formula.append(current_binding_potential_energy_per_atom * num_atoms_in_formula)
            pressure_array.append(pressure_scale)

    return volume_per_atom, volume_per_formula, binding_potential_energy_per_atom, binding_potential_energy_per_formula, pressure_array
    