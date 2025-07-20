'''
import csv
import os
import matplotlib.pyplot as plt

import csv
import matplotlib.pyplot as plt
import os
import numpy as np
from math import sqrt

# Function to read data from the CSV file
def read_csv_file(file_path):
    volume_array = []
    binding_energy_array = []
    pressure = []
    with open(file_path, 'r') as csvfile:
        reader = csv.DictReader(csvfile)  # Use DictReader to access columns by header names
        for row in reader:
            volume_array.append(float(row['Volume per Atom']))  # Replace with the column header for volume
            binding_energy_array.append(float(row['Binding Potential Energy per Atom']))  # Replace with the column header for binding energy
            pressure.append(float(row['Pressure Array']))
    return volume_array, binding_energy_array, pressure

# Equilibrium binding energy per atom in eV/atom is, which is -E(rws)
E_rwse = {
    'Al': 3.365,
    'Fe': 4.29,
    'Mo': 10.853,
    'Pd': 3.91
}

# Bulk modulus of each element in eV/Å^3
B = {
    'Al': 0.4955,
    'Fe': 1.0760,
    'Mo': 1.6559,
    'Pd': 1.2202
}

# Polynomial function for E* vs a*
def polynomial_a_star(x):
    return (-1 - x - 0.05 * (x**3)) * np.exp(-x)

# Main function
def main():
    # List of CSV files and corresponding V0 values
    csv_files = {
        'Al': {
            'file': '/home/openkim/test-drivers/BindingEnergyvsHydrostaticPressure/output/extracted_arrays_Al.csv',
            'V0': 16.5424
        },
        'Fe': {
            'file': '/home/openkim/test-drivers/BindingEnergyvsHydrostaticPressure/output/extracted_arrays_Fe.csv',
            'V0': 11.7410
        },
        'Mo': {
            'file': '/home/openkim/test-drivers/BindingEnergyvsHydrostaticPressure/output/extracted_arrays_Mo.csv',
            'V0': 15.9380
        },
        'Pd': {
            'file': '/home/openkim/test-drivers/BindingEnergyvsHydrostaticPressure/output/extracted_arrays_Pd.csv',
            'V0': 14.7057
        }
    }
    
    # Plot the data for each file
    plt.figure()
    for element, info in csv_files.items():
        file_path = info['file']
        V0 = info['V0']
        
        if not os.path.exists(file_path):
            print(f"File not found: {file_path}")
            continue
        
        # Read the data
        volume_array, binding_energy_array, pressure = read_csv_file(file_path)
        
        # Compute normalized volumes
        normalized_volumes = [v / V0 for v in volume_array]
        
        # Dynamically compute C for the current element
        C = 3 * sqrt((V0 * B[element]) / E_rwse[element])

    import numpy as np
    #scale the pressure
    # Convert to NumPy array
    Al_pressures = np.array(Al_pressures)
    Fe_pressures = np.array(Fe_pressures)
    Mo_pressures = np.array(Mo_pressures)
    Pd_pressures = np.array(Pd_pressures)
    Al_pressures_scaled = Al_pressures / B_Al
    Fe_pressures_scaled = Fe_pressures / B_Fe
    Mo_pressures_scaled = Mo_pressures / B_Mo
    Pd_pressures_scaled = Pd_pressures / B_Pd
    # Plot all pressure vs. normalized volume data together in one plot
    plt.figure()
    plt.plot(Al_normalized_volumes, Al_pressures_scaled, 'o', label='Al', markersize=4)
    #plt.plot(Fe_normalized_volumes, Fe_pressures_scaled, 'o', label='Fe', markersize=4)
    plt.plot(Mo_normalized_volumes, Mo_pressures_scaled, 'o', label='Mo', markersize=4)
    #plt.plot(Pd_normalized_volumes, Pd_pressures_scaled, 'o', label='Pd', markersize=4)
    plt.xlabel('Normalized Volume (V/Vo)')
    plt.ylabel('Pressure (eV/Angstrom^3)/B')
    plt.title('Pressure vs. Normalized Volume')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()  # Adjust layout to prevent clipping
    plt.savefig('pressure_vs_normalized_volume_combined.png')
    plt.show()

    # Define the pressure array p/B
    p_over_B_Al = [p / B_Al for p in Al_pressures]
    p_over_B_Fe = [p / B_Fe for p in Fe_pressures]
    p_over_B_Mo = [p / B_Mo for p in Mo_pressures]
    p_over_B_Pd = [p / B_Pd for p in Pd_pressures]

    #Volume at zero pressure
    V_Al = 16.5424
    V_Fe = 11.7410
    V_Mo = 15.9380
    V_Pd = 14.7057

    # Equilibrium binding energy per atom in eV/atom is, which is -E(rws)
    E_rwse_Al = 3.365
    E_rwse_Fe = 4.29
    E_rwse_Mo = 10.853
    E_rwse_Pd = 3.91

    from math import sqrt
    #find C for all elements
    C_Al = 3 * sqrt((V_Al * B_Al)/(E_rwse_Al))
    C_Fe = 3 * sqrt((V_Fe * B_Fe)/(E_rwse_Fe))        
    C_Mo = 3 * sqrt((V_Mo * B_Mo)/(E_rwse_Mo))
    C_Pd = 3 * sqrt((V_Pd * B_Pd)/(E_rwse_Pd))
    print(f"C_Al: {C_Al}")
    print(f"C_Fe: {C_Fe}")
    print(f"C_Mo: {C_Mo}")
    print(f"C_Pd: {C_Pd}")
    # scale V/V0 values

    from math import exp
    x_Al = [(3*((1/v)**(1/3))*((1/v)**(1/3)-1)*(1-0.15*C_Al*((v**(1/3))-1)+0.05*(C_Al**2)*((v**(1/3))-1)**2)*exp(-C_Al*((v**(1/3))-1))) for v in Al_normalized_volumes]
    x_Fe = [(3*((1/v)**(1/3))*((1/v)**(1/3)-1)*(1-0.15*C_Fe*((v**(1/3))-1)+0.05*(C_Fe**2)*((v**(1/3))-1)**2)*exp(-C_Fe*((v**(1/3))-1))) for v in Fe_normalized_volumes]
    x_Mo = [(3*((1/v)**(1/3))*((1/v)**(1/3)-1)*(1-0.15*C_Mo*((v**(1/3))-1)+0.05*(C_Mo**2)*((v**(1/3))-1)**2)*exp(-C_Mo*((v**(1/3))-1))) for v in Mo_normalized_volumes]
    x_Pd = [(3*((1/v)**(1/3))*((1/v)**(1/3)-1)*(1-0.15*C_Pd*((v**(1/3))-1)+0.05*(C_Pd**2)*((v**(1/3))-1)**2)*exp(-C_Pd*((v**(1/3))-1))) for v in Pd_normalized_volumes]
    
    import numpy as np
    def y(x, C):
        #C_Al=15
        x_cbrt = np.cbrt(x)
        x_cbrt = np.cbrt(x)  # Compute cube root safely
        numerator = 3 * ((1 / x_cbrt)**2-(1 / x_cbrt))
        #denominator = x**(2/3)  # This is fine for positive x
        factor = (1 - 0.15 * C * (x_cbrt - 1) + 0.05 * (C**2) * (x_cbrt - 1)**2)
        exponent = np.exp(- C * (x_cbrt - 1))
        y_value = (numerator) * factor * exponent
        y_secondary = numerator   # Second return value
        
        return y_value, y_secondary
    #generate x values
    x = np.linspace(0.9, 1.4, 100)
    print(x)
   # Compute y values (unpacking the tuple)
    y1_with_C_Fe, y2_with_C_Fe = y(x, C_Fe)
    y1_with_C_Al, y2_with_C_Al = y(x, C_Al)
    y1_with_C_Mo, y2_with_C_Mo = y(x, C_Mo)
    y1_with_C_Pd, y2_with_C_Pd = y(x, C_Pd)

    # Plot the scaled p/B vs x for all elements in one plot
    plt.figure()
    #plt.plot(x_Al, p_over_B_Al, label='Al')
    #plt.plot(x_Fe, p_over_B_Fe, label='Fe')
    #plt.plot(x_Mo, p_over_B_Mo, label='Mo')
    #plt.plot(x_Pd, p_over_B_Pd, label='Pd')
    plt.plot(Al_normalized_volumes, Al_pressures, 'o', label='Al', markersize=4)
    #plt.plot(Fe_normalized_volumes, Fe_pressures, 'o', label='Fe', markersize=4)
    plt.plot(Mo_normalized_volumes, Mo_pressures, 'o', label='Mo', markersize=4)
    #plt.plot(Pd_normalized_volumes, Pd_pressures, 'o', label='Pd', markersize=4)
    #plt.plot(x, y1_with_C_Fe, label='analytical')
    #plt.plot(x, y2_with_C_Fe, label='analytical without exponential')
    plt.plot(x, y1_with_C_Al, label='analytical')
    #plt.plot(x, y2_with_C_Al, label='analytical without exponential')
    plt.plot(x, y1_with_C_Mo, label='analytical', color='red', linestyle='--')
    #plt.plot(x, y2_with_C_Mo, label='analytical without exponential')
    #plt.plot(x, y1_with_C_Pd, label='analytical')
    #plt.plot(x, y2_with_C_Pd, label='analytical without exponential')

    plt.xlabel('x')
    plt.ylabel('p/B')
    plt.title('Scaled p/B vs x')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()  # Adjust layout to prevent clipping
    plt.savefig('scaled_p_over_B_vs_x_combined.png')
    plt.show()

if __name__ == '__main__':
    main()

'''
import csv
import os
import matplotlib.pyplot as plt
import numpy as np
from math import sqrt, exp

# Function to read data from the CSV file
def read_csv_file(file_path):
    volume_array = []
    binding_energy_array = []
    pressure = []
    with open(file_path, 'r') as csvfile:
        reader = csv.DictReader(csvfile)  # Use DictReader to access columns by header names
        for row in reader:
            volume_array.append(float(row['Volume per Atom']))  # Replace with the column header for volume
            binding_energy_array.append(float(row['Binding Potential Energy per Atom']))  # Replace with the column header for binding energy
            pressure.append(float(row['Pressure Array']))
    return volume_array, binding_energy_array, pressure

# Equilibrium binding energy per atom in eV/atom is, which is -E(rws)
E_rwse = {
    'Al': 3.365,
    'Mo': 10.853,
    'Pd': 3.91,
    'Fe': 4.29
}

# Bulk modulus of each element in eV/Å^3
B = {
    'Al': 0.4955,
    'Mo': 1.6559,
    'Pd': 1.2202,
    'Fe': 1.0760
}

def y(x, C):
    x_cbrt = np.cbrt(x)  # Compute cube root safely
    numerator = 3 * ((1 / x_cbrt)**2 - (1 / x_cbrt))
    factor = (1 - 0.15 * C * (x_cbrt - 1) + 0.05 * (C**2) * (x_cbrt - 1)**2)
    exponent = np.exp(-C * (x_cbrt - 1))
    y_value = numerator * factor * exponent
    return y_value

def main():
    # File and V0 values for Al and Mo
    csv_files = {
        'Al': {
            'file': '/home/openkim/test-drivers/BindingEnergyvsHydrostaticPressure/output/extracted_arrays_Al.csv',
            'V0': 16.5424
        },
        'Mo': {
            'file': '/home/openkim/test-drivers/BindingEnergyvsHydrostaticPressure/output/extracted_arrays_Mo.csv',
            'V0': 15.9380
        },
        'Pd': {
            'file': '/home/openkim/test-drivers/BindingEnergyvsHydrostaticPressure/output/extracted_arrays_Pd.csv',
            'V0': 14.7057 
        },
        'Fe': {
            'file': '/home/openkim/test-drivers/BindingEnergyvsHydrostaticPressure/output/extracted_arrays_Fe.csv',
            'V0': 11.7410
        }
    }

    # Compute C values for Al and Mo
    C_values = {}
    for element, info in csv_files.items():
        V0 = info['V0']
        C_values[element] = 3 * sqrt((V0 * B[element]) / E_rwse[element])
        print(f"C_{element}: {C_values[element]}")

    # Generate x values for analytical curves
    x = np.linspace(0.9, 1.4, 100)

    # Plot the data for Al and Mo
    plt.figure()
    for element, info in csv_files.items():
        file_path = info['file']
        V0 = info['V0']

        if not os.path.exists(file_path):
            print(f"File not found: {file_path}")
            continue

        # Read the data
        volume_array, binding_energy_array, pressure_array = read_csv_file(file_path)

        # Compute normalized volumes
        normalized_volumes = [v / V0 for v in volume_array]

        # Compute p/B
        p_over_B = [p / B[element] for p in pressure_array]

        # Plot normalized volumes vs p/B
        plt.plot(normalized_volumes, p_over_B, 'o', label=f'{element}', markersize=4)

        # Compute and plot analytical curve
        y_values = y(x, C_values[element])
        plt.plot(x, y_values, label=f'{element} (analytical)', linestyle='--')

    # Add labels, title, and legend
    plt.xlabel('Normalized Volume (V/V0)')
    plt.ylabel('p/B')
    plt.title('Scaled p/B vs Normalized Volume for Al and Mo')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()  # Adjust layout to prevent clipping
    plt.savefig('scaled_p_over_B_vs_normalized_volume_Al_Mo_Pb_Fe.png')
    plt.show()

if __name__ == '__main__':
    main()
