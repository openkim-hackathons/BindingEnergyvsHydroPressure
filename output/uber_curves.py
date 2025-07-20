import csv
import matplotlib.pyplot as plt
import os
import numpy as np
from math import sqrt

# Function to read data from the CSV file
def read_csv_file(file_path):
    volume_array = []
    binding_energy_array = []
    with open(file_path, 'r') as csvfile:
        reader = csv.DictReader(csvfile)  # Use DictReader to access columns by header names
        for row in reader:
            volume_array.append(float(row['Volume per Atom']))  # Replace with the column header for volume
            binding_energy_array.append(float(row['Binding Potential Energy per Atom']))  # Replace with the column header for binding energy
    return volume_array, binding_energy_array

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
        volume_array, binding_energy_array = read_csv_file(file_path)
        
        # Compute normalized volumes
        normalized_volumes = [v / V0 for v in volume_array]
        
        # Dynamically compute C for the current element
        C = 3 * sqrt((V0 * B[element]) / E_rwse[element])
        
        # Compute a_star for the current element
        a_star = [C * (np.cbrt(v / V0) - 1) for v in volume_array]

        # Compute E_star i.e., scaled binding energy
        E_star = [Erw / (E_rwse[element]) for Erw in binding_energy_array]
        # Plot the data
        plt.plot(a_star, E_star, 'o', label=f'{element}', markersize=4)

    # Add the polynomial curve to the same plot
    x_values = np.linspace(-0.8, 1.5, 100)  # Define the range of a_star values
    y_values = polynomial_a_star(x_values)  # Compute the polynomial values
    plt.plot(x_values, y_values, label='Polynomial Fit', color='red', linestyle='--')  # Add the polynomial curve


    # Add labels, title, and legend for a_star plot
    plt.xlabel('a* (Scaled Volume)')
    plt.ylabel('Binding Energy per Atom (eV/atom)')
    plt.title('Binding Energy per Atom vs Scaled Volume (a*)')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()  # Adjust layout to prevent clipping
    plt.savefig('binding_energy_vs_scaled_volume.png')  # Save the plot as a PNG file
    plt.close()  # Close the figure to avoid overlapping plots

    # Plot binding energy vs volume_array
    plt.figure()
    for element, info in csv_files.items():
        file_path = info['file']
        V0 = info['V0']
        
        if not os.path.exists(file_path):
            print(f"File not found: {file_path}")
            continue
        
        # Read the data
        volume_array, binding_energy_array = read_csv_file(file_path)
        
        # Plot binding energy vs volume_array
        plt.plot(volume_array, binding_energy_array, 'o', label=f'{element}')

    # Add labels, title, and legend for volume_array plot
    plt.xlabel('Volume per Atom (Angstrom^3)')
    plt.ylabel('Binding Energy per Atom (eV/atom)')
    plt.title('Binding Energy per Atom vs Volume per Atom')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()  # Adjust layout to prevent clipping
    plt.savefig('binding_energy_vs_volume.png')  # Save the plot as a PNG file
    plt.close()  # Close the figure to avoid overlapping plots

if __name__ == '__main__':
    main()

