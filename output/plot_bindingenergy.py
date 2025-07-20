import edn_format
import matplotlib.pyplot as plt 
import os
import csv

#function to read and parse the edn data
def read_edn_file(file_path):
    with open(file_path, 'r') as file:
        data = file.read()
    return edn_format.loads(data)

#function to extract data from the parsed edn 
def extract_data(data):
    volume_per_atom = data.get("volume-per-atom", {}).get("source-value", [])
    volume_per_formula = data.get("volume-per-formula", {}).get("source-value", [])
    binding_potential_energy_per_atom = data.get("binding-potential-energy-per-atom", {}).get("source-value", [])
    binding_potential_energy_per_formula = data.get("binding-potential-energy-per-formula", {}).get("source-value", [])
    pressure_array = data.get("pressure", {}).get("source-value", [])
    return volume_per_atom, volume_per_formula, binding_potential_energy_per_atom, binding_potential_energy_per_formula, pressure_array

# Function to plot the data
def plot_data(X, Y, xlabel, ylabel, title, filename):
    plt.scatter(X, Y)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(True)
    plt.tight_layout()  # Adjust layout to prevent clipping
    plt.show()
    # Save the plot to a file
    plt.savefig(filename)
    plt.show()
    plt.clf() # Clear the plot
    plt.close() # Close the plot

# Function to save data to a CSV file
def save_data_to_csv(filename, headers, data):
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(headers)
        writer.writerows(zip(*data))

# Main function
def main():
    # Use an absolute path
    file_path = '/home/openkim/test-drivers/BindingEnergyvsHydrostaticPressure/output/results.edn'
    
    # Check if the file exists
    if not os.path.exists(file_path):
        print(f"File not found: {file_path}")
        return
    
    # Read the EDN file
    data = read_edn_file(file_path)
    
    # Extract the data
    volume_per_atom, volume_per_formula, binding_potential_energy_per_atom, binding_potential_energy_per_formula, pressure_array = extract_data(data)
    
    #find the index of the pressure exactly equal to zero
    try:
        zero_pressure_index = pressure_array.index(0.0)
        print(f"Index of zero pressure: {zero_pressure_index}")
    except ValueError:
        print("Zero pressure not found in the pressure array")
    
    #find the volume of the cell at the zero pressure
    volume_Vo = volume_per_atom[zero_pressure_index]
    print(f"Volume at zero pressure (V0): {volume_Vo}")
    # Divide each element of the volume_per_atom array by volume_Vo
    normalized_volume_per_atom = [v / volume_Vo for v in volume_per_atom]

    # save the pressure array and normalized volume array to a csv file
    #save_data_to_csv(filename='pressure_vs_normalized_volume_Mo.csv', headers=['Pressure (eV/Angstrom^3)', 'Normalized Volume (V/Vo)'], data=[pressure_array, normalized_volume_per_atom])
    # Plot the data
    plot_data(
        X=pressure_array,
        Y=binding_potential_energy_per_atom,
        xlabel='Pressure (eV/Angstrom^3)',
        ylabel='Binding Energy per Atom (eV/atom)',
        title='Binding Energy per Atom vs. Pressure',
        filename='binding_energy_vs_pressure.png'
    )
    # Plot the data for volume vs binding energy per formula
    plot_data(
        X=normalized_volume_per_atom,
        Y=pressure_array,
        xlabel='V/Vo',
        ylabel='Pressure (eV/Angstrom^3)',
        title='Pressure vs. V/Vo',
        filename='volume_vs_pressure.png'
    )
if __name__ == '__main__':
    main()

