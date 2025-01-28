import edn_format
import matplotlib.pyplot as plt 
import os

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
    print('volume_per_atom:', volume_per_atom)
    # Plot the data
    plot_data(
        X=pressure_array,
        Y=binding_potential_energy_per_atom,
        xlabel='Pressure (GPa)',
        ylabel='Binding Energy per Atom (eV)',
        title='Binding Energy per Atom vs. Pressure',
        filename='binding_energy_vs_pressure.png'
    )
    # Plot the data for volume vs binding energy per formula
    plot_data(
        X=pressure_array,
        Y=volume_per_atom,
        xlabel='Pressure (GPa)',
        ylabel='Volume per Atom (Angstrom^3)',
        title='Volume per Atom vs. Pressure',
        filename='volume_vs_pressure.png'
    )
if __name__ == '__main__':
    main()
