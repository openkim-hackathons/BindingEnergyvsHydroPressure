import edn_format
import os
import csv

# Function to read and parse the EDN file
def read_edn_file(file_path):
    with open(file_path, 'r') as file:
        data = file.read()
    return edn_format.loads(data)

# Function to extract arrays from the parsed EDN data
def extract_arrays(data):
    volume_per_atom = data.get("volume-per-atom", {}).get("source-value", [])
    volume_per_formula = data.get("volume-per-formula", {}).get("source-value", [])
    binding_potential_energy_per_atom = data.get("binding-potential-energy-per-atom", {}).get("source-value", [])
    binding_potential_energy_per_formula = data.get("binding-potential-energy-per-formula", {}).get("source-value", [])
    pressure_array = data.get("pressure", {}).get("source-value", [])
    return volume_per_atom, volume_per_formula, binding_potential_energy_per_atom, binding_potential_energy_per_formula, pressure_array

# Function to save arrays to a CSV file
def save_arrays_to_csv(filename, headers, data):
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(headers)  # Write the headers
        writer.writerows(zip(*data))  # Write the data row by row

# Main function
def main():
    # Path to the EDN file
    file_path = '/home/openkim/test-drivers/BindingEnergyvsHydrostaticPressure/output/results.edn'
    
    # Check if the file exists
    if not os.path.exists(file_path):
        print(f"File not found: {file_path}")
        return
    
    # Read the EDN file
    data = read_edn_file(file_path)
    
    # Extract arrays
    volume_per_atom, volume_per_formula, binding_potential_energy_per_atom, binding_potential_energy_per_formula, pressure_array = extract_arrays(data)
    
    # Print the extracted arrays
    print("Volume per Atom:", volume_per_atom)
    print("Volume per Formula:", volume_per_formula)
    print("Binding Potential Energy per Atom:", binding_potential_energy_per_atom)
    print("Binding Potential Energy per Formula:", binding_potential_energy_per_formula)
    print("Pressure Array:", pressure_array)
    
    # Save the arrays to a CSV file
    output_csv = '/home/openkim/test-drivers/BindingEnergyvsHydrostaticPressure/output/extracted_arrays.csv'
    headers = ['Volume per Atom', 'Volume per Formula', 'Binding Potential Energy per Atom', 'Binding Potential Energy per Formula', 'Pressure Array']
    data_to_save = [volume_per_atom, volume_per_formula, binding_potential_energy_per_atom, binding_potential_energy_per_formula, pressure_array]
    save_arrays_to_csv(output_csv, headers, data_to_save)
    print(f"Arrays saved to {output_csv}")

if __name__ == '__main__':
    main()