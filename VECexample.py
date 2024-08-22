from smact.properties import valence_electron_count

# Define the compound and the path to the valence data file
compound = "Ba5In4Bi5"
valence_file = "element_valence_modifiedry.csv"

# Calculate the Valence Electron Count (VEC)
vec = valence_electron_count(compound, valence_file)

# Print the result
print(f"The Valence Electron Count (VEC) for {compound} is: {vec:.2f}")