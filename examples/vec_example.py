from smact.properties import valence_electron_count

# Define the compound
compound = "Ba5In4Bi5"

# Calculate the Valence Electron Count (VEC)
vec = valence_electron_count(compound)

# Print the result
print(f"The Valence Electron Count (VEC) for {compound} is: {vec:.2f}")
