import smact
from smact.properties.Band_gap_simple import band_gap_simple
import copy

# Get correct directory if calling from elsewhere
# (if statement needed in case of empty response)

# Generate a dictionary elements, form the dataset oxidationstates.data
# Dictionary contains elements and their oxidation states
# Reduce the regions of the periodic table to visit, by using search_space
search_space = ('Li', 'Be', 'Na', 'Mg', 'K', 'Ca', 'Rb', 'Sr', 'Cs', 'Ba',
    'Al', 'Si', 'Ga', 'Ge', 'As', 'In', 'Sn', 'Sb', 'Te', 'Tl', 'Pb', 'Bi',
    'Po', 'At', 'S', 'O', 'Se', 'F', 'Cl', 'Br', 'Zn', 'Cu', 'I')

# Get the list of possible constituent elements

# Generate list of binary compositions which satisfy charge neutrality
neutral_species = []
for i, ele_a in enumerate(search_space):
    for ox_a in smact.Element(ele_a).oxidation_states:
        for j, ele_b in enumerate(search_space):
            if j > i:
                for ox_b in smact.Element(ele_b).oxidation_states:
                    neutral_exists, neutral_ratios = smact.charge_neutrality([ox_a,ox_b],threshold = 4)
                    if neutral_exists:
                        combo = [smact.Element(ele_a).symbol, smact.Element(ele_b).symbol]
                        neutral_species.append([combo,neutral_ratios])
print len(neutral_species)

# Calculate Band Gap and Mulliken electronegativity for each composition and write to file
out = open("output.txt","w")
print "no. ", "band gap ", "Mulliken e-neg ", "internuclear dist"
out.write("no. " +", " + "band gap " +", " + "Mulliken e-neg " +", " + "internuclear dist" + "\n")

element_compositions = copy.deepcopy(neutral_species)

for combo in neutral_species:
    dist = smact.Element(combo[0][0]).covalent_radius + smact.Element(combo[0][1]).covalent_radius
    bandgap = band_gap_simple(False, combo[0][0], combo[0][1], dist)

    print combo[0][0], combo[0][1], bandgap
