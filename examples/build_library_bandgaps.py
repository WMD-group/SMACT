import smact
from smact.properties.Band_gap_simple import band_gap_simple
from smact.properties.compound_electroneg import compound_electroneg
import itertools

# Examine a subset of periodic table
search_space = ('Li', 'Be', 'Na', 'Mg', 'K', 'Ca', 'Rb', 'Sr', 'Cs', 'Ba',
    'Al', 'Si', 'Ga', 'Ge', 'As', 'In', 'Sn', 'Sb', 'Te', 'Tl', 'Pb', 'Bi',
    'Po', 'At', 'S', 'O', 'Se', 'F', 'Cl', 'Br', 'Zn', 'Cu', 'I')

elements = smact.element_dictionary(search_space)

# Set up a generator for charge-neutral binary compounds
def binary_generator(search_space, elements):
    # Iterate over all binary combinations (e.g. (Li, Be), (Li, Na)...)    
    for el_a, el_b in itertools.combinations(search_space, 2):
        # Read possible oxidation states from element dictionary
        for ox_combo in itertools.product(elements[el_a].oxidation_states,
                                          elements[el_b].oxidation_states):
            # When a legal combination is found, yield it
            # and move onto next binary combination
            success, stoichs = smact.charge_neutrality(ox_combo)
            if success:
                yield (el_a, el_b, stoichs[0])
                break

print "no. \tband gap \tMulliken e-neg \tinternuclear dist"

# Run the generator and calculate properties
for i, (el_a, el_b, stoich) in enumerate(binary_generator(search_space, elements)):
    distance = elements[el_a].covalent_radius + elements[el_b].covalent_radius
    band_gap = band_gap_simple(False, el_a, el_b, distance,
                               elements_dict=elements)
    electroneg = compound_electroneg(elements=(el_a, el_b), stoichs=stoich,
                                     elements_dict=elements)

    print "{0}  \t{1:5.2f}  \t\t{2:5.3f}  \t\t{3:6.2f}".format(i, band_gap,
                                                electroneg, distance)
