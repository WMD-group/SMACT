### Use SMACT to generate a list of sensible chemical compositions ###
### and save in some useful format.                                ###

# Imports
from smact.screening import smact_test
from smact import Element, element_dictionary
import itertools
import multiprocessing
from functools import partial
from pymatgen import Composition
from datetime import datetime
import pickle

# Get a dictionary of all Element objects keyed by symbol
all_el = element_dictionary()
all_symbols = [k for k,i in all_el.items()]

### === Set up parameters here === ###
# Either enter symbols manually as str, or use the element dictionary
# that has been imported if there are lots.
symbols = all_symbols[:20]
symbols_to_ignore = ['O']
symbols = [x for x in symbols if x not in symbols_to_ignore]

# Enter elements that must be in every compound, e.g. ['O'] if you
# are only interested in oxides. Else leave as None or remove.
# These elements should also appear in symbols_to_ignore above.
always_include = ['O']

# Set stoichiometry threshold
threshold = 4

# Set order of combinations, e.g 2 for binaries, 3 for ternaries.
# MIND: If you are using always_include, order = elements per composition - 1
# e.g. to explore ternary oxides, order = 2.
order = 3

# Set name of files for saving the list:
filekey = 'Test_quaternary_oxides'

# Convert to pretty formulas using Pymatgen?
pretty_formulas_export = True

### === You don't need to edit anything below here === ###

elements = [all_el[x] for x in symbols]
if always_include:
    include = [all_el[x] for x in always_include]

# Multiprocessing-friendly function
def MapFunction(els):
    return smact_test(els, threshold = threshold, include = include)

# Function to convert into pretty formulas if desired
def comp_maker(comp):
    form = []
    for el, ammt in zip(comp[0], comp[1]):
        form.append(el)
        form.append(ammt)
    form = ''.join(str(e) for e in form)
    pmg_form = Composition(form).reduced_formula
    return pmg_form

if __name__ == "__main__":
    # Set up iterator
    all_combos = itertools.combinations(elements, order)

    print("Elements to consider: {0}".format(symbols))
    print("Elements to include in every composition: {0}".format(always_include))

    # Do the smact tests using multiprocessing
    print("Running SMACT tests...")
    start_time = datetime.now()
    p = multiprocessing.Pool()
    result = p.map(MapFunction, all_combos)
    print("Time to complete SMACT tests: {0}".format(datetime.now() - start_time))

    # Tidy output into one list and do pickling
    print("Tidying everything into one list...")
    flat_list = [item for sublist in result for item in sublist]

    print("Pickling list of compositions to {0}_compositions.pkl...".format(filekey))
    with open('{0}_compositions.pkl'.format(filekey), 'wb') as f:
        pickle.dump(flat_list, f)

    if pretty_formulas_export:
        print("Converting to a list of unique pretty formulas... ")
        pretty_formulas = p.map(comp_maker, flat_list)
        pretty_formulas = list(set(pretty_formulas))
        print("Pickling list of pretty formulas to {0}_prettyform.pkl...".format(filekey))
        with open('{0}prettyform.pkl'.format(filekey), 'wb') as f:
            pickle.dump(pretty_formulas, f)

    # Print summary of results
    print("SUMMARY:")
    print("Number of compositions: {0}".format(len(flat_list)))
    if pretty_formulas_export:
        print("Number of unique pretty formulas: {0}".format(len(pretty_formulas)))
    print("Total time: {0}".format(datetime.now() - start_time))
    print("Thanks for playing.")
