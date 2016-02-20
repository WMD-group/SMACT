"""
 Example script using distorter

 Generate all possible (symmetry inequivalent) substitutions of Sr on Ba
 sites; single and double substitutions.

"""

import smact.builder as builder
import smact.distorter as distort


def pretty_print_atoms(atoms, linewrap=15):
    entries = ["{0:5.3f} {1:5.3f} {2:5.3f} {symbol}".format(*position,
                                                            symbol=symbol)
               for position, symbol in zip(atoms.get_positions(),
                                           atoms.get_chemical_symbols())]

    for output_i in range(linewrap):
        line = ''
        for i, entry in enumerate(entries):
            if (output_i == i % linewrap):
                line = line + entry + '\t'
        print line

# Build the input
smact_lattice, test_case = builder.cubic_perovskite(['Ba', 'Ti', 'O'],
                                                    repetitions=[2, 2, 2])

hlinewidth = 68


print ('-' * hlinewidth)
print "Original coordinates: "
pretty_print_atoms(test_case)
print ('-' * hlinewidth)

# Do the single substitution first, it is trivial as all Ba
#  sites are equivalent we will choose the first Ba
subs_site = [0.0, 0.0, 0.0]
single_substitution = distort.make_substitution(test_case, subs_site, "Sr")
print "Single: "
pretty_print_atoms(single_substitution)

# Build a sub-lattice you wish to disorder [test case do the Ba sub-lattice]
sub_lattice = distort.build_sub_lattice(single_substitution, "Ba")

# Enumerate the inequivalent sites
inequivalent_sites = distort.get_inequivalent_sites(sub_lattice,
                                                    single_substitution)

# Replace Ba at inequivalent sites with Sr
for i, inequivalent_site in enumerate(inequivalent_sites):
    print ('-' * hlinewidth)
    print "Substituted coordinates {0}".format(i)
    # print test_case,inequivalent_site
    distorted = distort.make_substitution(single_substitution,
                                          inequivalent_site,
                                          "Sr")
    pretty_print_atoms(distorted)

print('='*hlinewidth)
