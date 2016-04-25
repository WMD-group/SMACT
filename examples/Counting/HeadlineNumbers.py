#! /usr/bin/env python

import itertools
from fractions import gcd

def sum_iter(iterator):
    return sum(1 for _ in iterator)


def _gcd_recursive(*args):
    """
    Get the greatest common denominator among any number of ints
    """
    if len(args) == 2:
        return gcd(*args)
    else:
        return gcd(args[0], _gcd_recursive(*args[1:]))


def _is_irreducible(stoichs):
    for i in range(1, len(stoichs)):
        if gcd(stoichs[i-1], stoichs[i]) == 1:
            return True
    else:
        return False


# Set to 103 for combinations of elements
# or 403 for combinations of element in oxidation states 
search_space = 103

            
def main():

    compound_names = {2:'binary', 3: 'ternary', 
                      4: 'quaternary', 5:'quinternary'}

    print "In a search space of {0} elements:".format(search_space)
    print ""

    def count_combinations(n):
        c = sum_iter(itertools.combinations(range(search_space), n))
        print "Number of unique combinations of {0} elements: {1}".format(n, c)
        return c

    combinations = {n: count_combinations(n) for n in range(2,5)}


    print ""

    # Count all the combinations of inequivalent ratios
    # e.g. for ternaries up to x=2: (1,1,1), (1,1,2), (1,2,1), 
    #                               (1,2,2), (2,1,1), (2,1,2), (2,1,2)
    def ineq_ratios_with_coeff(n_species, max_coeff):
        print "Number of inequivalent {0} ratios with max coefficient {1}: ".format(
                   compound_names[n_species], max_coeff),
        n_ratios = sum_iter(itertools.ifilter(
            _is_irreducible,   # Filter operation cuts equivalent
                               # stoichs e.g. accept (1,1), reject (2,2)
            itertools.product(*( # All stoich combinations, e.g.
                                  # Product of (1,2), (1,2) is 
                                  # (1,1), (2,1), (2,1), (2,2)
                                  #
                n_species * (tuple(range(1,max_coeff+1)),)
                                    # Multiply nested tuples to form
                                    # initial set of possibilities
                                    # (1,2,..), (1,2,..),...
                                ))))

        print "{0}".format(n_ratios)
        return n_ratios

    ratios = {n: {x: ineq_ratios_with_coeff(n, x) for x in range(2,9,2)} for n in range(2,5)}
    print ""


    def unique_compounds(n_species, max_coeff):
        print "Unique {0} compounds (no stability constraints) with max coefficient {1}: ".format(
        compound_names[n_species], max_coeff),
        n_compounds = combinations[n_species] * ratios[n_species][max_coeff]
        print "{0:5.3e}".format(n_compounds)
        return n_compounds

    compounds = {n: {x: unique_compounds(n, x) for x in range(4,9,2)} for n in range(2,5)}

if __name__ == '__main__':
    main()
