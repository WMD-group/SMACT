#! /usr/bin/env python

import smact
import itertools

elements = smact.element_dictionary()

def main():
    oxidation_states = [element.oxidation_states for _, element in elements.iteritems()]

    oxidation_states = set([state for sublist in oxidation_states for state in sublist])

    oxidation_state_combinations = {n: list(
        itertools.combinations_with_replacement(oxidation_states, n))
                                        for n in range(2,5)}

    # Sort combinations and cast to tuples for unambiguous keys
    oxidation_state_combinations = {n: [tuple(sorted(x)) for x in l]
                    for n, l in oxidation_state_combinations.items()}

    print "Combinations of known oxidation states:"
    print "Binary:     {0}".format(len(oxidation_state_combinations[2]))
    print "Ternary:    {0}".format(len(oxidation_state_combinations[3]))
    print "Quaternary: {0}".format(len(oxidation_state_combinations[4]))

    neutral_stoichs = {}
    for n in range(2,5):
        neutral_stoichs.update({n: {
            ox_states: tuple(smact.charge_neutrality_iter(ox_states, threshold=8)) 
                for ox_states in oxidation_state_combinations[n]
            }
        })

    print "Number of charge-neutral stoichiometries for combinations of 2 elements"
    print "(using known oxidation states, not including zero): ",
    print sum(neutral_stoichs_count(combo, neutral_stoichs[2]) for combo in element_ox_combos(2))

    print "Number of charge-neutral stoichiometries for combinations of 3 elements"
    print "(using known oxidation states, not including zero): ",
    print sum(neutral_stoichs_count(combo, neutral_stoichs[3]) for combo in element_ox_combos(3))

    print "Number of charge-neutral stoichiometries for combinations of 4 elements"
    print "(using known oxidation states, not including zero): ",
    print sum(neutral_stoichs_count(combo, neutral_stoichs[4]) for combo in element_ox_combos(4))


def neutral_stoichs_count(ox_combo, neutral_stoichs_dict):
    return len(neutral_stoichs_dict[tuple(sorted(ox_combo))])

def element_ox_combos(n, elements=elements):
    for element_combo in element_combinations(n):
        for ox_combo in itertools.product(*[el.oxidation_states 
                                        for el in element_combo]):
            yield ox_combo

def element_combinations(n, elements=elements):
    symbols = elements.keys()
    for el_combo in itertools.combinations(symbols, n):
        yield [elements[sym] for sym in el_combo]

if __name__ == '__main__':
    main()
