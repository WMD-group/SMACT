from itertools import izip, combinations


def pauling_test(oxidation_states, electronegativities,
                 symbols=[], repeat_anions=True,
                 repeat_cations=True, threshold=0.):
    """ Check if a combination of ions makes chemical sense,
        (i.e. positive ions should be of lower electronegativity)

    Args:
        ox (list):  oxidation states of the compound
        paul (list): the corresponding  Pauling electronegativities
            of the elements in the compound
        symbols (list) :  chemical symbols of each site.
        threshold (float): a tolerance for the allowed deviation from
            the Pauling criterion
        repeat_anions : boolean, allow an anion to repeat in different
            oxidation states in the same compound.
        repeat_cations : as above, but for cations.

    Returns:
        bool:
            True if positive ions have lower
            electronegativity than negative ions

"""

    if repeat_anions and repeat_cations and threshold == 0.:
        return eneg_states_test(oxidation_states,
                                electronegativities)

    elif repeat_anions and repeat_cations:
        return eneg_states_test_threshold(oxidation_states,
                                          electronegativities,
                                          threshold=threshold)

    else:
        if _no_repeats(oxidation_states, symbols,
                       repeat_anions=repeat_anions,
                       repeat_cations=repeat_cations):
            if threshold == 0.:
                return eneg_states_test(oxidation_states,
                                        electronegativities)
            else:
                return eneg_states_test_threshold(oxidation_states,
                                                  electronegativities,
                                                  threshold=threshold)
        else:
            return False


def _no_repeats(oxidation_states, symbols,
                repeat_anions=False, repeat_cations=False):
    """
    Check if any anion or cation appears twice.

    Args:
        oxidation_states (list): oxidation states of species
        symbols (list): chemical symbols corresponding to oxidation
            states
        repeat_anions (bool): If True, anions may be repeated (e.g. O
            in -1 and -2 states)
        repeat_cations (bool): if True, cations may be repeated (e.g.
            Cu in +1 and +2 states)

    Returns: bool
    """
    if repeat_anions is False and repeat_cations is False:
        return (len(symbols) == len(set(symbols)))
    else:

        anions, cations = [], []
        for state, symbol in izip(oxidation_states, symbols):
            if state > 0:
                cations.append(symbol)
            else:
                anions.append(symbol)
        if not repeat_anions and len(anions) != len(set(anions)):
            return False
        elif not repeat_cations and len(cations) != len(set(cations)):
            return False
        else:
            return True


def pauling_test_old(ox, paul, symbols, repeat_anions=True,
                     repeat_cations=True, threshold=0.5):
    """ Check if a combination of ions makes chemical sense,
        (i.e. positive ions should be of lower Pauling electronegativity)

    Args:
        ox (list):  oxidation states of the compound
        paul (list): the corresponding  Pauling electronegativities
            of the elements in the compound
        symbols (list) :  chemical symbols of each site.
        threshold (float): a tolerance for the allowed deviation from
            the Pauling criterion
        repeat_anions : boolean, allow an anion to repeat in different
            oxidation states in the same compound.
        repeat_cations : as above, but for cations.

    Returns:
        makes_sense (bool):
            True if positive ions have lower
            electronegativity than negative ions

    """
    if None in paul:
        return False
    positive = []
    negative = []
    pos_ele = []
    neg_ele = []
    for i, state in enumerate(ox):
        if state > 0:
            if repeat_cations:
                positive.append(paul[i])
            elif symbols[i] in pos_ele:     # Reject materials where the same
                                            # cation occupies two sites.
                return False
            else:
                positive.append(paul[i])
                pos_ele.append(symbols[i])

        if state < 0:
            if repeat_anions:
                negative.append(paul[i])
            elif symbols[i] in neg_ele:     # Reject materials where the same
                                            # anion occupies two sites.
                return False
            else:
                negative.append(paul[i])
                neg_ele.append(symbols[i])

    if len(positive) == 0 or len(negative) == 0:
        return False
    if max(positive) == min(negative):
        return False
    if max(positive) - min(negative) > threshold:
        return False
    else:
        return True


def eneg_states_test(ox_states, enegs):
    """
    Internal function for checking electronegativity criterion

    This implementation is fast as it 'short-circuits' as soon as it
    finds an invalid combination. However it may be that in some cases
    redundant comparisons are made. Performance is very close between
    this method and eneg_states_test_alternate.

    Args:
        ox_states (list): oxidation states corresponding to species
            in compound
        enegs (list): Electronegativities corresponding to species in
            compound

    Returns:
        bool : True if cations have higher electronegativity than
            anions, otherwise False

    """
    for ((ox1, eneg1), (ox2, eneg2)) in combinations(zip(ox_states, enegs), 2):
        if (ox1 > 0) and (ox2 < 0) and (eneg1 >= eneg2):
            return False
        elif (ox1 < 0) and (ox2 > 0) and (eneg1 <= eneg2):
            return False
        elif eneg1 is None or eneg2 is None:
            return False
    else:
        return True


def eneg_states_test_threshold(ox_states, enegs, threshold=0):
    """Internal function for checking electronegativity criterion

    This implementation is fast as it 'short-circuits' as soon as it
    finds an invalid combination. However it may be that in some cases
    redundant comparisons are made. Performance is very close between
    this method and eneg_states_test_alternate.

    A 'threshold' option is added so that this constraint may be
    relaxed somewhat.

    Args:
        ox_states (list): oxidation states corresponding to species
            in compound
        enegs (list): Electronegativities corresponding to species in
            compound
        threshold (Option(float)): a tolerance for the allowed deviation from
            the Pauling criterion

    Returns:
        bool : True if cations have higher electronegativity than
            anions, otherwise False

    """
    for ((ox1, eneg1), (ox2, eneg2)) in combinations(zip(ox_states, enegs), 2):
        if (ox1 > 0) and (ox2 < 0) and ((eneg1 - eneg2) > threshold):
            return False
        elif (ox1 < 0) and (ox2 > 0) and (eneg2 - eneg1) > threshold:
            return False
    else:
        return True



def eneg_states_test_alternate(ox_states, enegs):
    """Internal function for checking electronegativity criterion

    This implementation appears to be slightly slower than
    eneg_states_test, but further testing is needed.

    Args:
        ox_states (list): oxidation states corresponding to species
            in compound
        enegs (list): Electronegativities corresponding to species in
            compound

    Returns:
        bool : True if cations have higher electronegativity than
            anions, otherwise False

    """
    min_cation_eneg, max_anion_eneg = 10, 0
    for ox_state, eneg in izip(ox_states, enegs):
        if ox_state < 1:
            min_cation_eneg = min(eneg, min_cation_eneg)
        else:
            max_anion_eneg = max(eneg, max_anion_eneg)
    return min_cation_eneg > max_anion_eneg
