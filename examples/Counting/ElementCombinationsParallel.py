#! /usr/bin/env python

import itertools
import time
from multiprocessing import Pool
from sys import stdout

import smact
from smact.screening import eneg_states_test

element_list = smact.ordered_elements(1, 103)
elements = smact.element_dictionary(element_list)

max_n = 4
neutral_stoichiometries_threshold = 8
include_pauling_test = True
pauling_test_threshold = 0.0

# Number of times to report progress during the counting loop.

count_progress_interval = 100

# Parameters for the threaded version of the code.

mp_use = True
mp_processes = 4
mp_chunk_size = 10


def print_status(text):
    "Refresh progress meter on one line"
    stdout.write(text + "\r")
    stdout.flush()


def count_iter(iterator):
    """Efficiently count number of items in an iterator"""
    return sum(1 for _ in iterator)


def count_element_combination(args):
    """Count the oxidation state combinations for a given set of elements.

    The arguments are grouped together into a single tuple. This allows use
    with the mapping functions of the multiprocessing.Pool class.

    The format of the tuple should be:
    (element_combination, n, neutral_stoichs_lookup)
    where element_combination is a tuple of smact.Element objects,
    n is the number of species to be included in a compound (e.g. n=3 for
    AxByCz materials), neutral_stoichs_lookup is a dictionary of
    the number of permitted neutral stoichiometry combinations given a sorted
    list of oxidation states.

    """
    element_combination, n, neutral_stoichs_lookup = args

    count = 0

    # Create a list of all the unique ionic species,
    # each represented as (symbol, state, eneg)
    oxidation_states_list = [
        (element.symbol, oxidation_state, element.pauling_eneg)
        for element in element_combination
        for oxidation_state in element.oxidation_states
    ]

    # Screen on charge neutrality
    # and on electronegativity ordering if include_pauling_test = True
    # (eneg_states_test is a special case of pauling_test where threshold=0,
    # repeat_anions=True, repeat_cations=True)

    for state_combination in itertools.combinations(oxidation_states_list, n):
        oxidation_states = [x[1] for x in state_combination]
        pauling_electronegativities = [x[2] for x in state_combination]

        sorted_states = tuple(sorted(oxidation_states))

        if include_pauling_test:
            if neutral_stoichs_lookup[sorted_states] and eneg_states_test(
                oxidation_states, pauling_electronegativities
            ):
                count += neutral_stoichs_lookup[sorted_states]
        else:
            count += neutral_stoichs_lookup[sorted_states]
    return count


def main():
    """
    Count element combinations which satisfy charge neutrality and Pauling
    electronegativity constraints.
    """

    # Obtain the unique oxidation states across all the elements considered.
    oxidation_states = {
        oxidation_state
        for element in list(elements.values())
        for oxidation_state in element.oxidation_states
    }

    # List unique oxidation-state combinations for each set of
    # n-element combinations
    oxidation_state_combinations = {
        n: list(itertools.combinations_with_replacement(oxidation_states, n))
        for n in range(2, max_n + 1)
    }

    # Sort combinations and cast to tuples --
    # these will be used as keys to a lookup table below.
    oxidation_state_combinations = {
        key: [tuple(sorted(item)) for item in value]
        for key, value in list(oxidation_state_combinations.items())
    }

    # Print the number of combinations of oxidation states for each value
    # of n to be analysed. Note that this will include ridiculous combinations
    # such as (+1, +2, +1)
    print("Combinations of known oxidation states:")

    for i in range(2, max_n + 1):
        print(f"m = {i}: {len(oxidation_state_combinations[i])}")

    print("")

    # Initialise a thread pool for multi-threaded calculations if required
    thread_pool = Pool(processes=mp_processes) if mp_use else None

    for n in range(2, max_n + 1):
        # Build lookup table with the number of charge-neutral stoichiometries
        # for each combination of oxidation-states.
        # Uses the tuples in oxidation_state_combinations as keys.
        # For combinations where all the states are positive or negative, the
        # number will be zero.

        # This code is rather ugly, consisting of nested functions inside a
        #  list comprehension. This is done in order to minimise memory access
        # but is equivalent to:
        #
        # neutral_stoichiometries = {}
        # for oxidation_states in oxidation_state_combinations[n]:
        #     count = n_neutral_ratios(oxidation_states)
        #     neutral_stoichiometries.update({oxidation_states: count})

        def n_neutral_ratios(oxidation_states, threshold=8):
            return len(smact.neutral_ratios(oxidation_states, threshold=threshold)[1])

        neutral_stoichiometries = {
            oxidation_states: n_neutral_ratios(
                oxidation_states, threshold=neutral_stoichiometries_threshold
            )
            for oxidation_states in oxidation_state_combinations[n]
        }

        start_time = time.time()

        # Count the number of element combinations - this is needed for the
        # progress indicator.

        combination_count = sum(
            count_iter(itertools.combinations(element_list, i)) for i in range(2, n + 1)
        )

        print("Counting ({} element combinations)" "...".format(combination_count))

        # Combinations are counted in chunks set by count_progress_interval.
        # In Python 2.7 the // symbol is "integer division" which rounds
        # downwards to the nearest integer.

        batch_size = combination_count // count_progress_interval
        data_pointer = 0

        # Set up a generator expression to return combinations of elements.
        # Generators are evaluated lazily in Python, which avoids keeping all
        # the combinations in memory - doing so can require quite a lot for
        # large n (~300 Mb for quaternaries).  To count the combinations in
        # chunks, a second generator expression is used in the while loop
        # (below) which calls next() on the generator the required number of
        # times.  As well as keeping memory usage fairly constant, this seems
        # to carry little, if any, CPU overhead compared to creating and
        # slicing a list.  So, it's perhaps a bit "filthy", but...

        combination_generator = (
            [elements[symbol] for symbol in element_combination]
            for i in range(2, n + 1)
            for element_combination in itertools.combinations(element_list, i)
        )

        # Count and sum charge-neutral combinations of oxidation states for
        # each of the combinations in element_combinations.

        count = 0

        while data_pointer < combination_count:
            iterations = min(batch_size, combination_count - data_pointer)

            imap_arg_generator = (
                (next(combination_generator), n, neutral_stoichiometries)
                for i in range(0, iterations)
            )

            if mp_use:
                # Threaded code path -- iteration over element combinations is
                # done using the multiprocessing.Pool.imap_unordered()
                # function.

                count = count + sum(
                    thread_pool.imap_unordered(
                        count_element_combination,
                        imap_arg_generator,
                        chunksize=mp_chunk_size,
                    )
                )
            else:
                # Serial code path -- iteration over element combinations is
                # done using the itertools.imap() function.

                count = count + sum(map(count_element_combination, imap_arg_generator))

            # After each chunk, report the % progress, elapsed time and an
            # estimate of the remaining time.  The smact.pauling_test() calls
            # take a variable amount of time, so the estimate of the remaining
            # time is not always accurate, particularly for the first few
            # cycles.

            data_pointer = data_pointer + iterations

            percent_done = 100 * float(data_pointer) / float(combination_count)

            time_elapsed = time.time() - start_time
            time_remaining = (
                combination_count * (time_elapsed / data_pointer) - time_elapsed
            )

            print_status(
                "  -> {}/{} done ({:.2f} %); {:.2f} s elapsed, "
                "~{:.2f} s remaining".format(
                    data_pointer,
                    combination_count,
                    percent_done,
                    time_elapsed,
                    time_remaining,
                )
            )

        print("")

        total_time = time.time() - start_time

        # Print results and total time for counting.

        print(
            "Number of charge-neutral stoichiometries for combinations "
            "of {} elements".format(n)
        )
        print("(using known oxidation states, not including zero): " "{}".format(count))
        print("")

        print(f"Total time for counting: {total_time:.3f} sec")
        print("")


if __name__ == "__main__":
    main()
