"""SMACT benchmarking."""

from __future__ import annotations

from smact.structure_prediction.mutation import CationMutator

from .utilities import timeit


class MutatorBenchmarker:
    """Benchmarking tests for CationMutator."""

    @timeit
    def run_tests(self):
        """Initialize Mutator and perform tests."""
        self.__cm_setup()
        self.__pair_corr()

    @timeit
    def __cm_setup(self):
        """Create a CationMutator."""
        self.cm = CationMutator.from_json()

    @timeit
    def __pair_corr(self):
        """Get pair correlations."""
        self.cm.complete_pair_corrs()


@timeit(delim=True, n=100)
def mutator_test_run():
    """Run benchmark tests for CationMutator."""
    MutatorBenchmarker().run_tests()
