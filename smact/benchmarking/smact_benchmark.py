"""SMACT benchmarking."""

from .utilities import timeit
from ..structure_prediction.mutation import CationMutator


class MutatorBenchmarker:
    """Benchmarking tests for CationMutator."""

    @timeit
    def run_tests(self):
        """Initialize Mutator and perform tests."""
        self.__cm_setup()
        self.__pair_corr()

    @timeit
    def __cm_setup(self) -> CationMutator:
        """Create a CationMutator."""
        self.cm = CationMutator.from_json()

    @timeit
    def __pair_corr(self) -> pd.DataFrame:
        """Get pair correlations."""
        return self.cm.complete_pair_corrs()


if __name__ == "__main__":
    MutatorBenchmarker().run_tests()
