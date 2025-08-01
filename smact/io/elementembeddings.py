"""Interface to ElementEmbeddings. See https://github.com/WMD-group/ElementEmbeddings."""

from __future__ import annotations

from enum import Enum
from typing import TYPE_CHECKING

from elementembeddings.composition import composition_featuriser as ee_composition_featuriser
from elementembeddings.composition import species_composition_featuriser as ee_species_composition_featuriser

try:
    from enum import StrEnum
except ImportError:

    class StrEnum(str, Enum):
        """Backport of Python 3.11's StrEnum for Python 3.10."""

        def __str__(self):
            return str(self.value)


if TYPE_CHECKING:
    import pandas as pd
    from elementembeddings.composition import CompositionalEmbedding
    from elementembeddings.core import Embedding


# Should be moved to element embeddings codebase
class AllowedElementEmbeddings(StrEnum):
    """ElementEmbeddings implemented in ElementEmbeddings."""

    magpie = "magpie"
    mat2vec = "mat2vec"
    skipatom = "skipatom"
    cgnf = "cgnf"
    xenonpy = "xenonpy"
    random = "random"
    oliynyk = "oliynyk"
    matscholar = "matscholar"
    crystallm = "crystallm"
    megnet = "megnet"


class AllowedSpeciesEmbeddings(StrEnum):
    """Allowed Species Embeddings."""

    skipspecies = "skipspecies"


# Should be moved to element embeddings codebase
class PoolingStats(StrEnum):
    """Pooling statistical operations."""

    mean = "mean"
    variance = "variance"
    minpool = "minpool"
    maxpool = "maxpool"
    range = "range"
    sum = "sum"
    geometric_mean = "geometric_mean"
    harmonic_mean = "harmonic_mean"


def composition_featuriser(
    composition_data: pd.DataFrame | pd.Series | CompositionalEmbedding | list,
    formula_column: str = "formula",
    embedding: Embedding | AllowedElementEmbeddings = AllowedElementEmbeddings.magpie,
    stats: PoolingStats | list[PoolingStats] = PoolingStats.mean,
    inplace: bool = False,
) -> pd.DataFrame:
    """Wrapper to `composition_featuriser` in ElementEmbeddings."""
    return ee_composition_featuriser(
        data=composition_data,
        formula_column=formula_column,
        embedding=embedding,
        stats=stats,
        inplace=inplace,
    )


def species_composition_featuriser(
    composition_data: list[dict[str, float]],
    embedding: AllowedSpeciesEmbeddings | str = AllowedSpeciesEmbeddings.skipspecies,
    stats: PoolingStats | list[PoolingStats] = PoolingStats.mean,
    to_dataframe: bool = False,
) -> list | pd.DataFrame:
    """Compute a feature vector for a composition.

    The feature vector is based on the statistics specified
    in the stats argument.

    Args:
    ----
        composition_data: list[dict[str, float]]:
            a list of composition dictionaries
        embedding (Union[AllowedSpeciesEmbeddings, str], optional): An AllowedSpeciesEmbeddings class
            or a string
        stats (Union[str, list], optional): A list of statistics to be computed.
            The default is ['mean'].
        to_dataframe (bool, optional): Whether to return the feature vectors
            as a DataFrame. The default is False.

    Returns:
    -------
        Union[pd.DataFrame,list]: A pandas DataFrame containing the feature vector,
        or a list of feature vectors is returned
    """
    return ee_species_composition_featuriser(
        data=composition_data,
        embedding=embedding,
        stats=stats,
        to_dataframe=to_dataframe,
    )
