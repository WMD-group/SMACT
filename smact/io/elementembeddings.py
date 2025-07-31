"""Interface to ElementEmbeddings. See https://github.com/WMD-group/ElementEmbeddings."""

from __future__ import annotations

from enum import StrEnum, auto
from typing import TYPE_CHECKING

from elementembeddings.composition import composition_featuriser as ee_composition_featuriser
from elementembeddings.composition import species_composition_featuriser as ee_species_composition_featuriser

if TYPE_CHECKING:
    import pandas as pd
    from elementembeddings.composition import CompositionalEmbedding
    from elementembeddings.core import Embedding


# Should be moved to element embeddings codebase
class AllowedElementEmbeddings(StrEnum):
    """ElementEmbeddings implemented in ElementEmbeddings."""

    magpie = auto()
    mat2vec = auto()
    skipatom = auto()
    cgnf = auto()
    xenonpy = auto()
    random = auto()
    oliynyk = auto()
    matscholar = auto()
    crystallm = auto()
    megnet = auto()


class AllowedSpeciesEmbeddings(StrEnum):
    """Allowed Species Embeddings."""

    skipspecies = auto()


# Should be moved to element embeddings codebase
class PoolingStats(StrEnum):
    """Pooling statistical operations."""

    mean = auto()
    variance = auto()
    minpool = auto()
    maxpool = auto()
    range = auto()
    sum = auto()
    geometric_mean = auto()
    harmonic_mean = auto()


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
