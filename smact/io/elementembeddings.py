"""Interface to ElementEmbeddings. See https://github.com/WMD-group/ElementEmbeddings."""

from __future__ import annotations

from typing import TYPE_CHECKING

from smact.utils.compat import StrEnum

if TYPE_CHECKING:
    import pandas as pd
    from elementembeddings.composition import CompositionalEmbedding  # type: ignore[import-untyped]
    from elementembeddings.core import Embedding  # type: ignore[import-untyped]

try:
    from elementembeddings.composition import (  # type: ignore[import-untyped]
        composition_featuriser as ee_composition_featuriser,
    )
    from elementembeddings.composition import (  # type: ignore[import-untyped]
        species_composition_featuriser as ee_species_composition_featuriser,
    )

    HAS_ELEMENTEMBEDDINGS = True
except ImportError:
    HAS_ELEMENTEMBEDDINGS = False

    def ee_composition_featuriser(*args, **kwargs):  # type: ignore[misc]
        """Stub — never called due to HAS_ELEMENTEMBEDDINGS guard."""
        raise ImportError("ElementEmbeddings not installed")

    def ee_species_composition_featuriser(*args, **kwargs):  # type: ignore[misc]
        """Stub — never called due to HAS_ELEMENTEMBEDDINGS guard."""
        raise ImportError("ElementEmbeddings not installed")


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
    """Compute feature vectors for compositions using element embeddings.

    Wrapper around `composition_featuriser` in ElementEmbeddings.

    Args:
    ----
        composition_data: Input compositions as a DataFrame, Series,
            CompositionalEmbedding, or list of formula strings.
        formula_column (str): Column name containing formulas if input
            is a DataFrame. Default is "formula".
        embedding: Element embedding to use for featurisation.
            Default is MAGPIE.
        stats: Pooling statistic(s) to compute. Default is mean.
        inplace (bool): Whether to modify the input DataFrame in place.
            Default is False.

    Returns:
    -------
        pd.DataFrame: DataFrame containing the computed feature vectors.
    """
    if not HAS_ELEMENTEMBEDDINGS:
        raise ImportError(
            "The 'ElementEmbeddings' package is required for this function. "
            "Install it with: pip install ElementEmbeddings"
        )

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
    if not HAS_ELEMENTEMBEDDINGS:
        raise ImportError(
            "The 'ElementEmbeddings' package is required for this function. "
            "Install it with: pip install ElementEmbeddings"
        )

    return ee_species_composition_featuriser(
        data=composition_data,
        embedding=embedding,
        stats=stats,
        to_dataframe=to_dataframe,
    )
