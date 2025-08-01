import numpy as np
import pydantic
import pandas as pd
import pathlib as pl
from typing import TypeVar
from elementembeddings.core import Embedding


T = TypeVar('T') #used internally for duck typing

#included for reference
def _PCA_transormation(coords: np.typing.NDArray[np.float64], n_pca_dims: int) -> np.typing.NDArray[np.float64]:
    '''transforms input array to principle component axes without sklearn. Each row represents one datapoint'''
    centered_coords = coords - np.mean(coords, axis=0, keepdims=True)
    variances, basis_vectors = np.linalg.eigh(np.cov(centered_coords, rowvar=False))  # PCA step
    # firts principle axis should encode the most positional variance
    sorted_basis_vectors = basis_vectors[:, np.argsort(variances)[::-1]]  # each col in basis_vectors in a vector
    transformed_coords = (centered_coords @ sorted_basis_vectors).T  # basis transformation
    return transformed_coords[:, :n_pca_dims]


class DatasetTransformations(pydantic.BaseModel):
    '''
    Helper class that facilitates running Principal Component Analysis (PCA) dimenstionality reduction and
    Kernel Density Estimator (KDE) similarity scoring on a materials dataset using a given featuriser. 

    Args:
        dataset_name (str): Name of the materials database, including version number and featuriser used.
            Should follow XXX format.
        PCA_basis_vectors (array[float]): [dataset_vector_len, PCA_vector_len] array with each column representing
            a principle component vector.
        PCA_dimension_means (array[float]): [dataset_vector_len] array where the nth elemnt corresponds to the mean
            of the nth column in the featurised dataset. Used to centre data during PCA.
        KDE_lookup_table (array[float]): array with PCA_vector_len dimensions, each of size X. the first element of a
            given dimension would map to a value of -5 and the last to 5. Works as PCA encoding could be standardised?
            Maybe would be better to just include a kernel and bandwidth, recalculating the KDE only once on class
            instantiation. Not sure what the most safe and portable approach would be.
        
    '''
    database: str
    PCA_basis_vectors: np.typing.NDArray[np.float64]
    PCA_dimension_means: np.typing.NDArray[np.float64]
    KDE_lookup_table: np.typing.NDArray[np.float64]


    @classmethod
    def _from_dataset(cls: T, path: str, featuriser: Embedding) -> T:
        '''
        Loads in materials compositional database from given path. Featurises data, calculates principle component
        analysis and kernel density estimation parameters and saves them into a standardised path.

        Args:
            path (str): file location of materials composition database. Must follow XXX strucutre with XXX Columns.
            featuriser (Embedding): function from elementembeddings library that transforms composition into chemically
                descriptive representation.
        '''
        raise NotImplementedError

    @classmethod
    def _from_pretrained(cls: T, dataset_name: str) -> T:
        '''
        Searches internally stored DatasetTransformations object for given name.

        Args:
            dataset_name (str): Name of internally stored DatasetTransformations object. Must follow XXX format.
        '''
        data_folder = pl.Path(__file__).parent / 'tranformation_files'
        available_datasets = data_folder.glob(pattern='*.json')
        if dataset_name not in available_datasets:
            raise NameError(f'{dataset_name} not found. {available_datasets=}')
        raise NotImplementedError


    def get_similarities(self, compositions: list[str]) -> list[float]:
        '''
        Featurises given compositions, calculates the Principle Component Analysis transform, and finds the Kernel
        Densitry Estimation similarity using the internally stored parameteres trained on the class' dataset.
        '''
        raise NotImplementedError
