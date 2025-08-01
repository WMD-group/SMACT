from sklearn.decomposition import PCA
import numpy as np
import pydantic
import pandas as pd

class metadata(pydantic.BaseModel):
    '''
    Data class that reads in stored data required to run Principal Component Analysis (PCA)
    dimenstionality reduction and Kernel Density Estimator (KDE) similarity scorer pretrained on
    a given dataset. 
    '''
    dataset: str
    PCA_transform_matrix: np.typing.NDArray[np.float64]
    PCA_dims: int
    KDE_Transform_matrix: np.typing.NDArray[np.float64]
    SMACT_version: str    
    
    def PCA_transform(input_data: pd.DataFrame) -> None:
        '''
        Recalculates the PCA transformation using the given data and the number of
        dimensions in the transformed basis. Stores resulting PCA transfomration
        matrix internally.

        Args:
            input_data (pd.DataFrame): The data whose dimensions are to be reduced.
            n_dims (int): The second parameter.

        Returns:
            bool: The return value. True for success, False otherwise.

        
        '''

	def transform_data(data: pd.DataFrame) -> None;
