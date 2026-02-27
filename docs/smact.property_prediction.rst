Property Prediction Module
===========================

The property prediction module enables prediction of material properties
directly from chemical composition, without requiring crystal structure
information. It uses pre-trained `ROOST <https://doi.org/10.1038/s41467-020-19964-7>`_
(Representation Learning from Stoichiometry) models to provide fast
inference from stoichiometry alone.

The typical workflow is:

#. Use the convenience function :func:`~smact.property_prediction.convenience.predict_band_gap`
   for quick predictions:

   .. code-block:: python

      from smact.property_prediction import predict_band_gap

      predictions = predict_band_gap(["NaCl", "TiO2", "GaN"])

#. Or use the :class:`~smact.property_prediction.roost.predictor.RoostPropertyPredictor`
   class directly for more control:

   .. code-block:: python

      from smact.property_prediction import RoostPropertyPredictor

      predictor = RoostPropertyPredictor(property_name="band_gap")
      result = predictor.predict(["GaN", "ZnO"], return_uncertainty=True)

#. Query the model registry for available properties and models:

   .. code-block:: python

      from smact.property_prediction import get_supported_properties, get_available_models

      print(get_supported_properties())  # ['band_gap']
      print(get_available_models())  # ['Roost-MP-2024.12.0-band_gap']

Installation
------------

The property prediction module requires additional dependencies that are
not installed by default. Install them with:

.. code-block:: bash

   pip install smact[property_prediction]

This installs `aviary-models <https://pypi.org/project/aviary-models/>`_ (>=1.2.0)
and `PyTorch <https://pytorch.org/>`_ (>=2.0.0).

Architecture
------------

The module follows an extensible abstract base class pattern, making it
straightforward to add new predictor backends in the future.

.. code-block:: text

   BasePropertyPredictor (ABC)
   ├── predict()             # Abstract — must be implemented
   ├── supported_properties  # Abstract property
   └── _validate_compositions()  # Shared composition validation
        │
        └── RoostPropertyPredictor  # Concrete implementation
             ├── predict()          # ROOST inference with uncertainty
             └── from_checkpoint()  # Load from saved checkpoint

**PredictionResult** is a dataclass returned when ``return_uncertainty=True``.
It carries ``predictions``, ``uncertainties``, ``compositions``, and ``unit``.

The design is inspired by `MatGL <https://github.com/materialsvirtuallab/matgl>`_'s
pretrained model management, following AntObi's recommendation for model versioning
with per-model directories and READMEs.

How ROOST Works
~~~~~~~~~~~~~~~

ROOST (Representation Learning from Stoichiometry) predicts material properties
from chemical composition alone, without crystal structure information.

.. code-block:: text

   Composition "NaCl"
       │
       ▼
   ┌─────────────────────────────────────────┐
   │ CompositionData                         │
   │ Parse formula → elements + weights      │
   │ Na: 0.5, Cl: 0.5 → element embeddings  │
   └─────────────────────────────────────────┘
       │
       ▼
   ┌─────────────────────────────────────────┐
   │ Message Passing (3 layers)              │
   │ Elements exchange information via       │
   │ attention-weighted graph neural network │
   └─────────────────────────────────────────┘
       │
       ▼
   ┌─────────────────────────────────────────┐
   │ Trunk + Output Networks                 │
   │ [1024, 512] → [256, 128, 64] → output  │
   │ If robust: also outputs log(std)        │
   └─────────────────────────────────────────┘
       │
       ▼
   Prediction (+ uncertainty if robust=True)

Key components:

- **Element embeddings**: 200-dimensional MatScholar embeddings capture chemical similarity.
- **Message passing**: Elements in a composition exchange information through 3 graph neural network layers with multi-head attention (3 heads).
- **Attention pooling**: Weighted aggregation of element features into a single composition vector.
- **Robust loss**: Heteroscedastic loss enables per-prediction uncertainty estimates.

Available Models
----------------

.. list-table::
   :header-rows: 1

   * - Model
     - Property
     - Dataset
     - Test MAE
   * - Roost-MP-2024.12.0-band_gap
     - Band gap (eV)
     - Materials Project 103K
     - 0.283 eV

Roost-MP-2024.12.0-band_gap
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The default band gap model was trained on 103,644 DFT (GGA/GGA+U)
band gap calculations from the Materials Project.

**Training details:**

.. list-table::
   :widths: 30 70

   * - Dataset
     - Materials Project (103,644 compositions)
   * - Split
     - 82,914 train / 10,365 val / 10,365 test (80/10/10)
   * - Epochs trained
     - 232 (early stopping)
   * - Learning rate
     - 3e-4 with cosine annealing
   * - Loss
     - Heteroscedastic (robust) L1

**Performance:**

.. list-table::
   :header-rows: 1

   * - Metric
     - Value
   * - Test MAE
     - 0.283 eV
   * - Test RMSE
     - 0.599 eV
   * - Val MAE
     - 0.282 eV

**Baseline comparison** (MP test set):

.. list-table::
   :header-rows: 1

   * - Model
     - MAE (eV)
   * - Mean baseline
     - 1.361
   * - Median baseline
     - 1.251
   * - Linear + Magpie
     - 0.862
   * - Random Forest + Magpie
     - 0.339
   * - **ROOST**
     - **0.283**

ROOST beats RF+Magpie by 17% (0.339 → 0.283 eV), consistent with
improvements reported in the original paper.

Known Limitations
~~~~~~~~~~~~~~~~~

These are fundamental limitations of composition-only prediction,
not implementation issues:

- **Polymorph ambiguity**: The same composition can be metallic in one
  crystal structure and semiconducting in another (e.g. Fe₃O₄ magnetite
  vs maghemite). ROOST predicts an "average" behaviour across polymorphs.
- **Metals at zero gap**: The MP dataset contains 42% metals (gap=0).
  Materials with metastable metallic phases may be predicted with a
  non-zero gap if the semiconducting polymorph is more common in the
  training data.
- **Extremes**: Both very low (<1 eV) and very high (>7 eV) band gap
  materials show more scatter, likely due to fewer training examples.

Model Versioning and Fidelity
-----------------------------

The module supports model versioning following the naming convention:

.. code-block:: text

   Roost-<dataset>-<version>-<property>[-<fidelity>]

For example: ``Roost-MP-2024.12.0-band_gap``, or in future
``Roost-MP-2025.6.0-band_gap-hse06``.

The config supports multiple fidelities per property:

.. code-block:: python

   # Current (single model):
   DEFAULT_MODELS = {
       "band_gap": {
           "default": "Roost-MP-2024.12.0-band_gap",
       },
   }

   # Future expansion example:
   DEFAULT_MODELS = {
       "band_gap": {
           "default": "pbe",
           "pbe": "Roost-MP-2025.6.0-band_gap-pbe",
           "hse06": "Roost-MP-2025.6.0-band_gap-hse06",
       },
   }

Fidelity can be selected at construction time:

.. code-block:: python

   # Use default fidelity
   predictor = RoostPropertyPredictor(property_name="band_gap")

   # Select specific fidelity (when multiple are available)
   predictor = RoostPropertyPredictor(property_name="band_gap", fidelity="hse06")

   # Or specify a model directly
   predictor = RoostPropertyPredictor(model_name="Roost-MP-2024.12.0-band_gap")

Model Resolution Order
----------------------

When loading a model by name, SMACT checks the following locations in order:

#. **Local directory path** -- if an explicit path is provided.
#. **pretrained_models/** -- models bundled in the repository (source installs only).
#. **Cache** ``~/.cache/smact/models/`` -- previously downloaded models.
#. **Remote download** -- fetched from GitHub releases on first use.

For pip installs (where ``pretrained_models/`` is not included), models are
downloaded automatically on first use and cached locally.

Model Management
----------------

.. code-block:: python

   from smact.property_prediction import (
       get_supported_properties,
       get_property_unit,
       get_available_models,
       list_cached_models,
       clear_cache,
   )

   # List available properties
   get_supported_properties()  # ['band_gap']

   # Get property unit
   get_property_unit("band_gap")  # 'eV'

   # List all discoverable models
   get_available_models()  # ['Roost-MP-2024.12.0-band_gap']

   # List models in local cache
   list_cached_models()

   # Clear the download cache
   clear_cache()

Training New Models
-------------------

The module includes a full training pipeline for ROOST models.
Training requires a CSV file with a ``composition`` column and a
numeric target column.

**Using the Python API:**

.. code-block:: python

   from smact.property_prediction.roost.train import train_roost_model

   model_path = train_roost_model(
       data_path="band_gap_data.csv",
       target_column="band_gap",
       property_name="band_gap",
       output_dir="./trained_models",
       epochs=100,
       batch_size=128,
       learning_rate=3e-4,
       robust=True,  # Enable uncertainty estimation
       device="cuda",  # Use "cpu" if no GPU
   )

**Using the command line:**

.. code-block:: bash

   python -m smact.property_prediction.roost.train \
       --data-path band_gap_data.csv \
       --target band_gap \
       --property band_gap \
       --output-dir ./trained_models \
       --epochs 100 \
       --batch-size 128 \
       --robust \
       --device cuda

**Loading a custom trained model:**

.. code-block:: python

   from smact.property_prediction import RoostPropertyPredictor

   predictor = RoostPropertyPredictor(
       property_name="band_gap",
       model_path="./trained_models/Roost-Custom-band_gap",
   )
   predictions = predictor.predict(["NaCl", "GaN", "TiO2"])

References
----------

- Goodall, R.E.A. & Lee, A.A. "Predicting materials properties without crystal structure: deep representation learning from stoichiometry." *Nature Communications* 11, 6280 (2020). `doi:10.1038/s41467-020-19964-7 <https://doi.org/10.1038/s41467-020-19964-7>`_

- Jain, A. et al. "Commentary: The Materials Project: A materials genome approach to accelerating materials innovation." *APL Materials* 1, 011002 (2013). `doi:10.1063/1.4812323 <https://doi.org/10.1063/1.4812323>`_

Submodules
----------

.. toctree::

    smact.property_prediction.base_predictor
    smact.property_prediction.convenience
    smact.property_prediction.config
    smact.property_prediction.io
    smact.property_prediction.registry
    smact.property_prediction.roost
