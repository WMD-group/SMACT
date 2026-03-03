# Roost-MP-2024.12.0-band_gap

## Aim

Predict the DFT-computed electronic band gap of inorganic materials from composition alone, using the ROOST (Representation Learning from Stoichiometry) architecture.

## Training Dataset

- **Source**: Materials Project (MP), accessed December 2024
- **Size**: 103,644 compositions (82,914 train / 10,365 val / 10,365 test)
- **Target**: PBE band gap (eV)
- **Split**: Random 80/10/10

## Model Details

| Parameter                           | Value                                    |
| ----------------------------------- | ---------------------------------------- |
| Architecture                        | ROOST (composition-only message passing) |
| Element embedding                   | matscholar200 (200-dim)                  |
| Element feature length              | 64                                       |
| Graph layers                        | 3                                        |
| Attention heads (element)           | 3                                        |
| Gate / message hidden (element)     | [256] / [256]                            |
| Attention heads (composition)       | 3                                        |
| Gate / message hidden (composition) | [256] / [256]                            |
| Trunk hidden layers                 | [1024, 512]                              |
| Output hidden layers                | [256, 128, 64]                           |
| Loss                                | Heteroscedastic (robust)                 |
| Epochs trained                      | 232                                      |

## Performance

| Metric    | Value    |
| --------- | -------- |
| Test MAE  | 0.283 eV |
| Test RMSE | 0.599 eV |
| Val MAE   | 0.282 eV |

## Usage

```python
from smact.property_prediction import predict_band_gap

result = predict_band_gap("NaCl")
print(f"Band gap: {result[0]:.2f} eV")
```

## References

- Goodall, R.E.A. & Lee, A.A. "Predicting materials properties without crystal structure: deep representation learning from stoichiometry." _Nature Communications_ 11, 6280 (2020). [doi:10.1038/s41467-020-19964-7](https://doi.org/10.1038/s41467-020-19964-7)
- Jain, A. et al. "Commentary: The Materials Project: A materials genome approach to accelerating materials innovation." _APL Materials_ 1, 011002 (2013). [doi:10.1063/1.4812323](https://doi.org/10.1063/1.4812323)
