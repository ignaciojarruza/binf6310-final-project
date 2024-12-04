# binf6310-final-project

## Pipeline Recreation

This is a recreation of the Two-sample MR causality analysis performed in the following research publication: https://doi.org/10.1038/s41467-024-51824-6
The computational methods performed in this recreation are:

- MR-Egger
- Weighted Median
- Inverse variance weighted
- Weighted mode

## Requirements

To install dependencies, run the following command:

```
pip3 install -r requirements.txt
```

## How to run

To run bash script that runs all cancer types with the appropriate inputs, run the following command:

```
bash run.sh
```

To run a specific cancer type, run the following command:

```
python3 main.py psoriasis [cancer type dataset filepath]
```

For example:

```
python3 main.py psoriasis datasets/breast_cancer.tsv
```

## Datasets

The data for this pipeline recreation is all included in the datasets and psoriasis folders.

## Database References

References:

- Burdett, T. (2024a). GWAS Catalog. Ebi.ac.uk. https://www.ebi.ac.uk/gwas/efotraits/MONDO_0004647
- Burdett, T. (2024b). GWAS Catalog. Ebi.ac.uk. https://www.ebi.ac.uk/gwas/efotraits/MONDO_0008903
- Burdett, T. (2024c). GWAS Catalog. Ebi.ac.uk. https://www.ebi.ac.uk/gwas/efotraits/EFO_0009260
- Burdett, T. (2024d). GWAS Catalog. Ebi.ac.uk. https://www.ebi.ac.uk/gwas/efotraits/MONDO_0007254
- Burdett, T. (2024e). GWAS Catalog. Ebi.ac.uk. https://www.ebi.ac.uk/gwas/efotraits/MONDO_0002367
- Burdett, T. (2024f). GWAS Catalog. Ebi.ac.uk. https://www.ebi.ac.uk/gwas/efotraits/MONDO_0018906
