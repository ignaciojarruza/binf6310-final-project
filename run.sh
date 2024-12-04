# !/bin/bash
# https://www.ebi.ac.uk/gwas/efotraits/MONDO_0004647
echo "Anal canal cancer"
echo "___________________"
python3 main.py psoriasis datasets/in_situ_carcinoma.tsv

# https://www.ebi.ac.uk/gwas/efotraits/MONDO_0008903
echo "Lung Cancer"
echo "___________________"
python3 main.py psoriasis datasets/lung_cancer.tsv

# https://www.ebi.ac.uk/gwas/efotraits/EFO_0009260
echo "Non-Melanoma skin cancer"
echo "___________________"
python3 main.py psoriasis datasets/non_melanoma_skin_cancer.tsv

# https://www.ebi.ac.uk/gwas/efotraits/MONDO_0007254
echo "Breast Cancer Stats"
echo "___________________"
python3 main.py psoriasis datasets/breast_cancer.tsv

# https://www.ebi.ac.uk/gwas/efotraits/MONDO_0002367
echo "Kidney Cancer Stats"
echo "___________________"
python3 main.py psoriasis datasets/kidney_cancer.tsv

# https://www.ebi.ac.uk/gwas/efotraits/MONDO_0018906
echo "follicular non-hodgkins lymphoma"
echo "___________________"
python3 main.py psoriasis datasets/follicular_non_hodgkins.tsv