# !/bin/bash
echo "Lung Cancer Stats"
python3 main.py psoriasis.tsv lung_cancer_NIH.tsv
echo "Breast Cancer Stats"
python3 main.py psoriasis.tsv breast_cancer.tsv
echo "Kidney Cancer Stats"
python3 main.py psoriasis.tsv kidney_second.tsv