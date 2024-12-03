# !/bin/bash
echo "Lung Cancer Stats"
python3 main.py psoriasis.tsv lung_cancer_NIH.tsv
echo "Kidney Cancer Stats"
python3 main.py psoriasis.tsv kidney_second.tsv