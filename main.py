from gwas import *
import sys

# python3 main.py psoriasis.tsv lung_cancer_NIH.csv
exposure_data, outcome_data = handleExposureOutcomeFiles(sys.argv[1], sys.argv[2])

merged_data = findCommonSNP(exposure_data, outcome_data)

inverseVarianceWeightedEstimante(merged_data)
mrEggerRegression(merged_data)

wald_ratios = merged_data['Wald_ratio'].values
weights = 1 / (merged_data['SE_Wald_ratio'] ** 2)
median_causal = weighted_median(wald_ratios, weights)

weightedMedianOddsRatioPValue(wald_ratios, median_causal)
weightedMode(wald_ratios, weights)

