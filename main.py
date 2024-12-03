from gwas import *
import sys

# python3 main.py psoriasis.tsv lung_cancer_NIH.csv
exposure_data, outcome_data = handleExposureOutcomeFiles(sys.argv[1], sys.argv[2])

merged_data = findCommonSNP(exposure_data, outcome_data)

odds_ratio, p_value = inverseVarianceWeightedEstimante(merged_data)
print(f"Odds Ratio (IVW): {odds_ratio}")
print(f"P-Value (IVW): {p_value}")
print("-------------------------------------------------")
egger_odds_ratio, egger_p_value = mrEggerRegression(merged_data)
print(f"MR-Egger Odds Ratio: {egger_odds_ratio}")
print(f"MR-Egger P-Value: {egger_p_value}")
print("------------------------------------------------")

wald_ratios = merged_data['Wald_ratio'].values
weights = 1 / (merged_data['SE_Wald_ratio'] ** 2)
median_causal = weighted_median(wald_ratios, weights)

median_odds_ratio, median_p_value = weightedMedianOddsRatioPValue(wald_ratios, median_causal)
print(f"Weighted Median Odds Ratio: {median_odds_ratio}")
print(f"Weighted Median P-Value: {median_p_value}")
print("---------------------------------------------------")

mode_odds_ratio, mode_p_value = weightedMode(wald_ratios, weights)
print(f"Weighted Mode Odds Ratio: {mode_odds_ratio}")
print(f"Weighted Mode P-Value: {mode_p_value}")
