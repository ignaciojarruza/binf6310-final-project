from gwas import *
import sys
from tabulate import tabulate

# Input Handling
exposure_data, outcome_data = handleExposureOutcomeFiles(sys.argv[1], sys.argv[2])
merged_data = findCommonSNP(exposure_data, outcome_data)
merged_data = calculateWaldRatios(merged_data)

# Method Calculations
odds_ratio, p_value = inverseVarianceWeightedEstimate(merged_data)
egger_odds_ratio, egger_p_value = mrEggerRegression(merged_data)

wald_ratios = merged_data['Wald_ratio'].values
weights = 1 / (merged_data['SE_Wald_ratio'] ** 2)
median_causal = weighted_median(wald_ratios, weights)

median_odds_ratio, median_p_value = weightedMedianOddsRatioPValue(wald_ratios, median_causal)
mode_odds_ratio, mode_p_value = weightedMode(wald_ratios, weights)

# Tabulate Output Printing
methods = ["MR-Egger", "Weighted Median", "Inverse variance weighted", "Weighted Mode"]
OR = [egger_odds_ratio, median_odds_ratio, odds_ratio, mode_odds_ratio]
p_values = [egger_p_value, median_p_value, p_value, mode_p_value]
table = list(zip(methods, OR, p_values))
print(tabulate(table, headers=["Methods", "OR (95% CI)", "p value"], tablefmt="pretty"))