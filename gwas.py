import pandas as pd
import sys
import numpy as np
from scipy.stats import norm
from sklearn.linear_model import LinearRegression
from scipy.stats import gaussian_kde

def handleExposureOutcomeFiles(exposure_file: str, outcome_file: str):
    try:
        # Load exposure and outcome datasets
        exposure_data = pd.read_csv(exposure_file, sep="\t")
        outcome_data = pd.read_csv(outcome_file, sep="\t")
        return exposure_data, outcome_data
    except FileNotFoundError:
        print(f"Exposure or Outcome file not found.")
        sys.exit(1)

def findCommonSNP(exposure_data, outcome_data):
    try:
        # Convert SNPs to string type for reliable comparison
        exposure_data['SNPS'] = exposure_data['SNPS'].astype(str)
        outcome_data['SNPS'] = outcome_data['SNPS'].astype(str)

        # Find common SNPs
        common_snps = set(exposure_data['SNPS']).intersection(set(outcome_data['SNPS']))
        print(f"Number of common SNPs: {len(common_snps)}")

        # Filter data for common SNPs
        filtered_exposure = exposure_data[exposure_data['SNPS'].isin(common_snps)]
        filtered_outcome = outcome_data[outcome_data['SNPS'].isin(common_snps)]
    
        # Merge datasets on SNPs
        merged_data = pd.merge(
            filtered_exposure[['SNPS', 'OR or BETA', 'P-VALUE']],
            filtered_outcome[['SNPS', 'OR or BETA', 'P-VALUE']],
            on='SNPS',
            suffixes=('_exposure', '_outcome')
        )

        # Rename columns for clarity
        merged_data.rename(columns={
        'OR or BETA_exposure': 'Beta_exposure',
        'OR or BETA_outcome': 'Beta_outcome',
        'P-VALUE_exposure': 'P_value_exposure',
        'P-VALUE_outcome': 'P_value_outcome'
        }, inplace=True)

         # Drop rows with missing values in the required columns
        merged_data.dropna(subset=['Beta_exposure', 'Beta_outcome'], inplace=True)

        return merged_data
    except KeyError:
        print(f"SNP data not available in input data.")
        sys.exit(1)

def inverseVarianceWeightedEstimante(merged_data):
    # Calculate Wald ratios
    merged_data['Wald_ratio'] = merged_data['Beta_outcome'] / merged_data['Beta_exposure']

    # Calculate standard error of Wald ratios
    merged_data['SE_Wald_ratio'] = np.sqrt(
        (merged_data['P_value_exposure'] ** -1) + (merged_data['P_value_outcome'] ** -1)
    )

    # Inverse Variance Weighted (IVW) Method
    weights = 1 / (merged_data['SE_Wald_ratio'] ** 2)
    ivw_estimate = np.sum(weights * merged_data['Wald_ratio']) / np.sum(weights)

    print(f"Inverse Variance Weighted Estimate: {ivw_estimate}")
    ivw_se = np.sqrt(1 / np.sum(weights))
    odds_ratio = np.exp(ivw_estimate)
    z_score = ivw_estimate / ivw_se
    p_value = 2 * (1 - norm.cdf(abs(z_score)))

    print(f"Odds Ratio (IVW): {odds_ratio}")
    print(f"P-Value (IVW): {p_value}")
    print("-------------------------------------------------")
    return odds_ratio, p_value

def mrEggerRegression(merged_data):
    # MR-Egger Regression
    X = merged_data['Beta_exposure'].values.reshape(-1, 1)
    y = merged_data['Beta_outcome'].values
    weights = 1 / (merged_data['SE_Wald_ratio'] ** 2)

    reg = LinearRegression()
    reg.fit(X, y, sample_weight=weights)

    egger_causal = reg.coef_[0]
    egger_intercept = reg.intercept_

    egger_odds_ratio = np.exp(egger_causal)
    egger_se = np.sqrt(1 / np.sum(weights))
    egger_z_score = egger_causal / egger_se
    egger_p_value = 2 * (1 - norm.cdf(abs(egger_z_score)))

    print(f"MR-Egger Odds Ratio: {egger_odds_ratio}")
    print(f"MR-Egger P-Value: {egger_p_value}")
    print("------------------------------------------------")
    return egger_odds_ratio, egger_p_value

# Weighted Median
def weighted_median(values, weights):
    sorted_indices = np.argsort(values)
    values = np.array(values)[sorted_indices]
    weights = np.array(weights)[sorted_indices]
    cumulative_weight = np.cumsum(weights) - 0.5 * weights
    total_weight = np.sum(weights)
    median_idx = np.where(cumulative_weight >= 0.5 * total_weight)[0][0]
    return values[median_idx]

def weightedMedianOddsRatioPValue(wald_ratios, median_causal):
    median_odds_ratio = np.exp(median_causal)
    median_se = np.std(wald_ratios)
    median_z_score = median_causal / median_se if median_se != 0 else 0
    median_p_value = 2 * (1 - norm.cdf(abs(median_z_score)))

    print(f"Weighted Median Odds Ratio: {median_odds_ratio}")
    print(f"Weighted Median P-Value: {median_p_value}")
    print("---------------------------------------------------")
    return median_odds_ratio, median_p_value

def weightedMode(wald_ratios, weights):
    # Weighted Mode
    if len(wald_ratios) > 1:
        kde = gaussian_kde(wald_ratios, weights=weights)
        x_grid = np.linspace(min(wald_ratios), max(wald_ratios), 1000)
        kde_values = kde(x_grid)
        mode_causal = x_grid[np.argmax(kde_values)]
    else:
        print("Error: Weighted Mode calculation requires multiple Wald ratios. Using median as fallback.")
        mode_causal = np.median(wald_ratios) if len(wald_ratios) > 0 else np.nan

    weights_array = np.array(weights)
    mode_odds_ratio = np.exp(mode_causal) if not np.isnan(mode_causal) else np.nan
    mode_se = np.std(wald_ratios) if len(wald_ratios) > 1 else (1 / weights_array[0] if len(wald_ratios) == 1 else np.nan)
    mode_z_score = mode_causal / mode_se if mode_se != 0 else 0
    mode_p_value = 2 * (1 - norm.cdf(abs(mode_z_score))) if not np.isnan(mode_causal) else np.nan

    print(f"Weighted Mode Odds Ratio: {mode_odds_ratio}")
    print(f"Weighted Mode P-Value: {mode_p_value}")
    return mode_odds_ratio, mode_p_value




    

   