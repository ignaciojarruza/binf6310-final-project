import pandas as pd
import sys

def handleExposureOutcomeFiles(exposure_file: str, outcome_file: str):
    try:
        # Load exposure and outcome datasets
        exposure_data = pd.read_csv(exposure_file, sep="\t")
        outcome_data = pd.read_csv(outcome_file, sep=",")
        return exposure_data, outcome_data
    except FileNotFoundError:
        print(f"Exposure or Outcome file not found.")
        sys.exit(1)

def findCommonSNP(exposure_data, outcome_data):
    try:
        # Convert SNPs to string type for reliable comparison
        exposure_data['SNPS'] = exposure_data['SNPS'].astype(str)
        outcome_data['snp'] = outcome_data['snp'].astype(str)

        # Find common SNPs
        common_snps = set(exposure_data['SNPS']).intersection(set(outcome_data['snp']))
        print(f"Number of common SNPs: {len(common_snps)}")

        # Filter data for common SNPs
        filtered_exposure = exposure_data[exposure_data['SNPS'].isin(common_snps)]
        filtered_outcome = outcome_data[outcome_data['snp'].isin(common_snps)]
    
        # Merge datasets on SNPs
        merged_data = pd.merge(
            filtered_exposure[['SNPS', 'OR or BETA', 'P-VALUE']],
            filtered_outcome[['snp', 'beta', 'p_value']],
            left_on='SNPS',
            right_on='snp',
            suffixes=('_exposure', '_outcome')
        )

        # Rename columns for clarity
        merged_data.rename(columns={
            'OR or BETA': 'Beta_exposure',
            'beta': 'Beta_outcome',
            'P-VALUE': 'P_value_exposure',
            'p_value': 'P_value_outcome'
        }, inplace=True)
        
         # Drop rows with missing values in the required columns
        merged_data.dropna(subset=['Beta_exposure', 'Beta_outcome'], inplace=True)

        return merged_data
    except KeyError:
        print(f"SNP data not available in input data.")
        sys.exit(1)



    

   