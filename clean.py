import pandas as pd
import sys

# Read the TSV file into a DataFrame
df = pd.read_csv(sys.argv[1], sep="\t")

# Rename the columns
column_mapping = {
    'snp': 'SNPS',
    'beta': 'OR or BETA',
    'p_value': 'P-VALUE'
}
df.rename(columns=column_mapping, inplace=True)

# Save the updated DataFrame back to a TSV file
df.to_csv(sys.argv[1], sep="\t", index=False)

print(f"Columns renamed and saved to {sys.argv[1]}")
