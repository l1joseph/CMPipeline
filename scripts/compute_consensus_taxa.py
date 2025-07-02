import pandas as pd
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text
import sys

# Argument parser
parser = argparse.ArgumentParser(
    description="Compute consensus microbial taxa from MetaPhlAn and Bracken outputs, "
                "filter common taxa based on samples with at least 100,000 reads, and generate a visualization."
)

parser.add_argument("--metaphlan", required=True, help="Path to the MetaPhlAn merged abundance table (genus level).")
parser.add_argument("--bracken_genus", required=True, help="Path to the Bracken genus MPA report.")
parser.add_argument("--bracken_species", required=True, help="Path to the Bracken species MPA report.")
parser.add_argument("--output", required=True, help="Path to save the output visualization (PDF).")
parser.add_argument("--output_common_genus",required=True, help="Path to save filtered common genus report (TXT).")
parser.add_argument("--output_common_species", required=True, help="Path to save filtered common species report (TXT).")

args = parser.parse_args()


# Load MetaPhlAn data
metaphlan_df = pd.read_csv(args.metaphlan, sep="\t",header=1)
metaphlan_df['clade_name'] = metaphlan_df['clade_name'].str.split('g__', expand=True)[1]
metaphlan_df = metaphlan_df.loc[~metaphlan_df['clade_name'].str.contains('GGB|_unclassified', na=False)]
metaphlan_df['tool_type'] = 'metaphlan'


# Load Bracken genus data
bracken_genus_df = pd.read_csv(args.bracken_genus, sep="\t")
bracken_genus_df.rename(columns={'#Classification': 'clade_name'}, inplace=True)
bracken_genus_df = bracken_genus_df[~bracken_genus_df['clade_name'].str.contains('Viruses|\|s__', na=False)]
bracken_genus_df.columns = bracken_genus_df.columns.str.replace(r'\.bracken\.G\.krakenreport\.txt$', '', regex=True)
sample_names  = list(bracken_genus_df.columns.drop('clade_name'))
bracken_genus_df = bracken_genus_df.dropna(subset=['clade_name'])
bracken_genus_df['genus_name'] = bracken_genus_df['clade_name'].str.split('\|g__', expand=True)[1]
bracken_genus_df = bracken_genus_df.dropna(subset=['genus_name'])
bracken_genus_df = bracken_genus_df.loc[bracken_genus_df['genus_name'] != 'Cutibacterium']
bracken_genus_df['tool_type'] = 'bracken'
sample_sums = bracken_genus_df[sample_names].sum()
samples_above_threshold = sample_sums[sample_sums >= 100000].index.tolist()

if len(samples_above_threshold) < 2:
    print("Consensus not possible: fewer than 2 samples have at least 100,000 reads.")
    sys.exit(0)


# Load Bracken species data
bracken_species_df = pd.read_csv(args.bracken_species, sep="\t")
bracken_species_df.rename(columns={'#Classification': 'clade_name'}, inplace=True)
bracken_species_df = bracken_species_df[~bracken_species_df['clade_name'].str.contains('Viruses', na=False)]
bracken_species_df.columns = bracken_species_df.columns.str.replace(r'\.bracken\.S\.krakenreport\.txt$', '', regex=True)
bracken_species_df = bracken_species_df.dropna(subset=['clade_name'])
bracken_species_df['genus_name'] = bracken_species_df['clade_name'].str.split('\|g__', expand=True)[1].str.split('\|s__', expand=True)[0]
bracken_species_df['species_name'] = bracken_species_df['clade_name'].str.split('\|g__', expand=True)[1].str.split('\|s__', expand=True)[1]
bracken_species_df = bracken_species_df.dropna(subset=['species_name'])
bracken_species_df = bracken_species_df.loc[bracken_species_df['genus_name'] != 'Cutibacterium']
bracken_species_df['tool_type'] = 'bracken'

# Filter Bracken genus/species to only samples above threshold for consensus calculation
bracken_genus_filtered = bracken_genus_df[['clade_name','genus_name', 'tool_type'] + samples_above_threshold]
bracken_species_filtered = bracken_species_df[['clade_name','genus_name', 'tool_type'] +samples_above_threshold]

# Merge for consensus
common_columns = list(set(bracken_genus_filtered.columns).intersection(metaphlan_df.columns))
merged_df = pd.merge(bracken_genus_filtered, metaphlan_df, on=common_columns, how='outer').fillna(0.0)


# Calculate Proportions
proportions = []
total_samples = len(samples_above_threshold)

for clade_name, group in merged_df.loc[:,samples_above_threshold+['clade_name', 'tool_type']].groupby('clade_name'):
    proportion_tool1 = group[group['tool_type'] == 'bracken'][samples_above_threshold].apply(lambda x: (int(x) > 0).sum() / len(x), axis=1).mean()
    proportion_tool2 = group[group['tool_type'] == 'metaplan'][samples_above_threshold].apply(lambda x: (int(x) > 0).sum() / len(x), axis=1).mean()
    proportions.append([clade_name, proportion_tool1, proportion_tool2])

# Create a DataFrame from proportions list
proportions_df = pd.DataFrame(proportions, columns=['clade_name', 'Bracken', 'Metaphlan'])
proportions_df.fillna(0, inplace=True)
common_taxa = proportions_df[(proportions_df['Bracken'] >= 0.01) & (proportions_df['Metaphlan'] >= 0.01)]


# Filter all samples based on common clades (retain all samples, not just filtered ones)
common_clades = common_taxa['clade_name'].tolist()
bracken_genus_df = bracken_genus_df[bracken_genus_df['genus_name'].isin(common_clades)]
bracken_species_df = bracken_species_df[bracken_species_df['genus_name'].isin(common_clades)]

# Save filtered Bracken genus/species data
bracken_genus_df.to_csv(args.output_common_genus, sep="\t", index=False)
bracken_species_df.to_csv(args.output_common_species, sep="\t", index=False)

# Plot Results
plt.figure(figsize=(6, 6))
sns.scatterplot(data=proportions_df, x='Bracken', y='Metaphlan', s=60, alpha=0.7, color='lightgray', edgecolor='k', linewidth=0.5)
sns.scatterplot(data=common_taxa, x='Bracken', y='Metaphlan', s=60, color='firebrick', edgecolor='k', linewidth=0.5)

texts = [plt.text(row['Bracken'], row['Metaphlan'], row['clade_name'], color='firebrick', fontsize=7) 
         for _, row in common_taxa.iterrows()]
adjust_text(texts, arrowprops=dict(arrowstyle="->", color='k', lw=1))

plt.axvline(x=0.01, color='gray', linestyle='--', linewidth=0.7)
plt.axhline(y=0.01, color='gray', linestyle='--', linewidth=0.7)
plt.xlabel('Proportion of samples\n Bracken', fontsize=12)
plt.ylabel('Proportion of samples\n MetaPhlan', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.4)
plt.tight_layout()
plt.savefig(args.output, dpi=300)
plt.show()
