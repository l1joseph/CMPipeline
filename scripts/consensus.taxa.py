import pandas as pd
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text

# Argument parser
parser = argparse.ArgumentParser(
    description="Compute consensus microbial taxa from MetaPhlAn and Bracken outputs, "
                "filter common taxa, and generate a visualization."
)

parser.add_argument(
    "--metaphlan", required=True, help="Path to the MetaPhlAn merged abundance table (genus level)."
)
parser.add_argument(
    "--bracken_genus", required=True, help="Path to the Bracken genus MPA report."
)
parser.add_argument(
    "--bracken_species", required=True, help="Path to the Bracken species MPA report."
)
parser.add_argument(
    "--output", required=True, help="Path to save the output visualization (PDF)."
)
parser.add_argument(
    "--output_common_genus", required=True,
    help="Path to save the filtered common genus-level taxa report (TXT)."
)
parser.add_argument(
    "--output_common_species", required=True,
    help="Path to save the filtered common species-level taxa report (TXT)."
)

args = parser.parse_args()

# Load MetaPhlAn data
metaphlan_df = pd.read_csv(args.metaphlan, sep="\t")
metaphlan_df['clade_name'] = metaphlan_df['clade_name'].str.split('g__', expand=True)[1]
metaphlan_df = metaphlan_df.loc[~metaphlan_df['clade_name'].str.contains('GGB|_unclassified', na=False)]
metaphlan_df['tool_type'] = 'metaphlan'

# Load Bracken genus data
bracken_genus_df = pd.read_csv(args.bracken_genus, sep="\t")
bracken_genus_df.rename(columns={'#Classification': 'clade_name'}, inplace=True)
bracken_genus_df = bracken_genus_df[~bracken_genus_df['clade_name'].str.contains('Viruses|\|s__', na=False)]
bracken_genus_df = bracken_genus_df.dropna(subset=['clade_name'])
bracken_genus_df['genus_name'] = bracken_genus_df['clade_name'].str.split('\|g__', expand=True)[1]
bracken_genus_df = bracken_genus_df.dropna(subset=['genus_name'])
bracken_genus_df = bracken_genus_df.loc[bracken_genus_df['genus_name'] != 'Cutibacterium']
bracken_genus_df['tool_type'] = 'bracken'

# Load Bracken species data
bracken_species_df = pd.read_csv(args.bracken_species, sep="\t")
bracken_species_df.rename(columns={'#Classification': 'clade_name'}, inplace=True)
bracken_species_df = bracken_species_df[~bracken_species_df['clade_name'].str.contains('Viruses', na=False)]
bracken_species_df = bracken_species_df.dropna(subset=['clade_name'])
bracken_species_df['genus_name'] = bracken_species_df['clade_name'].str.split('\|g__', expand=True)[1].str.split('\|s__', expand=True)[0]
bracken_species_df['species_name'] = bracken_species_df['clade_name'].str.split('\|g__', expand=True)[1].str.split('\|s__', expand=True)[1]
bracken_species_df = bracken_species_df.dropna(subset=['species_name'])
bracken_species_df = bracken_species_df.loc[bracken_species_df['genus_name'] != 'Cutibacterium']
bracken_species_df['tool_type'] = 'bracken'

# Merge Data
common_columns = list(set(bracken_genus_df.columns).intersection(metaphlan_df.columns))
merged_df = pd.merge(bracken_genus_df, metaphlan_df, on=common_columns, how='outer').fillna(0.0)

# Calculate Proportions
proportions = []

for clade_name, group in merged_df.groupby('clade_name'):
    proportion_bracken = group[group['tool_type'] == 'bracken'].gt(0).sum(axis=1).mean()
    proportion_metaphlan = group[group['tool_type'] == 'metaplan'].gt(0).sum(axis=1).mean()
    proportions.append([clade_name, proportion_bracken, proportion_metaphlan])

proportions_df = pd.DataFrame(proportions, columns=['clade_name', 'Bracken(tumor)', 'Metaphlan(tumor)'])
taxa_common = proportions_df[(proportions_df['Bracken(tumor)'] >= 0.01) & (proportions_df['Metaphlan(tumor)'] >= 0.01)]

# Filter Bracken genus/species using common clades
common_clades = taxa_common['clade_name'].tolist()
bracken_genus_df = bracken_genus_df[bracken_genus_df['genus_name'].isin(common_clades)]
bracken_species_df = bracken_species_df[bracken_species_df['genus_name'].isin(common_clades)]

# Save filtered Bracken genus/species data
bracken_genus_df.to_csv(args.output_common_genus, sep="\t", index=False)
bracken_species_df.to_csv(args.output_common_species, sep="\t", index=False)

# Plot Results
plt.figure(figsize=(6, 6))
sns.scatterplot(data=proportions_df, x='Bracken(tumor)', y='Metaphlan(tumor)', s=60, alpha=0.7, color='lightgray', edgecolor='k', linewidth=0.5)
sns.scatterplot(data=taxa_common, x='Bracken(tumor)', y='Metaphlan(tumor)', s=60, color='firebrick', edgecolor='k', linewidth=0.5)

texts = [plt.text(row['Bracken(tumor)'], row['Metaphlan(tumor)'], row['clade_name'], color='firebrick', fontsize=7) 
         for _, row in taxa_common.iterrows()]
adjust_text(texts, arrowprops=dict(arrowstyle="->", color='k', lw=1))

plt.axvline(x=0.01, color='gray', linestyle='--', linewidth=0.7)
plt.axhline(y=0.01, color='gray', linestyle='--', linewidth=0.7)
plt.xlabel('Proportion of samples\n Bracken(tumor)', fontsize=12)
plt.ylabel('Proportion of samples\n MetaPhlan(tumor)', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.4)
plt.tight_layout()
plt.savefig(args.output, dpi=300)
plt.show()
