import pandas as pd
import os
import matplotlib.pyplot as plt

file_path_detailed_data = "Data/gx_details_refseq.20230416.tsv"
file_path_summary_data = "Data/gx_summary_refseq.20230416.tsv"

# Read the TSV file
raw_detailed_data = pd.read_csv(file_path_detailed_data, sep='\t')
raw_summary_data = pd.read_csv(file_path_summary_data, sep='\t')
print(print(f"Rows: {raw_detailed_data.shape[0]}, Columns: {raw_detailed_data.shape[1]}"))

raw_detailed_data = pd.merge(raw_detailed_data , raw_summary_data, on='assembly_accession', how='left')
raw_detailed_data.to_csv(f"{file_path_detailed_data[:-4]}_with_summary_info.csv")
# Filter rows where 'action' column is either 'TRIM' or 'EXCLUDE'
data_without_manual_review = raw_detailed_data[raw_detailed_data['action'].isin(['TRIM', 'EXCLUDE'])]
print(print(f"Rows: {data_without_manual_review.shape[0]}, Columns: {data_without_manual_review.shape[1]}"))

# Further filter rows where 'contam_type' starts with 'prok:'
data_prokaryotes = data_without_manual_review[data_without_manual_review['contam_type'].str.startswith('prok:')]
print(print(f"Rows: {data_prokaryotes.shape[0]}, Columns: {data_prokaryotes.shape[1]}"))

data_prokaryotes['contam_type_clean'] = data_prokaryotes['contam_type'].str.replace('prok:', '')
data_prokaryotes['taxon'] = data_prokaryotes['contam_details'].str.split().str[0]


# Filtering out very short contaminations

data_prokaryotes_long_contaminations = data_prokaryotes[data_prokaryotes['length_contamination'] > 5000]

data_prokaryotes_long_contaminations.to_csv(f"{file_path_detailed_data[:-4]}_long_contaminations_from_prokaryotes.csv")


# Plotting

# Capitalizing labels for the bar chart
contam_counts = data_prokaryotes['contam_type_clean'].value_counts()
contam_counts.index = contam_counts.index.map(lambda x: x.capitalize())

contam_counts_taxon = data_prokaryotes['taxon'].value_counts()
contam_counts_taxon.index = contam_counts_taxon.index.map(lambda x: x.capitalize())
# Bar Chart
plt.figure(figsize=(10, 6))
bar_plot = contam_counts.plot(kind='bar', color='skyblue', edgecolor='black')
plt.title('Bar Chart of Contaminant Type Distribution')
plt.xlabel('Contaminant Types')
plt.ylabel('Count')
plt.savefig('contaminant_type_distribution_bar_chart.png')  # Saving the bar chart
plt.show()

# Histogram
data_prokaryotes['length_contamination'] = data_prokaryotes['end_pos'] - data_prokaryotes['start_pos']
plt.figure(figsize=(10, 6))
plt.hist(data_prokaryotes['length_contamination'], bins=50, color='blue', edgecolor='black', log=True)
plt.title('Histogram of Length of Contaminations (Log Scale)')
plt.xlabel('Length of Contamination')
plt.ylabel('Frequency (Log Scale)')
plt.savefig('length_of_contaminations_histogram_log_scale.png')  # Saving the histogram
plt.show()

# Box Plot
unique_phyla = data_prokaryotes['contam_type_clean'].unique()
capitalized_phyla = [phylum.capitalize() for phylum in unique_phyla]
boxplot_data = [data_prokaryotes[data_prokaryotes['contam_type_clean'] == phylum]['length_contamination'] for phylum in unique_phyla]
data_prokaryotes['length_contamination'] = data_prokaryotes['end_pos'] - data_prokaryotes['start_pos']
plt.figure(figsize=(12, 8))
plt.boxplot(boxplot_data, labels=capitalized_phyla, vert=False)
plt.title('Box Plot of Length of Contamination Across Phyla (Log Scale)')
plt.xscale('log')  # Applying logarithmic scale
plt.xlabel('Length of Contamination (Log Scale)')
plt.ylabel('Phylum')
plt.xticks(rotation=45)  # Rotating x-axis labels to prevent overlap
plt.savefig('length_of_contamination_across_phyla_box_plot.png')  # Saving the box plot
plt.show()

# Bar Chart
plt.figure(figsize=(10, 6))
bar_plot = contam_counts_taxon.plot(kind='bar', color='skyblue', edgecolor='black')
plt.title('Bar Chart of Contaminant Type Distribution')
plt.xlabel('Contaminant Types')
plt.ylabel('Count')
plt.savefig('contaminant_type_distribution_bar_chart.png')  # Saving the bar chart
plt.show()

# Box Plot
unique_phyla = data_prokaryotes['taxon'].unique()
print(unique_phyla)
capitalized_phyla = [phylum.capitalize() for phylum in unique_phyla]
boxplot_data = [data_prokaryotes[data_prokaryotes['taxon'] == phylum]['taxon'] for phylum in unique_phyla]
plt.figure(figsize=(12, 8))
plt.boxplot(boxplot_data, labels=capitalized_phyla, vert=False)
plt.title('Box Plot of Length of Contamination Across Taxa (Log Scale)')
plt.xscale('log')  # Applying logarithmic scale
plt.xlabel('Length of Contamination (Log Scale)')
plt.ylabel('Taxon')
plt.xticks(rotation=45)  # Rotating x-axis labels to prevent overlap
plt.savefig('length_of_contamination_across_taxa_box_plot.png')  # Saving the box plot
plt.show()
