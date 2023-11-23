
import matplotlib.pyplot as plt
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
 
