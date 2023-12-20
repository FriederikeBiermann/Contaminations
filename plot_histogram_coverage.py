import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

data = pd.read_csv("Data/gx_details_genbank.20230416_long_contaminations_from_prokaryotes.csv")

# Function to calculate the number of entries above each threshold
def calculate_threshold_data(data, column_name):
    thresholds = np.linspace(data[column_name].min(), data[column_name].max(), 100)
    counts = [data[data[column_name] >= t].shape[0] for t in thresholds]
    return pd.DataFrame({'threshold': thresholds, 'counts': counts})

# Calculate threshold data
threshold_data = calculate_threshold_data(data, 'coverage')

fig, axs = plt.subplots(1, 2, figsize=(15, 6))

# Histogram
axs[0].hist(data['coverage'], bins=30, color='#6982B5')
axs[0].set_title('Coverage of contamination')
axs[0].set_xlabel('Coverage')
axs[0].set_ylabel('Frequency')

# Threshold curve with the new color #A96073
axs[1].plot(threshold_data['threshold'], threshold_data['counts'], color='#A96073')
axs[1].set_title('Entries above Threshold')
axs[1].set_xlabel('Coverage Threshold')
axs[1].set_ylabel('Number of Entries')

plt.tight_layout()
plt.savefig("gx_details_genbank.20230416_long_contaminations_from_prokaryotes_coverage_histogram_coverage.png")


