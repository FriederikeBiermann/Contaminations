import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.ticker as ticker

# Assuming the data is loaded correctly from the CSV file
data = pd.read_csv(
    "Data/gx_details_genbank.20230416_long_contaminations_from_prokaryotes_coverage_less_100.csv"
)


# Function to calculate the number of entries above each threshold
def calculate_threshold_data(data, column_name):
    thresholds = np.linspace(data[column_name].min(), data[column_name].max(), 100)
    counts = [data[data[column_name] >= t].shape[0] for t in thresholds]
    return pd.DataFrame({"threshold": thresholds, "counts": counts})


data["length_contamination_kb"] = data["length_contamination"] / 1000

# Calculate threshold data for the converted column
threshold_data = calculate_threshold_data(data, "length_contamination_kb")

fig, axs = plt.subplots(1, 2, figsize=(15, 6))

# Histogram in kilobases
axs[0].hist(data["length_contamination_kb"], bins=30, color="#6982B5")
axs[0].set_title("Histogram of Coverage")
axs[0].set_xlabel("Length (kb)")
axs[0].set_yscale("log")
axs[0].set_ylabel("Frequency")

# Add more major and minor ticks to the first subplot (histogram)
axs[0].xaxis.set_major_locator(ticker.AutoLocator())
axs[0].xaxis.set_minor_locator(ticker.AutoMinorLocator())

# Threshold curve in kilobases
axs[1].plot(threshold_data["threshold"], threshold_data["counts"], color="#A96073")
axs[1].set_title("Entries above Threshold")
axs[1].set_xlabel("Length Threshold (kb)")
axs[1].set_yscale("log")
axs[1].set_xscale("log")
axs[1].set_ylabel("Number of Entries")

# Add more major and minor ticks to the second subplot (threshold curve)
axs[1].xaxis.set_major_locator(
    ticker.LogLocator(base=10.0, numticks=15)
)  # Base 10 logarithm, adjust numticks as needed
axs[1].xaxis.set_minor_locator(
    ticker.LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1, numticks=15)
)  # Minor ticks


plt.tight_layout()
plt.savefig(
    "gx_details_genbank.20230416_long_contaminations_from_prokaryotes_coverage_less_100_histogram_length_kb.png"
)
