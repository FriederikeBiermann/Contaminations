import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


def calculate_statistics_per_coverage(data, coverage_column, length_column):
    """
    Calculate mean, median, and mode of sequence length for each coverage value.

    :param data: DataFrame containing the data
    :param coverage_column: Name of the column containing coverage values
    :param length_column: Name of the column containing sequence length values
    :return: DataFrame with mean, median, and mode of sequence lengths per coverage
    """
    grouped_data = data.groupby(coverage_column)[length_column].agg(
        ["mean", "median", pd.Series.mode]
    )
    # Handling multi-modal data
    grouped_data["mode"] = grouped_data["mode"].apply(
        lambda x: x[0] if isinstance(x, np.ndarray) else x
    )
    return grouped_data.reset_index()


def plot_statistics_per_coverage(
    statistics_data,
    coverage_column,
    color_mean="#A96073",
    color_median="#6982B5",
    color_mode="#615E6E",
):
    """
    Plots the mean, median, and mode of sequence length per coverage.

    :param statistics_data: DataFrame containing the statistics of sequence lengths per coverage
    :param coverage_column: Name of the column containing coverage values
    :param color_mean: Color of the mean plot line
    :param color_median: Color of the median plot line
    :param color_mode: Color of the mode plot line
    """
    print(statistics_data["mode"])
    plt.figure(figsize=(12, 8))
    plt.plot(
        statistics_data[coverage_column],
        statistics_data["mean"],
        color=color_mean,
        label="Mean",
    )
    plt.plot(
        statistics_data[coverage_column],
        statistics_data["median"],
        color=color_median,
        label="Median",
    )
    plt.plot(
        statistics_data[coverage_column],
        statistics_data["mode"],
        color=color_mode,
        label="Mode",
        linestyle="--",
    )
    plt.yscale("log")
    plt.title("Sequence Length Statistics per Coverage")
    plt.xlabel("Coverage")
    plt.ylabel("Sequence Length (Log Scale)")
    plt.legend()
    plt.grid(True)
    plt.show()
    plt.savefig(
        "gx_details_genbank.20230416_long_contaminations_from_prokaryotes_coverage_less_100_length_coverage_median_modus.png"
    )


data = pd.read_csv(
    "Data/gx_details_genbank.20230416_long_contaminations_from_prokaryotes_coverage_less_100.csv"
)

# Calculate the average sequence length per coverage
average_length_data = calculate_statistics_per_coverage(
    data, "coverage", "length_contamination"
)

# Plot the result
plot_statistics_per_coverage(average_length_data, "coverage")
