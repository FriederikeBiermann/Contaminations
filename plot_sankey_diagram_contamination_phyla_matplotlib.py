import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Assuming df is your DataFrame
filename = "Data/gx_details_refseq.20230416_contaminations_from_prokaryotes.csv"
df = pd.read_csv(filename)


# Assuming df is your DataFrame
# Use Host_Kingdom if available, otherwise use Host_Superkingdom
df["Host_Kingdom_Group"] = df["Host_Kingdom"].fillna(df["Host_Superkingdom"])

# Prepare the data
grouped_data = (
    df.groupby(["Host_Kingdom_Group", "Host_Phylum", "Contamination_Phylum"])
    .size()
    .reset_index(name="counts")
)

# Create mappings for the axes
host_kingdom_groups = grouped_data["Host_Kingdom_Group"].unique()
host_phyla = grouped_data["Host_Phylum"].unique()
contam_phyla = grouped_data["Contamination_Phylum"].unique()

# Create figure and axis
fig, ax = plt.subplots(figsize=(12, 8))

# Custom colors
colors = ["red", "green", "blue", "orange", "purple", "cyan", "magenta"]

# Draw the lines
for i, row in grouped_data.iterrows():
    host_kingdom_index = np.where(host_kingdom_groups == row["Host_Kingdom_Group"])[0][
        0
    ]
    host_phylum_index = np.where(host_phyla == row["Host_Phylum"])[0][0] + len(
        host_kingdom_groups
    )
    contam_phylum_index = (
        np.where(contam_phyla == row["Contamination_Phylum"])[0][0]
        + len(host_kingdom_groups)
        + len(host_phyla)
    )

    # Draw line from host kingdom group to host phylum
    ax.plot(
        [0, 1], [host_kingdom_index, host_phylum_index], color=colors[i % len(colors)]
    )

    # Draw line from host phylum to contamination phylum
    ax.plot(
        [1, 2], [host_phylum_index, contam_phylum_index], color=colors[i % len(colors)]
    )

# Set labels for host kingdoms, host phyla, and contam phyla
ax.set_yticks(
    list(range(len(host_kingdom_groups) + len(host_phyla) + len(contam_phyla)))
)
ax.set_yticklabels(list(host_kingdom_groups) + list(host_phyla) + list(contam_phyla))

# Set x-axis labels
ax.set_xticks([0, 1, 2])
ax.set_xticklabels(["Host Kingdom/Group", "Host Phylum", "Contamination Phylum"])

# Customizing the plot
ax.set_title("Host to Contamination Phylum Flow")
plt.tight_layout()

# Save the figure
plt.savefig("custom_sankey_diagram.png")
plt.show()
