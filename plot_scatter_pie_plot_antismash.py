import matplotlib.pyplot as plt
import json
import pandas as pd
import numpy as np
import matplotlib.colors as mcolors
import seaborn as sns
import matplotlib.patches as mpatches

from matplotlib.lines import Line2D

minsize = 30
refseq_threshold = 0.3
mibig_threshold = 0.3
data = "Data/gx_details_genbank.20230416_contaminations_from_prokaryotes_with_antismash_80_coverage_only_hits.csv"


def parse_refseq_details(row):
    try:
        details = json.loads(row["RefSeq_details"].replace("'", '"'))
        return details
    except json.JSONDecodeError:
        return {}  # Return an empty dict if there's an error


def parse_mibig_details(row):
    try:
        details = json.loads(row["MIBiG_details"].replace("'", '"'))
        return details
    except json.JSONDecodeError:
        return {}  # Return an empty dict if there's an error


def adjust_color_transparency(color, alpha=0.5):
    """Adjust the alpha value of a given color."""
    return mcolors.to_rgba(color, alpha)


def aggregate_by_group(row):
    product_class = row["Product_class"].lower()  # Convert to lowercase
    for group, classes in product_classes.items():
        # Convert each class in the dictionary to lowercase before comparison
        if product_class in [c.lower() for c in classes]:
            return group
    return "Other"


product_classes = {
    "RiPP": [
        "RiPP-like",
        "thiopeptide",
        "redox-cofactor",
        "RRE-containing",
        "triceptide",
        "thioamitides",
        "lassopeptide",
        "lanthipeptide-class-v",
        "lanthipeptide-class-i",
        "RaS-RiPP",
        "ranthipeptide",
        "fungal-cd RiPP-like",
        "lanthipeptide-class-ii",
        "lanthipeptide-class-iii",
        "microviridin",
        "lanthipeptide-class-iv",
        "sactipeptide",
        "proteusin",
        "darobactin",
        "crocagin",
        "methanobactin",
        "bottromycin",
        "fungal-RiPP",
        "epipeptide",
        "guanidinotides",
        "achaeal-RiPP",
        "spliceotide",
        "lipolanthine",
    ],
    "NRP": [
        "NRP-metallophore",
        "NRPS-like",
        "NRPS",
        "thioamide-NRP",
        "isocyanide-nrp",
        "aminocoumarin",
    ],
    "PK": [
        "hserlactone",
        "T3PKS",
        "arylpolyene",
        "T1PKS",
        "transAT-PKS",
        "PKS-like",
        "hglE-KS",
        "HR-T2PKS",
        "T2PKS",
        "transAT-PKS-like",
    ],
    "Terpene": ["terpene"],
    "Alkaloid": ["ectoine", "prodigiosin", "indole", "pyrrolidine"],
    "Other": [
        "phosphonate",
        "betalactone",
        "NI-siderophore",
        "NAPAA",
        "NAGGN",
        "cyclic-lactone-autoinducer",
        "hydroxycyanide",
        "opine-like-metallophore",
        "CDPS",
        "amglyccycl",
        "furan",
        "resorcinol",
        "other",
        "butyrolactone",
        "phenazine",
        "LAP",
        "acyl_amino_acids",
        "blactam",
        "linaridin",
        "oligosaccharide",
        "melanin",
        "glycocin",
        "phosphonate-like",
        "nucleoside",
        "2dos",
        "PUFA",
        "tropodithietic-acid",
        "aminopolycarboxylic-acid",
        "mycosporine-like",
        "phosphoglycolipid",
        "PBDE",
        "ladderane",
        "RCDPS",
    ],
}

category_colors = {
    "NRP": "#516DA4",
    "PK": "#8B765A",
    "Alkaloid": "#658F65",
    "Terpene": "#924F60",
    "RiPP": "#73668F",
    "Other": "#44424D",
}

# Reverse mapping from specific product class to broader category
category_map = {
    specific: broad
    for broad, specifics in product_classes.items()
    for specific in specifics
}

df = pd.read_csv(data)
df = df.dropna(subset=["Host_Phylum", "Contamination_Order"])

# filter only ones with BGCs
df = df[df["Product_class"].notna()]
print(df["Product_class"].unique())
df["Host_Kingdom_Group"] = df["Host_Kingdom"].fillna(df["Host_Superkingdom"])
df = df.sort_values(by=["Host_Kingdom_Group"])

# Aggregate data
df["Group"] = df.apply(aggregate_by_group, axis=1)
bgc_counts = df.groupby(["Contamination_Order", "Host_Phylum"]).size()


# Calculate the distribution of 'Group' within each host/contamination combination
group_distribution = df.groupby(["Contamination_Order", "Host_Phylum", "Group"]).size()
group_distribution_df = group_distribution.reset_index(name="count")
bgc_counts_df = bgc_counts.reset_index(name="total_count")

# Merge to align the indices
merged_df = pd.merge(
    group_distribution_df,
    bgc_counts_df,
    on=["Contamination_Order", "Host_Phylum"],
    how="left",
)

# Calculate percentages
merged_df["percentage"] = merged_df["count"] / merged_df["total_count"]
group_counts = merged_df.set_index(["Contamination_Order", "Host_Phylum", "Group"])[
    "count"
]

df["RefSeq_details_dict"] = df.apply(parse_refseq_details, axis=1)

# Now filter based on the 'Percentage_genes_significant_hits' value
RefSeq_df = df[
    df["RefSeq_details_dict"].apply(
        lambda x: x.get("Percentage_genes_significant_hits", 0) > refseq_threshold
    )
]

# Calculate refseq percentages for each group
refseq_group_distribution = (
    RefSeq_df[RefSeq_df["RefSeq"].notna()]
    .groupby(["Contamination_Order", "Host_Phylum", "Group"])
    .size()
)

df["MIBiG_details_dict"] = df.apply(parse_mibig_details, axis=1)

# Now filter based on the 'Percentage_genes_significant_hits' value
MIBiG_df = df[
    df["MIBiG_details_dict"].apply(
        lambda x: x.get("Percentage_genes_significant_hits", 0) > mibig_threshold
    )
]


# Calculate MIBiG percentages for each group
mibig_group_distribution = (
    MIBiG_df[MIBiG_df["MIBiG_ID"].notna()]
    .groupby(["Contamination_Order", "Host_Phylum", "Group"])
    .size()
)


# Unique categories for x and y axes
unique_contaminations = df["Contamination_Order"].unique()
unique_hosts = df["Host_Phylum"].unique()
unique_product_classes = df["Product_class"].dropna().unique()

print(unique_product_classes)
# Mapping of categories to numeric values
contamination_to_num = {
    contamination: i for i, contamination in enumerate(unique_contaminations)
}
host_to_num = {host: i for i, host in enumerate(unique_hosts)}

# Create scatter plot
plt.figure(figsize=(30, 30))

for (contamination, host), count in bgc_counts.items():
    y = contamination_to_num[contamination]
    x = host_to_num[host]
    pie_size = minsize + count * 10

    if (contamination, host) in group_counts.index:
        group_sizes = group_counts.loc[(contamination, host)].values
        group_labels = group_counts.loc[(contamination, host)].index
        total_size = sum(group_sizes)
        group_colors = [category_colors[group] for group in group_labels]
        # Main pie chart segment
        plt.pie(
            x=group_sizes,
            colors=[
                adjust_color_transparency(color, alpha=0.40) for color in group_colors
            ],
            center=(x, y),
            radius=pie_size * 0.01,
            counterclock=False,
            wedgeprops={"width": pie_size * 0.01, "edgecolor": "w"},
        )

        # nested chart
        inx_group = 0
        for group_size, group in zip(group_sizes, group_labels):
            # Calculate the angle for each segment
            segment_angle = (group_size / total_size) * 360
            # Nested pie chart for MIBiG entries
            mibig_size = mibig_group_distribution.get((contamination, host, group), 0)
            if mibig_size > 0:
                inner_radius = (mibig_size / group_size) * pie_size

                patches, texts = plt.pie(
                    x=group_sizes,
                    colors=group_colors,
                    center=(x, y),
                    radius=inner_radius * 0.01,
                    counterclock=False,
                    wedgeprops={"width": inner_radius * 0.01, "edgecolor": "w"},
                )
                # Iterate over slices
                for idx, wedge in enumerate(patches):
                    # Hide all slices except the Age Group 35-54
                    if idx != inx_group:
                        wedge.set_visible(False)
                        texts[idx].set_visible(False)
                inx_group += 1
            refseq_size = refseq_group_distribution.get((contamination, host, group), 0)
            # print(
            #     group, mibig_size, refseq_size, group_size, (refseq_size / group_size)
            # )
            if refseq_size > 0:
                inner_radius = (refseq_size / group_size) * pie_size

                patches, texts = plt.pie(
                    x=group_sizes,
                    colors=[
                        adjust_color_transparency(color, alpha=0.50)
                        for color in group_colors
                    ],
                    center=(x, y),
                    radius=inner_radius * 0.01,
                    counterclock=False,
                    wedgeprops={"width": inner_radius * 0.01, "edgecolor": "w"},
                )
                # Iterate over slices
                for idx, wedge in enumerate(patches):
                    # Hide all slices except the Age Group 35-54
                    if idx != inx_group:
                        wedge.set_visible(False)
                        texts[idx].set_visible(False)
                inx_group += 1


plt.yticks(range(len(unique_contaminations)), unique_contaminations, fontsize=20)
plt.xticks(range(len(unique_hosts)), unique_hosts, rotation=90, fontsize=20)
plt.ylabel("Contamination Phylum")
plt.xlabel("Host Phylum")
# Adjust axis label positions
ax = plt.gca()  # Get the current axes instance

# Move the x-axis label to the left and down
# The first number moves it left (negative) or right (positive)
# The second number moves it down (negative) or up (positive)
ax.xaxis.set_label_coords(5, 5)

# Move the y-axis label down and to the left
ax.yaxis.set_label_coords(-5, 5)

# Create a list of patches to use as handles in the legend
legend_handles = [
    mpatches.Patch(color=color, label=category)
    for category, color in category_colors.items()
]
# Create custom legend handles for pie chart sizes
size_legend_handles = [
    Line2D(
        [],
        [],
        marker="o",
        color="w",
        label=f"{factor} counts",
        markerfacecolor="#44424D",
        markersize=(minsize + factor * 10) / np.pi,
        linestyle="None"
    )
    for factor in [1, 5, 10]  # The factors are square-rooted to scale the marker size
]

# Plot the size legend below the existing legend
size_legend = ax.legend(
    handles=size_legend_handles,
    title="Pie Chart Sizes",
    loc="lower left",
    bbox_to_anchor=(1.05, 0),
    borderaxespad=0.0,
)
ax.add_artist(size_legend)

# Add the legend to the plot
plt.legend(
    handles=legend_handles,
    title="Categories",
    bbox_to_anchor=(1.05, 1),
    loc="upper left",
)
plt.savefig("scatter_plot_antismash_coverage_under_80_genbank_contamination_order.svg")
