import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib.colors as mcolors
import seaborn as sns
import matplotlib.patches as mpatches
import json

minsize = 50
refseq_threshold = 0.5
mibig_threshold = 0.5
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
    "PK": "#44424D",
    "Alkaloid": "#658F65",
    "Terpene": "#924F60",
    "RiPP": "#73668F",
    "Other": "#727272",
}

# Reverse mapping from specific product class to broader category
category_map = {
    specific: broad
    for broad, specifics in product_classes.items()
    for specific in specifics
}

df = pd.read_csv(data)
df = df.dropna(subset=["Host_Phylum", "Contamination_Family"])

# filter only ones with BGCs
df = df[df["Product_class"].notna()]
print(df["Product_class"].unique())

# Aggregate data
df["Group"] = df.apply(aggregate_by_group, axis=1)
bgc_counts = df.groupby(["Contamination_Family", "Host_Phylum"]).size()

# Calculate the distribution of 'Group' within each host/contamination combination
group_distribution = df.groupby(["Contamination_Family", "Host_Phylum", "Group"]).size()
group_distribution_df = group_distribution.reset_index(name="count")
bgc_counts_df = bgc_counts.reset_index(name="total_count")

# Merge to align the indices
merged_df = pd.merge(
    group_distribution_df,
    bgc_counts_df,
    on=["Contamination_Family", "Host_Phylum"],
    how="left",
)

# Calculate percentages
merged_df["percentage"] = merged_df["count"] / merged_df["total_count"]
group_counts = merged_df.set_index(["Contamination_Family", "Host_Phylum", "Group"])[
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
    .groupby(["Contamination_Family", "Host_Phylum", "Group"])
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
    .groupby(["Contamination_Family", "Host_Phylum", "Group"])
    .size()
)



# Unique categories for x and y axes
unique_contaminations = df["Contamination_Family"].unique()
unique_hosts = df["Host_Phylum"].unique()
unique_product_classes = df["Product_class"].dropna().unique()


# Mapping of categories to numeric values
contamination_to_num = {
    contamination: i for i, contamination in enumerate(unique_contaminations)
}
host_to_num = {host: i for i, host in enumerate(unique_hosts)}

# Create scatter plot
plt.figure(figsize=(50, 50))

for (contamination, host), count in bgc_counts.items():
    x = contamination_to_num[contamination]
    y = host_to_num[host]
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
                adjust_color_transparency(color, alpha=0.50) for color in group_colors
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
            # Nested pie chart for MIBiG entries
            mibig_size = mibig_group_distribution.get((contamination, host, group), 0)
            if mibig_size != group_size:
                inner_radius = (1 - mibig_size / group_size) * pie_size

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
            refseq_size = refseq_group_distribution.get((contamination, host, group), 0)
            print(
                group, mibig_size, refseq_size, group_size, (refseq_size / group_size)
            )
            if refseq_size != group_size:
                inner_radius = (1 - refseq_size / group_size) * pie_size

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


plt.xticks(range(len(unique_contaminations)), unique_contaminations, rotation=90)
plt.yticks(range(len(unique_hosts)), unique_hosts)
plt.xlabel("Contamination Phylum")
plt.ylabel("Host Phylum")
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

# Add the legend to the plot
plt.legend(
    handles=legend_handles,
    title="Categories",
    bbox_to_anchor=(1.05, 1),
    loc="upper left",
)
plt.savefig(
    "scatter_plot_antismash_coverage_under_80_genbank_contamination_family_reverse.svg"
)
plt.show()
