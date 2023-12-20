import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import json

minsize = 50
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

df = pd.read_csv(data)
df = df[df["Product_class"].notna()]
# Filter data for specific Host_Phylum and Contamination_Order
# Filter the DataFrame
filtered_df = df[
    (df["Host_Phylum"] == "Pseudomonadota")
    & (df["Contamination_Order"] == "Desulfobacterales")
]


# Apply the functions to create new columns
filtered_df["Category"] = filtered_df.apply(aggregate_by_group, axis=1)
filtered_df["RefSeq_details_dict"] = filtered_df.apply(parse_refseq_details, axis=1)
filtered_df["MIBiG_details_dict"] = filtered_df.apply(parse_mibig_details, axis=1)
filtered_df.to_csv(
    "Data/gx_details_genbank.20230416_contaminations_from_prokaryotes_with_antismash_80_coverage_only_hits_Pseudomonadota_Desulfobacterales.csv"
)
filtered_df = filtered_df.sort_values("Category")
# Compute the counts for each product class
# Compute the counts for each product class
counts = filtered_df.groupby(["Category", "Product_class"]).size().rename("counts")

mibig_counts = (
    filtered_df[
        filtered_df["MIBiG_ID"].notna()
        & (
            filtered_df["MIBiG_details_dict"].apply(
                lambda x: x.get("Percentage_genes_significant_hits", 0)
                > mibig_threshold
            )
        )
    ]
    .groupby(["Category", "Product_class"])
    .size()
    .rename("mibig_counts")
)

refseq_counts = (
    filtered_df[
        filtered_df["RefSeq"].notna()
        & (
            filtered_df["RefSeq_details_dict"].apply(
                lambda x: x.get("Percentage_genes_significant_hits", 0)
                > refseq_threshold
            )
        )
    ]
    .groupby(["Category", "Product_class"])
    .size()
    .rename("refseq_counts")
)

# Merge the counts into one DataFrame and fill missing values with 0
grouped = (
    counts.to_frame()
    .join(mibig_counts, how="outer")
    .join(refseq_counts, how="outer")
    .fillna(0)
)


grouped = grouped.sort_values(by=["Category", "Product_class"])
# Create a figure and axis
fig, ax = plt.subplots(figsize=(10, 6))

# Iterate over each product class
for i, (product_class, row) in enumerate(grouped.iterrows()):
    # Get the category color
    print(row)
    category = row.name[0]
    print(category)

    color = category_colors.get(category, "gray")

    # Plot each bar with different transparency
    ax.bar(i, row["counts"], color=color, label="Total", alpha=0.4, edgecolor="black")
    ax.bar(
        i, row["mibig_counts"], color=color, label="MiBiG", alpha=1, edgecolor="black"
    )
    ax.bar(
        i,
        row["refseq_counts"],
        color=color,
        label="RefSeq",
        alpha=0.7,
        edgecolor="black",
    )

# Set the x-ticks and labels
ax.set_xticks(range(len(grouped)))
ax.set_xticklabels([index[1] for index in grouped.index], rotation=90)

# Set labels and title
plt.xlabel("Product Class")
plt.ylabel("Counts")
plt.title("Stacked Bar Graph for Pseudomonadota and Desulfobacterales")

# Create a legend
handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles))
plt.legend(by_label.values(), by_label.keys())

plt.savefig("Desulfobacterales_Pseudomonadota_product_classes.svg")
