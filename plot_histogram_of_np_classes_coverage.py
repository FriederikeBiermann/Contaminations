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

# Load the data
df = pd.read_csv(data)
df = df[df["Product_class"].notna()]

# Apply the functions to create new columns
df["Category"] = df.apply(aggregate_by_group, axis=1)
df["RefSeq_details_dict"] = df.apply(parse_refseq_details, axis=1)
df["MIBiG_details_dict"] = df.apply(parse_mibig_details, axis=1)

# Group the data by 'Category' and 'Product_class', and count the occurrences
grouped = df.groupby(["Category", "Product_class"]).size().rename("counts")


# Determine coverage bins
coverage_bins = np.linspace(
    df["coverage"].min(), df["coverage"].max(), 30
)  # Adjusted the number of bins

# Create a figure and axis
fig, ax = plt.subplots(figsize=(10, 6))

# Prepare data for each category
hist_data = []
colors = []
labels = []
for category in df["Category"].unique():
    category_data = df[df["Category"] == category]["coverage"]
    hist_data.append(category_data)
    colors.append(category_colors.get(category, "gray"))
    labels.append(category)

# Plot a stacked histogram for all categories
ax.hist(
    hist_data, bins=coverage_bins, label=labels, color=colors, alpha=0.8, stacked=True, edgecolor="black",
)

# Set labels and title
plt.xlabel("Coverage")
plt.ylabel("Counts")
plt.title("Stacked Histogram of Coverage by Category")

# Create a legend
plt.legend()

# Save the plot
plt.savefig("coverage_histogram_by_category.svg")

plt.show()
