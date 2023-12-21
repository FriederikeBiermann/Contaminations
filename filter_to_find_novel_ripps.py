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

RefSeq_df.to_csv(f"{data[:-4]}_without_refseq.csv")
RiPP_RefSeq = RefSeq_df[RefSeq_df["Group"] == "RiPP"]
RiPP_RefSeq.to_csv(f"{data[:-4]}_without_refseq_ripps.csv")
print(RefSeq_df)
print(RiPP_RefSeq)
