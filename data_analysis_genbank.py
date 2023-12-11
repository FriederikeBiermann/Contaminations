import pandas as pd
import os

file_path_detailed_data = "Data/gx_details_genbank.20230416.tsv"
file_path_summary_data = "Data/gx_summary_genbank.20230416.tsv"

# # Read the TSV file
# raw_detailed_data = pd.read_csv(file_path_detailed_data, sep="\t")
# raw_summary_data = pd.read_csv(file_path_summary_data, sep="\t", encoding="latin1")
# print(f"Rows: {raw_detailed_data.shape[0]}, Columns: {raw_detailed_data.shape[1]}")

# raw_detailed_data = pd.merge(
#     raw_detailed_data, raw_summary_data, on="assembly_accession", how="left"
# )
# raw_detailed_data.to_csv(f"{file_path_detailed_data[:-4]}_with_summary_info.csv")
raw_detailed_data = pd.read_csv(f"{file_path_detailed_data[:-4]}_with_summary_info.csv")
# Filter rows where 'action' column is either 'TRIM' or 'EXCLUDE'
data_without_manual_review = raw_detailed_data[
    raw_detailed_data["action"].isin(["TRIM", "EXCLUDE"])
]
print(
    f"Rows: {data_without_manual_review.shape[0]}, Columns: {data_without_manual_review.shape[1]}"
)

# Further filter rows where 'contam_type' starts with 'prok:'
data_prokaryotes = data_without_manual_review[
    data_without_manual_review["contam_type"].str.startswith("prok:")
]
print(f"Rows: {data_prokaryotes.shape[0]}, Columns: {data_prokaryotes.shape[1]}")

data_prokaryotes.to_csv(
    f"{file_path_detailed_data[:-4]}_contaminations_from_prokaryotes.csv"
)

# merge with lineages
file_path_prokaryotes_lineage = "Data/gx_details_genbank.20230416_contaminations_from_prokaryotes_lineage.tsv"
file_path_prokaryotes_lineage_host = "Data/gx_details_genbank.20230416_contaminations_from_prokaryotes_lineage_host.tsv"
# Load lineage data
data_lineage = pd.read_csv(file_path_prokaryotes_lineage, sep="\t")
# Merge with lineage data on specified columns
merged_data = pd.merge(data_prokaryotes, data_lineage, left_on="contam_details", right_on="Original_Name", how="left")

# Load lineage host data
data_lineage_host = pd.read_csv(file_path_prokaryotes_lineage_host, sep="\t")
# Merge with lineage host data on specified columns
data_prokaryotes = pd.merge(merged_data, data_lineage_host, left_on="taxid", right_on="Original_Name", how="left")

# Filtering out very short contaminations

data_prokaryotes["length_contamination"] = (
    data_prokaryotes["end_pos"] - data_prokaryotes["start_pos"]
)
data_prokaryotes_long_contaminations = data_prokaryotes[
    data_prokaryotes["length_contamination"] > 5000
]

print(
    f"Rows: {data_prokaryotes_long_contaminations.shape[0]}, Columns: {data_prokaryotes_long_contaminations.shape[1]}"
)
data_prokaryotes_long_contaminations.to_csv(
    f"{file_path_detailed_data[:-4]}_long_contaminations_from_prokaryotes.csv"
)



data_prokaryotes_long_contaminations_without_known_hit = data_prokaryotes_long_contaminations[
        data_prokaryotes_long_contaminations["coverage"] < 100
    ]

print(
    f"Rows: {data_prokaryotes_long_contaminations_without_known_hit[0]}, Columns: {data_prokaryotes_long_contaminations_without_known_hit[1]}"
)

data_prokaryotes_long_contaminations_without_known_hit.to_csv(
    f"{file_path_detailed_data[:-4]}_long_contaminations_from_prokaryotes_coverage_less_100.csv"
)
