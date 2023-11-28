import pandas as pd
import os

file_path_detailed_data = "Data/gx_details_refseq.20230416.tsv"
file_path_summary_data = "Data/gx_summary_refseq.20230416.tsv"

# Read the TSV file
raw_detailed_data = pd.read_csv(file_path_detailed_data, sep='\t')
raw_summary_data = pd.read_csv(file_path_summary_data, sep='\t')
print(print(f"Rows: {raw_detailed_data.shape[0]}, Columns: {raw_detailed_data.shape[1]}"))

raw_detailed_data = pd.merge(raw_detailed_data , raw_summary_data, on='assembly_accession', how='left')
raw_detailed_data.to_csv(f"{file_path_detailed_data[:-4]}_with_summary_info.csv")
# Filter rows where 'action' column is either 'TRIM' or 'EXCLUDE'
data_without_manual_review = raw_detailed_data[raw_detailed_data['action'].isin(['TRIM', 'EXCLUDE'])]
print(print(f"Rows: {data_without_manual_review.shape[0]}, Columns: {data_without_manual_review.shape[1]}"))

# Further filter rows where 'contam_type' starts with 'prok:'
data_prokaryotes = data_without_manual_review[data_without_manual_review['contam_type'].str.startswith('prok:')]
print(print(f"Rows: {data_prokaryotes.shape[0]}, Columns: {data_prokaryotes.shape[1]}"))

data_prokaryotes.to_csv(f"{file_path_detailed_data[:-4]}_contaminations_from_prokaryotes.csv")

# Filtering out very short contaminations

data_prokaryotes['length_contamination'] = data_prokaryotes['end_pos'] - data_prokaryotes['start_pos']
data_prokaryotes_long_contaminations = data_prokaryotes[data_prokaryotes['length_contamination'] > 5000]

data_prokaryotes_long_contaminations.to_csv(f"{file_path_detailed_data[:-4]}_long_contaminations_from_prokaryotes.csv")


