import pandas as pd
import numpy as np
import concurrent.futures
from Bio import Entrez

# Set your email here for NCBI Entrez
Entrez.email = "f.biermann@bio.uni-frankfurt.de"


def get_lineage(taxid):
    """Get the full lineage of a species from its TaxID."""
    try:
        with Entrez.efetch(db="taxonomy", id=str(taxid), retmode="xml") as handle:
            records = Entrez.read(handle)
            if records:
                return records[0]["LineageEx"]
            else:
                print(f"No lineage found for TaxID: {taxid}")
                return []
    except Exception as e:
        print(f"Error fetching lineage for TaxID {taxid}: {e}")
        return []

def process_species(taxid):
    if taxid:
        lineage_info = {}
        print(taxid)
        lineage_ex = get_lineage(taxid)
        if lineage_ex:
            for entry in lineage_ex:
                if entry["Rank"] != "no rank":
                    rank = entry["Rank"].capitalize()
                    lineage_info[f'Host_{rank}'] = entry["ScientificName"]
                    lineage_info[f'Host_{rank}_TaxID'] = entry["TaxId"]
            lineage_info['Original_Name'] = taxid
            return lineage_info
    return {'Original_Name': taxid}

# Read input and output datasets
input_csv = "Data/gx_details_genbank.20230416_contaminations_from_prokaryotes.csv"
output_csv = "Data/gx_details_genbank.20230416_contaminations_from_prokaryotes_lineage_host.tsv"
input_df = pd.read_csv(input_csv, sep=',', low_memory=False)
input_df = input_df[np.isfinite(input_df['taxid'])]
input_df['taxid'] = input_df['taxid'].fillna(input_df['species_taxid'])
input_df['taxid'] = input_df['taxid'].dropna()
input_df['taxid'] = input_df['taxid'].astype(int)
unique_taxids = set(input_df['taxid'].unique())

# Read the existing output file, if it exists, to get already processed taxids
try:
    output_df = pd.read_csv(output_csv, sep='\t')
    processed_taxids = set(output_df[output_df['Host_Superkingdom'].notna()]['Original_Name'])
    output_df = output_df[output_df['Host_Superkingdom'].notna()]
except FileNotFoundError:
    processed_taxids = set()

# Filter out the already processed taxids
taxids_to_process = unique_taxids - processed_taxids
print(f"Processing {len(taxids_to_process)} taxids.")

# Use ThreadPoolExecutor to parallelize the process
with concurrent.futures.ThreadPoolExecutor(max_workers=1) as executor:
    results = executor.map(process_species, taxids_to_process)

# Convert the results to a DataFrame and append to the existing DataFrame
new_lineage_data = pd.DataFrame(list(results))
if len(processed_taxids)>0:
    lineage_data = pd.concat([output_df, new_lineage_data], ignore_index=True)
else:
    lineage_data = new_lineage_data

# Write to new CSV
lineage_data.to_csv(output_csv, sep='\t', index=False)

