import pandas as pd
import concurrent.futures
from Bio import Entrez

# Set your email here for NCBI Entrez
Entrez.email = "friederike@biermann-erfurt.de"

def get_taxid(species):
    """Get the TaxID for a species from its name."""
    try:
        with Entrez.esearch(db="taxonomy", term=species) as search:
            result = Entrez.read(search)
            if result["IdList"]:
                taxid = result["IdList"][0]
                return taxid
            else:
                print(f"No TaxID found for species: {species}")
                return None
    except Exception as e:
        print(f"Error fetching TaxID for {species}: {e}")
        return None

def get_lineage(taxid):
    """Get the full lineage of a species from its TaxID."""
    try:
        with Entrez.efetch(db="taxonomy", id=taxid, retmode="xml") as handle:
            records = Entrez.read(handle)
            if records:
                return records[0]["LineageEx"]
            else:
                print(f"No lineage found for TaxID: {taxid}")
                return []
    except Exception as e:
        print(f"Error fetching lineage for TaxID {taxid}: {e}")
        return []

def process_species(name):
    taxid = get_taxid(name)
    print(name)
    if taxid:
        lineage_ex = get_lineage(taxid)
        lineage_info = {}
        if lineage_ex:
            for entry in lineage_ex:
                if entry["Rank"] != "no rank":
                    rank = entry["Rank"].capitalize()
                    lineage_info[f'Contamination_{rank}'] = entry["ScientificName"]
                    lineage_info[f'Contamination_{rank}_TaxID'] = entry["TaxId"]
            lineage_info['Original_Name'] = name
            return lineage_info
    return {'Original_Name': name}

# Read input and output datasets
input_csv = "Data/gx_details_genbank.20230416_contaminations_from_prokaryotes.csv"
output_csv = "Data/gx_details_genbank.20230416_contaminations_from_prokaryotes_lineage.tsv"
input_df = pd.read_csv(input_csv, sep=',')
unique_names = set(input_df['contam_details'].unique())

# Check which names have been fully processed
try:
    output_df = pd.read_csv(output_csv, sep='\t')
    # Assuming that fully processed entries have non-empty values in certain columns
    # Modify this condition based on your actual data structure
    fully_processed_names = set(output_df[output_df['Contamination_Superkingdom'].notna()]['Original_Name'])
    output_df = output_df[output_df['Contamination_Superkingdom'].notna()]
except FileNotFoundError:
    fully_processed_names = set()

# Filter out the fully processed names
names_to_process = unique_names - fully_processed_names
print(f"Processing {len(names_to_process)} names.")

# Use ThreadPoolExecutor to parallelize the process
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as executor:
    results = executor.map(process_species, names_to_process)

# Convert the results to a DataFrame and append to the existing DataFrame
new_lineage_data = pd.DataFrame(list(results))
if len(fully_processed_names)>0:
    lineage_data = pd.concat([output_df, new_lineage_data], ignore_index=True)
else:
    lineage_data = new_lineage_data

# Write to new CSV
lineage_data.to_csv(output_csv, sep='\t', index=False)

