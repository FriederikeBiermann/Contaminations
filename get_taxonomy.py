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
        if lineage_ex:
            for entry in lineage_ex:
                if entry["Rank"] != "no rank":
                    rank = entry["Rank"].capitalize()
                    lineage_info[f'Contamination_{rank}'] = entry["ScientificName"]
                    lineage_info[f'Contamination_{rank}_TaxID'] = entry["TaxId"]
            lineage_info['Original_Name'] = name
            return lineage_info
    return {'Original_Name': name}

# Read dataset
input_csv = "Data/gx_details_refseq.20230416_contaminations_from_prokaryotes.csv"  # Replace with your dataset path
output_csv = "Data/gx_details_refseq.20230416_contaminations_from_prokaryotes_lineage.tsv"  # Replace with your desired output path
df = pd.read_csv(input_csv, sep=',')

# Extract unique names
unique_names = df['contam_details'].unique()
print(len(unique_names))

# Use ThreadPoolExecutor to parallelize the process
with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
    results = executor.map(process_species, unique_names)

# Convert the results to a DataFrame
lineage_data = pd.DataFrame(list(results))

# Write to new CSV
lineage_data.to_csv(output_csv, index=False)
