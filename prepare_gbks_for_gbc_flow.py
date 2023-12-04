import os
import csv
from Bio import SeqIO
import argparse

# Set up command line arguments
parser = argparse.ArgumentParser(description="Update GenBank files with detailed lineage information and create a CSV file.")
parser.add_argument("gbk_dir", help="Directory containing GenBank files")
parser.add_argument("csv_file", help="CSV file containing organism information")
parser.add_argument("lineage_file", help="TSV file containing detailed lineage information")
parser.add_argument("output_dir", help="Directory to output updated GenBank files")
args = parser.parse_args()
get_string_before_second_to_last_underscore = lambda s: '_'.join(s.split('_')[:-2]) if s.count('_') >= 2 else None

# Directories and files from command line arguments
gbk_directory = args.gbk_dir
csv_file = args.csv_file
lineage_file = args.lineage_file
output_directory = args.output_dir

# Ensure the output directory exists
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

# CSV file setup for output
samples_csv_file = 'samples.csv'
csv_headers = ['genome_id', 'source', 'organism', 'genus', 'species', 'strain', 'closest_placement_reference']
csv_data = []

# Load organism information from the CSV file into a dictionary
organism_info = {}
contam_to_seq_id = {}

# Load organism information from the CSV file into a dictionary
with open(csv_file, "r") as csv_file:
    reader = csv.DictReader(csv_file)
    for row in reader:
        seq_id = row["#seq_id"]
        contam_details = row["contam_details"]
        organism_info[seq_id] = {
            'contam_details': contam_details,
            'genus': None,  # Initialize genus as None
            'species': None,  # Initialize species as None
            'strain': None,  # Initialize strain as None
        }
        contam_to_seq_id[contam_details] = seq_id  # Map contam_details to seq_id

# Load detailed lineage information from the lineage TSV file
with open(lineage_file, "r") as lineage_file:
    reader = csv.DictReader(lineage_file, delimiter='\t')
    for row in reader:
        original_name = row["Original_Name"]
        # Directly access the seq_id using the contam_details
        seq_id = contam_to_seq_id.get(original_name)
        if seq_id:
            organism_info[seq_id]['genus'] = row["Contamination_Genus"]
            organism_info[seq_id]['species'] = row["Contamination_Species"]
            organism_info[seq_id]['strain'] = row["Contamination_Strain"]

# Iterate over the gbk files and update them
for filename in os.listdir(gbk_directory):
    if filename.endswith('.gbk'):
        filepath = os.path.join(gbk_directory, filename)
        genome_id = get_string_before_second_to_last_underscore(filename)  # Extract genome ID from file name

        # Extract and update the organism info
        organism_info_entry = organism_info.get(genome_id, {})
        organism_name = organism_info_entry.get('contam_details', '') or 'Unknown Organism'

        # Using organism_name.split(" ")[0] as the default for genus if it's not available or empty
        genus = organism_info_entry.get('genus', '') or organism_name.split(" ")[0]

        # Using organism_name as the default for species if it's not available or empty
        species = organism_info_entry.get('species', '') or organism_name

        # Using a static default for strain if it's not available or empty
        strain = organism_info_entry.get('strain', '') or ''

        # Parse the GenBank file and update fields
        records = list(SeqIO.parse(filepath, "genbank"))
        for record in records:
            record.annotations['organism'] = organism_name
            record.annotations['source'] = organism_name
            record.annotations['genus'] = genus
            record.annotations['species'] = species
            record.annotations['strain'] = strain

        # Write the updated records to the new output directory
        output_filepath = os.path.join(output_directory, filename)
        with open(output_filepath, "w") as gbk_file:
            SeqIO.write(records, gbk_file, "genbank")

        # Extract information for CSV file
        csv_data.append([
            genome_id,
            'custom',
            organism_name,
            genus,
            species,
            strain,
            ''  # closest_placement_reference left empty
        ])

# Write the data to CSV
with open(samples_csv_file, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(csv_headers)
    writer.writerows(csv_data)

print("GenBank files updated and samples CSV file created.")

