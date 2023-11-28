import os
import argparse
from Bio import SeqIO
import csv

# Set up command line arguments
parser = argparse.ArgumentParser(description="Convert GenBank files to FASTA format and create a samples CSV file.")
parser.add_argument("gbk_dir", help="Directory containing GenBank files")
parser.add_argument("fasta_dir", help="Directory to output FASTA files")
args = parser.parse_args()

# Directories from command line arguments
gbk_directory = args.gbk_dir
fasta_directory = args.fasta_dir

# Ensure the FASTA output directory exists
if not os.path.exists(fasta_directory):
    os.makedirs(fasta_directory)

# CSV file setup
csv_file = 'samples.csv'
csv_headers = ['genome_id', 'source', 'organism', 'genus', 'species', 'strain', 'closest_placement_reference']

# Initialize list to store CSV data
csv_data = []

# Iterate over the gbk files
for filename in os.listdir(gbk_directory):
    if filename.endswith('.gbk'):
        filepath = os.path.join(gbk_directory, filename)
        genome_id = os.path.splitext(filename)[0]  # Extract genome ID from file name

        # Parse the GenBank file
        for record in SeqIO.parse(filepath, "genbank"):
            # Extract data from the record as per your requirements
            source = 'custom' 
            organism = record.annotations.get('organism', '')
            organism_parts = organism.split(' ')
            genus = organism_parts[0] if organism else ''
            species = ' '.join(organism_parts[1:-1]) if len(organism_parts) > 1 else ''
            strain = organism_parts[-1] if len(organism_parts) > 1 else ''
            closest_placement_reference = ''  # Modify as needed

            # Add to CSV data
            csv_data.append([genome_id, source, organism, genus, species, strain, closest_placement_reference])

            # Write to FASTA file
            fasta_filepath = os.path.join(fasta_directory, f"{genome_id}.fna")
            with open(fasta_filepath, "w") as fasta_file:
                SeqIO.write(record, fasta_file, "fasta")

# Write the data to CSV
with open(csv_file, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(csv_headers)
    writer.writerows(csv_data)

print("Conversion complete.")

