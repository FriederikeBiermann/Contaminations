import os
import shutil
from Bio import SeqIO

# Define the source and destination directories
source_dir = '/projects/p450/NCBI_contaminations/Contaminations/Data/Genbank_genbank_over_5000_coverage_smaller_80'  # Replace with the path to your source directory
destination_dir = '/projects/p450/NCBI_contaminations/Contaminations/Data/Genbank_genbank_over_5000_coverage_smaller_80_with_gene_annotations'  # Replace with the path to your destination directory

# Check and create destination directory if it doesn't exist
if not os.path.exists(destination_dir):
    os.makedirs(destination_dir)

# Function to check if a file contains at least one gene annotation
def has_gene_annotation(file_path):
    with open(file_path, 'r') as file:
        for record in SeqIO.parse(file, 'genbank'):
            for feature in record.features:
                if feature.type == 'gene':
                    return True
    return False

# Iterate over all files in the source directory
for filename in os.listdir(source_dir):
    if filename.endswith('.gb'):
        file_path = os.path.join(source_dir, filename)
        if has_gene_annotation(file_path):
            shutil.copy(file_path, destination_dir)
            print(f"Copied {filename}")

