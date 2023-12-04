#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
from Bio import Entrez, SeqIO
import os

# Set your email here for NCBI Entrez
Entrez.email = "friederike@biermann-erfurt.de"

def fetch_partial_genome(accession, start, end):
    with Entrez.efetch(db="nucleotide", id=accession, rettype='gb', retmode="text", seq_start=start, seq_stop=end) as handle:
        return SeqIO.read(handle, "genbank")

# Import CSV file and output folder
csv_file = "Data/gx_details_genbank.20230416_long_contaminations_from_prokaryotes_coverage_less_80.csv"
output_folder = "/projects/p450/NCBI_contaminations/Contaminations/Data/Genbank_genbank_over_5000_coverage_smaller_80/"
input_data = pd.read_csv(csv_file)

# Get a list of already processed files
existing_files = os.listdir(output_folder)
processed_accessions = set(file for file in existing_files)
# Process each row in the CSV
for index, row in input_data.iterrows():
    accession = row["#seq_id"]

    start = int(row["start_pos"]) - 1  # Adjusting to 0-based index
    end = int(row["end_pos"])
    filename_raw = f"{accession}_{start}_{end}.gb"
    output_file = f"{output_folder}{filename_raw}"
    print(filename_raw)
    if filename_raw  in processed_accessions:
        print(f"Already processed {accession}, skipping.")
        continue

    try:
        genbank_record = fetch_partial_genome(accession, start, end)
        SeqIO.write(genbank_record, output_file, "genbank")
    except Exception as e:
        print(f"Error processing {accession}: {e}")

