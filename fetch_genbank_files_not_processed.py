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
csv_file = "Data/gx_details_refseq.20230416_long_contaminations_from_prokaryotes.csv"
output_folder = "Data/Refseq_genbank_over_5000/"
input_data = pd.read_csv(csv_file)

# Get a list of already processed files
existing_files = os.listdir(output_folder)
processed_accessions = set('_'.join(file.split('_')[:-2]) for file in existing_files)

# Process each row in the CSV
for index, row in input_data.iterrows():
    accession = row["#seq_id"]
    if accession in processed_accessions:
        print(f"Already processed {accession}, skipping.")
        continue

    start = int(row["start_pos"]) - 1  # Adjusting to 0-based index
    end = int(row["end_pos"])
    output_file = f"{output_folder}{accession}_{start}_{end}.gb"
    print(accession)

    try:
        genbank_record = fetch_partial_genome(accession, start, end)
        SeqIO.write(genbank_record, output_file, "genbank")
    except Exception as e:
        print(f"Error processing {accession}: {e}")

