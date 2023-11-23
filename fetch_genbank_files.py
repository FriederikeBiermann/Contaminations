#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import Seq

# Set your email here for NCBI Entrez
Entrez.email = "friederike@biermann-erfurt.de"

def fetch_partial_genome(accession, start, end):
    with Entrez.efetch(db="nucleotide", id=accession, rettype='gb', retmode="text", seq_start=start, seq_stop=end) as handle:
        return SeqIO.read(handle, "genbank")

# Import CSV file and output folder
csv_file = "Data/gx_details_refseq.20230416_long_contaminations_from_prokaryotes.csv"
output_folder = "Data/Refseq_genbank_over_5000/"
input_data = pd.read_csv(csv_file)

# Process each row in the CSV
for index, row in input_data.iterrows():
    accession = row["#seq_id"]
    start = int(row["start_pos"]) - 1  # Adjusting to 0-based index
    end = int(row["end_pos"])
    print(accession)
    try:
        genbank_record = fetch_partial_genome(accession, start, end)
        output_file = f"{output_folder}{accession}_{start}_{end}.gb"
        SeqIO.write(genbank_record, output_file, "genbank")
    except Exception as e:
        print(f"Error processing {accession}: {e}")

 
