from Bio import SeqIO
import os
from concurrent.futures import ProcessPoolExecutor

def remove_invalid_translations(record):
    for feature in record.features:
        if 'translation' in feature.qualifiers:
            location_length = len(feature.location)
            expected_translation_length = location_length // 3

            actual_translation = feature.qualifiers['translation'][0]
            if len(actual_translation) > expected_translation_length:
                del feature.qualifiers['translation']
    return record

def process_file(input_path, output_path):
    records = list(SeqIO.parse(input_path, "genbank"))
    adjusted_records = [remove_invalid_translations(record) for record in records]
    SeqIO.write(adjusted_records, output_path, "genbank")

def process_genbank_files(input_directory, output_directory, num_workers=8):
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        for filename in os.listdir(input_directory):
            if filename.endswith(".gbk") or filename.endswith(".gb"):
                input_path = os.path.join(input_directory, filename)
                output_path = os.path.join(output_directory, filename)
                executor.submit(process_file, input_path, output_path)

input_directory = "Data/Genbank_genbank_over_5000_coverage_smaller_80_with_gene_annotations"  # Replace with your input directory path
output_directory = "Data/Genbank_genbank_over_5000_coverage_smaller_80_with_gene_annotations_without_nonfunctional_records_without_incorrect_translations"  # Replace with your output directory path

process_genbank_files(input_directory, output_directory)

