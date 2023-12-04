#!/bin/bash 
#SBATCH --job-name=antismash
#SBATCH --ntasks=1
#SBATCH -c 32
#SBATCH -o /beegfs/projects/p450/out_files/Julia_antismash_%A_%a.out  
#SBATCH -t 1-00:00  
#SBATCH -p batch 
#SBATCH --array=0-10%100  # Execute 100 files per node
#SBATCH -e /beegfs/projects/p450/error_files/Julia_antismash_error_%A_%a.txt



source ~/.bashrc
conda activate /beegfs/home/fbiermann/miniconda3_supernew/envs/antismash7


# Define the input and output directories
input_dirs=(
   "/projects/p450/NCBI_contaminations/Contaminations/Data/Genbank_genbank_over_5000_coverage_smaller_40_with_gene_annotations/genbank/"
)  # Update with your input directories
output_parent_dir="/projects/p450/NCBI_contaminations/Contaminations/Data/Genbank_genbank_over_5000_coverage_smaller_40_with_gene_annotations_antismash_out/"  # Replace with your output parent directory

# Get the index of the input directory for this job based on the SLURM_ARRAY_TASK_ID
# Calculate start and end index for files to process
start_index=$(($SLURM_ARRAY_TASK_ID * 500))
end_index=$(($start_index + 499))

# Process 100 files per job
for input_dir in "${input_dirs[@]}"; do
    # Create the output directory for this input directory
    output_dir="${output_parent_dir}/$(basename "${input_dir}")"
    mkdir -p "$output_dir"
    
    # Get the list of files in the directory
    input_files=($(find "${input_dir}" -maxdepth 2 -type f -name "*.gb"))
    
    # Process each file in the current job's range
    for index in $(seq $start_index $end_index); do
        if [ $index -lt ${#input_files[@]} ]; then
            input_file="${input_files[$index]}"
            filename=$(basename "$input_file" .gb)
            
            # Execute the antismash script for the input file
            antismash --fullhmmer  --tigrfam --cb-knownclusters --cb-general --rre  "$input_file" --output-dir "$output_dir/$filename" 
        fi
    done
done

