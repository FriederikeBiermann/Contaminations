#!/bin/bash
#SBATCH --job-name=convert_gbk
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH -o /beegfs/projects/p450/out_files/%A__NCBI_contaminations_genbank.out
#SBATCH --time=24:00:00
#SBATCH -p batch
#SBATCH --cpus-per-task=1
#SBATCH -e /beegfs/projects/p450/error_files/%A_NCBI_contaminations_genbank.txt


source ~/.bashrc
conda activate /beegfs/home/fbiermann/miniconda3_supernew/envs/Noemi

python3 prepare_gbks_for_gbc_flow.py Data/Refseq_genbank_over_5000/ Data/Refseq_fasta_over_5000/


