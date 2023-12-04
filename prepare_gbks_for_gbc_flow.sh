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

python3 prepare_gbks_for_gbc_flow.py /projects/p450/NCBI_contaminations/Contaminations/Data/test_for_bgcflow/ /projects/p450/NCBI_contaminations/Contaminations/Data/gx_details_refseq.20230416_long_contaminations_from_prokaryotes.csv /projects/p450/NCBI_contaminations/Contaminations/Data/gx_details_refseq.20230416_contaminations_from_prokaryotes_lineage.tsv /projects/p450/NCBI_contaminations/Contaminations/Data/test_for_bgcflow_correct_species/


