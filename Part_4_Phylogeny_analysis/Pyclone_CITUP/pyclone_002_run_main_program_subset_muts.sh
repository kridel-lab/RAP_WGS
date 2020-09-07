#!/bin/bash
#SBATCH -p long
#SBATCH --mem=30720M
#SBATCH -J pyclone_subset
#SBATCH -c 10
#SBATCH -t 21-00:00 # Runtime in D-HH:MM

cd /cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/Pyclone

module load python2
module load pyclone/0.13.1

#patient
PyClone run_analysis_pipeline --tumour_contents 0.789 0.7359 0.6651 0.96266 \
0.91772 0.7864 \
0.96735 0.824 0.94581 0.9119 0.8574 0.804 0.96301 \
0.93695 0.8643 0.95118 0.95338 0.8487 0.95055 0.95822 \
--in_files FFPE_left_breast_15\:S12966E_subset_muts_pyclone_input.tsv \
FFPE_left_axilla_LN_15:S12966C_subset_muts_pyclone_input.tsv \
FFPE_right_neck_LN_15\:S12966A_subset_muts_pyclone_input.tsv \
FT_Adrenal\ gland\,\ NOS_341502_subset_muts_pyclone_input.tsv \
FT_Stomach\,\ NOS_341480_subset_muts_pyclone_input.tsv \
FT_Bladder\,\ NOS_341468_subset_muts_pyclone_input.tsv \
FT_Spleen_341462_subset_muts_pyclone_input.tsv \
FT_Pancreas\,\ NOS_341373_subset_muts_pyclone_input.tsv \
FT_Parotid\ gland_341369_subset_muts_pyclone_input.tsv \
FT_Kidney\,\ NOS_341364_subset_muts_pyclone_input.tsv \
FT_Kidney\,\ NOS_341360_subset_muts_pyclone_input.tsv \
FT_Mediastinum\,\ NOS_341355_subset_muts_pyclone_input.tsv \
FT_Cecum_341166_subset_muts_pyclone_input.tsv \
FT_Retroperitoneum\,\ NOS_340898_subset_muts_pyclone_input.tsv \
FT_Inguinal\ region\,\ NOS_340893_subset_muts_pyclone_input.tsv \
FT_Abdomen\,\ NOS_340889_subset_muts_pyclone_input.tsv \
FT_Omentum_340885_subset_muts_pyclone_input.tsv \
FT_Shoulder\,\ NOS_340881_subset_muts_pyclone_input.tsv \
FT_Cervical\ lymph\ node_340867_subset_muts_pyclone_input.tsv \
FT_Axilla\,\ NOS_340852_subset_muts_pyclone_input.tsv \
--init_method connected --seed 123 --min_cluster_size 5 --working_dir Pyclone_subset_muts_Sept2020

PyClone build_table --config_file Pyclone_subset_muts_Sept2020/config.yaml --out_file Pyclone_subset_muts_Sept2020_table_file.txt --table_type loci
PyClone build_table --config_file Pyclone_subset_muts_Sept2020/config.yaml --out_file Pyclone_subset_muts_Sept2020_table_file_cluster.txt --table_type cluster
