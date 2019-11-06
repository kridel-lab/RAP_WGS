#!/bin/bash
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -J pyclone_RAP
#SBATCH -c 8
#SBATCH -t 5-00:00 # Runtime in D-HH:MM

pwd
#/cluster/projects/kridelgroup/RAP_ANALYSIS/pyclone/rap_pyclone/rap_pyclone

module load pyclone/0.13.1

#patient
PyClone run_analysis_pipeline --tumour_contents 0.96028 0.95311 0.84850 0.95339 0.95379 0.86410 0.93685 0.96303 0.80420 0.85730 0.90871 \
0.94584 0.82400 0.96622 0.78610 0.91774 0.96269 0.66560 0.73560 0.51490 --in_files LY_RAP_0003_Aut_FzT_01_pyclone_input.tsv \
LY_RAP_0003_Aut_FzT_02_pyclone_input.tsv LY_RAP_0003_Aut_FzT_03_pyclone_input.tsv LY_RAP_0003_Aut_FzT_04_pyclone_input.tsv \
LY_RAP_0003_Aut_FzT_05_pyclone_input.tsv LY_RAP_0003_Aut_FzT_06_pyclone_input.tsv LY_RAP_0003_Aut_FzT_07_pyclone_input.tsv LY_RAP_0003_Aut_FzT_09_pyclone_input.tsv \
LY_RAP_0003_Aut_FzT_10_pyclone_input.tsv LY_RAP_0003_Aut_FzT_11_pyclone_input.tsv LY_RAP_0003_Aut_FzT_12_pyclone_input.tsv LY_RAP_0003_Aut_FzT_13_pyclone_input.tsv \
LY_RAP_0003_Aut_FzT_14_pyclone_input.tsv LY_RAP_0003_Aut_FzT_15_pyclone_input.tsv LY_RAP_0003_Aut_FzT_16_pyclone_input.tsv LY_RAP_0003_Aut_FzT_17_pyclone_input.tsv \
LY_RAP_0003_Aut_FzT_18_pyclone_input.tsv LY_RAP_0003_Dia_FoT_01_pyclone_input.tsv LY_RAP_0003_Dia_FoT_03_pyclone_input.tsv LY_RAP_0003_Dia_FoT_05_pyclone_input.tsv \
 --working_dir RAP_WGS_pyclone

PyClone build_table --config_file RAP_WGS_pyclone/config.yaml --out_file RAP_WGS_pyclone_table_file.txt --table_type loci
PyClone build_table --config_file RAP_WGS_pyclone/config.yaml --out_file RAP_WGS_pyclone_table_file_cluster.txt --table_type cluster



