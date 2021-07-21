#One time make symlinks for raw data uploads
#manually renamed files using this conversion provided by Michael Hong

#KLCS_0084=LY_RAP_0001
#KLCS_0085=LY_RAP_0002
#KLCS_0086=LY_RAP_0003

cd /cluster/projects/kridelgroup/RAP_ANALYSIS

#Batch 1 data (May 2021)

input_folder=/cluster/projects/kridelgroup/RAP_ANALYSIS/RAW_DATA/ctDNA/210419_A00469_0169_BH53HYDRXY.KLCS.fastqs
output_folder=/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher/April29OICRupload

ln -s ${input_folder}/KLCS_0084_Pl_n_PE_480_TS_210419_A00469_0169_BH53HYDRXY_1_CTGATCGT-GCGCATAT_R1.fastq.gz ${output_folder}/LY_0001_R1.fastq.gz
ln -s ${input_folder}/KLCS_0084_Pl_n_PE_480_TS_210419_A00469_0169_BH53HYDRXY_1_CTGATCGT-GCGCATAT_R2.fastq.gz ${output_folder}/LY_0001_R2.fastq.gz

ln -s ${input_folder}/KLCS_0085_Pl_n_PE_401_TS_210419_A00469_0169_BH53HYDRXY_1_ACTCTCGA-CTGTACCA_R1.fastq.gz ${output_folder}/LY_0002_R1.fastq.gz
ln -s ${input_folder}/KLCS_0085_Pl_n_PE_401_TS_210419_A00469_0169_BH53HYDRXY_1_ACTCTCGA-CTGTACCA_R2.fastq.gz ${output_folder}/LY_0002_R2.fastq.gz

ln -s ${input_folder}/KLCS_0086_Pl_n_PE_439_TS_210419_A00469_0169_BH53HYDRXY_1_TGAGCTAG-GAACGGTT_R1.fastq.gz ${output_folder}/LY_0003_R1.fastq.gz
ln -s ${input_folder}/KLCS_0086_Pl_n_PE_439_TS_210419_A00469_0169_BH53HYDRXY_1_TGAGCTAG-GAACGGTT_R2.fastq.gz ${output_folder}/LY_0003_R2.fastq.gz

#Batch 2 data (June 30 2021)

input_folder=/cluster/projects/kridelgroup/RAP_ANALYSIS/RAW_DATA/ctDNA/210617_A00469_0183_AHCNJCDRXY.KLCS.fastqs
output_folder=/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher/June30OICRupload

ln -s ${input_folder}/KLCS_0084_Pl_n_PE_415_TS_210617_A00469_0183_AHCNJCDRXY_1_CTGATCGT-GCGCATAT_R1.fastq.gz ${output_folder}/LY_0001_R1.fastq.gz
ln -s ${input_folder}/KLCS_0084_Pl_n_PE_415_TS_210617_A00469_0183_AHCNJCDRXY_1_CTGATCGT-GCGCATAT_R2.fastq.gz ${output_folder}/LY_0001_R2.fastq.gz

ln -s ${input_folder}/KLCS_0085_Pl_n_PE_350_TS_210617_A00469_0183_AHCNJCDRXY_1_ACTCTCGA-CTGTACCA_R1.fastq.gz ${output_folder}/LY_0002_R1.fastq.gz
ln -s ${input_folder}/KLCS_0085_Pl_n_PE_350_TS_210617_A00469_0183_AHCNJCDRXY_1_ACTCTCGA-CTGTACCA_R2.fastq.gz ${output_folder}/LY_0002_R2.fastq.gz

ln -s ${input_folder}/KLCS_0086_Pl_n_PE_381_TS_210617_A00469_0183_AHCNJCDRXY_1_TGAGCTAG-GAACGGTT_R1.fastq.gz ${output_folder}/LY_0003_R1.fastq.gz
ln -s ${input_folder}/KLCS_0086_Pl_n_PE_381_TS_210617_A00469_0183_AHCNJCDRXY_1_TGAGCTAG-GAACGGTT_R2.fastq.gz ${output_folder}/LY_0003_R2.fastq.gz
