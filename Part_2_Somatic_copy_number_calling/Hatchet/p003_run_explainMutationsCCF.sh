conda activate hatchet

cd /Users/kisaev/Documents/Hatchet_analysis

snvs=LY_RAP_0003.csv
seg_file=/Users/kisaev/Documents/Hatchet_analysis/p001/results/chosen.diploid.seg.ucn

python /Users/kisaev/github/RAP_WGS/Part_2_Somatic_copy_number_calling/Hatchet/explainMutationsCCF.py $snvs -s $seg_file> output.txt
