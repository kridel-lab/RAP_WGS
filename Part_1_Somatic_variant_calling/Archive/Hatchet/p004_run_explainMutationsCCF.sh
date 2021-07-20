conda activate hatchet

cd /Users/kisaev/Documents/Hatchet_analysis

#LY_RAP_0001

snvs=/Users/kisaev/Documents/Hatchet_analysis/p001/LY_RAP_0001.csv
awk '$4="LY_RAP_0001-"$4' /Users/kisaev/Documents/Hatchet_analysis/p001/results/best.seg.ucn > /Users/kisaev/Documents/Hatchet_analysis/p001/cns_file.seg.ucn
seg_file=/Users/kisaev/Documents/Hatchet_analysis/p001/cns_file.seg.ucn
#nano $seg_file and remove patient name from SAMPLE column 
python /Users/kisaev/github/RAP_WGS/Part_2_Somatic_copy_number_calling/Hatchet/explainMutationsCCF.py $snvs -s $seg_file > p001_SNVs_explained_by_hatchet_output.txt

#LY_RAP_0002

snvs=/Users/kisaev/Documents/Hatchet_analysis/p002/LY_RAP_0002.csv
awk '$4="LY_RAP_0002-"$4' /Users/kisaev/Documents/Hatchet_analysis/p002/results/best.seg.ucn > /Users/kisaev/Documents/Hatchet_analysis/p002/cns_file.seg.ucn
seg_file=/Users/kisaev/Documents/Hatchet_analysis/p002/cns_file.seg.ucn
python /Users/kisaev/github/RAP_WGS/Part_2_Somatic_copy_number_calling/Hatchet/explainMutationsCCF.py $snvs -s $seg_file> p002_SNVs_explained_by_hatchet_output.txt

#LY_RAP_0003

snvs=/Users/kisaev/Documents/Hatchet_analysis/p003/LY_RAP_0003.csv
awk '$4="LY_RAP_0003-"$4' /Users/kisaev/Documents/Hatchet_analysis/p003/results/best.seg.ucn > /Users/kisaev/Documents/Hatchet_analysis/p003/cns_file.seg.ucn
seg_file=/Users/kisaev/Documents/Hatchet_analysis/p003/cns_file.seg.ucn
python /Users/kisaev/github/RAP_WGS/Part_2_Somatic_copy_number_calling/Hatchet/explainMutationsCCF.py $snvs -s $seg_file> p003_SNVs_explained_by_hatchet_output.txt
