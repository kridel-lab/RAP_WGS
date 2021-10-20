#!/bin/bash
#SBATCH -p long
#SBATCH --mem=30720M
#SBATCH -J pairtree
#SBATCH -c 10
#SBATCH -t 21-00:00 # Runtime in D-HH:MM

#echo ". /cluster/tools/software/python/3.6.5/etc/profile.d/conda.sh" >> ~/.bashrc

#conda deactivate
#conda activate pairtree

source /cluster/home/kisaev/.bashrc
source activate pairtree

#set up input files
input_files=/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/Pairtree/2021-10-18_input_files
cd $input_files

#don't forget to manually edit input *params.json files to remove the "" in the garbage []

#set-up pairtree
PTDIR=$HOME/pairtree
cd min100_muts

mkdir final_chosen_tree

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#P001+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Plot best tree results in an HTML file.
$PTDIR/bin/plottree --runid p001 $input_files/p001_ssm_input.ssm $input_files/p001_input.params.json p001.results.npz final_chosen_tree/p001.results.html --reorder-subclones --tree-json final_chosen_tree/p001_solution.json --tree-index 1 --remove-normal

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#P002+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Plot best tree results in an HTML file.
#$PTDIR/bin/plottree --runid p002 $input_files/p002_ssm_input.ssm $input_files/p002_input.params.json p002.results.npz final_chosen_tree/p002.results.html --reorder-subclones --tree-json final_chosen_tree/p002_solution.json --tree-index 0 --remove-normal

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#P003+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Plot best tree results in an HTML file.
#$PTDIR/bin/plottree --runid p003 $input_files/p003_ssm_input.ssm $input_files/p003_input.params.json p003.results.npz final_chosen_tree/p003.results.html --reorder-subclones --tree-json final_chosen_tree/p003_solution.json --tree-index 7 --remove-normal
