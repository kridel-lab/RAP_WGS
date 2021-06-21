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
input_files=/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/Pairtree/2021-06-21_input_files
cd $input_files

#don't forget to manually edit input *params.json files to remove the "" in the garbage []

#set-up pairtree
PTDIR=$HOME/pairtree
mkdir min150_muts
cd min150_muts

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#P001+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#nano  $input_files/p001_input.params.json to remove "" from garabge variable

# Run Pairtree.
$PTDIR/bin/pairtree --params $input_files/p001_input.params.json $input_files/p001_ssm_input.ssm p001.results.npz --seed 123 --trees-per-chain 4000

# Plot best tree results in an HTML file.
$PTDIR/bin/plottree --runid p001 $input_files/p001_ssm_input.ssm $input_files/p001_input.params.json p001.results.npz p001.results.html --reorder-subclones

#plot all posteriors
$PTDIR/bin/summposterior --runid p001 $input_files/p001_ssm_input.ssm $input_files/p001_input.params.json p001.results.npz p001.posterior_plots.html

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#P002+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#nano  $input_files/p002_input.params.json to remove "" from garabge variable

# Run Pairtree.
$PTDIR/bin/pairtree --params $input_files/p002_input.params.json $input_files/p002_ssm_input.ssm p002.results.npz --seed 123 --trees-per-chain 4000

# Plot best tree results in an HTML file.
$PTDIR/bin/plottree --runid p002 $input_files/p002_ssm_input.ssm $input_files/p002_input.params.json p002.results.npz p002.results.html --reorder-subclones

#plot all posteriors
$PTDIR/bin/summposterior --runid p002 $input_files/p002_ssm_input.ssm $input_files/p002_input.params.json p002.results.npz p002.posterior_plots.html

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#P003+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#nano  $input_files/p003_input.params.json to remove "" from garabge variable

# Run Pairtree.
$PTDIR/bin/pairtree --params $input_files/p003_input.params.json $input_files/p003_ssm_input.ssm p003.results.npz --seed 123 --trees-per-chain 4000

# Plot best tree results in an HTML file.
$PTDIR/bin/plottree --runid p003 $input_files/p003_ssm_input.ssm $input_files/p003_input.params.json p003.results.npz p003.results.html --reorder-subclones

#plot all posteriors
$PTDIR/bin/summposterior --runid p003 $input_files/p003_ssm_input.ssm $input_files/p003_input.params.json p003.results.npz p003.posterior_plots.html
