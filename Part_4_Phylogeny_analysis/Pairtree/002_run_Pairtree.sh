#!/bin/bash
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J run pairtree

conda deactivate
conda activate pairtree

#set up input files
input_files=/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/Pairtree/2021-05-12_input_files
cd $input_files

#set-up pairtree
PTDIR=$HOME/pairtree
#mkdir results
cd results

#nano  $input_files/p001_input.params.json to remove "" from garabge variable

# Run Pairtree.
$PTDIR/bin/pairtree --params $input_files/p001_input.params.json $input_files/p001_ssm_input.ssm p001.results.npz

# Plot results in an HTML file.
$PTDIR/bin/plottree --runid p001 $input_files/p001_ssm_input.ssm $input_files/p001_input.params.json p001.results.npz p001.results.html

# View the HTML file.
#firefox example.results.html
