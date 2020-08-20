#!/bin/bash
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -J citup
#SBATCH -c 8
#SBATCH -t 5-00:00 # Runtime in D-HH:MM

cd /cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/Pyclone

#running citup
module load citup
module load python3

run_citup_qip.py 2020-08-20_rap_citup_input_freqs.txt \
2020-08-20_rap_citup_input_clusters.txt p003_results.h5 --min_nodes 2 --max_nodes 5 \
--submit local --max_children_per_node 5

#usage: run_citup_iter.py [-h] [--tmpdir TMPDIR] [--pipelinedir PIPELINEDIR]
#                         [--loglevel {CRITICAL,ERROR,WARNING,INFO,DEBUG}]
#                         [--submit SUBMIT] [--submit_config SUBMIT_CONFIG]
#                         [--nativespec NATIVESPEC] [--storage STORAGE]
#                         [--storage_config STORAGE_CONFIG] [--maxjobs MAXJOBS]
#                         [--repopulate] [--rerun] [--nocleanup]
#                         [--interactive] [--sentinel_only]
#                         [--context_config CONTEXT_CONFIG]
#                         [--min_nodes MIN_NODES] [--max_nodes MAX_NODES]
#                         [--max_children_per_node MAX_CHILDREN_PER_NODE]
#                         input_freqs output_results

#--min_nodes MIN_NODES
#                      Output For All Trees (default: 1)
#--max_nodes MAX_NODES
#                      Output For All Trees (default: 8)
#--max_children_per_node MAX_CHILDREN_PER_NODE
#                      Output For All Trees (default: 100)
