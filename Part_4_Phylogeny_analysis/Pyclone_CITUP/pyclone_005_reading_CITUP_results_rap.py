 #Within the hdf5 results, the /results/optimal can be used to
 #identify the id of the optimal tree solution.

 #The clonal phylogeny as an adjacency list is then
 #the /trees/{tree_solution}/adjacency_list
 #entry and the clone frequencies are the /trees/{tree_solution}/clone_freq entry in the hdf5 file.

 #The adjacency list can be written as a TSV with the column names source, target to be input into E-scape,

 #and the clone frequencies should be reshaped such that each row represents
 #a clonal frequency in a specific sample for a specific clone, with the columns
 #representing the time or space ID, the clone ID, and the clonal prevalence.

import numpy as np
import pandas as pd

#patient 003
hdf = pd.HDFStore('p003_results.h5', mode='r')
hdf.groups()

optimal = hdf.get('/results/optimal')
optimal

adjacency = hdf.get('/trees/7/adjacency_list') #source and target A-->B, B--C, A-->D for example (edges)
clone = hdf.get('/trees/7/clone_freq')

adjacency.to_csv("pat003_adjacency.csv", header=False)
clone.to_csv("pat003_clone.csv", header=False)
