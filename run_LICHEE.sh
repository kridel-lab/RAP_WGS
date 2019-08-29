#this tool: 
#was run locally in order to simplify visualization 
#Karin Isaev
#August 23rd, 2019

#need to run tool from /release directory

cd /Users/kisaev/lichee/LICHeE/release

lichee_input=/Users/kisaev/Documents/RAP_ANALYSIS/2019-08-23_liche_input_somatic_muts.txt	
less $lichee_input | head -6000 > test_lichee.txt

#run 
./lichee -build -i $lichee_input -n 0 -showTree 1 -color -dot -sampleProfile 
#test
./lichee -build -i test_lichee.txt -n 0 -showTree 1 -color -dot -sampleProfile 

#on the cluster 
lichee_input=/cluster/projects/kridelgroup/RAP_ANALYSIS/chr/vcfs_final/2019-08-23_liche_input_somatic_muts.txt	
./lichee -build -i $lichee_input -n 0 -showTree 1 -color -dot -sampleProfile 
