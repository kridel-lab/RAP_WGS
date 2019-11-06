	#setting_up_symlinks.sh
	#author: Karin Isaev
	#date started: June 7th , 2019
	#June 21st, now more files have been downloaded and they are spread out across multiple folders...

	#generate symlinks for all tsv files since they are all in different folders
	#1 normal control sample from previous sequencing run (March 2019)
	#6 samples from recent June sequencing upload in TCAG (June 6 2019)
	#14 samples from most recent June sequencing upload in TCGA (June 18 2019)

	mkdir RAP_ANALYSIS
	cd RAP_ANALYSIS

	#RAP = 1 control sample from March 2019
	ln -s /cluster/projects/kridelgroup/RAP/BAK8973/LY_RAP_0003_Ctl_FzG_01_files/ /cluster/projects/kridelgroup/RAP_ANALYSIS/

	#VO18AHH = 6 samples 
	#/cluster/projects/kridelgroup/VO18AHH/BAK8973_topup_analysis
	ln -s /cluster/projects/kridelgroup/VO18AHH/BAK8973_topup_analysis/* /cluster/projects/kridelgroup/RAP_ANALYSIS

	#IA6CWXM = 14 samples 
	#/cluster/projects/kridelgroup/IA6CWXM/BAK10122

	ln -s /cluster/projects/kridelgroup/IA6CWXM/BAK10122/* /cluster/projects/kridelgroup/RAP_ANALYSIS

	tree

	#21 files total 
	