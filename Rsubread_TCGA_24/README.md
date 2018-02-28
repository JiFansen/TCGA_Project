TCGA_Project
============
Rsubread_TCGA_24:
-----------------
	In this directory, we will use the Rsubread Data to do the analysis.

-----------------------------------------------------------------------------
The DataSet is in this location: 
	/Share/home/lanxun5/Data/TCGA/FPKM_9000/Rsubread_TCGA_24

Run the following command:
	1. bash sampleID_2_Cancertype.sh
	2. run the R scripts.
	
① sample.sh: 
	This script will extract the samples header i the expression matrix, and then display according to the column. 

② match.pl: 
	This script will match the TCGA samples to their corresponding cancer types.

③ sampleID_2_Cancertype.sh : 
	This script will organize the functions above.

④ Rsubread_No_Correction_Plot.R: 
	This R script will use all the genes or certain gene lists to plot the cluster, no correction has been done, normalize only.

⑤ Rsubread_lymphocyte_Correction_Plot.R: 
	This R script will use the slides data to correct the immune gene expression, then plot the cluster.
=============================================================================