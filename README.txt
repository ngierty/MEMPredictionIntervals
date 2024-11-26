README: 

MEM Prediction Intervals GitHub Repository Outline:

	functions.R contains all of our functions that are needed to generate data, estimate the regression relationship, and calculate the non-conformity scores


	sims_01_calc_ncss.R calculates non-conformity scores for a specific scenario of simulated data described in the manuscript
	sims_01_plots.R generates the plots used in the manuscript that were based on the simulated data
	run_sims.sh runs the simulation for a variety of assumptions on the measurement errors

	

	MR_01_prep.R cleans the data obtained from the Plantery Systems table from the NASA Exoplanet Archive
	MR_02_calc_ncss.R calculates non-conformity scores for the mass-radius relation known in the literature
	MR_03_output.R generates the plots used in the manuscript that were based on the real data
	run_mr.sh runs MR_02_calc_ncss.R to produce the non-conformity scores needed as input for MR_03_output.R




RAW folder:

	NASA_MRF_confirmed.csv: data from Ma (2021) and used as a check against the actual data used
	NASA_exoplanet_archive.csv: data from NASA Exoplanet Archive that was used for the real data analysis




INTER folder:

	inter/sim_data: contains the intermediate output from all the sims* code described above
	inter/mr_data: contains the intermediate output from all the MR* code described above



OUTPUT folder:

	output/sim_data: contains the final output from all the sims* code described above
	output/mr_data: contains the final output from all the MR* code described above