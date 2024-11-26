#!/bin/tcsh
#BSUB -J sims[1-300]  #job name 
#BSUB -n 8          #number of threads
#BSUB -W 04:00          #walltime limit: hh:mm
#BSUB -R span[hosts=1]
#BSUB -o /share/hmmrs/ngierty/meas_error/cluster_out/Output_%J_%I.out 
#BSUB -e /share/hmmrs/ngierty/meas_error/cluster_out/Error_%J_%I.err  #error - %J is the job-id %I is the job-array index 

setenv XDG_RUNTIME_DIR $TMP

conda activate /usr/local/usrapps/hmmrs/ngierty/meas_err_env


Rscript sims_01_calc_ncss.R --mc $LSB_JOBINDEX --beta '1,5,4,7' --tau 'log(3),0' --zetas '10,25' --nobs 100 --dist 'norm'
Rscript sims_01_calc_ncss.R --mc $LSB_JOBINDEX --beta '1,5,4,7' --tau 'log(3),0.2' --zetas '10,25' --nobs 100 --dist 'norm'



Rscript sims_01_calc_ncss.R --mc $LSB_JOBINDEX --beta '1,5,4,7' --tau 'log(3),0' --zetas '10,25' --nobs 100 --dist 'gnorm' --dist_param 0.365
Rscript sims_01_calc_ncss.R --mc $LSB_JOBINDEX --beta '1,5,4,7' --tau 'log(3),0.2' --zetas '10,25' --nobs 100 --dist 'gnorm' --dist_param 0.365


Rscript sims_01_calc_ncss.R --mc $LSB_JOBINDEX --beta '1,5,4,7' --tau 'log(3),0' --zetas '10,25' --nobs 100 --dist 't' --dist_param 3
Rscript sims_01_calc_ncss.R --mc $LSB_JOBINDEX --beta '1,5,4,7' --tau 'log(3),0.2' --zetas '10,25' --nobs 100 --dist 't' --dist_param 3


conda deactivate
