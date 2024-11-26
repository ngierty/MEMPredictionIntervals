#!/bin/tcsh
#BSUB -J mr[1-300]  #job name 
#BSUB -n 8          #number of threads
#BSUB -W 04:00          #walltime limit: hh:mm
#BSUB -R span[hosts=1]
#BSUB -o /share/hmmrs/ngierty/meas_error/cluster_out/Output_%J_%I.out 
#BSUB -e /share/hmmrs/ngierty/meas_error/cluster_out/Error_%J_%I.err  #error - %J is the job-id %I is the job-array index 

setenv XDG_RUNTIME_DIR $TMP

conda activate /usr/local/usrapps/hmmrs/ngierty/meas_err_env


Rscript MR_02_calc_ncss.R --mc $LSB_JOBINDEX


conda deactivate
