#!/bin/bash

#SBATCH -N 1
#SBATCH -n 50
#SBATCH --qos=normal
#SBATCH --job-name=amygdala_theta
#SBATCH --output=amygdala_batch.out
#SBATCH --time 0-24:00

START=$(date)
mpiexec nrniv -mpi -quiet -python run_network.py simulation_configECP_vpsi_homogenous.json
#mpiexec ./components_homogenous/mechanisms/x86_64/special -mpi run_network.py simulation_configECP_vpsi_homogenous.json
END=$(date)

printf "Start: $START \nEnd:   $END\n"
python analysis_hom.py --save-plots 




