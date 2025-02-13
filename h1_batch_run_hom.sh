#!/bin/bash

#SBATCH -N 1
#SBATCH -n 24
#SBATCH --qos=normal
#SBATCH --job-name=amygdala_theta
#SBATCH --output=amygdala_batch_%j.out
#SBATCH --time 0-12:00

START=$(date)
mpiexec nrniv -mpi -python run_network.py simulation_configECP_base_homogenous.json
#mpiexec ./components_homogenous/mechanisms/x86_64/special -mpi run_network.py simulation_configECP_base_homogenous.json
END=$(date)

printf "Start: $START \nEnd:   $END\n" 
python analysis_hom.py --save-plots

echo "Done running model at $(date)"
