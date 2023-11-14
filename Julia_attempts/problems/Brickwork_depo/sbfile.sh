#!/bin/bash -l

# Number of CPUs per task
#SBATCH --cpus-per-task=1
# Walltime (job duration)
#SBATCH --time=1:00:00
# Assign to particular partition
#SBATCH --partition=standard
#SBATCH --mem=4G

# Then finally, our code we want to execute. 
julia pt_mpo.jl