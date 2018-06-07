#!/bin/bash
#SBATCH -n 1
#SBATCH -t 00:30:00
#SBATCH -p main
#SBATCH --export=ALL
#SBATCH -D /home/ag1349/Roberto/latencyPorj/model 
#SBATCH -o matlab_ks-%j.out
#SBATCH -e matlab_ks-%j.err

module load MATLAB/R2017a 

matlab -nodisplay < latency_kstest.m 
