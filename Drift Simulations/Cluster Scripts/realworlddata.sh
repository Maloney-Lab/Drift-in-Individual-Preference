#!/bin/bash
#SBATCH -n 24 # Number of cores requested
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 220 # Runtime in minutes
#SBATCH -p shared # Partition to submit to
#SBATCH --mem-per-cpu=1536 # Memory per node in MB (see also --mem-per-cpu)
#SBATCH --open-mode=append # Append when writing files
#SBATCH --output=Logs/hostname_%j_array_%A-%a.out    # Standard output and error log
#SBATCH --array=30-40               # Array range
#SBATCH 
module load Anaconda3/2019.10
source activate drift
python realworlddata.py