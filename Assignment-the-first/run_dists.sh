#!/bin/bash

#SBATCH -A bgmp
#SBATCH -p bgmp
#SBATCH --cpus-per-task=4
#SBATCH --job-name=dists_%j
#SBATCH --output=dists_%j.out
#SBATCH --error=dists_%j.out

python /projects/bgmp/imre/bioinfo/bi622/Demultiplex/Assignment-the-first/dists.py