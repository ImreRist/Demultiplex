#!/bin/bash

#SBATCH -A bgmp
#SBATCH -p bgmp
#SBATCH --job-name=demultiplex_%j
#SBATCH --output=demultiplex_%j.out
#SBATCH --error=demultiplex_%j.out

in_dir='/projects/bgmp/shared/2017_sequencing/'
out_dir='/scratch/bgmp/imre/demux/'

/usr/bin/time -v python /projects/bgmp/imre/bioinfo/bi622/Demultiplex/Assignment-the-third/demultiplex.py -i $in_dir -o $out_dir -t /projects/bgmp/shared/2017_sequencing/indexes.txt -q 25