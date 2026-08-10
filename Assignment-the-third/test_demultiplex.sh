#!/bin/bash

#SBATCH -A bgmp
#SBATCH -p bgmp
#SBATCH --job-name=test_demultiplex_%j
#SBATCH --output=test_demultiplex_%j.out
#SBATCH --error=test_demultiplex_%j.out

in_dir='/projects/bgmp/imre/bioinfo/bi622/Demultiplex/TEST-input_FASTQ/'
out_dir='/projects/bgmp/imre/bioinfo/bi622/Demultiplex/TEST-program_output_FASTQ/'

/usr/bin/time -v python /projects/bgmp/imre/bioinfo/bi622/Demultiplex/Assignment-the-third/demultiplex.py -i $in_dir -o $out_dir -t /projects/bgmp/shared/2017_sequencing/indexes.txt -q 25