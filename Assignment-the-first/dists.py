#!/usr/bin/env python

import bioinfo
import numpy as np
import matplotlib.pyplot as plt
import gzip
import multiprocessing


def average_quality_bar(args: tuple):
    '''Creates bar chart of average quality score per base.
    Args is a tuple of the following parameters:
    in_file is input FASTQ, out_file is output image name,
    prefix is what to start the chart title with (e.g. "Read 1").'''

    in_file, out_file, prefix = args
    print(f"Function starting with {prefix}")

    with gzip.open(in_file, 'rt') as fh: # add gzip
        first = True
        count = 0
        for i, line in enumerate(fh):
            if i == 1000:
                print(f"Line {i} in {prefix}")
            if i % 30000000 == 0:
                print(f"Line {i} in {prefix}")

            if i % 4 == 3: # Score line
                line = line.strip()
                if first:
                    scores_sum = np.zeros(len(line))
                    first = False
                scores = []
                for char in line:
                    scores.append(bioinfo.convert_phred(char)) # type: ignore
                scores_sum = scores_sum + np.array(scores)
                count += 1

    plt.figure(figsize=((20, 5) if len(scores_sum) > 50 else (6.4, 4.8)))
    plt.bar(range(len(scores_sum)), scores_sum/count)
    plt.title(f"{prefix} Average Quality Score Per Base")
    plt.xlabel("Base Index")
    plt.ylabel("Average Quality Score")
    plt.savefig(out_file)

    print(f"Finished {prefix}")
    
    return

# Test files
# read1 = '../TEST-input_FASTQ/R1_test.fastq'
# index1 = '../TEST-input_FASTQ/R2_test.fastq'
# index2 = '../TEST-input_FASTQ/R3_test.fastq'
# read2 = '../TEST-input_FASTQ/R4_test.fastq'

# Paths are hardcoded for initial exploration, I know it's not good practice generally
read_path = '/projects/bgmp/shared/2017_sequencing/'
read1 = read_path + '1294_S1_L008_R1_001.fastq.gz'
index1 = read_path + '1294_S1_L008_R2_001.fastq.gz'
index2 = read_path + '1294_S1_L008_R3_001.fastq.gz'
read2 = read_path + '1294_S1_L008_R4_001.fastq.gz'

print(f'Running with {multiprocessing.cpu_count()} CPUs')
out_path = '/projects/bgmp/imre/bioinfo/bi622/Demultiplex/Assignment-the-first/'
if __name__ == '__main__':
    print('Main')
    args_list = [(read1, out_path + 'read1_hist.png', 'Read 1'), (index1, out_path + 'index1_hist.png', 'Index 1'), (index2, 'index2_hist.png', 'Index 2'),(read2, 'read2_hist.png', 'Read 2')]
    pool = multiprocessing.Pool()
    print('Pool created')
    pool.map(average_quality_bar, args_list)