#!/usr/bin/env python

import gzip
import bioinfo
import argparse
from itertools import product
import glob

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', help='Input directory with 4 fastq.gz files', type=str)
    parser.add_argument('-o', help='Output directory', type=str)
    parser.add_argument('-t', help='TSV file with a list of index sequences (5th column) and their associated samples (4th column)', type=str)
    parser.add_argument('-q', help='Quality score cutoff', type=int)
    return parser.parse_args()

args = get_args()

# Add trailing slash to directories if not present:
in_dir = args.i
if in_dir[-1] != '/':
    in_dir = in_dir + '/'

out_dir = args.o
if out_dir[-1] != '/':
    out_dir = out_dir + '/'

tsv = args.t
quality_score_cutoff = args.q

indices = set()
with open(tsv) as fh: # Populate indices using the TSV file
    for i, line in enumerate(fh):
        if i != 0: # Skip header
            index = line.strip().split('\t')[4] # Get the 4th item
            indices.add(index)

R1 = glob.glob(f'{in_dir}*R1*.fastq.gz')[0]
R2 = glob.glob(f'{in_dir}*R2*.fastq.gz')[0]
R3 = glob.glob(f'{in_dir}*R3*.fastq.gz')[0]
R4 = glob.glob(f'{in_dir}*R4*.fastq.gz')[0]


def correct_seq(seq: str, index_set: set) -> str:
    '''Takes a sequence string and a dict_keys object (of indices) as input.
    Replaces up to 1 N in the sequence with a nucleotide, and if the resulting sequence matches an index, returns that index.
    Otherwise returns an empty string.'''

    if seq.count('N') > 1:
        return ''

    for rep in (seq.replace('N', 'A'), seq.replace('N', 'T'), seq.replace('N', 'C'), seq.replace('N', 'G')):
        if rep in index_set:
            return rep

    return ''


def demultiplex(R1: str, R2: str, R3: str, R4: str, out_dir: str, index_set: set, qual_cutoff: int) -> tuple:
    '''Takes 4 input FASTQ file names (R1 = read 1, R2 = index 1, R3 = index 2, R4 = read 2),
    a list of indices, and a quality score cutoff value. Sorts read records based on indices
    into files labelled with sample names (separate files for read 1 and read 2).
    If either index 1 or the reverse complement of index 2 is not in index_dict,
    or if either does not meet the quality score cutoff, the records will be written to a
    pair of files for unknown reads. If both indices are in index_dict but do not match each other,
    the records will be written to a pair of files for index-hopped reads. Returns counts for how many records were dual matched,
    index-hopped, and unknown, and a dictionary with counts for all possible pairs of indices in index_dict.'''

    index_pairs = {} # Dictionary to keep track of how many times each possible index pair shows up
    for pair in product(index_set, index_set):
        index_pairs[pair[0] + '-' + pair[1]] = 0

    paired_count, hopped_count, unknown_count = 0, 0, 0

    print('Opening files...')
    with gzip.open(R1, 'rt') as r1, gzip.open(R2, 'rt') as r2, gzip.open(R3, 'rt') as r3, gzip.open(R4, 'rt') as r4, open(f'{out_dir}Unknown_R1.fastq', 'w') as u1,  open(f'{out_dir}Unknown_R2.fastq', 'w') as u2,  open(f'{out_dir}Hopped_R1.fastq', 'w') as h1, open(f'{out_dir}Hopped_R2.fastq', 'w') as h2:
        index_out_files = {}
        for index in index_set:
            index_out_files[index] = (open(f'{out_dir + index}_R1.fastq', 'w'), open(f'{out_dir + index}_R2.fastq', 'w'))
        print('Files opened')

        counter = 1
        while True:
            if counter % 36000000 == 0:
                print(f'Processing record {counter}')

            # if counter > 1: # DELETE AFTER TEST
            #     break


            counter += 1

            record_1 = ''
            record_2 = ''

            for i in range(4):
                record_1 += r1.readline()
                record_2 += r4.readline()

                if i == 1: # Index reads
                    index_1 = r2.readline().strip()
                    index_2 = bioinfo.reverse_complement(r3.readline().strip())
                elif i == 3: # Index quality scores
                    qscore_1 = r2.readline().strip()
                    qscore_2 = r3.readline().strip()
                else: # Skip
                    r2.readline()
                    r3.readline()

            if record_1 == '': # Reached the end of the files
                print('End of files reached')
                break

            # If the indices are in the set, set the corrected indices to be identical. Otherwise, call correct_seq().
            index_1_corrected = index_1 if index_1 in index_set else correct_seq(index_1, index_set)
            index_2_corrected = index_2 if index_2 in index_set else correct_seq(index_2, index_set)

            # print(f'Index 1: {index_1}')
            # print(f'Index 1 Corrected: {index_1_corrected}')
            # print(f'Index 2: {index_2}')
            # print(f'Index 2 Corrected: {index_1_corrected}')
            # print(f'Quality 1: {bioinfo.qual_score(qscore_1)}')
            # print(f'Quality 2: {bioinfo.qual_score(qscore_2)}')
            # print(f'Score 1: {qscore_1}')
            # print(f'Score 2: {qscore_1}')
            if index_1_corrected == '' or index_2_corrected == '' or bioinfo.qual_score(qscore_1) < qual_cutoff or bioinfo.qual_score(qscore_2) < qual_cutoff:
                # If correction didn't work or one of the indices has average quality below the cutoff
                # Add indices to headers and write records to unknown files:
                # print('Unknown')
                u1.write(record_1.replace('\n', f' {index_1}-{index_2}\n', 1))
                u2.write(record_2.replace('\n', f' {index_1}-{index_2}\n', 1))
                unknown_count += 1
                continue

            # String to append to header (I'm pretty sure when I asked Leslie she told me to use the corrected indices here,
                # but Hannah K said she was told to use the uncorrected ones.)
            # Also the second index is the reverse complement of the actual text in the R3 file,
                # so it will match an actual index
            pair = f'{index_1_corrected}-{index_2_corrected}'
            index_pairs[pair] += 1

            # Add indices to headers:
            record_1 = record_1.replace('\n', f' {pair}\n', 1)
            record_2 = record_2.replace('\n', f' {pair}\n', 1)

            if index_1_corrected == index_2_corrected:
                # Look up file handles for the index and write records to them:
                files = index_out_files[index_1_corrected]
                files[0].write(record_1)
                files[1].write(record_2)
                paired_count += 1
                continue
            else:
                # Write records to hopped files:
                h1.write(record_1)
                h2.write(record_2)
                hopped_count += 1
                continue

        # Close files that weren't opened by 'with open() as':
        for file in index_out_files.values():
            file[0].close()
            file[1].close()

    return paired_count, hopped_count, unknown_count, index_pairs


paired_count, hopped_count, unknown_count, index_pairs = demultiplex(R1, R2, R3, R4, out_dir, indices, quality_score_cutoff)

with open(out_dir + 'summary.tsv', 'w') as fh:
    fh.write(f'Paired Reads\t{paired_count}\n')
    fh.write(f'Hopped Reads\t{hopped_count}\n')
    fh.write(f'Unknown Reads\t{unknown_count}\n')

    for pair in index_pairs:
        fh.write(f'{pair}\t{index_pairs[pair]}\n')