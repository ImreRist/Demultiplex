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

def correct_seq(seq: str, barcodes: list) -> str:
    '''Takes a sequence string and a dict_keys object (of indices) as input.
    If the sequence unambiguosly matches one of the indices, return that index.
    If the sequence matches zero or more than one index, return an empty string.'''

    # Needs to be adjusted to accomodate the set, try not to use lists(?)
    # Maybe have an explicit cutoff so if there are more than 2 Ns, return an empty string
    n_indices = [i for i,  base in enumerate(seq) if base == 'N']
    barcode_list = list(barcodes)
    for i in n_indices:
        for n, barcode in enumerate(barcode_list):
            barcode_list[n] = barcode[:i] + 'N' + barcode[i+1:]

    if barcode_list.count(seq) == 1:
        return barcodes[barcode_list.index(seq)]
    else:
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

    index_pairs = {}
    for pair in product(index_set, index_set):
        index_pairs[pair[0] + '-' + pair[1]] = 0

    paired_count, hopped_count, unknown_count = 0, 0, 0

    with open(R1, 'r') as r1, open(R2, 'r') as r2, open(R3, 'r') as r3, open(R4, 'r') as r4, open(f'{out_dir}Unknown_R1.fq', 'w') as u1,  open(f'{out_dir}Unknown_R2.fq', 'w') as u2,  open(f'{out_dir}Hopped_R1.fq', 'w') as h1, open(f'{out_dir}Hopped_R2.fq', 'w') as h2:
        index_out_files = {}
        for index in index_set:
            index_out_files[index] = (open('{index}_R1.fq', 'w'), open('{index}_R2.fq', 'w'))


        while True:

            record_1 = ''
            record_2 = ''

            for i in range(4):
                record_1 += r1.readline()
                record_2 += r4.readline()

                if i == 1:
                    index_1 = r2.readline().strip()
                    index_2 = bioinfo.reverse_complement(r3.readline().strip())
                elif i == 2:
                    qscore_1 = r2.readline().strip()
                    qscore_2 = r3.readline().strip()
                else:
                    r2.readline()
                    r3.readline()

            if record_1 == '':
                break

            index_1_corrected = index_1
            index_2_corrected = index_2

            if index_1 not in index_set:
                index_1_corrected = correct_seq(index_1, index_set)

            if index_2 not in index_set:
                index_2_corrected = correct_seq(index_2, index_set)

            if index_1_corrected == '' or index_2_corrected == '' or (QSCORE CONDITION):
                # Add indices to headers and write records to unknown files: 
                u1.write(record_1.replace('\n', f' {index_1}-{index_2}\n', 1))
                u2.write(record_2.replace('\n', f' {index_1}-{index_2}\n', 1))
                unknown_count += 1
                continue

            pair = f'{index_1_corrected}-{index_2_corrected}'
            index_pairs[pair] += 1

            record_1 = record_1.replace('\n', f' {pair}\n', 1)
            record_2 = record_2.replace('\n', f' {pair}\n', 1)

            if index_1_corrected == index_2_corrected:
                files = index_out_files[index_1_corrected]
                files[0].write(record_1)
                files[1].write(record_2)
                paired_count += 1
                continue
            else:
                h1.write(record_1)
                h2.write(record_2)
                hopped_count += 1
                continue

        for file in index_out_files.values():
            file[0].close()
            file[1].close()


    return paired_count, hopped_count, unknown_count, index_pairs


print(demultiplex(R1, R2, R3, R4, out_dir, indices, 30))

# Write report values to a file