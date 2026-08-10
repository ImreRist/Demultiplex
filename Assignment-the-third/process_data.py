#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt
import numpy as np

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', help='Summary TSV input file', type=str)
    return parser.parse_args()

in_file = get_args().i

pairs = {}
with open(in_file, 'r') as fh:
    for i, line in enumerate(fh):
        if i == 0:
            paired_reads = int(line.strip().split('\t')[1])
        elif i == 1:
            hopped_reads = int(line.strip().split('\t')[1])
        elif i == 2:
            unknown_reads = int(line.strip().split('\t')[1])
        else:
            indices, count = line.strip().split('\t')
            pairs[indices] = int(count)

total = paired_reads + hopped_reads + unknown_reads # type: ignore
index_counts = {}
hopped_counts = {}
for index_pair in pairs: # Isolate the counts with matching indices
    i1, i2 = index_pair.split('-')
    if i1 == i2:
        index_counts[i1] = pairs[index_pair]
    else:
        hopped_counts[index_pair] = pairs[index_pair]

sorted_indices = sorted(index_counts, key=index_counts.get, reverse=True) # type: ignore
sorted_index_counts = sorted(index_counts.values(), reverse=True)
index_pcts = np.array(sorted_index_counts) / total * 100

sorted_hopped_indices = sorted(hopped_counts, key=hopped_counts.get, reverse=True) # type: ignore

plt.bar(sorted_indices, index_pcts)
plt.title('Percentage of Total Records From Each Sample')
plt.xlabel('Sample Index')
plt.ylabel('Percent of Total Records')
plt.xticks(rotation=-50, ha='left')
plt.savefig('index_pcts.png', bbox_inches='tight')

with open('results.md', 'w') as fh:
    fh.write('## Results\n')
    fh.write('### General Stats\n')
    fh.write('|Record Type|Count|Percent|\n')
    fh.write('|-----------|-----|-------|\n')
    fh.write(f'|Paired|{paired_reads}|{round(paired_reads/total * 100, 2)}|\n')
    fh.write(f'|Hopped|{hopped_reads}|{round(hopped_reads/total * 100, 2)}|\n')
    fh.write(f'|Unknown|{unknown_reads}|{round(unknown_reads/total * 100, 2)}|\n')
    fh.write(f'|Total|{total}|100|\n')
    fh.write('### Samples\n')
    fh.write('|Index|Count|Percent|\n')
    fh.write('|-----|-----|-------|\n')
    for index in sorted_indices:
        fh.write(f'|{index}|{index_counts[index]}|{round(index_counts[index]/total * 100, 2)}|\n')
    fh.write('\n<img src="index_pcts.png">\n\n')
    fh.write('### Hopped Records\n\n')
    fh.write('|Index Pair|Count|\n')
    fh.write('|-----|-----|\n')
    for pair in sorted_hopped_indices:
        fh.write(f'|{pair}|{hopped_counts[pair]}|\n')