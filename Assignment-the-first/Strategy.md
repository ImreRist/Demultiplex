### Demultiplex Outline

The algorithm needs to demultiplex Illumina reads, sorting them into separate sample files based on paired indices.

The input is 4 FASTQ files (Read 1, Index 1, Read 2, Index 2) and a TSV file containing an index/sample key.

The files are located in `/projects/bgmp/shared/2017_sequencing` on Talapas.

The algorithm should write each record to a new FASTQ file based on the indices. Paired reads with matching indices that match one of the samples will be written to R1 and R2 FASTQ files for that sample. Paired reads with with indices that match the samples but don't match each other are considered index hopped, and will be written to R1 and R2 hopped files. Paired reads with at least one index that does not match one of the samples will be written to R1 and R2 unknown files.

There will be 2 FASTQ files for each matching index pair, plus 2 hopped and 2 unknown FASTQ files. With 24 indices, there will be a total of 52 output files (24*2 + 2 + 2).

### Unit Tests

I want to include at least the following records in the unit tests:
- 2 different sets of matching indices that are in the list (should go to their respective sample files)
- 1 set of indices that match and are in the list, but have Ns (make the N positions different to test the correct_seq funciton more robustly). Should go to the sample files (make it the same as one of the above to keep outpus files low)
- 1 set of indices that are in the list but do not match (should go to hopped files)
- 1 set of indices that are in the list but do not match, with some Ns (should go to hopped files)
- 1 set of indices that match but are both not in the list (unknown)
- 1 set of indices that do not match, with one being in the list and one not (unknown)
- 1 set of indices that do not match, with neither being in the list (unknown)
- 1 set with the index 'NTAGCNNA' to test how correct_seq handles ambiguity (unknown)

Should end up with 4 unknown, 3 matched (between 2 samples), and 2 hopped. 8 output files total.

Indices used:
TAGCCATG / TAGCCATG - C9
GCTACTCT / GCTACTCT - B2
TAGNCANG / NANCCATG - C9
CGATCGAT / TATGGCAC - Hopped (A5 / B7)
CGATCGAT / TANGNCAC - Hopped (A5 / B7)
ATGATCGA / ATGATCGA - Unknown (X)
TCGAGAGT / ATGATCGA - Unknown (A10 / X)
ATGATCGA / TGATGATC - Unknown (X / X)
NTAGCNNA / NTAGCNNA - Unknown (B1/A11 / B1/A11)

### Psuedocode

```
import gzip
import bioinfo

Define file paths (strings) R1, R2, R3, R4, and indexes (either explicitly or using argparse)
Define quality score cutoff (int) (either explicitly or using argparse)

define make_index_dict function, input is TSV path:
    Create empty dicitonary
    For each line/sample, add an entry to the dictionary with the index as the key and the sample name as the value
    Return the dictionary

Define correct_seq function:
    Input is a sequence (string) and a dict_keys object (containing all indices that correspond to our samples)
    If the sequence unambiguously matches one of the indices (i.e. treating Ns as any base), return that index
    If it doesn't match any or matches more than one, return an empty string

    Note: The minimum hamming distance between any 2 indices in the list is 3, and ~90% of pairs have a hamming distance of at least 5 (of 8 total bases), so I think it's fairly safe to do this. The only way a match will be ambiguous is if at least 3/8 bases are N, at which point hopefully the quality score check will flag it.


Define demultiplex function:
    Takes 4 input file names, an index dictionary, and a quality score cutoff as parameters
    
    Initialize index_pairs dictionary with each possible pair of indices as keys and counts as values, all initialized to 0 (itertools.product() may be useful). Keys will be the strings which will be appended to the header lines.
    
    Initialize matched_count, hopped_count, and unknown_count to 0

    Open 52 output FASTQ files:
        Open 4 input FASTQ files:
            while True
                initialize 2 empty strings to store the records from each read file
                for i in range(4)
                    readline for R1 and R4, append to the appropriate string (i.e. record += line)
                    if i == 1:
                        readline for R2 and R3, store to index1 and index2 variables (set index2 to the reverse complement of the R3 index)
                    else:
                        readline for R2 and R3, don't store anywhere

                If the R1 string is empty, close files (if necessary) and return (we've reached the end of the file)
                    The return value will be the 3 count variables and the index_pairs dictionary

                Initialize index1_corrected and index2_corrected to index1 and index2

                If index1 isn't in index_dict:
                    Call correct_seq on index1, save to index1_corrected
                ^ Do the same for index2

                If index1_corrected or index2_corrected is an empty string OR either index's quality score doesn't meet the cutoff:
                    Note: not sure if it should assess average quality score or the individual scores. Will use functions from bioinfo

                    Insert the indices (uncorrected) at the end of the header and write the read records to the unknown files
                    Add 1 to unknown_count
                    Continue (iterate the while loop)

                Add 1 to the index pair's count in the index_pairs dictionary

                Insert the indices at the end of the headers
                    Use str.replace() with count=1 to replace the first newline character with the thing I want to append plus a newline
                
                If index1_corrected == index2_corrected:
                    Write the read records to the files for that sample (lookup the sample using index_dict)
                    Add 1 to matched_count
                    Continue
                Else:
                    Write the read records to the hopped files
                    Add 1 to hopped_count
                    Continue

Call demultiplex function with input file names, dictionary from make_index_dict (using TSV file), and quality cutoff as arguments
    Report counts for matched/hopped/unknown as well as counts for each index pair (print or write to summary file)
```

### Functions

```python
def make_index_dict(file: str) -> dict:
    '''Takes a TSV file path (string) as input and returns a dictionary with indices as keys and sample names as values, taken from the 5th and 4th columns of the TSV file respectively, excluding the header line'''
    return index_dict
# Input: '/projects/bgmp/shared/2017_sequencing/indexes.txt'
# Expected output: {'GTAGCGTA': 'B1', 'CGATCGAT': 'A5', 'GATCAAGG': C1, ... , 'AGGATAGC': 'A8'}

def reverse_complement(seq: str) -> str: # Put in bioinfo.py
    '''Takes a sequence string, returns the reverse complement'''
    return reverse_complement
# Input: 'ATCNTGAC'
# Expected output: 'GTCANGAT'

def correct_seq(seq: str, indices: <class 'dict_keys'>) -> str:
    '''Takes a sequence string and a dict_keys object (of indices) as input. If the sequence unambiguosly matches one of the indices, return that index. If the sequence matches zero or more than one index, return an empty string.'''
    return corrected_seq
# Input: 'CTNTGGAN', dict_keys(['CTCTGGAT'])
# Expected output: 'CTCTGGAT'
# Input: 'NTAGCNNA', dict_keys(['CTAGCTCA', 'GTAGCGTA'])
# Expected output: ''

def demultiplex(R1: str, R2: str, R3: str, R4: str, index_dict: dict, qual_cutoff: int) -> tuple:
    '''Takes 4 input FASTQ file names (R1 = read 1, R2 = index 1, R3 = index 2, R4 = read 2), a dictionary with indices as keys and sample names as values, and a quality score cutoff value. Sorts read records based on indices into files labelled with sample names (separate files for read 1 and read 2). If either index 1 or the reverse complement of index 2 is not in index_dict, or if either does not meet the quality score cutoff, the records will be written to a pair of files for unknown reads. If both indices are in index_dict but do not match each other, the records will be written to a pair of files for index-hopped reads. Returns counts for how many records were dual matched, index-hopped, and unknown, and a dictionary with counts for all possible pairs of indices in index_dict.'''
    return matched_count, hopped_count, unknown_count, index_pairs
# Input: '/projects/bgmp/imre/bioinfo/bi622/Demultiplex/TEST-input_FASTQ/R1_test.fastq', '/projects/bgmp/imre/bioinfo/bi622/Demultiplex/TEST-input_FASTQ/R2_test.fastq', '/projects/bgmp/imre/bioinfo/bi622/Demultiplex/TEST-input_FASTQ/R3_test.fastq', '/projects/bgmp/imre/bioinfo/bi622/Demultiplex/TEST-input_FASTQ/R4_test.fastq', make_index_dict('/projects/bgmp/shared/2017_sequencing/indexes.txt')
# Expected output: (3, 2, 4, {'TAGCCATG-TAGCCATG': 2, 'GCTACTCT-GCTACTCT': 1, 'CGATCGAT-TATGGCAC': 2, 'GTAGCGTA-GTAGCGTA': 0, ... (576 pairs in total, all others should have a value of 0) ..., 'AGGATAGC-AGGATAGC': 0})
```