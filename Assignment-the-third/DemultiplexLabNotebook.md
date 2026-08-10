## Demultiplex Lab Notebook

### Data Exploration (Part 1)

FASTQ files and index info table are in `/projects/bgmp/shared/2017_sequencing/` on Talapas.

Bash commands:
```bash
cd /projects/bgmp/shared/2017_sequencing/
ls
zcat 1294_S1_L008_R1_001.fastq.gz | less
zcat 1294_S1_L008_R2_001.fastq.gz | less
zcat 1294_S1_L008_R3_001.fastq.gz | less
zcat 1294_S1_L008_R4_001.fastq.gz | less
```

It looks like R1 and R4 contain reads, while R2 and R3 contain indices. The Illumina notes suggest that R2 will correspond to R1 and R3 to R4.

|File|Data|
|----|----|
|1294_S1_L008_R1_001.fastq.gz|read1|
|1294_S1_L008_R2_001.fastq.gz|index1|
|1294_S1_L008_R3_001.fastq.gz|index2|
|1294_S1_L008_R4_001.fastq.gz|read2|

```bash
zcat 1294_S1_L008_R1_001.fastq.gz | head -2 | tail -1 | wc -c
# 102
zcat 1294_S1_L008_R4_001.fastq.gz | head -2 | tail -1 | wc -c
# 102
zcat 1294_S1_L008_R1_001.fastq.gz | head -2 | tail -1 | wc -c
# 9
zcat 1294_S1_L008_R1_001.fastq.gz | head -2 | tail -1 | wc -c
# 9
```

The read length is 101 for the read files, and 8 for the indices (the count above includes the newline at the end)

```bash
zcat 1294_S1_L008_R1_001.fastq.gz | head -4 | tail -1 | wc -c
```

The `#` character in the score line corresponds with the `N` character in the read, which should have a score of 2. This is consistent with Phred+33 encoding'

Created `dists.py` to create a histogram of average quality score per base for each file. Experimented with the `multiprocessing` module, which ran but the output file implies that the processes didn't run in parallel, so I don't know if I'll try multithreading/multiprocessing on the main script.

The minimum hamming distance between any 2 indices is 3. I'm planning on doing error correction (I think just correcting up to 2 Ns), so I can't just throw out all reads with an N in the index. I think a good cutoff would be if an index has an average quality score below 25. At 25, the probability of error is ~0.00316. The probability that 3 or more bases in a sequence of 8 would be misread is ~1.75x10^-6. There are ~360 million records, so even with these odds we'd expect a fair number of such events, but only 18 of 576 index combinations have a hamming distance of 3, so the odds that the base call errors would occur in the right index, in the right locations of the index, and be changed to a base that matches another index seem astronomically low. Even if 2 Ns have been corrected, the odds that a third base has had an error that makes it match another index seem very low. I also want to balance this consideration with wanting to keep as much data as possible. If there are 2 Ns in a sequence of 8, the maximum average quality score possible is 30.5 (if the other 6 have scores of 40), so a cutoff of 30 would make it fairly unlikely that any indices with more than 1 N are saved, since we would need all the other bases to be effectively perfect. With all this in mind, I think 25 is a reasonable cutoff to avoid losing too much data while also keeping the likelihood of error low.

Bash commands for counting how many indices have Ns:
```bash
zcat 1294_S1_L008_R2_001.fastq.gz | grep -A 1 -E '^@' | grep -v -E '^@' | grep -c 'N'
# 3976613
zcat 1294_S1_L008_R3_001.fastq.gz | grep -A 1 -E '^@' | grep -v -E '^@' | grep -c 'N'
# 3328051

zcat 1294_S1_L008_R2_001.fastq.gz | grep -A 1 -E '^@' | grep -v -E '^@' | grep -c -E 'N.*N'
```

Copied above info to `Answers.md`



Created [demultiplex.py](demultiplex.py) in accordance with pseudocode outlined in `Strategy.md` in `Assignment-the-first`. Output is a simple tsv file called `summary.tsv` where the first 3 lines are the counts for paired, unknown, and hopped records, and the rest of the lines are the individual counts for each index pair.

Tested `demultiplex.py` using test files. Initially failed because the R3 test files had the forward indices, not the reverse complements. One of the records that was supposed to be hopped was unknown because its quality score was below the threshold, which I didn't account for when I made the tests. Other than that the records went where they were supposed to.

Ran `demultiplex.py` with `run_demultiplex.sh` using a quality cutoff of 25. Output is in `/scratch/bgmp/imre/demux` on Talapas. Finished in ~37 minutes, used 82% of one CPU, and had a maximum RAM usage of ~250 MB. Copied `summary.tsv` to `Assignment-the-third/`.