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

Created `dists.py` to create a histogram of average quality score per base for each file. Experimented with multithreading using the `multiprocessing` module, which seemed to work.

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