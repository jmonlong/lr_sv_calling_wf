# Tiny test dataset

Unrealistic dataset, just to check that the workflow run locally.

Prepare dataset with:

```sh
python3 sim_tiny_test_data.py
gzip -f reads.fastq
samtools import reads.fastq.gz -o reads.bam -O BAM
```
