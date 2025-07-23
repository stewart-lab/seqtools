# seqtools
dockerized tools for tasks related to RNA-Seq

## To run STAR aligner
1. In rnatools/build folder, modify run_data.sh
```
-v /w5home/bmoore/sequencing:/mnt/sequencing # <your data directory>:<data directory on docker> (mounts volume dir to docker)
--fastq_dir /mnt/sequencing/01.RawData # where your fastqs are located
--reference_genome_dir /mnt/sequencing/Homo_sapiens_CRCh38 # where your reference genome is located: need gtf and primary fasta file
--output_dir /mnt/sequencing/02.alignment # directory to put the alignment output
--resume # flag to resume if stopped before alignment is finished
```
2. Run run_data.sh
```
source run_data.sh
```

## To run DEseq2
1. In deseq2 folder, modify run_deseq.sh file with options:
```
-v /w5home/bmoore/data/:/mnt/data # <your data directory>:<data directory on docker> (mounts volume dir to docker)
--design /mnt/data/design.csv # design matrix with Sample and Condition columns
--counts /mnt/data/counts.csv # counts file from mapping with genes as rows and samples as columns
--contrasts /mnt/data/contrasts.csv # what to compare with columns Treatment and Control (should be same values as what is in the design matrix Condition column
--output_dir /mnt/data/ # where you want output to go: same as data directory on docker
```
2. Run run_deseq.sh
```
source run_deseq.sh
```

## Running GO enrichment
NOTE: Due to some library conflicts, we have not been able to get this to run on Docker, so you need to install packages
