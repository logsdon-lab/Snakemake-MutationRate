# Mutation Rate Analysis
See https://www.nature.com/articles/s41586-024-07278-3#Sec6.

# Setup

> TODO

Extend coordinates 1.5 Mbp out from central HOR array.
```bash
bedtools slop -i /project/logsdon_shared/projects/CHM13_CDR-Finder/regions.bed -g data/T2T-CHM13v2.fasta.fai -b 1000000 > data/regions.bed 
```

Modify configfile.

Run.
```bash
snakemake -p --configfile config.yaml -j 40 --workflow-profile ~/profiles/lpc/ --sdm conda apptainer -n
```