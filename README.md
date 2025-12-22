# Mutation Rate Analysis

Calculate mutation rate of query regions against a reference region.

![](docs/images/chr1@6.png)

## Getting Started
Clone the repo and setup Snakemake.
```bash
git clone https://github.com/logsdon-lab/Snakemake-MutationRate.git --recursive
cd Snakemake-MutationRate

# Load or install apptainer
# module load apptainer
which apptainer
conda env create --name smk-mut-rate -f env.yaml
conda activate smk-mut-rate
```

## Configuration
```yaml
# Alignment memory.
mem: 100GB
# Alignment threads.
threads: 8
# Directories
output_dir: "results"
log_dir: "logs"
benchmarks: "benchmarks"
# Reference config.
reference:
  - name: reference
    # Sample fasta.
    path: reference.fasta
    # Chromosome specific bed files of windows to evaluate. Filename is region name.
    # Each window must have a chromosome name
    # Can be created from a region bedfile via, bedtools makewindows -w "${window_size}" -b "${bedfile}"
    bed:
      - reference_chr1.bed
    # Calculate rate for this region relative to everything else not in bed.
    bed_comparison: reference.bed
    # Reference region cenplot format configuration
    # The mutation rate plot is added as the final track with legend elements for each track added in relative order.
    plot_format: annotation.toml
# Regular expression patterns within fasta headers to find haplotype and matching chromosome.
# Used to group and filter alignments by exact match.
regex_sm_hap: "'(mat|pat|haplotype1|haplotype2|hap1|hap2|h1|h2)'"
regex_sm_chrom: "'_(hsa(.*?)|chr(.*?)|cen(.*?))[:_]'"
regex_ref_chrom: "(?:chr|hsa|cen)([0-9XY]{1,2})"
# Samples to align against each reference.
samples:
  - name: query
    # Query fasta
    path: query.fasta
    # Query bedfile. If none, entire region used.
    bed: query.bed
    # Minimap2 parameters. Important.
    # Preset of asm20 recommmended if highly divergent region.
    mm2_opts: -x asm20 --secondary=no -K 8G -t 8
    # Assumed divergence time from reference.
    divergence_time: 0
```

## Input
The main inputs are:
* Reference fasta (`reference.path`)
* Reference BED coordinates (`reference.bed`)
* Query fasta (`samples[n].path`)
* Query BED coordinates (`samples[n].bed`)

> [!NOTE]
> Should include chromosome and haplotype information in fasta names to filter alignments with `config.regex_*`.
> ex. `HG00096_chr1_haplotype1-00000X`

## Output
The main output files have the following prefix: `results/{clade}/mutation_rate/{contig_name}/{clade}`

|file|desc|
|-|-|
|`{clade}_annot.bedpe`|BEDPE-like file with paired intervals and their mutation rate. Includes header with columns: `ref_chrom,ref_st,ref_end,mu,qry_chrom,qry_st,qry_end,ref_sm,qry_sm,category` |
|`{clade}_mut.tsv`|TSV file with summary of mutation rates with format: `ref_chrom,clade,cmp_mutation_rate,non_cmp_mutation-rate,fold_mutation_rate,`|
|`{clade}.pdf`|PDF of mutation rate plot.|
|`{clade}.png`|PNG of mutation rate plot.|

## Usage
```bash
snakemake -p --configfile config/config.yaml -j 40
```

## Test
```bash
# Rerun mutation rate estimation.
# https://www.nature.com/articles/s41586-024-07278-3#Sec6
snakemake -p --configfile test/config/config_cen_var_glogsdon.yaml -c 40 -n
```

## Sources
* https://www.nature.com/articles/s41586-024-07278-3#Sec6.
