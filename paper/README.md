
Requires the following files:
```python
# See generate_config.smk:L4

# TSV with columns: clade, sm, chrom, contig.
CLADE_KEY = "clade_key.tsv"
# Reference fasta directory
# {REF_FA_DIR}/{sm}/{sm}_{chrom}_{contig}.fa
REF_FA_DIR = "/project/logsdon_shared/projects/HGSVC3/HGSVC_centromere_annotation/buildTree/HGSVC3_ref/"
# Reference bed directory
# Contains columns: chrom, st, end, time, sm, chrom_name, clade_n, full_chrom
# {REF_BED_DIR}/{clade}.bed
REF_BED_DIR = "/project/logsdon_shared/projects/HGSVC3/Snakemake-MutationRate/SG_working/tree_based_mutation_rate/ref_clades"
# Query fasta directory
# {QRY_FA_DIR}/{clade}/{qry}.fa
QRY_FA_DIR = "/project/logsdon_shared/projects/HGSVC3/Snakemake-MutationRate/SG_working/tree_based_mutation_rate/others"
# Non-human primate directory
# {NHP_FA_DIR}/{species}/{chrom}_{hap}_{hchrom}.fa
NHP_FA_DIR = "/project/logsdon_shared/projects/HGSVC3/Snakemake-MutationRate/data/cens_primates"
# Base config template per clade.
TEMPLATE = "config_template.yaml"
# Sample colors
SAMPLE_COLOR_CODE = "sample_colors.tsv"
# Plot format.
# If annotations change, must replace path with name abbreviation used in generate_config.py:L16. ex. cdr
PLOT_FORMAT = "base.toml"

# Annotations where ref in BED chrom column
# See generate_config.py:L16
DATA_SOURCES = {
    # BED4
    "hor_arr": lambda chrom_name: "/project/logsdon_shared/projects/HGSVC3/HGSVC_centromere_annotation/all_array_length.bed",
    # BED9
    "rm": lambda chrom_name: f"/project/logsdon_shared/projects/HGSVC3/HGSVC_centromere_annotation/RM/all_cens_{chrom_name}.annotation.fa.out",
    # BED9
    "hor": lambda chrom_name: f"/project/logsdon_shared/projects/HGSVC3/HGSVC_centromere_annotation/HOR/{chrom_name}_AS-HOR_stv_row.all.bed",
    # BED3
    "cdr": lambda chrom_name: "/project/logsdon_shared/projects/HGSVC3/HGSVC_centromere_annotation/all_cdrs.bed",
    # BED4
    "mdp": lambda chrom_name: "/project/logsdon_shared/projects/HGSVC3/HGSVC_centromere_annotation/1DModplot.bed",
}
```

Generate config.
```bash
snakemake -p -s paper/generate_config.smk -c 12 -d paper/
```

Run mutation rate workflow per clade.
```bash
snakemake -np -s paper/run_multiple_mut.smk -c 12 --config input_dir=paper/config
```
