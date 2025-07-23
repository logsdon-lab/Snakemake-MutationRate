import os


CLADE_KEY = "clade_key.tsv"
# /project/logsdon_shared/projects/HGSVC3/HGSVC_centromere_annotation/buildTree/HGSVC3_ref/{sm}/{sm}_{chrom}_{contig}.fa
REF_FA_DIR = "/project/logsdon_shared/projects/HGSVC3/HGSVC_centromere_annotation/buildTree/HGSVC3_ref/"
# /project/logsdon_shared/projects/HGSVC3/Snakemake-MutationRate/SG_working/tree_based_mutation_rate/ref_clades/{clade}.bed
REF_BED_DIR = "/project/logsdon_shared/projects/HGSVC3/Snakemake-MutationRate/SG_working/tree_based_mutation_rate/ref_clades"
# /project/logsdon_shared/projects/HGSVC3/Snakemake-MutationRate/SG_working/tree_based_mutation_rate/others/{clade}/{qry}.fa
QRY_FA_DIR = "/project/logsdon_shared/projects/HGSVC3/Snakemake-MutationRate/SG_working/tree_based_mutation_rate/others"
# /project/logsdon_shared/projects/HGSVC3/Snakemake-MutationRate/data/cens_primates/{species}/{chrom}_{hap}_{hchrom}.fa
NHP_FA_DIR = (
    "/project/logsdon_shared/projects/HGSVC3/Snakemake-MutationRate/data/cens_primates"
)
# Base config template per clade.
TEMPLATE = "config_template.yaml"
# Sample colors
SAMPLE_COLOR_CODE = "sample_colors.tsv"
# Plot format.
# If annotations change, must replace path with name abbreviation used in generate_config.py:L16. ex. cdr
PLOT_FORMAT = "base.toml"

with open(CLADE_KEY) as fh:
    # clade sm chrom contig
    REF_FA = {}
    for line in fh:
        clade, sm, chrom, contig = line.strip().split()
        REF_FA[clade] = os.path.join(REF_FA_DIR, sm, f"{sm}_{chrom}_{contig}.fa")


rule generate_config:
    input:
        script="generate_config.py",
        template=TEMPLATE,
        color_code=SAMPLE_COLOR_CODE,
        plot_format=PLOT_FORMAT,
        ref_fa=lambda wc: REF_FA[wc.clade],
        # Expects format:
        # chrom, st, end, time, sm, chrom_name, clade_n, full_chrom
        # HG02953_chr1_haplotype1-0000013	16773102	25725928	0.07	HG02953	chr1	3	HG02953_chr1_haplotype1-0000013:17273102-25225928
        ref_bed=os.path.join(REF_BED_DIR, "{clade}.bed"),
        qry_fa_dir=os.path.join(QRY_FA_DIR, "{clade}"),
        nhp_fa_dir=NHP_FA_DIR,
    output:
        config="config/config_{clade}.yaml",
    shell:
        """
        python {input.script} \
        {input.template} \
        {input.ref_fa} \
        {input.ref_bed} \
        {input.qry_fa_dir} \
        {input.nhp_fa_dir} \
        {input.color_code} \
        {input.plot_format} > {output}
        """


rule all:
    input:
        expand(rules.generate_config.output, clade=REF_FA.keys()),
    default_target: True
