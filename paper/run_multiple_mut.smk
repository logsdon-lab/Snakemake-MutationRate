import os


INPUT_CONFIG_DIR = config["input_dir"]
CFG_GLOB = os.path.join(INPUT_CONFIG_DIR, "config_{clade}.yaml")
WCS = glob_wildcards(CFG_GLOB)


rule run_mutation_rate:
    input:
        cfg=CFG_GLOB,
    output:
        touch("checkpoints/{clade}.done"),
    resources:
        jobs=30,
    log:
        "logs_ref_clade/{clade}.log",
    threads: 1
    shell:
        """
        snakemake -pk --configfile {input.cfg} -j {resources.jobs} --workflow-profile workflow/profiles/lpc &> {log}
        """


rule all:
    input:
        expand(rules.run_mutation_rate.output, clade=WCS.clade),
    default_target: True
